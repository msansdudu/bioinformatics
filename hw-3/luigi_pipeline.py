import luigi
import subprocess
import os
import re

class IndexReference(luigi.Task):
    ref_fasta = luigi.Parameter(default="data/ecoli.fa")

    def output(self):
        return luigi.LocalTarget(f"{self.ref_fasta}.bwt")

    def run(self):
        subprocess.run(["bwa", "index", self.ref_fasta], check=True)


class RunFastQC(luigi.Task):
    fastq1 = luigi.Parameter(default="data/SRR390728_1.fastq.gz")
    fastq2 = luigi.Parameter(default="data/SRR390728_2.fastq.gz")

    def output(self):
        out1 = self.fastq1.replace(".fastq.gz", "_fastqc.html")
        out2 = self.fastq2.replace(".fastq.gz", "_fastqc.html")
        return [luigi.LocalTarget(out1), luigi.LocalTarget(out2)]

    def run(self):
        subprocess.run(["fastqc", self.fastq1, self.fastq2], check=True)


class MapReads(luigi.Task):
    ref_fasta = luigi.Parameter(default="data/ecoli.fa")
    fastq1 = luigi.Parameter(default="data/SRR390728_1.fastq.gz")
    fastq2 = luigi.Parameter(default="data/SRR390728_2.fastq.gz")

    def requires(self):
        return IndexReference(ref_fasta=self.ref_fasta)

    def output(self):
        return luigi.LocalTarget("data/sample_luigi.sam")

    def run(self):
        with self.output().open('w') as f:
            subprocess.run(["bwa", "mem", self.ref_fasta, self.fastq1, self.fastq2], stdout=f, check=True)


class ViewBam(luigi.Task):
    def requires(self):
        return MapReads()

    def output(self):
        return luigi.LocalTarget("data/sample_luigi.bam")

    def run(self):
        with self.output().open('w') as f:
            subprocess.run(["samtools", "view", "-S", "-b", self.input().path], stdout=f, check=True)


class RunFlagstat(luigi.Task):
    def requires(self):
        return ViewBam()

    def output(self):
        return luigi.LocalTarget("data/sample_luigi_flagstat.txt")

    def run(self):
        with self.output().open('w') as f:
            subprocess.run(["samtools", "flagstat", self.input().path], stdout=f, check=True)


class EvaluateMapping(luigi.Task):
    def requires(self):
        return RunFlagstat()

    def output(self):
        return luigi.LocalTarget("data/evaluation.txt")

    def run(self):
        with self.input().open('r') as f:
            content = f.read()
        
        match = re.search(r'mapped \(([\d\.]+)%', content)
        if match:
            percent = float(match.group(1))
            status = "OK" if percent >= 0.1 else "not OK"
            with self.output().open('w') as f:
                f.write(f"Mapped: {percent}%\nStatus: {status}\n")
            if status != "OK":
                raise ValueError("Mapping percentage is too low!")
        else:
            raise ValueError("Could not parse mapped percentage.")


class SortBam(luigi.Task):
    def requires(self):
        return {"bam": ViewBam(), "eval": EvaluateMapping()}

    def output(self):
        return luigi.LocalTarget("data/sample_luigi.sorted.bam")

    def run(self):
        subprocess.run(["samtools", "sort", self.input()["bam"].path, "-o", self.output().path], check=True)


class CallVariants(luigi.Task):
    ref_fasta = luigi.Parameter(default="data/ecoli.fa")

    def requires(self):
        return SortBam()

    def output(self):
        return luigi.LocalTarget("data/sample.vcf")

    def run(self):
        with self.output().open('w') as f:
            subprocess.run(["freebayes", "-f", self.ref_fasta, self.input().path], stdout=f, check=True)


class RunVariantCalling(luigi.WrapperTask):
    """
    This is the main entry point to run the entire pipeline.
    It requires the final variant calling output and FastQC reports.
    """
    def requires(self):
        yield CallVariants()
        yield RunFastQC()

if __name__ == '__main__':
    luigi.run()
