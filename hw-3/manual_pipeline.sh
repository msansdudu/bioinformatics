#!/bin/bash
# manual_pipeline.sh
# Bash script to run the variant calling algorithm steps manually

echo "1. Indexing reference..."
bwa index data/ecoli.fa

echo "2. Running FastQC..."
fastqc data/SRR390728_1.fastq.gz data/SRR390728_2.fastq.gz

echo "3. Mapping reads..."
bwa mem data/ecoli.fa data/SRR390728_1.fastq.gz data/SRR390728_2.fastq.gz > data/sample.sam

echo "4. Converting SAM to BAM..."
samtools view -S -b data/sample.sam > data/sample.bam

echo "5. Generating Flagstat..."
samtools flagstat data/sample.bam > data/sample_flagstat.txt

echo "6. Parsing mapped percentage..."
mapped_percent=$(python3 parse_flagstat.py data/sample_flagstat.txt)
echo "Mapped percentage is: $mapped_percent%"

echo "7. Evaluating quality..."
if (( $(echo "$mapped_percent >= 0.1" | bc -l) )); then
    echo "OK"
    echo "8. Sorting BAM..."
    samtools sort data/sample.bam -o data/sample.sorted.bam
    echo "9. Calling variants..."
    freebayes -f data/ecoli.fa data/sample.sorted.bam > data/sample.vcf
    echo "Finished"
else
    echo "not OK"
fi
