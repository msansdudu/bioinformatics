import luigi

class HelloWorldTask(luigi.Task):
    """A minimal 'Hello world' Luigi task."""

    def output(self):
        return luigi.LocalTarget("hello_world.txt")

    def run(self):
        with self.output().open('w') as f:
            f.write("Hello world from Luigi!\n")

if __name__ == '__main__':
    luigi.run()
