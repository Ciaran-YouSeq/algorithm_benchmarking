#!/usr/bin/python3

import os
import argparse
import yaml
import base_benchmark
import time
import subprocess
import sys

class alignment_quality(base_benchmark.base_benchmark):
    def __init__(self):
        base_benchmark.base_benchmark.__init__(self)
        self.name = __class__.__name__

    def run_benchmark(self):
        print("running " + self.get_name())
        print(self.settings)
        with open("benchmarks/params_file.txt",'r') as params_file:
            params = params_file.read()
        with open("nextflow.config", 'w') as config:
            config.write(params)
            config.write('params.aligned_bam = "'+self.settings['alignment_file']+'"')
#        command = "./nextflow benchmarks/quality_measurement_pipe.nf"
        command = "./nextflow benchmarks/full_pipe.nf"
        with subprocess.Popen(command, stdout=sys.stderr, stderr=subprocess.PIPE, universal_newlines=True, shell=True) as process:
                _, stderr = process.communicate()
                if process.returncode != 0:
                    print(stderr.splitlines()[:-1], file=sys.stderr)
                    print("job failed.")
