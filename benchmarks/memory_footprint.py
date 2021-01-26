#!/usr/bin/python3

import argparse
import yaml
import base_benchmark
import subprocess
import sys

class memory_footprint(base_benchmark.base_benchmark):
    def __init__(self):
        base_benchmark.base_benchmark.__init__(self)
        self.name = __class__.__name__

    def run_command(self, command):
        with subprocess.Popen(command, stdout=sys.stderr, stderr=subprocess.PIPE, universal_newlines=True, shell=True) as process:
                stdout, stderr = process.communicate()
                if process.returncode != 0:
                    print(command)
                    exit()

    def run_benchmark(self):
        print("running " + self.get_name())
        print(self.settings)
        if "threads" in self.settings.keys():
            threads = self.settings['threads']
        else:
            threads = 32
        command = "valgrind --tool=massif --massif-out-file=massif.out "+self.settings['command'].format(threads=threads, output=self.settings['output'])
        print(command)
        self.run_command(command)
        command = "/snap/bin/ms_print massif.out > massif.txt"
        self.run_command(command)
        with open("massif.txt", "r") as massif_output:
            massif_output = massif_output.readlines()
            data_extract = []
            for line in massif_output:
                if not line.startswith("|") and not line.startswith("->") and not line.endswith("etc.\n"):
                    data_extract.append(line)
        return data_extract
        
