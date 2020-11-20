#!/usr/bin/python3

import argparse
import yaml
import base_benchmark
import time
import subprocess
import sys

class thread_scaling(base_benchmark.base_benchmark):
    def __init__(self):
        base_benchmark.base_benchmark.__init__(self)
        self.name = __class__.__name__

    def run_benchmark(self):
        print("running " + self.get_name())
        print(self.settings)
        command = self.settings['command']
        threads_list = self.settings['threads']
        results = {}
        for threads in threads_list:
            command.format(threads=threads)
            print(command)
            start = time.time()
            with subprocess.Popen(command, stdout=sys.stderr, stderr=subprocess.PIPE, universal_newlines=True, shell=True) as process:
                _, stderr = process.communicate()
                if process.returncode != 0:
                    print(stderr.splitlines()[:-1], file=sys.stderr)
                    print("fucked it.")
            end = time.time() - start
            print(end)
            results[threads] = end
        return end
