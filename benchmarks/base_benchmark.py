#!/usr/bin/python3

import argparse
import yaml

class base_benchmark:
    def __init__(self, settings=None):
        print("HERE at init")

    def get_name(self):
        return self.name

    def parse_settings(self, settings):
        print("here at parse")
        self.settings = settings
        print(self.settings)

    def run(self, run_benchmark_function):
        results = {}
        repeats = self.settings['repeats']
        for repeat in range(repeats):
            results["repeat_{repeat}".format(repeat=repeat+1)] = eval("self.run_benchmark()")
        return results
