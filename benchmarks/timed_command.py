#!/usr/bin/python3

import argparse
import yaml
import base_benchmark

class timed_command(base_benchmark.base_benchmark):
    def __init__(self):
        base_benchmark.base_benchmark.__init__(self)
        self.name = __class__.__name__

    def run(self):
        return "running" + get_name(self)
        
