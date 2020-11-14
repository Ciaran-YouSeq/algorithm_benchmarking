#!/usr/bin/python3

import argparse
import yaml
import base_benchmark

class settings_optimizer(base_benchmark.base_benchmark):
    def __init__(self):
        base_benchmark.base_benchmark.__init__(self)

    def run(self):
        return "running" + get_name(self)
        
