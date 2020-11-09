#!/usr/bin/python3

import argparse
import yaml
import sys
sys.path.insert(1, './benchmarks/')

def get_args():
    """Parse the input parameters for benchmarking"""

    description =  "proofreader for yaml, checking that all neccessary settings are set"
    yml_file_help = "YML file name"

    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "config_file",
        type=str,
        help="""The config file""")

    args = argparser.parse_args()

    with open(args.config_file, "r") as yamlconfig:
        config = yaml.safe_load(yamlconfig)
    return config

config = get_args()

def parse_config(config):
    benchmarks = []
    for benchmark_settings in config:
        benchmark_settings_section = list(benchmark_settings.keys())[0]
        for benchmark in "general_settings", "timed_command", "thread_scaling", "memory_footprint", "settings_optimizer":
            if benchmark_settings_section == benchmark:
                exec('import {}'.format(benchmark))
                benchmarks.append(eval(benchmark).parse_settings(benchmark_settings))
    #if list item starts with X, send this section to benchmark script for X
    return benchmarks

print(parse_config(config))


# for thing in thing:
#run thing with settings_dict
#read config

#run tests if specified in config

#
