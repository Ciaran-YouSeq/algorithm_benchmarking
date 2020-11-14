#!/usr/bin/python3

import argparse
import yaml
import os
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

def get_benchmarks():
    """find all benchmarks"""
    benchmarks = []
    for benchmark in os.listdir("./benchmarks/"):
        if ".py" in benchmark:
            benchmarks.append(benchmark.split(".")[0])
    return benchmarks


def load_benchmarks(config, benchmarks):
    """parse config file, and load benchmarks from ./benchmarks/"""
    loaded_benchmarks = []
    for benchmark_settings in config:
        benchmark_settings_section = list(benchmark_settings.keys())[0]
        print(benchmark_settings)
        for benchmark in benchmarks:
            if benchmark_settings_section == benchmark:
                exec('import {}'.format(benchmark))
                #loaded_benchmarks.append(eval(benchmark).parse_settings(benchmark_settings))
                loaded_benchmarks.append(exec("{benchmark}.{benchmark}".format(benchmark=benchmark)))
                #loaded_benchmarks.append(exec("{benchmark}.{benchmark}.parse_settings({benchmark_settings})".format(benchmark=benchmark, benchmark_settings=benchmark_settings)))
    #if list item starts with X, send this section to benchmark script for X
    return loaded_benchmarks

def parse_settings(config, loaded_benchmarks):
    parsed_settings = []
    print(config)
    for benchmark_settings in config:
        benchmark_settings_section = list(benchmark_settings.keys())[0]
    for benchmark in loaded_benchmarks:
        parsed_settings.append(exec("{benchmark}.parse_settings({benchmark_settings})".format(benchmark=benchmark, benchmark_settings=benchmark_settings)))


def run_benchmarks(loaded_benchmarks):
    print(loaded_benchmarks)
    for benchmark in loaded_benchmarks:
        print(benchmark)
        exec("{}.run(settings)".format(benchmark.keys()[0]))

#{'general_settings': {'algorithm_name': 'bwa', 'algorithm_version': '1.0.0'}}

print("get benchmarks")
benchmarks = get_benchmarks()
print("loaded benchmarks")
loaded_benchmarks = load_benchmarks(config, benchmarks)
print("parsed benchmarks")
parse_settings = parse_settings(config, loaded_benchmarks)
print("run benchmarks")
run_benchmarks(loaded_benchmarks)

# for thing in thing:
#run thing with settings_dict
#read config

#run tests if specified in config

#
