#!/usr/bin/python3

import argparse
import yaml
import os
import sys
sys.path.insert(1, './benchmarks/')

def get_args():
    """Parse the input yaml, and store settings as a nested dictionary"""

    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "config_file",
        type=str,
        help="""The config file""")

    args = argparser.parse_args()
    
    config = {}

    with open(args.config_file, "r") as yamlconfig:
        yamlconfig = yaml.load(yamlconfig)
        for setting in yamlconfig:
            #config[list(setting.keys())[0]]['command'] = list(yamlconfig)[0]['general_settings']['command']
            config[list(setting.keys())[0]] = list(setting.values())[0]
            config[list(setting.keys())[0]]['command'] = list(yamlconfig)[0]['general_settings']['command']
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
    for benchmark_settings_section in config.keys():
        for benchmark in benchmarks:
            if benchmark_settings_section == benchmark:
                exec('import {}'.format(benchmark))
                loaded_benchmarks.append(eval("{benchmark}.{benchmark}()".format(benchmark=benchmark)))
    return loaded_benchmarks


def parse_settings(config, loaded_benchmarks):
    print(config)
    for benchmark in loaded_benchmarks:
        parsed_settings_benchmark_objects = benchmark.parse_settings(config[benchmark.get_name()])
    return parsed_settings_benchmark_objects

def run_benchmarks(loaded_benchmarks):
    for benchmark in loaded_benchmarks:
        benchmark.run()

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
