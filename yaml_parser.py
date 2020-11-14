#!/usr/bin/python3

import argparse
import yaml
import os
import json
from datetime import datetime
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
    """load benchmarks listed in config from ./benchmarks/"""
    loaded_benchmarks = []
    for benchmark_settings_section in config.keys():
        for benchmark in benchmarks:
            if benchmark_settings_section == benchmark:
                exec('import {}'.format(benchmark))
                loaded_benchmarks.append(eval("{benchmark}.{benchmark}()".format(benchmark=benchmark)))
    return loaded_benchmarks


def parse_settings(config, loaded_benchmarks):
    """parse settings from config, and send them to the relevant loaded benchmark class objects"""
    print(config)
    for benchmark in loaded_benchmarks:
        parsed_settings_benchmark_objects = benchmark.parse_settings(config[benchmark.get_name()])
    return parsed_settings_benchmark_objects

def run_benchmarks(loaded_benchmarks):
    """run all loaded benchmark objects and output results to a dictionary"""
    results = {}
    for benchmark in loaded_benchmarks:
        results[benchmark.get_name()] = benchmark.run(benchmark.run(benchmark.run_benchmark()))
    return results

def write_results(config, results):
    """json dump results dictionary to date+time stamped file tagged with program name and version"""
    #join config and results dicts together, then json dump
    date_and_time_stamp = datetime.now().strftime("%d_%m_%Y_%H_%M_%S")
    algorithm_name = config['general_settings']['algorithm_name']
    file_name = algorithm_name+"_"+date_and_time_stamp+".json"
    with open(file_name, 'w') as json_file:
        json.dump(results, json_file)
    


if __name__ == "__main__":
    print("get benchmarks")
    benchmarks = get_benchmarks()
    print("loaded benchmarks")
    loaded_benchmarks = load_benchmarks(config, benchmarks)
    print("parsed benchmarks")
    parse_settings = parse_settings(config, loaded_benchmarks)
    print("run benchmarks")
    results = run_benchmarks(loaded_benchmarks)
    print("write results")
    write_results(config, results)
# for thing in thing:
#run thing with settings_dict
#read config

#run tests if specified in config

#
