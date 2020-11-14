#!/usr/bin/python3

import argparse
import yaml

class base_benchmark:
    def __init__(self, settings=None):
        print("HERE at init")

    def get_name(self):
        return __class__.__name__

    def parse_settings(self, settings):
        print("here at parse")
        print(settings)


