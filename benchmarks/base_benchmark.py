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
        self.settings = list(settings)[0]
        print(self.settings)

