#!/usr/bin/python3

import argparse
import yaml

class general_settings:
    def __init__(self, settings=None):
        print("HERE at init")
        pass

    def get_name(self):
        print(__class__.__name__)

    def parse_settings(self, settings):
        print("here at parse")
        print(settings)
