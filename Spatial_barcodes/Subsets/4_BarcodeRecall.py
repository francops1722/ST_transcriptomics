#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 09:57:40 2022

@author: fpomasot
"""

# DECODING Step 3a (for simulations done within this notebook).
# If there is an answers file, compute precision and recall of the decode.
import argparse
import yaml
from numpy import *

parser = argparse.ArgumentParser()
parser.add_argument('--config',required=True, help='path to yaml config file')
yaml_file = parser.parse_args()

with open(yaml_file.config, "r") as fp:
        args = yaml.safe_load(fp)

Qq = args["Q"]
best_path = args["best_path"]
answersfilename = args["answersfilename"]#user set to output file path
answers = genfromtxt(answersfilename, dtype=int)
best = loadtxt(best_path)
ngood = sum(logical_and(answers[:Qq] == best[:Qq], best[:Qq] >= 0))
nbad = sum(logical_and(answers[:Qq] != best[:Qq], best[:Qq] >= 0))
nerase = sum(best[:Qq] < 0)
recall = 1. - nerase / Qq
precision = ngood/(ngood + nbad + 1.e-30)

print(f"In {Qq} decodes, {ngood} were correct, {nbad} were wrong, {nerase} were erasures.")
print(f"precision = {precision:.4f}, recall = {recall:.4f}.")

