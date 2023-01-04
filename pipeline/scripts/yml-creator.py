#!/usr/bin/env python
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-d','--default', dest='defyaml', default=None)
parser.add_argument('-A','--pathA', dest='path1', default=None)
parser.add_argument('-B','--pathB', dest='path2', default=None)
args = parser.parse_args()

if '[' in args.path1:
  A = re.sub(r'\]', '', re.sub(r'\[', '', args.path1))
  B = re.sub(r'\]', '', re.sub(r'\[', '', args.path2))
else:
  A = args.path1
  B = args.path2

## defualts with {} for filling
with open(args.defyaml, 'r') as file:
	content = file.read().format(A, B)

print(content)
