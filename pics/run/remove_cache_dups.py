#!/usr/bin/python

#This script will remove any duplicate entries in the cache file
import sys
import gzip

pos1="a0"
pos2="a1"

for line_c1 in sys.stdin:
        line_c=line_c1.rstrip().split('\t')
        if line_c[0]!=pos1 or line_c[1]!=pos2:
                sys.stdout.write('\t'.join(line_c[0:5])+"\n")
        pos1=line_c[0]
        pos2=line_c[1]
