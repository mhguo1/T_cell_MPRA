#!/usr/bin/python

#This file will force the variant with the lower chromosome position to always be the first variant. The purpose of this script is to help remove duplicate entries in the cache file
import sys
import gzip

for line_c1 in sys.stdin:
        line_c=line_c1.rstrip().split('\t')
        pos1=int(line_c[0].split(":")[1])
        pos2=int(line_c[1].split(":")[1])
        if pos1<=pos2:
                sys.stdout.write('\t'.join(line_c)+"\t"+str(pos1)+"\t"+str(pos2)+"\n")
        else:
                a1=line_c[4].split(",")[0][0]
                a2=line_c[4].split(",")[1][0]
                sys.stdout.write(line_c[1]+"\t"+line_c[0]+"\t"+line_c[2]+"\t"+str(line_c[3])+"\t"+a2+","+a1+"\t"+str(pos2)+"\t"+str(pos1)+"\n")
