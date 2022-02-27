#!/usr/bin/python

import sys
fasta=open("mpra_input_snps_ref.fa", "r")
outfile=open("mpra_input_snps_alt.fa", "w")

for line_f1 in fasta:
        line_f=line_f1.strip()
        if line_f[0]==">":
                name=line_f.replace(">", "")
                ref=name.split(":")[2]
                alt=name.split(":")[3]
                outfile.write(line_f1)
        else:
                str_out=line_f[:15]+alt+line_f[16:]
                outfile.write(str_out+"\n")
outfile.close()
fasta.close()
