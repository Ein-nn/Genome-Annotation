# -*- coding: utf-8 -*-
# @Author: EIN
# @Date:   2020-09-11 10:44:40
# @Last Modified by:   SSOJK
# @Last Modified time: 2020-09-13 18:22:48
import os
from FilterGff3 import filter

p = 100
times = 1
flag = "F"
fi = filter()
gff = "augustus_out_0-0.gff3"

# 从blast结果gff中提取序列并预测、筛选，并判断预测结果是否完整
os.system("Rscript SeqFetch_0.R 0 0 " + gff)
gff = "augustus_out_" + str(times) + "-0.gff3"
# augustus --gff3=on --outfile=Sc_augustus_out.gff3 --species=saccharomyces_cerevisiae_S288C S288C.fna
os.system("augustus --gff3=on --species=saccharomyces_cerevisiae_S288C --outfile=" + gff + " extend_0bp.fasta")
[flag, gff, incomplete_num] = fi.filter_au(gff, times, 0)
with open("incomplete_num.txt", "a") as n:
    n.write("extend\tincomplete_num\n0bp\t" + str(incomplete_num) + "\n")

while flag == "T" or times <= 80:
    # call R script to extract sequences
    os.system("Rscript SeqFetch.R 100 " + str(times) + gff)
    # run to predict
    times += 1
    gff = "augustus_out_" + str(times) + "-0.gff3"
    fasta = " extend_" + str(times*p) + "bp.fasta"
    os.system("augustus --gff3=on --species=saccharomyces_cerevisiae_S288C --outfile=" + gff + fasta)
    # filter incomplete
    [flag, gff, incomplete_num] = fi.filter_au(gff)
    with open("incomplete_num.txt", "a") as n:
        n.write(str(times*p) + "bp\t" + str(incomplete_num) + "\n")
