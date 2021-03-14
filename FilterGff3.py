# -*- coding: utf-8 -*-
# @Author: EIN
# @Date:   2020-09-10 11:14:33
# @Last Modified by:   SSOJK
# @Last Modified time: 2020-09-14 18:36:38

import re
# for Augustus
# filename = augustus_out_1-0.gff3
class filter():
    """docstring for filter
    """
    def __init__(self):
        self.flag = "F"
        self.out_gff = ""
        self.pass_gff = ""
        self.gene_dict = {}
        self.incomplete_num = 0

    def whether_complete(self, gene, extend_len):
        items = gene.split("\n")
        initial_start = int(items[1].split("\t")[0].split(";")[2])
        [now_start, now_end] = [int(items[1].split("\t")[3]), int(items[1].split("\t")[4])]
        new_start = initial_start - extend_len + now_start - 1
        new_end = initial_start - extend_len + now_end - 1
        gene = re.sub(r"\t[0-9]+\t[0-9]+\t", "\t"+str(new_start)+"\t"+str(new_end)+"\t", gene)
        if "start_codon" in gene and "stop_codon" in gene:
            with open(self.pass_gff, "a") as p:
                p.write(re.sub(r";gene.+;[\+\-]","",gene))
        else:
            self.flag = "T"
            self.incomplete_num += 1
            with open(self.out_gff, "a") as o:
                o.write(re.findall(r"g[0-9]+\n(.+\n?)NC_", gene)[0])

    def write_error(self, gff, times):
        pre_dict = {}
        if times == 1:
            with open("extend_0bp.fasta", "r") as f:
                for line in f.readlines():
                    if ">" in line:
                        id = line.strip()[1:]
                        if id not in pre_dict:
                            pre_dict[id] = 1
        else:
            pre_num = int(re.findall(r"out_(.+)\-", gff)[0])
            pre_gff = re.sub(r"[0-9]+\-0", str(pre_num - 1)+"-2", gff)
            with open(pre_gff, "r") as p:
                for line in p.readlines():
                    if "NC_" in line:
                        id = line.split("\t")[0]
                        if id not in pre_dict:
                            pre_dict[id] = 1
        for key in pre_dict.keys():
            if key not in self.gene_dict:
                with open("predicate_error.txt", "a") as f:
                    f.write(key.replace(";", "\t") + "\tpredicate_null\n")

    def max_overlap(self, genes, pre_dat, strand, extend_len):
        max_index = -1
        max_overlap = 0
        native_interval = range(extend_len + 1, int(pre_dat[3]) - int(pre_dat[2]) + 2 + extend_len)
        # a = set(native_interval)
        for i in range(len(genes)):
            if pre_dat[4] == strand:
                new_interval = range(int(genes[i].split("\t")[3]), int(genes[i].split("\t")[4])+1)
                new_overlap = len(list(set(native_interval).intersection(new_interval)))
                if new_overlap > max_overlap:
                    max_overlap = new_overlap
                    max_index = i
        if max_overlap/len(native_interval) < 0.6:
            max_index = -1
            with open("predicate_error.txt", "a") as f:
                f.write("\t".join(pre_dat)+"\tdislocation\n")
        return(max_index)

    def filter_au(self, gff, times, p):
        id = ""
        new_gene = ""
        old_gene = []
        self.pass_gff = gff.replace("-0.gff3","-1.gff3")
        self.out_gff = gff.replace("-0.gff3","-2.gff3")
        with open(gff, "r") as f:
            for line in f.readlines():
                if "# start gene" in line:
                    new_gene = ""
                if "# protein sequence" in line:
                    gene_head = re.findall(r"g[0-9]+\n(.+\n?)NC_", new_gene)[0].split("\t")
                    strand = gene_head[6]
                    pre_dat = gene_head[0].split(";")
                    if gene_head[0] == id:
                        old_gene.append(new_gene)
                    else:
                        if len(old_gene) > 0:  
                            ind = self.max_overlap(old_gene, pre_dat, strand, (times-1)*p)
                            if ind >= 0:
                                self.whether_complete(old_gene[ind], (times-1)*p)
                                self.gene_dict[id] = 1
                        id = gene_head[0]
                        old_gene = [new_gene]
                    new_gene = ""
                elif "# command line" in line:
                    ind = self.max_overlap(old_gene, pre_dat, strand, (times-1)*p)
                    if ind >= 0:
                        self.whether_complete(old_gene[ind], (times-1)*p)
                        gene_head = re.findall(r"g[0-9]+\n(.+\n?)NC_", old_gene[ind])[0].split("\t")
                        id = gene_head[0]
                        self.gene_dict[id] = 1
                else:
                    new_gene += line
        self.write_error(gff, times)
        n = self.incomplete_num
        self.incomplete_num = 0
        return([self.flag, self.out_gff, n])

fi = filter()
fi.filter_au("augustus_out_24-0.gff3", 24, 100)
