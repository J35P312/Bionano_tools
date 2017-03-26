import sys
from operator import itemgetter
import argparse
import numpy
import math
def load_fasta(fafile):
    sequence={}
    chromosome_order=[]
    #read the fast file
    with open(fafile, 'r+') as f:
        reference = f.read()
    split_reference=reference.split(">")
    del reference
    del split_reference[0]
    #store the reference as a dictionary
    chromosome_len={}
    chromosomes=[]
    simulated_bases=0
    for chromosome in split_reference:
        content=chromosome.split("\n",1)
        sequence[content[0].split()[0]]=content[1].replace("\n","")
        chromosome_order.append(content[0].split()[0])
    del split_reference
    return(sequence,chromosome_order)

def load_xmap(xmapfile,contig_conversion):
    query_contigs={}
    for line in open(xmapfile):
        if "#" == line[0]:
            continue
        content=line.split("\t")
        query_id=content[1]
        if not query_id in query_contigs:
            query_contigs[content[1]] = {}
    
        ref_id=content[2]
        if ref_id in contig_coversion:
            ref_id = contig_coversion[ref_id]
        query_pos=[int(float(content[3])),int(float(content[4]))]
        ref_pos=[int(float(content[5])),int(float(content[6]))]

        query_contigs[content[1]][content[0]]={"ref_contig_id":ref_id,"ref_start":min(ref_pos),"ref_end":max(ref_pos),"query_size":int(float(content[10])) ,"query_start":min(query_pos),"query_end":max(query_pos),"confidence":int(float(content[8])),"orientation":content[7]}

    return query_contigs



parser = argparse.ArgumentParser("""xmap2vcf: a script used to convert xmap to fasta""")
parser.add_argument('--conversion_table', type=str, help="The contig id conversion file, used to convert the xmap id to those of the fasta file")
parser.add_argument('--fa', type=str, help="The reference fasta",required=True)
parser.add_argument('--xmap', type=str, help="The xmap file",required=True)
parser.add_argument('--bin_size', type=int,default=1000, help="bin size (default = 1000) ")


args= parser.parse_args()

#read the reference
reference_dictionary,chromosomes=load_fasta(args.fa)
coverage_structure={}
for chromosome in chromosomes:
    chromosome_length=math.ceil( len(reference_dictionary[chromosome])/float(args.bin_size) )
    coverage_structure[chromosome]=numpy.zeros(chromosome_length)
    del reference_dictionary[chromosome]

#read the reference contig conversion file
contig_coversion={}
if args.conversion_table:
    for line in open(args.conversion_table):
        content=line.strip().split("=")
        contig_coversion[content[0]]=content[1]


#read the xmap file and store all the alignments
query_contigs=load_xmap(args.xmap,contig_coversion)
for query_contig in query_contigs:
    for alignment in query_contigs[query_contig]:
        chromosome=query_contigs[query_contig][alignment]["ref_contig_id"]
        pos=query_contigs[query_contig][alignment]["ref_start"]
        length=query_contigs[query_contig][alignment]["ref_end"]-query_contigs[query_contig][alignment]["ref_start"]
        while length > 0:
            index=math.floor(pos/args.bin_size)
            if length > args.bin_size:
                coverage_structure[chromosome][index] += 1
                pos += args.bin_size
                length += -args.bin_size
            else:
                coverage_structure[chromosome][index] +=  length/float(args.bin_size)
                break

for chromosome in chromosomes:
    pos=0
    for bin in coverage_structure[chromosome]:
        print"{}\t{}\t{}\t{}".format(chromosome,pos,pos+args.bin_size,bin)
        pos+=args.bin_size






