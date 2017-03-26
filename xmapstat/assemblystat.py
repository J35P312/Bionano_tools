import sys
from operator import itemgetter
import argparse

def load_xmap(xmapfile):
    query_contigs={}
    for line in open(xmapfile):
        if "#" == line[0]:
            continue
        content=line.split("\t")
        query_id=content[1]
        if not query_id in query_contigs:
            query_contigs[content[1]] = {}
    
        ref_id=content[2]
        query_pos=[int(float(content[3])),int(float(content[4]))]
        ref_pos=[int(float(content[5])),int(float(content[6]))]

        query_contigs[content[1]][content[0]]={"ref_contig_id":ref_id,"ref_start":min(ref_pos),"ref_end":max(ref_pos),"query_size":int(float(content[10])) ,"query_start":min(query_pos),"query_end":max(query_pos),"confidence":int(float(content[8])),"orientation":content[7]}

    return query_contigs



parser = argparse.ArgumentParser("""xmap2vcf: a script used to convert xmap to fasta""")
parser.add_argument('--xmap', type=str, help="The xmap file",required=True)

args= parser.parse_args()

#read the xmap file and store all the alignments
query_contigs=load_xmap(args.xmap)


contig_sizes=[]
for contig in query_contigs:
    contig_size=0
    for alignment in query_contigs[contig]:
        contig_size += query_contigs[contig][alignment]["ref_end"]-query_contigs[contig][alignment]["ref_start"]
    contig_sizes.append(contig_size)


assembly_size=sum(contig_sizes)
n=0
N=0
for contig in sorted(contig_sizes):
    n+= contig
    if n >= assembly_size/2: 
       N= contig
       break

print "\t".join(["file","N50","number_of_contigs","assembly_size","longest_contig"])
print "{}\t{}\t{}\t{}\t{}".format(args.xmap,N,len(contig_sizes),assembly_size,max(contig_sizes))







