import sys
from operator import itemgetter
import argparse

#return complementary sequence
def reverse_comp(sequence):
    complementary=[]
    reverse_dictionary={"A":"T","a":"t","T":"A","t":"a","G":"C","g":"c","C":"G","c":"g","N":"N","n":"n"}
    for i in range(0,len(sequence)):
        complementary.append( reverse_dictionary[ sequence[len(sequence)-1-i] ] )

    return "".join(complementary)

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
parser.add_argument('--min_aligned_len', type=int,default=10000, help="alignments shorter than this length will not be printed(default = 10000 bp) ")
parser.add_argument('--max_overlap', type=int,default=100000, help="the maximum overlap of split queries")
parser.add_argument('--print_unsplit'       , help="print unsplit alignments to fasta",action="store_true")

args= parser.parse_args()

#read the reference
reference_dictionary,chromosomes=load_fasta(args.fa)

#read the reference contig conversion file
contig_coversion={}
if args.conversion_table:
    for line in open(args.conversion_table):
        content=line.strip().split("=")
        contig_coversion[content[0]]=content[1]


#read the xmap file and store all the alignments
query_contigs=load_xmap(args.xmap,contig_coversion)


print_unsplit=args.print_unsplit
max_overlap=args.max_overlap
min_aligned_len=args.min_aligned_len
#print alignments in fasta format
for contig in query_contigs:
    if len(query_contigs[contig]) ==1 and print_unsplit:
        print ">{}".format(contig)
        query_id=query_contigs[contig].keys()[0]
        if query_contigs[contig][query_id]["orientation"] == "+":
            print reference_dictionary[query_contigs[contig][query_id]["ref_contig_id"]][query_contigs[contig][query_id]["ref_start"]-1:query_contigs[contig][query_id]["ref_end"]]
        else:
            print reverse_comp( reference_dictionary[query_contigs[contig][query_id]["ref_contig_id"]][query_contigs[contig][query_id]["ref_start"]-1:query_contigs[contig][query_id]["ref_end"]] )
        continue
    elif len(query_contigs[contig]) ==1:
        continue


    query_pos=[]
    for query_id in query_contigs[contig]:
        query_pos.append([query_contigs[contig][query_id]["query_start"],query_contigs[contig][query_id]["query_end"],query_id])

    sequence=""
    i=0
    for query in sorted(query_pos,key=lambda l:l[0], reverse=True):

        #quality control
        if i and abs(query_pos[i-1][1]-query_pos[i][0]) > max_overlap:
            continue
        elif min_aligned_len > (query_pos[i][1] - query_pos[i][0]):
            continue

        if query_contigs[contig][ query[-1] ]["orientation"] == "+":
            sequence += reference_dictionary[query_contigs[contig][query[-1]]["ref_contig_id"]][query_contigs[contig][query[-1]]["ref_start"]-1:query_contigs[contig][query[-1]]["ref_end"]]
        else:
            sequence += reverse_comp( reference_dictionary[query_contigs[contig][query[-1]]["ref_contig_id"]][query_contigs[contig][query[-1]]["ref_start"]-1:query_contigs[contig][query[-1]]["ref_end"]] )
        i += 1

    print ">{}".format(contig)
    print sequence






