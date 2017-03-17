import sys
from operator import itemgetter
#first argument: the xmap file
#second argument: the reference fasta
#third argument: reference contig id conversion file(visit the read me for more detail)

#return complementary sequence
def reverse_comp(sequence):
    complementary=""
    reverse_dictionary={"A":"T","a":"t","T":"A","t":"a","G":"C","g":"c","C":"G","c":"g","N":"N","n":"n"}
    for i in range(0,len(sequence)):
        complementary += reverse_dictionary[ sequence[len(sequence)-1-i] ]

    return complementary

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

#read the reference
reference_dictionary,chromosomes=load_fasta(sys.argv[2])

#read the reference contig conversion file
contig_coversion={}
if len(sys.argv) == 4:
    for line in open(sys.argv[3]):
        content=line.strip().split("=")
        contig_coversion[content[0]]=content[1]


#read the xmap file and store all the alignments
query_contigs=load_xmap(sys.argv[1],contig_coversion)


print_unsplit=False
max_overlap=100000
min_aligned_len=10000
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






