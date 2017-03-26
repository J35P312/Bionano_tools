import sys
from operator import itemgetter
import argparse
import glob

parser = argparse.ArgumentParser("""merge_xmap: a script that combines multiple xmap files into a single file""")
parser.add_argument('--pattern', type=str, help="load xmap files using a search pattern")
parser.add_argument('--files', type=str, help="list of xmap files")
parser.add_argument('--print_unsplit'       , help="print unsplit alignments",action="store_true")
args= parser.parse_args()

if args.pattern:
    args.files=[]
    for input_file in glob.glob(args.pattern):
        args.files.append(input_file)

if not args.pattern and not args.files:
    print "error: no xmap input, use --files or --pattern to load xmap files"
    quit()

first = True
j=1
i=1

for xmap in args.files:
    queries={}

    for line in open(xmap):
        if first and line[0] == "#":
            print line.strip()
            continue
        elif line[0] == "#":
            continue

        entry=line.strip().split("\t")
        if not entry[1] in queries:
            queries[ entry[1] ]=[]
        queries[entry[1]].append(entry)

    first = False
    for query in queries:
        if len(queries[query]) == 1 and not args.print_unsplit:
            continue

        for alignment in queries[query]:
            alignment[0]=str(i)
            alignment[1]=str(j)
            print "\t".join(alignment)
            i += 1
        j +=1

    queries.clear()
