# xmapstat
A collection of tools for validating the xmap file

compute the coverage of the xmap file
    python xmapcoverage.py --xmap xmapfile --fa reference.fasta --conversion_table contig_id --bin_size bin_size

compute the coverage across the entire genome, returns bed formated output to stoud

the contig_id file is used to translate the reference contig id of the xmap file to those of the reference. an example is given in contig_ids.txt, which is used to rename contig 23 and 24 of the xmap file to X and Y.

compute assembly statistics:
    python assemblystat.py --xmap xmapfile

prints assembly statistics to stdout, including N50, longest contig, assembly size, and number of contigs

