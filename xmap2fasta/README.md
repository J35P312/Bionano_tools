# xmap2fasta

Convert the Xmap file to fasta
    python xmap2fasta.py --xmap xmapfile --fa reference.fasta --conversion_table contig_id

the contig_id file is used to translate the reference contig id of the xmap file to those of the reference. an example is given in contig_ids.txt, which is used to rename contig 23 and 24 of the xmap file to X and Y.
