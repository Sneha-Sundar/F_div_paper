from Bio import SeqIO
import argparse
import sys
 
parser = argparse.ArgumentParser(prog='seq_len',
                    description='Print id and length of all sequences found in fasta record.')
parser.add_argument('-f','--input_fasta',help="Path to input multi-fasta file",required = True) 
args = parser.parse_args()
in_file = args.input_fasta

with open(in_file) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        sys.stdout.write(record.id+"\t"+str(len(record))+'\n')

# run:
# python seq_lengths.py query_tral.faa
#result: 
# NP_490564.1	103
# WP_000012097.1	103
# WP_000012100.1	103
# WP_000012106.1	103
# WP_000012113.1	103
# WP_000012116.1	103
# WP_000012119.1	103
# WP_000012129.1	103
# WP_000398849.1	101
# WP_000986844.1	104
# WP_001757027.1	98
# WP_004144424.1	101
# WP_012825259.1	103 