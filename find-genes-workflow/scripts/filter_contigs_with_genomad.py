import argparse 
from Bio import SeqIO
from helpers import extract_seqs

def main():
    parser = argparse.ArgumentParser(description="Filter sequences based on genomad classification scores")
    parser.add_argument("--genomic_unit", choices=["chromosome", "plasmid", "virus"], required=True, help="Either chromosome, plasmid or virus")
    parser.add_argument("--score_threshold",default = 0.7,type=float, required=True, help="Minimum score threshold for filtering sequences. Sequences with scores greater than or equal to this threshold will be retained.Note:genomad classificaton scores are not probabilities, but rather a measure of confidence in the classification. A higher score indicates a more confident classification, while a lower score indicates less confidence. The specific interpretation of the score may depend on the context of the analysis and the specific classification being performed.")
    parser.add_argument("--classification_file", required=True, help="Path to genomad aggregated classification TSV file")
    parser.add_argument("--genome", required=True, help="Path to genome fasta file")
    parser.add_argument("--outfile", required=True, help="Output path for filtered sequences")
  
    
    args = parser.parse_args()

    classification_file = args.classification_file
    genomic_unit = args.genomic_unit
    score_threshold = args.score_threshold
    genome_file = args.genome
    output_seq = args.outfile

    # Read the classification file and extract contig IDs that meet the score threshold for the specified genomic unit
    contigs_to_keep = list()
    with open(classification_file, 'r') as fin: 
        header = fin.readline().strip().split('\t')
        contig_id_index = header.index('seq_name')
        score_index = header.index(genomic_unit + '_score')

        for line in fin:
            fields = line.strip().split('\t')
            if float(fields[score_index]) >= score_threshold:
                contigs_to_keep.append(fields[contig_id_index])


    # Extract the sequences of the selected contigs from the genome fasta file
    extract_seqs(genome_file, contigs_to_keep, output_seq)

if __name__ == "__main__":
    main()
