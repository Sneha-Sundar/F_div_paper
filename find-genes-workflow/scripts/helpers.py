from Bio import SeqIO


def extract_seqs(genome_fna, sequence_ids, output_fasta):
    """
    Extract sequences of specified contigs from a genome fasta file and write them to an output fasta file.

    Parameters:
    genome_fna (str): Path to the genome fasta file.
    sequence_ids (list): List of sequence IDs to extract.
    output_fasta (str): Path to the output fasta file where extracted sequences will be saved.

    Writes the extracted sequences to the specified output fasta file.
    """
    contig_seq_to_write = list()
    with open(genome_fna, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in sequence_ids:
                contig_seq_to_write.append(record)
    with open(output_fasta, "w") as fout:
        SeqIO.write(contig_seq_to_write, fout, "fasta")


