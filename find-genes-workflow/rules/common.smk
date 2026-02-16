def get_input_samples():
    with open(config["samples"],'r') as fin:
        samples = [line.rstrip() for line in fin]
    return samples

def get_blast_param_combo():
    return "eval{}_codon.{}".format(config['tblastn']['eval'],config['tblastn']['codon_table'])

def get_tblastn_output():
    samples = get_input_samples()
    param_combo = get_blast_param_combo()

    return expand(config['out_folder']+'/blast_results/' + param_combo + '/{genome}.tblastn6out',genome = samples)

def get_bakta_output():
    samples = get_input_samples()
    return expand(config['out_folder'] + '/bakta_results' + '/{genome}.gff3',genome = samples)

def get_genomad_output():
    samples = get_input_samples()
    return expand(config['out_folder'] + '/genomad_results'+'/{genome}_annotate.log',genome = samples)

def get_contiglens_output():
    samples = get_input_samples()
    return expand(config['out_folder'] + '/contig_lengths' + '/{genome}.contig_len',genome = samples)

def get_plasmidfinder_output():
    samples = get_input_samples()
    return expand(config['out_folder'] + '/plasmidfinder_res' + '/{genome}/{genome}.tsv',genome = samples)