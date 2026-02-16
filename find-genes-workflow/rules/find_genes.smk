rule decompress_genome:
    input:
        genome = config['genomes_folder'] + '/{genome}.fna.gz'
    output: 
        genome = config['genomes_folder'] + '/{genome}.fna'
    shell:
        '''
        gzip -d {input.genome}
        '''

rule contig_lengths:
    input:
        genome = config['genomes_folder'] + '/{genome}.fna'
    output:
        contig_len = config['out_folder'] + '/contig_lengths' + '/{genome}.contig_len'
    conda:
        config['conda_env']+'/biotools'
    shell:
        '''
        python seq_lengths.py -f {input.genome} > {output.contig_len}
        '''


rule annotate_bakta:
    input:
        genome = config['genomes_folder'] + '/{genome}.fna',
        db = config['bakta']['bakta_db']
    output:
        gff = config['out_folder'] + '/bakta_results' + '/{genome}.gff3',
    conda:
        config['conda_env']+'/bakta',
    params:
        output_dir = config['out_folder'] + '/bakta_results'
    shell:
        '''
        bakta -f --db {input.db} -o {params.output_dir} --keep-contig-headers --threads {threads} {input.genome}
        '''


rule tblastn_tra_genes:
    input:
        genome = config['genomes_folder'] + '/{genome}.fna',
        reference = config['tblastn']['ref_seq_file']
    output:
        tblastn6out = config['out_folder']+'/blast_results/' + get_blast_param_combo() + '/{genome}.tblastn6out',
    params:
        query_name = config['tblastn']['ref_seq_file'].split('/')[-1],
        eval_threshold = config['tblastn']['eval'],
        codon_table = config['tblastn']['codon_table']
    envmodules:
        "gcc/8.2.0",
        "blast-plus/2.12.0"
    shell:
        '''
        scp {input.genome} $TMPDIR
        scp {input.reference} $TMPDIR
        cd $TMPDIR
        makeblastdb -in {wildcards.genome}.fna -dbtype nucl -parse_seqids -out {wildcards.genome} 
        tblastn -query {params.query_name} -db {wildcards.genome} -evalue {params.eval_threshold} -outfmt 6 -out {wildcards.genome}.tblastn6out -num_threads {threads} --db_gen_code {params.codon_table}
        scp {wildcards.genome}.tblastn6out {output.tblastn6out}
        '''


rule genomad:
    input:
        genome = config['genomes_folder'] + '/{genome}.fna',
        db = config['genomad']['genomad_db']
    output:
        result = config['out_folder'] + '/genomad_results/{genome}_annotate.log'
    conda:
        config['conda_env'] + '/genomad'
    params:
        outfolder = config['out_folder'] + '/genomad_results'
    shell:
        '''
        genomad end-to-end --cleanup --threads {threads} {input.genome} {params.outfolder} {input.db}
        '''

rule plasmidfinder:
    input:
        genome = config['genomes_folder'] + '/{genome}.fna',
        db = config['plasmidfinder']['plasmidfinder_db']
    output:
        config['out_folder'] + '/plasmidfinder_res/{genome}/{genome}.tsv'
    conda:
        config['conda_env'] + '/plasmidfinder'
    params:
        outdir = config['out_folder'] + '/plasmidfinder_res/{genome}'
    shell: 
        '''
        plasmidfinder.py -i {input.genome} -o {params.outdir} -p {input.db} -x -q
        sed '1s/^/genome\t/; 2,$s/^/{wildcards.genome}\t/' {params.outdir}/results_tab.tsv > {output}
        '''
