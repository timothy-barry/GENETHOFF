import re
import os
from glob import glob
import pandas as pd
from snakemake.utils import validate


if not os.path.isfile("guideSeq_GNT.yml"):
    raise SystemExit("\n  No config file found in current directory \n")

configfile: "guideSeq_GNT.yml"
#validate(config, "config.schema.yaml")


## Load the sample information files as a TSV file. 'sampleName' column is mandatory
samples = pd.read_table(config["sampleInfo_path"],sep=";").set_index("sampleName", drop=False)
validate(samples, "samples.schema.yaml")


#print(samples)
## RUN the pipeline in the project folder.
## snakemake -s ../00-pipeline/guideSeq_GNT.smk -k -j 24 --use-conda 

##########################################################
##########################################################

rule target:
    input: "05-Report/summary.tsv"
    #["05-Report/{sample}.tsv".format(sample=sample) for sample in samples["sampleName"]]


# Merge index1 and 2 in a new fastq file I3. Easier for demultiplexing.
rule merge_indexes:
    input: I1=config["I1"] , I2=config["I2"]
    output: I3="I3.fastq.gz"
    shell: """
        python ../00-pipeline/merge_I1_I2.py {input.I1} {input.I2} {output}
        """

# make a barcode file with I1 & I2 sequences for demultiplexing.
rule make_indexes_fasta:
    input: 
    output: temp("demultiplexing_barcodes.fa")
    run: 
        with open("demultiplexing_barcodes.fa", 'w') as fasta_file:
            for index, row in samples.iterrows():
                sequence = row['index1'] + row['index2']
                fasta_file.write('>' + row['sampleName'] + '\n')
                fasta_file.write(sequence + '\n')


# demultiplexe libraries based on barcodes in the sampleInfo file and I3 file generated before
rule demultiplex_library:
    input: R1=config["R1"], R2=config["R2"], I3=rules.merge_indexes.output.I3, barcodes=rules.make_indexes_fasta.output
    output: 
     R1=["00-demultiplexing/{sample}_R1.fastq.gz".format(sample=sample) for sample in samples["sampleName"]],
     R2=["00-demultiplexing/{sample}_R2.fastq.gz".format(sample=sample) for sample in samples["sampleName"]],
     I3=["00-demultiplexing/{sample}_I3.fastq.gz".format(sample=sample) for sample in samples["sampleName"]]
    conda: "../01-envs/env_tools.yml"
    log: R1= "00-demultiplexing/demultiplexing_R1.log", R2= "00-demultiplexing/demultiplexing_R2.log"
    threads: 12
    shell: """
        
        cutadapt -g ^file:{input.barcodes} -j {threads} -e 0  --action trim --no-indels -o 00-demultiplexing/{{name}}_I3.fastq.gz -p 00-demultiplexing/{{name}}_R1.fastq.gz {input.I3} {input.R1} > {log.R1}
        cutadapt -g ^file:{input.barcodes} -j {threads} -e 0  --action trim --no-indels -o 00-demultiplexing/{{name}}_I3.fastq.gz -p 00-demultiplexing/{{name}}_R2.fastq.gz {input.I3} {input.R2} > {log.R2}
        rm 00-demultiplexing/*unknown* 
        """
        
rule add_UMI:
    input: R1="00-demultiplexing/{sample}_R1.fastq.gz", R2="00-demultiplexing/{sample}_R2.fastq.gz",I3="00-demultiplexing/{sample}_I3.fastq.gz"
    output: R1="00-demultiplexing/{sample}_R1.UMI.fastq.gz",
     R2="00-demultiplexing/{sample}_R2.UMI.fastq.gz" ,
     I3="00-demultiplexing/{sample}_I3.UMI.fastq.gz"
    threads:12
    conda: "../01-envs/env_tools.yml"
    params: suffix_length={config["UMI_length_3prime"]} ## bp in 3' of index to considere as UMI
    shell: """
        cutadapt -j {threads} -u -{params.suffix_length} --rename='{{id}}_{{r1.cut_suffix}} {{comment}}' -o {output.I3} -p {output.R1} {input.I3} {input.R1}
        cutadapt -j {threads} -u -{params.suffix_length} --rename='{{id}}_{{r1.cut_suffix}} {{comment}}' -o {output.I3} -p {output.R2} {input.I3} {input.R2}
        """

# remove ODN and discard reads without ODN
rule trim_ODN:
    input: R1=rules.add_UMI.output.R1, R2=rules.add_UMI.output.R2
    output: R1="01-trimming/{sample}_R1.UMI.ODN.fastq.gz",
     R2="01-trimming/{sample}_R2.UMI.ODN.fastq.gz"
    threads: 12
    log:"01-trimming/{sample}.odn.log"
    conda: "../01-envs/env_tools.yml"
    message: "removing ODN sequence, discard reads without ODN sequence"
    params: R2_leading=lambda wildcards: config[samples["type"][wildcards.sample]][samples["orientation"][wildcards.sample]]["R2_leading"]
    shell: """
        cutadapt -j {threads} -G "^{params.R2_leading};max_error_rate=0" --discard-untrimmed  -o {output.R1} -p {output.R2} {input.R1} {input.R2} > {log}
        """


# remove leading and trailing ODN sequences
rule trim_reads:
    input: R1=rules.trim_ODN.output.R1, R2=rules.trim_ODN.output.R2
    output: R1="01-trimming/{sample}_R1.UMI.ODN.trimmed.fastq.gz", 
     R2="01-trimming/{sample}_R2.UMI.ODN.trimmed.fastq.gz"
    threads: 12
    log: "01-trimming/{sample}.trailing.log"
    conda: "../01-envs/env_tools.yml"
    message: "trimming ODN and adaptor sequences in reads"
    params: R2_leading=lambda wildcards: config[samples["type"][wildcards.sample]][samples["orientation"][wildcards.sample]]["R2_leading"], 
     R2_trailing=lambda wildcards: config[samples["type"][wildcards.sample]][samples["orientation"][wildcards.sample]]["R2_trailing"],
     R1_trailing=lambda wildcards: config[samples["type"][wildcards.sample]][samples["orientation"][wildcards.sample]]["R1_trailing"]
    shell: """
        cutadapt -j {threads} -A "{params.R2_trailing};min_overlap=6;max_error_rate=0.1" -a "{params.R1_trailing};min_overlap=6;max_error_rate=0.1" -o {output.R1} -p {output.R2} {input.R1} {input.R2} > {log}
        """



# make a size selection before mapping
rule filter_reads:
    input: R1=rules.trim_reads.output.R1, R2=rules.trim_reads.output.R2
    output:  R1="02-filtering/{sample}_R1.UMI.ODN.trimmed.filtered.fastq.gz", R2="02-filtering/{sample}_R2.UMI.ODN.trimmed.filtered.fastq.gz",
     R1short="02-filtering/{sample}_R1.UMI.ODN.trimmed.tooshort.fastq.gz", R2short="02-filtering/{sample}_R2.UMI.ODN.trimmed.tooshort.fastq.gz"
    threads: 12
    log: "02-filtering/{sample}.filter.log"
    params: length=config["minLength"]
    conda: "../01-envs/env_tools.yml"
    message: "remove pairs if one mate is shorter than x bp"
    shell: """
        cutadapt -j {threads} --pair-filter=any --minimum-length {params.length} --too-short-output {output.R1short} --too-short-paired-output  {output.R2short} --output {output.R1} --paired-output {output.R2} {input.R1} {input.R2} > {log}
        """

# map reads on the reference genome as pairs
if (config["aligner"]  == "bowtie2" or config["aligner"]  == "Bowtie2") :
    rule alignOnGenome:
        input: R1=rules.filter_reads.output.R1, R2=rules.filter_reads.output.R2
        output: sam=temp("03-align/{sample}.UMI.ODN.trimmed.filtered.sam")
        threads: 6
        log: "03-align/{sample}.UMI.ODN.trimmed.filtered.align.log"
        conda: "env_tools.yml"
        message: "Aligning PE reads on genome with Bowtie2"
        params: n=config["reportedAlignments"], index=config["genome"]["index"]
        shell: """
            bowtie2 -p {threads} --no-unal -X 1500 --dovetail --no-mixed --no-discordant --un-conc-gz 03-align/{wildcards.sample}_R%.UMI.ODN.trimmed.unmapped.fastq.gz   -x {params.index} -1 {input.R1} -2 {input.R2} -S {output.sam} 2> {log}
        """
elif (config["aligner"]  == "bwa" or config["aligner"]  == "Bwa") :
    rule alignOnGenome:
        input: R1=rules.filter_reads.output.R1, R2=rules.filter_reads.output.R2
        output: sam="03-align/{sample}.UMI.ODN.trimmed.filtered.sam", unfilterdsam="03-align/{sample}.UMI.ODN.trimmed.unfiltered.sam"
        threads: 6
        log: "03-align/{sample}.UMI.ODN.trimmed.filtered.align.log"
        conda: "tools"
        message: "Aligning PE reads on genome with BWA mem"
        params: n=config["reportedAlignments"], index=config["genome"]["index"]
        shell: """
            #bwa mem -t {threads} {params.index} {input.R1} {input.R2} > {output.sam} 2> {log}
            bwa mem -t {threads} {params.index} {input.R1} {input.R2} > {output.unfilterdsam} 2> {log}
            samtools view -F 0x4 -F 0x8 -F 0x100 -F 0x800 -f 0x2 -b {output.unfilterdsam} > {output.sam}
        """


# sort alignments by names (required for BEDPE conversion) and position (for viewing)
rule sort_aligned:
    input: sam=rules.alignOnGenome.output.sam
    output: bamPos="03-align/{sample}.UMI.ODN.trimmed.filtered.bam",  bamName="03-align/{sample}.UMI.ODN.trimmed.filtered.sortedName.bam"
    threads: 6
    log:
    conda: "../01-envs/env_tools.yml"
    message: "Sort reads by name"
    shell: """
        samtools sort -@ {threads} {input.sam} > {output.bamPos}
        samtools index {output.bamPos}
        samtools sort -n -@ {threads} {input.sam} > {output.bamName}
    """


# convert BAM file to BEDPE and identify insertion point, read length
# filter alignment with MAPQ score > threshold in config file
# aggregate reads per fragment size and genomic coordinates
# add a cluster ID to group close IS (distance defined in the config file) 
rule call_IS:
    input: rules.sort_aligned.output.bamName
    output: tmp="04-IScalling/{sample}.pebed", 
     frag="04-IScalling/{sample}.readsPerFragmentPerIS.bed",
     collapse="04-IScalling/{sample}.collapsefragPerISCluster.bed"
    threads: 1
    conda: "../01-envs/env_tools.yml"
    log:
    params: window=config["ISbinWindow"], minMAPQ=config["minMAPQ"],minReadsPerFragment=config["minReadsPerFrag"]
    shell: """
        bedtools bamtobed -bedpe -mate1 -i  {input} > {output.tmp}

        # read R2 contains the genome/ODN junction (use $10 of PE bed = mate 2)
        # count number of reads per UMI and per IS (pos and strand separated)

        #echo -e "#chromosome\tstart\tend\tUMI\tstrand\tNreads\tmedianMAPQ\tclusterID" > {output.frag}
        awk 'BEGIN{{OFS="\\t";FS="\\t"}} ($8>={params.minMAPQ}) && ($1 == $4) {{if($10=="+") print $4,$5,$5,$7,substr($7,length($7)-7,8),$8,$10,$3-$5; else print $4,$6-1,$6-1,$7,substr($7,length($7)-7,8),$8,$10,$6-$2}}' {output.tmp} |  sort -k1,1 -k2,3n -k5,5 -k7,7  | bedtools groupby -g 1,2,3,5,7 -c 4,6 -o count_distinct,median | bedtools cluster -d {params.window}  > {output.frag}
        
        # aggregate UMI per IS
        #echo -e "#chromosome\tstart\tend\tIS_ID\tMedianMAPQ\tstrand\tN_UMI\tNreads\tUMI_list\tReadPerUMI\tclusterID" > {output.collapse}
        awk '$6>{params.minReadsPerFragment}' {output.frag} | sort -k1,1 -k2,2n -k3,3n -k8,8n -k5,5 | bedtools groupby -g 1,2,3,5,8 -c 4,6,4,6,7  -o count_distinct,sum,collapse,collapse,median  | awk 'BEGIN{{OFS="\\t";FS="\\t"}} {{print $1,$2,$3,"ID_"NR,$10,$4,$6,$7,$8,$9,$5}}' > {output.collapse}
        """

rule get_chrom_length:
    input: config["genome"]["fasta"]
    output: config["genome"]["fasta"]+".fai"
    threads: 1
    conda: "../01-envs/env_tools.yml"
    params: index=config["genome"]["fasta"]
    shell: """
        samtools faidx {params.index}
        """

rule get_fasta_around_is:
    input: bed=rules.call_IS.output.collapse, chr_length=rules.get_chrom_length.output
    output: cluster="04-IScalling/{sample}.cluster_slop.bed", fa="04-IScalling/{sample}.cluster_slop.fa"
    threads: 1
    conda: "../01-envs/env_tools.yml"
    log:
    params: fasta=config["genome"]["fasta"], slop_size=config["slopSize"]
    shell: """
        #echo -e "#chromosome\tstart\tend\tclusterID\tN_orientations\tmedianMAPQ\tN_IS\tN_UMI\tN_reads" > {output.cluster}
        sort -k11n {input.bed} | bedtools groupby -g 11 -c 1,2,3,4,5,6,7,8 -o distinct,min,max,count_distinct,median,count_distinct,sum,sum |  awk 'BEGIN{{OFS="\\t"}}{{print $2,$3,$4,$1,$7,$6,$5,$8,$9}}' |  bedtools slop -i - -b {params.slop_size} -g {input.chr_length} > {output.cluster} 
        bedtools getfasta -name -fi {params.fasta} -bed {output.cluster} > {output.fa}
        """ 

rule get_stats_fq:
    input: rules.demultiplex_library.output.R1, rules.add_UMI.output.R1, rules.trim_ODN.output.R1, rules.trim_reads.output.R1, rules.filter_reads.output.R1
    output: "05-Report/{sample}.stat"
    conda: '../01-envs/env_tools.yml'
    threads: 6
    shell: """
        seqkit stat -j {threads} -a -T {input} > {output}
        """

rule report_data:
    input: fasta=rules.get_fasta_around_is.output.fa, cluster=rules.get_fasta_around_is.output.cluster, bed=rules.call_IS.output.collapse, statfq=rules.get_stats_fq.output,statal=rules.alignOnGenome.log
    output: "05-Report/{sample}.rdata"
    conda: "../01-envs/env_R4.3.2.yml"
    log:
    threads: 1
    params: gRNA_seq=lambda wildcards:samples["gRNA_sequence"][wildcards.sample], gRNA_name=lambda wildcards:samples["gRNA_name"][wildcards.sample]
    shell : """
        Rscript ../00-pipeline/multiple_alignments.R {input.fasta} {input.cluster} {input.bed} {params.gRNA_seq} {params.gRNA_name} {output}
        """

rule annotate_sites:
    input: ["05-Report/{sample}.rdata".format(sample=sample) for sample in samples["sampleName"]]
    output: "05-Report/summary.tsv"
    conda: "../01-envs/env_R4.3.2.yml"
    threads: 1
    params:
    shell: """
        Rscript ../00-pipeline/annotate_cuting_sites.R 
        """


onsuccess:
    print("Workflow finished, no error")
    shell("cowsay -e \♥♥  You rock !! ")

onerror:
    print("An error occurred")
    shell(" cowsay -e \XX  Houston, we have a problem !")
