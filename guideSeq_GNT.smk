import re
import os
from glob import glob

configfile: "guideSeq_GNT.yml"


## get the list of samples from the sampleInfo file first column
samples_list=[]
with open(config["sampleInfo_path"],'r') as file:
        next(file)
        for line in file:
                columns=line.strip().split(',')
                samples_list.append(columns[0])

print(samples_list)


##########################################################
##########################################################

rule target:
    input: expand("04-IScalling/{sample}_{orientation}.readsPerFragmentPerIS_annotated.bed", sample=samples_list, orientation=["POS","NEG"])


# Merge index1 and 2 in a new fastq file I3. Easier for demultiplexing.
rule merge_indexes:
    input: I1=config["I1"] , I2=config["I2"]
    output: I3=temp("I3.fastq.gz")
    shell: """
        python merge_I1_I2.py {input.I1} {input.I2} {output}
        """

# make a barcode file with I1 & I2 sequences for demultiplexing.
rule make_indexes_fasta:
    input: config["sampleInfo_path"]
    output: "demultiplexing_barcodes.fa"
    shell: """
        awk ' BEGIN{{FS=","}} NR> 1  {{print ">"$1"\\n"$3$4}}' {input} > {output}
        """


# demultiplexe libraries based on barcodes in the sampleInfo file and I3 file generated before
rule demultiplex_library:
    input: R1=config["R1"], R2=config["R2"], I3=rules.merge_indexes.output.I3, barcodes=rules.make_indexes_fasta.output
    output: temp(expand("00-demultiplexing/{sample}_R1.fastq.gz", sample=samples_list)),temp(expand("00-demultiplexing/{sample}_R2.fastq.gz", sample=samples_list))
    conda: "tools"
    log: R1= "00-demultiplexing/demultiplexing_R1.log", R2= "00-demultiplexing/demultiplexing_R2.log"
    threads: 12
    shell: """
        
        cutadapt -g ^file:{input.barcodes} -j {threads} -e 0.1  --action none --no-indels -o 00-demultiplexing/{{name}}_I3.fastq.gz -p 00-demultiplexing/{{name}}_R1.fastq.gz {input.I3} {input.R1} > {log.R1}
        cutadapt -g ^file:{input.barcodes} -j {threads} -e 0.1  --action none --no-indels -o 00-demultiplexing/{{name}}_I3.fastq.gz -p 00-demultiplexing/{{name}}_R2.fastq.gz {input.I3} {input.R2} > {log.R2}
        
        rm 00-demultiplexing/*I3* 
        rm 00-demultiplexing/*unknown*
        """


# make a fasta file with R2 leading sequence in order to recognize POS and NEG libraries
rule make_orientation_fasta:
    input:
    output: temp("00-demultiplexing/orientation.fa")
    params:POS=config[config["type"]]["POS"]["R2_leading"],NEG=config[config["type"]]["NEG"]["R2_leading"]
    shell: """
        echo -e ">POS\n"{params.POS}"\n>NEG\n"{params.NEG} > {output}
        """

# split libraries in PSO and NEG orientation (if they exist)
rule demultiplex_orientation:
    input: R1="00-demultiplexing/{NAME}_R1.fastq.gz", R2="00-demultiplexing/{NAME}_R2.fastq.gz", fa = rules.make_orientation_fasta.output
    output: "00-demultiplexing/{NAME}_orientation.log", 
     R1pos="00-demultiplexing/{NAME}_R1.POS.fastq.gz", R2neg="00-demultiplexing/{NAME}_R2.NEG.fastq.gz",
     R1neg="00-demultiplexing/{NAME}_R1.NEG.fastq.gz", R2pos="00-demultiplexing/{NAME}_R2.POS.fastq.gz"
    conda: "tools"
    message: "Spliting reads according to PCR direction"
    log:"00-demultiplexing/{NAME}_orientation.log"
    threads: 12
    params: 
    shell: """
        cutadapt  -j {threads} --action none -e 0.1 -g ^file:{input.fa} --no-indels -o 00-demultiplexing/{wildcards.NAME}_R2.{{name}}.fastq.gz -p 00-demultiplexing/{wildcards.NAME}_R1.{{name}}.fastq.gz {input.R2} {input.R1} > {log}
        rm 00-demultiplexing/*unknown*
    """

rule trim_reads:
    input: R1="00-demultiplexing/{NAME}_R1.{ORIENTATION}.fastq.gz", R2="00-demultiplexing/{NAME}_R2.{ORIENTATION}.fastq.gz"
    output: R1="01-trimming/{NAME}_R1.{ORIENTATION}.trimmed.fastq.gz", R2="01-trimming/{NAME}_R2.{ORIENTATION}.trimmed.fastq.gz"
    threads: 12
    log: "01-trimming/{NAME}_{ORIENTATION}.log"
    conda: "tools"
    message: "trimming ODN and adaptor sequences in reads"
    params: R2_leading=lambda wildcards: config[config["type"]][wildcards.ORIENTATION]["R2_leading"], 
     R2_trailing=lambda wildcards: config[config["type"]][wildcards.ORIENTATION]["R2_trailing"],
     R1_trailing=lambda wildcards: config[config["type"]][wildcards.ORIENTATION]["R1_trailing"]
    shell: """
        cutadapt -j {threads} -A "^{params.R2_leading};max_error_rate=0...{params.R2_trailing};min_overlap=6;max_error_rate=0.1" -a {params.R1_trailing} -o {output.R1} -p {output.R2} {input.R1} {input.R2} > {log}
        """


rule filter_reads:
    input: R1=rules.trim_reads.output.R1, R2=rules.trim_reads.output.R2
    output:  R1="02-filtering/{NAME}_R1.{ORIENTATION}.trimmed.filtered.fastq.gz", R2="02-filtering/{NAME}_R2.{ORIENTATION}.trimmed.filtered.fastq.gz",
     R1short="02-filtering/{NAME}_R1.{ORIENTATION}.trimmed.tooshort.fastq.gz", R2short="02-filtering/{NAME}_R2.{ORIENTATION}.trimmed.tooshort.fastq.gz"
    threads: 6
    log: "02-filtering/{NAME}_{ORIENTATION}.filter.log"
    params: length=config["minLength"]
    conda: "tools"
    message: "remove pairs if one mate is hsorter than x bp"
    shell: """
        cutadapt -j {threads} --pair-filter=any --minimum-length {params.length} --too-short-output {output.R1short} --too-short-paired-output  {output.R2short} --output {output.R1} --paired-output {output.R2} {input.R1} {input.R2} > {log}
        """


rule build_index:
    input: config["reference_genome"]+".fa"
    output:  index=config["reference_genome"]+".1.bt2", fai=config["reference_genome"]+"fa.fai"
    threads: 6
    shell: """
        bowtie-build --threads {threads} {input} {output.prefix}
        samtools faidx {input}
        cut -f1,2 {config[reference_genome]}.fa.fai
    """

rule alignOnGenome:
    input: R1=rules.filter_reads.output.R1, R2=rules.filter_reads.output.R2, index = rules.build_index.output.index
    output: sam=temp("03-align/{NAME}_{ORIENTATION}.sam")
    threads: 6
    log: "03-align/{NAME}_{ORIENTATION}.align.log"
    conda: "tools"
    message: "Aligning PE reads on genome"
    params: n=config["reportedAlignments"]
    shell: """
        bowtie2 -p {threads} --no-unal -X 1000  --no-mixed --no-discordant --un-conc-gz 03-align/{wildcards.NAME}_R%.{wildcards.ORIENTATION}.trimmed.unmapped.fastq.gz   -x {config[reference_genome]} -1 {input.R1} -2 {input.R2} -S {output.sam} 2> {log}
    """

rule sort_aligned:
    input: sam=rules.alignOnGenome.output.sam
    output: bamPos="03-align/{NAME}_{ORIENTATION}.bam",  bamName="03-align/{NAME}_{ORIENTATION}.sortedName.bam"
    threads: 6
    log:
    conda: "tools"
    message: "Sort reads by name"
    shell: """
        samtools sort -@ {threads} {input.sam} > {output.bamPos}
        samtools index {output.bamPos}
        samtools sort -n -@ {threads} {input.sam} > {output.bamName}
    """

rule call_IS:
    input: rules.sort_aligned.output.bamName
    output: tmp="04-IScalling/{NAME}_{ORIENTATION}.pebed", 
     frag="04-IScalling/{NAME}_{ORIENTATION}.readsPerFragmentPerIS.bed"
    threads: 1
    conda: "tools"
    log:
    params:window=config["ISbinWindow"],minMAPQ=config["minMAPQ"]
    shell: """
        bedtools bamtobed -bedpe -mate1 -i  {input} > {output.tmp}

        # count number of reads per fragment-size and per IS (pos and strand)
        awk 'BEGIN{{OFS="\\t";FS="\\t"}} $8>{params.minMAPQ} {{if($10=="+") print $1,$5,$5,$7,$8,$10,$3-$5; else print $1,$6-1,$6-1,$7,$8,$10,$6-$2}}' {output.tmp} | sort -k1,1 -k2,3n -k6,6 -k7,7n | bedtools groupby -g 1,2,3,6,7 -c 4 -o count_distinct | awk 'BEGIN{{OFS="\\t";FS="\\t"}} {{print $1,$2,$3,"ID_"NR,0,$4,$5,$6}}' | bedtools cluster -i - -d {params.window} -s  > {output.frag}
        """

rule makeAnnotationDB:
    input: config["annotation_gtf"]+".gtf"
    output: config["annotation_gtf"]+".genesOnly.gtf"
    threads: 1
    shell: """
         awk '$0~"#" || $3 == "gene" {{print $0}}' {input} > {output}
        """


rule annotate_cutSites:
    input: cutsites=rules.call_IS.output.frag, annotation=rules.makeAnnotationDB.output
    output: "04-IScalling/{NAME}_{ORIENTATION}.readsPerFragmentPerIS_annotated.bed"
    conda: "tools"
    params:
    message: "annotating cutting sites"
    threads: 1 
    shell: """
        bedtools intersect -a {input.cutsites} -b {input.annotation} -loj > {output}
        """

