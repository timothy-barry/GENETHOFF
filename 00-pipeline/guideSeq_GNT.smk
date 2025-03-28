import re
import os
from glob import glob
import pandas as pd
from snakemake.utils import validate
import sys
from collections import defaultdict


# check that the config file is present in current folder
if not os.path.isfile("guideSeq_GNT.yml"):
    raise SystemExit("\n  No config file found in current directory \n")

configfile: "guideSeq_GNT.yml"
config_file_path=os.getcwd()+'/'+workflow.configfiles[0]


## Load the sample information files as a TSV file.
# check that the config file is present in current folder
if not os.path.isfile(config["sampleInfo_path"]):
    raise SystemExit("\n  No Sample Dta Sheet file found in current directory \n")

samplesTable = pd.read_table(config["sampleInfo_path"],sep=";").set_index("sampleName", drop=False)

#check the validity of the sample Data Sheet (Path is relative to pipeline file, not current folder)
validate(samplesTable, "samples.schema.yaml")





## if multiple rows have the same sampleName, check that they have the same gRNA,PAM, Cas and genome
# The purpose of this code is to ensure data consistency by verifying that samples with the same name have identical features. 
# The check for hyphens in sample names is to prevent issues with delimiters in file names or identifiers.

RED= "\x1b[33m"
RESET= RESET = "\033[0m"

columns_to_check = ['sampleName','Genome','gRNA_name', 'gRNA_sequence','PAM_sequence','Cas','type',"Cut_Offset"]

identical_values=samplesTable.groupby(samplesTable.index).apply(lambda x: x[columns_to_check].nunique().eq(1).all())
    
if identical_values.all():
    samples = samplesTable[columns_to_check].drop_duplicates()
    print("################################################\n",samples,"\n################################################")
    
    my_list=samples[["Genome","gRNA_sequence","PAM_sequence"]].drop_duplicates()
    
    
    samples_unique=list(set(samples["sampleName"]))
    
    ## check there is no "-" in sample name
    strings_with_hyphen = [item for item in samples_unique if "-" in item]
    if strings_with_hyphen:
        print(f"{RED} !! '-' detected in the following sample name :", strings_with_hyphen,"\n Please modify sample name and replace '-' by another delimiter")
        sys.exit(1)
    else:
      
      genomes=my_list["Genome"].tolist()
      gRNAs=my_list["gRNA_sequence"].tolist()
      PAMs=my_list["PAM_sequence"].tolist()
      genomes_unique=list(set(genomes))
 
else:
    print(f"{RED}!! Samples with identical 'sampleName' value must share the same features : gRNA_name, gRNA_sequence, PAM_sequence, Genome, Cas\n---> Please correct metadata or use different samples name{RESET}")
    print(samplesTable[columns_to_check])
    sys.exit(1)










## RUN the pipeline in the project folder.
## snakemake -s ../00-pipeline/guideSeq_GNT.smk -k -j 12 --use-conda --conda-front-end mamba --conda-prefix ../01-envs -n

##########################################################
##########################################################

rule target:
    input: 
     expand("06-offPredict/{gen}_{grna}{pam}.csv", zip, grna=gRNAs, gen=genomes,pam=PAMs),
     ["05-Report/{sample}_summary.xlsx".format(sample=sample) for sample in samples_unique],
     "05-Report/report.html"
     

#rule make_genome_index:
#    input: lambda wildcards: config["genome"][wildcards.genome]["fasta"]
#    output: "../02-ressources/{genome}.log"
#    threads: 24
#    conda: "../01-envs/env_tools.yml"
#    params: idx=lambda wildcards: config["genome"][wildcards.genome]["index"]
#    shell: """
#        mkdir -p $(dirname {params.idx})
#        bowtie2-build --threads {threads} {input} {params.idx} > {params}.log
#        ln -s {params}.log {output}
#        """



rule get_chrom_length:
    input: lambda wildcards: config["genome"][wildcards.genome]["fasta"]
    output: "../02-ressources/{genome}.chrom.length"
    threads: 1
    conda: "../01-envs/env_tools.yml"
    params: lambda wildcards: config["genome"][wildcards.genome]["fasta"]+".fai"
    shell: """
        samtools faidx {input}
        ln -s  {params} {output}
        """





rule prepare_annotations:
    input: lambda wildcards: config["genome"][wildcards.genome]["annotation"]
    output: "../02-ressources/{genome}.rds"
    threads: 1
    conda: "../01-envs/env_R4.3.2.yml"
    params: config=os.path.abspath(workflow.configfiles[0])
    shell: """
        Rscript ../00-pipeline/prepare_annotations.R {wildcards.genome} {input}
    """     

     

# Merge index1 and 2 in a new fastq file I3. Easier for demultiplexing.
rule merge_indexes:
    input: I1=config["I1"] , I2=config["I2"]
    output: I3=("I3.fastq.gz")
    conda: "../01-envs/env_tools.yml"
    shell: """
        #python ../00-pipeline/merge_I1_I2.py {input.I1} {input.I2} {output}
        
        paste -d ";" <(gzip -dc {input.I1}) <(gzip -dc {input.I2}) | awk 'BEGIN{{FS=";";OFS=""}} {{if( NR % 4 ==1 || NR % 4 == 3) {{print $1}} else {{print $1,$2}}}}' > I3.fastq
        gzip I3.fastq

        """



# make a barcode file with I1 & I2 sequences for demultiplexing.
rule make_indexes_fasta:
    input: 
    output: ("demultiplexing_barcodes.fa")
    run: 
        sample_count = defaultdict(int)
        with open("demultiplexing_barcodes.fa", 'w') as fasta_file:
          for index, row in samplesTable.iterrows():
            sequence = row['index1'] + row['index2']
            sample_count[row['sampleName']] += 1
            sample_header = f">{row['sampleName']}-{sample_count[row['sampleName']]}"
            fasta_file.write(sample_header + '\n')
            fasta_file.write(sequence + '\n')


# demultiplexe libraries based on barcodes in the sampleInfo file and I3 file generated before
## 
# check whether an empty library will make the pip to crash (ie wrong index)
##
rule demultiplex_library:
    input: R1=config["R1"], R2=config["R2"], I1=config["I1"], I2=config["I2"],I3=rules.merge_indexes.output.I3, barcodes=rules.make_indexes_fasta.output
    output: 
     R1=temp(["00-demultiplexing/{sample}-1_R1.fastq.gz".format(sample=sample) for sample in samples["sampleName"]]),
     R2=temp(["00-demultiplexing/{sample}-1_R2.fastq.gz".format(sample=sample) for sample in samples["sampleName"]]),
     I3=temp(["00-demultiplexing/{sample}-1_I3.fastq.gz".format(sample=sample) for sample in samples["sampleName"]])
    conda: "../01-envs/env_tools.yml"
    log: R1= "00-demultiplexing/demultiplexing_R1.log", R2= "00-demultiplexing/demultiplexing_R2.log"
    threads: 6
    shell: """
        
        cutadapt -g ^file:{input.barcodes} -j {threads} -e 0  --action trim --no-indels -o 00-demultiplexing/{{name}}_I3.fastq.gz -p 00-demultiplexing/{{name}}_R1.fastq.gz {input.I3} {input.R1} > {log.R1}
        cutadapt -g ^file:{input.barcodes} -j {threads} -e 0  --action trim --no-indels -o 00-demultiplexing/{{name}}_I3.fastq.gz -p 00-demultiplexing/{{name}}_R2.fastq.gz {input.I3} {input.R2} > {log.R2}
        #cutadapt -g ^file:{input.barcodes} -j {threads} -e 0  --action trim --no-indels -o 00-demultiplexing/{{name}}_I3.fastq.gz -p 00-demultiplexing/{{name}}_I1.fastq.gz {input.I3} {input.I1} > {log.R1}
        #cutadapt -g ^file:{input.barcodes} -j {threads} -e 0  --action trim --no-indels -o 00-demultiplexing/{{name}}_I3.fastq.gz -p 00-demultiplexing/{{name}}_I2.fastq.gz {input.I3} {input.I2} > {log.R1}
        rm 00-demultiplexing/*unknown* 
        """
        


rule merge_sampleName:
    input: R1="00-demultiplexing/{sample}-1_R1.fastq.gz", R2="00-demultiplexing/{sample}-1_R2.fastq.gz",I3="00-demultiplexing/{sample}-1_I3.fastq.gz"
    output: R1="00-demultiplexing/{sample}_R1.fastq.gz", R2="00-demultiplexing/{sample}_R2.fastq.gz",I3="00-demultiplexing/{sample}_I3.fastq.gz"
    threads:1
    shell: """
        cat 00-demultiplexing/{wildcards.sample}-*_R1.fastq.gz > {output.R1}
        cat 00-demultiplexing/{wildcards.sample}-*_R2.fastq.gz > {output.R2}
        cat 00-demultiplexing/{wildcards.sample}-*_I3.fastq.gz > {output.I3}
    
    """



# remove ODN and discard reads without ODN
rule trim_ODN:
  input: R1="00-demultiplexing/{sample}_R1.fastq.gz", R2="00-demultiplexing/{sample}_R2.fastq.gz",I3="00-demultiplexing/{sample}_I3.fastq.gz"
  output: R1=temp("01-trimming/{sample}_R1.ODN.fastq.gz"),
   R2=temp("01-trimming/{sample}_R2.ODN.fastq.gz"),
   I3=temp("01-trimming/{sample}_I3.ODN.fastq.gz")
  threads: 6
  log:"01-trimming/{sample}.odn.log"
  conda: "../01-envs/env_tools.yml"
  message: "removing ODN sequence, discard reads without ODN sequence"
  params: ODN_pos=lambda wildcards: config[samples["type"][wildcards.sample]]["positive"]["R2_leading"],
   ODN_neg=lambda wildcards: config[samples["type"][wildcards.sample]]["negative"]["R2_leading"]
  shell: """
        cutadapt -j {threads} -G "negative={params.ODN_neg};max_error_rate=0;rightmost" -G "positive={params.ODN_pos};max_error_rate=0;rightmost" --discard-untrimmed  --rename='{{id}}_{{r2.adapter_name}} {{comment}}' -o {output.R1} -p {output.R2} {input.R1} {input.R2} > {log}
        
        cutadapt -j {threads} -G "negative={params.ODN_neg};max_error_rate=0;rightmost" -G "positive={params.ODN_pos};max_error_rate=0;rightmost" --discard-untrimmed  --rename='{{id}}_{{r2.adapter_name}} {{comment}}' -o {output.I3} -p {output.R2} {input.I3} {input.R2} > {log}
        
        """
        
        

  
  # Add UMI to read name
  
rule add_UMI:
    input: R1=rules.trim_ODN.output.R1, R2=rules.trim_ODN.output.R2,I3=rules.trim_ODN.output.I3
    output: R1=temp("01-trimming/{sample}_R1.ODN.UMI.fastq.gz"),
     R2=temp("01-trimming/{sample}_R2.ODN.UMI.fastq.gz"),
     I3=temp("01-trimming/{sample}_I3.ODN.UMI.fastq.gz")
    threads:6
    conda: "../01-envs/env_tools.yml"
    params: UMI=config["UMI_pattern"] ## bp in 3' of index to considere as UMI
    shell: """
        UMI_length=$(expr length {params.UMI})
        
        cutadapt -j {threads} -u -$UMI_length --rename='{{id}}_{{r1.cut_suffix}} {{comment}}' -o {output.I3} -p {output.R1} {input.I3} {input.R1}
        cutadapt -j {threads} -u -$UMI_length --rename='{{id}}_{{r1.cut_suffix}} {{comment}}' -o {output.I3} -p {output.R2} {input.I3} {input.R2}
        """


  # remove leading and trailing ODN sequences
rule trim_reads:
    input: R1=rules.add_UMI.output.R1, R2=rules.add_UMI.output.R2
    output: R1=temp("01-trimming/{sample}_R1.ODN.UMI.trimmed.fastq.gz"), 
     R2=temp("01-trimming/{sample}_R2.ODN.UMI.trimmed.fastq.gz")
    threads: 6
    log: "01-trimming/{sample}.trailing.log"
    conda: "../01-envs/env_tools.yml"
    message: "trimming ODN and adaptor sequences in reads"
    params: R2_trailing_pos=lambda wildcards: config[samples["type"][wildcards.sample]]["positive"]["R2_trailing"],
     R2_trailing_neg=lambda wildcards: config[samples["type"][wildcards.sample]]["negative"]["R2_trailing"],
     R1_trailing_pos=lambda wildcards: config[samples["type"][wildcards.sample]]["positive"]["R1_trailing"],
     R1_trailing_neg=lambda wildcards: config[samples["type"][wildcards.sample]]["negative"]["R1_trailing"]
    shell: """
        cutadapt -j {threads} -A "{params.R2_trailing_pos};min_overlap=6;max_error_rate=0.1" -A "{params.R2_trailing_neg};min_overlap=6;max_error_rate=0.1" -a "{params.R1_trailing_pos};min_overlap=6;max_error_rate=0.1" -a "{params.R1_trailing_neg};min_overlap=6;max_error_rate=0.1" -o {output.R1} -p {output.R2} {input.R1} {input.R2} > {log}
        """


  # make a size selection before mapping
rule filter_reads:
    input: R1=rules.trim_reads.output.R1, R2=rules.trim_reads.output.R2
    output:  R1="02-filtering/{sample}_R1.UMI.ODN.trimmed.filtered.fastq.gz", R2="02-filtering/{sample}_R2.UMI.ODN.trimmed.filtered.fastq.gz",
     R1short=temp("02-filtering/{sample}_R1.UMI.ODN.trimmed.tooshort.fastq.gz"), R2short=temp("02-filtering/{sample}_R2.UMI.ODN.trimmed.tooshort.fastq.gz")
    threads: 6
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
        conda: "../01-envs/env_tools.yml"
        message: "Aligning PE reads on genome with Bowtie2"
        params: index=lambda wildcards: config["genome"][samples["Genome"][wildcards.sample]]["index"],
         minFragLength=config["minFragLength"],
         maxFragLength=config["maxFragLength"]
        shell: """
                bowtie2 -p {threads} --no-unal -I {params.minFragLength} -X {params.maxFragLength} --dovetail --no-mixed --no-discordant --un-conc-gz 03-align/{wildcards.sample}_R%.UMI.ODN.trimmed.unmapped.fastq.gz   -x {params.index} -1 {input.R1} -2 {input.R2} -S {output.sam} 2> {log}
            """

elif (config["aligner"]  == "bwa" or config["aligner"]  == "Bwa") :
    rule alignOnGenome:
        input: R1=rules.filter_reads.output.R1, R2=rules.filter_reads.output.R2
        output: sam=temp("03-align/{sample}.UMI.ODN.trimmed.filtered.sam"), unfilterdsam=temp("03-align/{sample}.UMI.ODN.trimmed.unfiltered.sam")
        threads: 6
        log: "03-align/{sample}.UMI.ODN.trimmed.filtered.align.log"
        conda: "../01-envs/env_tools.yml"
        message: "Aligning PE reads on genome with BWA mem"
        params:  index=config["genome"][lambda wildcards:samples["Genome"][wildcards.sample]]["index"]
        shell: """
                #bwa mem -t {threads} {params.index} {input.R1} {input.R2} > {output.sam} 2> {log}
                bwa mem -t {threads} {params.index} {input.R1} {input.R2} > {output.unfilterdsam} 2> {log}
                samtools view -F 0x4 -F 0x8 -F 0x100 -F 0x800 -f 0x2 -b {output.unfilterdsam} > {output.sam}
            """
  
  
rule filter_alignments:
    input: sam=rules.alignOnGenome.output.sam
    output: list="03-align/{sample}_multi.txt",
    threads: 2
    conda: "../01-envs/env_tools.yml"
    message: "keep only alignments with high MAPQ or with equal score with best secondary alignment (multihits)"
    shell: """
        samtools view  {input} -e "((mapq < 30) && (mapq > 0) && ([AS] == [XS])) || (mapq >=20)" | cut -f1 | sort | uniq  > {output.list}
    """



# sort alignments by names (required for BEDPE conversion) and position (for viewing)
rule sort_aligned:
    input: sam=rules.alignOnGenome.output.sam, list = rules.filter_alignments.output.list
    output:  bam="03-align/{sample}.UMI.ODN.trimmed.filtered.sorted.filtered.bam", 
     bamPos="03-align/{sample}.UMI.ODN.trimmed.filtered.sorted.bam",
     bamName="03-align/{sample}.UMI.ODN.trimmed.filtered.sortedName.filtered.bam"
    threads: 6
    conda: "../01-envs/env_tools.yml"
    message: "Sort reads by name"
    shell: """
        samtools sort -@ {threads} {input.sam}  > {output.bamPos}

        samtools view -hb --qname-file {input.list}  {output.bamPos} > {output.bam}
        samtools index {output.bam}
        samtools sort -n -@ {threads} {output.bam} > {output.bamName}
    """



# convert BAM file to BEDPE and identify insertion point, read length
# filter alignment with MAPQ score > threshold in config file
# aggregate reads per fragment size and genomic coordinates
# add a cluster ID to group close IS (distance defined in the config file) 
rule call_IS:
    input: rules.sort_aligned.output.bamName
    output: tmp=temp("04-IScalling/{sample}.pebed"), 
     umi="04-IScalling/{sample}.reads_per_UMI_per_IS.bed"
    threads: 1
    conda: "../01-envs/env_tools.yml"
    log:
    params:  minMAPQ=config["minMAPQ"], UMI=config["UMI_pattern"]
    shell: """
    
        UMI_length=$(expr length {params.UMI})

        bedtools bamtobed -bedpe -mate1 -i {input} > {output.tmp}
        
        # read R2 contains the genome/ODN junction (use $10 of PE bed = mate 2)
        
        ##################################################################
        # count number of reads per UMI and per IS (pos and strand separated)
        ##################################################################

        awk 'BEGIN{{OFS="\\t";FS="\\t"}} ($8>={params.minMAPQ}) && ($1 == $4) {{split($7,a,"_"); if($10=="+") print $4,$5,$5,a[1],a[2],a[3],$8,$10,$3-$5; else print $4,$6-1,$6-1,a[1],a[2],a[3],$8,$10,$6-$2}}' {output.tmp} | sort -k1,1 -k2,3n -k6,6 -k5,5 -k8,8 | bedtools groupby -g 1,2,3,6,5,8 -c 4,7 -o count_distinct,median  > {output.umi}
        
        """


rule correct_UMI:
    input: bed=rules.call_IS.output.umi
    output: "04-IScalling/{sample}.reads_per_UMI_per_IS_corrected.bed"
    conda: "../01-envs/env_R4.3.2.yml"
    threads: 2
    params: hamming=config["UMI_hamming_distance"], filter = config["UMI_filter"],  UMI=config["UMI_pattern"], method = config["UMI_deduplication"]
    shell: """
        Rscript ../00-pipeline/correct_umi.R {input} {params.UMI} {params.filter} {params.hamming} {params.method} {output}
        """


rule Collapse_UMI_IS:
    input: rules.correct_UMI.output
    output: collapse="04-IScalling/{sample}.UMIs_per_IS_in_Cluster.bed"
    threads: 1
    conda: "../01-envs/env_tools.yml"
    log:
    params: window=config["ISbinWindow"], minMAPQ=config["minMAPQ"],minReadsPerUMI=config["minReadsPerUMI"],minUMIPerIS=config["minUMIPerIS"]
    shell: """
        
        ##################################################################
        # aggregate UMI per IS
        ##################################################################
        
        awk 'NR>1 && $7>{params.minReadsPerUMI}' {input} | sort -k1,1 -k2,2n -k3,3n -k5,5 -k 6,6 | bedtools groupby -g 1,2,3,5,6 -c 4,7,4,7,8  -o count_distinct,sum,collapse,collapse,median  | awk 'BEGIN{{OFS="\\t";FS="\\t"}} $6>{params.minUMIPerIS} {{print $1,$2,$3,"ID_"NR,$10,$5,$4,$6,$7,$8,$9}}' |   sort -k1,1 -k2,2n -k3,3n -k5,5 | bedtools cluster -d {params.window} > {output.collapse}
        """
    




    
rule get_fasta_around_is:
    input: bed=rules.Collapse_UMI_IS.output.collapse, length=["../02-ressources/{genome}.chrom.length".format(genome=genome) for genome in genomes_unique]
    output: cluster="04-IScalling/{sample}.cluster_slop.bed", fa="04-IScalling/{sample}.cluster_slop.fa"
    threads: 1
    conda: "../01-envs/env_tools.yml"
    params: fasta=lambda wildcards: config["genome"][samples["Genome"][wildcards.sample]]["fasta"], 
     chr_length=lambda wildcards: config["genome"][samples["Genome"][wildcards.sample]]["fasta"]+".fai",
     slop_size=config["slopSize"]
    shell: """
        # echo -e "#chromosome\tstart\tend\tclusterID\tN_orientations\tmedianMAPQ\tN_IS\tN_UMI\tN_reads" > {output.cluster}
        
        sort -k12,12n -k7 -k6 {input.bed} | bedtools groupby -g 12 -c 1,2,3,4,5,6,7,8,9 -o distinct,min,max,count_distinct,median,count_distinct,count_distinct,sum,sum |  awk 'BEGIN{{OFS="\\t"}}{{print $2,$3,$4,$1,$6,$5,$7,$8,$9,$10}}' | bedtools slop -i - -b {params.slop_size} -g {params.chr_length} > {output.cluster} 
        
        ## ORDER of columns
        ## chr start end clusterID MapQmedian nIS nOrientations nPCR nUMI, nReads
        
        
        bedtools getfasta -name -fi {params.fasta} -bed {output.cluster} > {output.fa}
        """ 




rule get_stats_fq:
    input: rules.merge_sampleName.output.R1, rules.trim_ODN.output.R1, rules.trim_reads.output.R1, rules.filter_reads.output.R1
    output: "05-Report/{sample}.stat"
    conda: '../01-envs/env_tools.yml'
    threads: 6
    shell: """
        seqkit stat -j {threads} -a -T {input} > {output}
        """



rule predict_offtarget_SWOffinder:
    input: 
    output: txt="06-offPredict/{gen}_{grna}{pam}.txt", csv="06-offPredict/{gen}_{grna}{pam}.csv"
    conda: '../01-envs/env_tools.yml'
    params: genome=lambda wildcards: config["genome"][wildcards.gen]["fasta"],
        gRNA=lambda wildcards:{wildcards.grna},
        maxE=config["max_edits_crRNA"],               #Max edits allowed (integer).
        maxM=config["max_edits_crRNA"],               #Max mismatches allowed without bulges (integer).
        maxMB=config["max_edits_crRNA"],             #Max mismatches allowed with bulges (integer).
        maxB=config["SWoffFinder"]["maxB"],               #Max bulges allowed (integer).
        window_size=config["SWoffFinder"]["window_size"], #The window size for choosing the best in a window
        bulges=config["tolerate_bulges"],
        PAM=lambda wildcards:wildcards.pam,
        SWoffFinder=config["SWoffFinder"]["path"]
    threads: 12
    shell: """
        echo {params.gRNA}{params.PAM} > {output.txt}
        
        if [ {params.bulges} = TRUE ] ; 
        then 
            bulge_size={params.maxB}; 
        else 
            bulge_size=0;
        fi

        
        java -cp {params.SWoffFinder}/bin SmithWatermanOffTarget.SmithWatermanOffTargetSearchAlign {params.genome} {output.txt} 06-offPredict/{wildcards.gen} {params.maxE} {params.maxM} {params.maxMB} $(echo $bulge_size) {threads} TRUE {params.window_size} {params.PAM} TRUE
        """
        
        
        


rule report_data:
    input: fasta=rules.get_fasta_around_is.output.fa, cluster=rules.get_fasta_around_is.output.cluster, bed=rules.Collapse_UMI_IS.output.collapse, statfq=rules.get_stats_fq.output,statal=rules.alignOnGenome.log
    output: "05-Report/{sample}.rdata"
    conda: "../01-envs/env_R4.3.2.yml"
    log:
    threads: 4
    params: gRNA_seq=lambda wildcards:samples["gRNA_sequence"][wildcards.sample], 
     gRNA_name=lambda wildcards:samples["gRNA_name"][wildcards.sample],
     PAM=lambda wildcards:samples["PAM_sequence"][wildcards.sample],
     offset=lambda wildcards:samples["Cut_Offset"][wildcards.sample],
     max_edits=config["max_edits_crRNA"],
     bulges=config["tolerate_bulges"]
    shell : """
        Rscript ../00-pipeline/multiple_alignments.R {input.fasta} {input.cluster} {input.bed} {params.gRNA_seq} {params.gRNA_name}  {params.PAM}  {params.offset} {params.max_edits} {params.bulges} {output}
        """        
        

rule annotate_sites:
    input: rdata=rules.report_data.output, annot=["../02-ressources/{genome}.rds".format(genome=genome) for genome in genomes_unique]
    output: "05-Report/{sample}_summary.xlsx"
    conda: "../01-envs/env_R4.3.2.yml"
    threads: 4
    params: species=lambda wildcards: samples["Genome"][wildcards.sample]
    shell: """
        Rscript ../00-pipeline/annotate_cuting_sites.R {params.species} {input.rdata} {output}
        """


rule report:
    input: summaries= ["05-Report/{sample}_summary.xlsx".format(sample=sample) for sample in samples_unique],
     predictions=expand("06-offPredict/{gen}_{grna}{pam}.csv", zip, grna=gRNAs, gen=genomes,pam=PAMs)
    output: report_rdata= "05-Report/report.rdata"
    conda: "../01-envs/env_R4.3.2.yml"
    threads: 1
    params: sampleInfo=config["sampleInfo_path"], 
     config=os.path.abspath(workflow.configfiles[0]), 
     minUMI_alignments_figure=config["minUMI_alignments_figure"],
     min_predicted_distance=config["min_predicted_distance"],
     max_clusters=config["max_clusters"]
    shell: """
        mkdir -p 05-Report/report-files/;
        Rscript ../00-pipeline/generate_tables.r "{input.summaries}" {params.sampleInfo} {params.config} "{input.predictions}" {params.max_clusters} {params.minUMI_alignments_figure} {params.min_predicted_distance}
        """
      
      
rule make_report:
    input:  rules.report.output.report_rdata
    output: report_html="05-Report/report.html"
    conda: "../01-envs/env_R4.3.2.yml"
    params: config_path=config_file_path
    threads: 1
    shell: """
        Rscript ../00-pipeline/publish_report.r {input} {params.config_path} {output}
    """

      
onsuccess:
    print("Workflow finished, no error")
    shell("cowsay -e \♥♥  You rock !! ")


onerror:
    print("An error occurred")
    shell(" cowsay -e \X\X  Houston, we have a problem !")


