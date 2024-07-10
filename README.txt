### GuideSeq pipeline
### 2024
### GENETHON

Copy the content of this folder to the remote server in a dedicated folder (ie ./guideseq)

Move multiplexed I1,I2,R1,R2 fastq.gz files to subfolder (ie ./guideseq/XXXXXX)

copy guideSeq_GNT.yml config file to (ie ./guideseq/XXXXXX) and edit :
- fastq files path
- genome reference fastq, annotation and index path
- alignment and filtering values 

copy sampleInfo.csv to (ie ./guideseq/XXXXXX) and edit with libraries informations


## FOLDER structure 
guideseq (or another name)
    |-00-pipeline
        |-guideSeq_GNT.smk
        |-merge_I1_I2.py
        |-multiple_alignments.R
        |-samples.schema.yaml
    |-01-envs
        |-env_R4.3.2.yml
        |-env_guideSeq.yml
    |-projetXXXXXX
        |-R1.fastq.gz
        |-R2.fastq.gz
        |-I1.fastq.gz
        |-I2.fastq.gz
    |-02-ressources
        |-guideSeq_GNT.yml
        |-sampleInfo.csv


run the pipeline:
snakemake -s ../00-pipeline/guideSeq_GNT.smk -k -j 12 --use-conda -n

snakemake -s ../00-pipeline/guideSeq_GNT.smk -k -j 12 --use-conda -n