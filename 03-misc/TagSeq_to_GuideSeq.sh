## convert TagSeq library to GuideSeq compatible library

# example
  # prefetch SRR13258532 SRR13258532 SRR13258533 SRR13258534
  # fasterq-dump --split-3 SRR13258531 SRR13258532 SRR13258533 SRR13258534


## requires cutadapt and sra-toolkit

for library in *_1.fastq; do ## change pattern if necessary
  lib="${library%_1.fastq}"
  echo "processing " $lib
  
  conda activate BI_tools #environment with cutadapt
  ## get reads with linker in R1
  cutadapt -j 24 -g "^NNNNNNNNNNNNNNNNACTACTAATACGACT"  \
    -O 31 \
    -e 0 \
    --no-indels \
    --action none \
    -o $lib"_R1temp.fastq.gz" \
    -p $lib"_R2.fastq.gz" \
    --discard-untrimmed \
    $lib"_1.fastq" $lib"_2.fastq" > $lib"_filter.log"
  
  ## get I2 from R1 leading sequence (8 +8 and reverse) 
  cutadapt --length 16 -j 12 $lib"_R1temp.fastq.gz" | seqkit seq -r -o $lib"_I2.fastq.gz"
  
  ## Trim UMI+BC+adaptor from R1
  cutadapt -j 12 -g "^NNNNNNNNNNNNNNNNACTACTAATACGACT" \
    -O 31 \
    -e 0 \
    --no-indels \
    -o $lib"_R1.fastq"  \
    $lib"_R1temp.fastq.gz"  > $lib"_trimming_R1.log"
  
  ## Make a fake I1 file as it doesn't exist in the SRA deposite
  # retreive the barcode from the SRA deposite
  conda activate sra-tool # environment with SRA toolkit
  index=$(sra-stat --quick $lib | awk 'BEGIN{FS="|"} {print $2}')
  echo $lib" barcode I1: "$index
  ~/projects/add_pattern_to_sequence.sh $index replace  $lib"_R1.fastq" > $lib"_I1.fastq"

  conda activate BI_tools #environment with cutadapt
  seqkit seq -j 12 -o $lib"_R1.fastq.gz" $lib"_R1.fastq"
  seqkit seq -j 12 -o $lib"_I1.fastq.gz" $lib"_I1.fastq"
done


rm *temp*

cat *_R1.fastq.gz > Undetermined_R1.fastq.gz
cat *_R2.fastq.gz > Undetermined_R2.fastq.gz
cat *_I1.fastq.gz > Undetermined_I1.fastq.gz
cat *_I2.fastq.gz > Undetermined_I2.fastq.gz

rm SRR*fastq*