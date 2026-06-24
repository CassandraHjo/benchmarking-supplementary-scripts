#!/bin/bash

#SBATCH --array=1-20                    
#SBATCH --nodes=1                            
#SBATCH --ntasks=1                           
#SBATCH --mem=2G                          
#SBATCH --time=01:00:00                      
#SBATCH --job-name=readMachine 
#SBATCH --output=readMachine_%j_%a.log

echo "I am running on the NODE $SLURM_NODELIST"
echo "I am running with $SLURM_CPUS_ON_NODE cpus"

echo "Starting $SLURM_JOB_ID at"
date


#################
### Settings
###
metadata_file="metadata.tsv"                                                               # Path to metadata file
amplicon_folder="amplicon_fasta/10000_species_10_dynamic"                                  # Path to folder containing per-community amplicon FASTA files                                   
out_folder="10000_species_10_dynamic/fastq"                                                # Path to output folder for simulated FASTQ files 
iss_app="apptainer exec /path/to/correct/software/insilicoseq:2.0.1--pyh7cba7a3_0.sif"     # Apptainer command to invoke the InSilicoSeq container  
tmp_folder="tmp"                                                                           # Path to temporary working folder for intermediate files
error_model="miseq"                                                                        # Sequencing error model passed to InSilicoSeq

###################
### Initiating
###
amplicon_file_column=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++){f[$i] = i}}{print $(f["amplicon_file"])}' $metadata_file)   # Extract the amplicon_file column from the metadata TSV
r1_file_column=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++){f[$i] = i}}{print $(f["R1_file"])}' $metadata_file)               # Extract the R1_file column from the metadata TSV
r2_file_column=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++){f[$i] = i}}{print $(f["R2_file"])}' $metadata_file)               # Extract the R2_file column from the metadata TSV
readpairs_column=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++){f[$i] = i}}{print $(f["readpairs"])}' $metadata_file)           # Extract the readpairs column from the metadata TSV

line=$(($SLURM_ARRAY_TASK_ID+1))                                                                                        # Convert array task ID to 1-indexed row number (offset by 1 to skip the header)
amplicon_file=$amplicon_folder/$(echo $amplicon_file_column | awk -vidx=$line '{print $idx}')                           # Retrieve the amplicon filename for this array job
r1_file=$(echo $r1_file_column | awk -vidx=$line '{print $idx}')                                                        # Retrieve the R1 output filename for this array job
r1_file=$(basename $r1_file .gz)                                                                                        # Strip the .gz extension; pigz will re-compress after simulation
r2_file=$(echo $r2_file_column | awk -vidx=$line '{print $idx}')                                                        # Retrieve the R2 output filename for this array job
r2_file=$(basename $r2_file .gz)                                                                                        # Strip the .gz extension; pigz will re-compress after simulation
readpairs=$(echo $readpairs_column | awk -vidx=$line '{print $idx}')                                                    # Retrieve the number of read pairs to simulate for this community

if [ ! -d $tmp_folder ]
then
  mkdir $tmp_folder         # Create the temporary folder if it does not already exist                                                         
fi
if [ ! -d $out_folder ]
then
  mkdir $out_folder         # Create the output folder if it does not already exist
fi

module load R/4.4.2         # Load the required R module


echo "*** Simulating $readpairs read pairs based on the amplicons in $amplicon_file"


#######################################
### Prepping the temporary files
###
echo "*** Creating temporary files with R..."
Rscript -e "suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(microseq))
args <- commandArgs(trailingOnly = T)
amplicon.file <- args[1]
tmp.folder <- args[2]
dataset <- args[3]
readpairs <- as.numeric(args[4])
cat('*** Reading fasta file with amplicons...\n')
amplicons.fa <- readFasta(amplicon.file)  |>                                                              # Read the amplicon FASTA file for this community
  mutate(abundance = str_extract(Header, 'read_abundance=.+'))  |>                                        # Extract the relative read abundance value embedded in each FASTA header
  mutate(abundance = str_remove(abundance, 'read_abundance='))  |>                                        # Remove the key prefix to isolate the numeric abundance value
  mutate(abundance = as.numeric(abundance)) |>                                                            # Convert the extracted abundance string to a numeric value
  mutate(Header = word(Header, 1))                                                                        # Retain only the first whitespace-delimited field of the header as the sequence ID
cat('*** Found ', nrow(amplicons.fa), 'amplicons...\n')
tt <- table(sample(amplicons.fa\$Header, size = readpairs, replace = T, prob = amplicons.fa\$abundance))  # Draw readpairs sequence IDs with replacement, weighted by abundance
ttb <- data.frame(Header = names(tt),
                  n_reads = 2 * as.numeric(tt))                                                           # Multiply sampled counts by 2 to account for paired-end reads (R1 + R2)
amplicons.fa <- amplicons.fa |>
  left_join(ttb, by = 'Header') |>                                                                        # Attach the sampled read counts to each amplicon
  filter(n_reads > 0)                                                                                     # Discard amplicons that received no reads in this simulation
cat('*** Amplicons with reads:', nrow(amplicons.fa), '\n')
amplicons.fa  |>  
  select(Header, n_reads)  |>  
  write.table(sep = '\t',
              file = file.path(tmp.folder, str_c('readcount_', dataset, '.txt')),
              quote = F,
              col.names = F,
              row.names = F)                                                                              # Write the read count table to a TSV file for use by InSilicoSeq
cat('*** Extracting the paired-ends of the amplicons...\n')
amplicons.fa |> 
  mutate(Sequence = str_sub(Sequence, 1, 302))  |>                                                        # Trim each amplicon to the first 302 bp to represent the forward (R1) read region
  writeFasta(out.file = file.path(tmp.folder, str_c('left_', dataset, '.fa')))                            # Write the forward read sequences to a temporary FASTA file
amplicons.fa  |>  
  mutate(Sequence = str_sub(Sequence, -302, -1))  |>                                                      # Extract the last 302 bp of each amplicon to represent the reverse (R2) read region
  mutate(Sequence = reverseComplement(Sequence))  |>                                                      # Reverse-complement the extracted region to match R2 orientation
  writeFasta(out.file = file.path(tmp.folder, str_c('right_', dataset, '.fa')))" $amplicon_file $tmp_folder $SLURM_ARRAY_TASK_ID $readpairs  # Write the reverse read sequences to a temporary FASTA file, passing runtime arguments


#################
### Simulations
###
echo "*** Simulating reads with InSilicoSeq"
echo "*** The R1 reads..."
$iss_app iss generate \
  --cpus 1 \
  --quiet \
  --model $error_model \
  --genomes $tmp_folder/left_$SLURM_ARRAY_TASK_ID.fa \               # Use the forward read FASTA as input
  --readcount_file $tmp_folder/readcount_$SLURM_ARRAY_TASK_ID.txt \  # Provide per-amplicon read counts
  --sequence_type amplicon \
  --output $tmp_folder/left_$SLURM_ARRAY_TASK_ID                     # Write simulated R1 reads to the temporary folder

  
echo "*** The R2 reads..."
$iss_app iss generate \
  --cpus 1 \
  --quiet \
  --model $error_model \
  --genomes $tmp_folder/right_$SLURM_ARRAY_TASK_ID.fa \              # Use the reverse read FASTA as input
  --readcount_file $tmp_folder/readcount_$SLURM_ARRAY_TASK_ID.txt \  # Provide the same per-amplicon read counts as for R1
  --sequence_type amplicon \
  --output $tmp_folder/right_$SLURM_ARRAY_TASK_ID                    # Write simulated R2 reads to the temporary folder



#######################
### Post-processing
###
echo "*** Post-processing..."
cp $tmp_folder/left_$SLURM_ARRAY_TASK_ID\_R1.fastq $out_folder/$r1_file   # Copy the simulated R1 reads to the output folder
cp $tmp_folder/right_$SLURM_ARRAY_TASK_ID\_R1.fastq $out_folder/$r2_file  # Copy the simulated R2 reads to the output folder (InSilicoSeq always names output _R1)
pigz $out_folder/$r1_file                                                 # Compress the R1 output file with pigz
pigz $out_folder/$r2_file                                                 # Compress the R2 output file with pigz

rm $tmp_folder/left_$SLURM_ARRAY_TASK_ID\_R1.fastq    # Remove temporary forward R1 file
rm $tmp_folder/left_$SLURM_ARRAY_TASK_ID\_R2.fastq    # Remove temporary forward R2 file (unused but generated by InSilicoSeq)
rm $tmp_folder/right_$SLURM_ARRAY_TASK_ID\_R1.fastq   # Remove temporary reverse R1 file (used as R2 output)
rm $tmp_folder/right_$SLURM_ARRAY_TASK_ID\_R2.fastq   # Remove temporary reverse R2 file (unused but generated by InSilicoSeq)
rm $tmp_folder/left_$SLURM_ARRAY_TASK_ID.fa            # Remove temporary forward amplicon FASTA
rm $tmp_folder/right_$SLURM_ARRAY_TASK_ID.fa           # Remove temporary reverse amplicon FASTA
rm $tmp_folder/readcount_$SLURM_ARRAY_TASK_ID.txt      # Remove temporary read count file
rm $tmp_folder/left_$SLURM_ARRAY_TASK_ID.*.vcf         # Remove VCF files generated as a side effect by InSilicoSeq
rm $tmp_folder/right_$SLURM_ARRAY_TASK_ID.*.vcf        # Remove VCF files generated as a side effect by InSilicoSeq
echo "done!"

echo "Ending $SLURM_JOB_ID at"
date