#!/bin/bash
#SBATCH --mem=150gb
#SBATCH --job-name=codeml_pipeline
#SBATCH --array=1-3336%10

start_time=$(date +%s)
echo "START TIME: $(date)"

source /home/slq4/anaconda3/bin/activate r-environment

cd selection # required files in selection directory
cp -n ../head.stats.table stats.table # for the stats table the headers are preset in this file
cp -n ../species.tree species.tree # add species tree to selection directory
#mkdir -p trees/ # create trees directory to save gene trees into tree directory for later comparison
#mkdir -p logs/ # create logs directory to save all output logs to log directory


LINE=$(head -"$SLURM_ARRAY_TASK_ID" filtered.protein.coding.transcript.IDs.txt | tail -1 ) # select line X corresponding to array of ID list

ID=$(echo $LINE | awk '{print $6}') #select human ID from LINE


echo "Parsing ${LINE}" #log the line bing parsed
#Run the R script with the current line as an argument
echo $LINE > ${ID}.output.txt
Rscript --save filter.fasta.R ${ID}.output.txt

cat ${ID}.human.fasta.txt > ${ID}.fasta.txt
cat ${ID}.sheep.fasta.txt >> ${ID}.fasta.txt
cat ${ID}.baboon.fasta.txt >> ${ID}.fasta.txt
cat ${ID}.mouse.fasta.txt >> ${ID}.fasta.txt
cat ${ID}.cow.fasta.txt >> ${ID}.fasta.txt
cat ${ID}.opossum.fasta.txt >> ${ID}.fasta.txt
cat ${ID}.pig.fasta.txt >> ${ID}.fasta.txt
cat ${ID}.rat.fasta.txt >> ${ID}.fasta.txt
cat ${ID}.chicken.fasta.txt >> ${ID}.fasta.txt
cat ${ID}.chimp.fasta.txt >> ${ID}.fasta.txt
cat ${ID}.gibbon.fasta.txt >> ${ID}.fasta.txt

rm ${ID}.human.fasta.txt
rm ${ID}.sheep.fasta.txt
rm ${ID}.baboon.fasta.txt
rm ${ID}.mouse.fasta.txt
rm ${ID}.cow.fasta.txt
rm ${ID}.opossum.fasta.txt
rm ${ID}.pig.fasta.txt
rm ${ID}.rat.fasta.txt
rm ${ID}.chicken.fasta.txt
rm ${ID}.chimp.fasta.txt
rm ${ID}.gibbon.fasta.txt

# Parse the first column from the current line
echo "Parsing ${ID} through pipeline script"
./pipeline.sh  ${ID}

rm ${ID}.output.txt
rm ${ID}.fasta.txt

cp -n ../species.tree trees/

end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "END TIME: $(date)"
echo "RUNTIME: ${runtime} seconds"

echo -e "${SLURM_ARRAY_TASK_ID}\t$(date)\t${runtime}" >> array_runtime_log.tsv
