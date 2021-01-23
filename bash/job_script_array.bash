#!/bin/bash
#-----------------------------Mail address-----------------------------rtffffy
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=a.l.vermeer2@students.uu.nl
#-----------------------------Required resources-----------------------
#SBATCH -N 1
#SBATCH --tasks-per-node 24
#SBATCH -t 12:00:00

#-----------------------------Environment, Operations and Job steps----
# load modules 
module load 2020
module load GDAL
module load R

cd $TMPDIR

input_scratch=$TMPDIR'/input_dir'
input_home=$HOME'/input_dir'

# Prepare input data (copy input data from home folder to scratch for quick access)
cp -r $input_home $TMPDIR

input_list=$(find $input_scratch -type f -name "*.grd")

echo $input_list

echo $SLURM_ARRAY_TASK_ID

input=($input_list)
infile=$(basename ${input[${SLURM_ARRAY_TASK_ID}-1]})

echo $infile 
	
outfile=${infile%.grd}_output.grd
echo $outfile 
	
rscript=$HOME/R_scripts/run_resIndSpatial_$infile'.R' 
echo $rscript
	
# Copy the default rscript and rename it to rscript
cp $HOME'/R_scripts/run_resIndSpatial.R' $rscript
	
sed -i "s|infile_rscript|'$infile'|g" $rscript
sed -i "s|outfile_rscript|'$outfile'|g" $rscript
	
# Run application 
R --no-save < $rscript 

# Remove input data from temporary directory 
rm -r $input_scratch 

# Copy data from scratch to home
cp -r $TMPDIR/* $HOME'/output_dir'
