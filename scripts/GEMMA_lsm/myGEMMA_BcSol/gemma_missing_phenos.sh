#!/bin/bash
#SBATCH -D /home/maktaylo/fibr_gwas
#SBATCH -o /home/maktaylo/fibr_gwas/slurm-logs/mark-stdout-%j.txt
#SBATCH -e /home/maktaylo/fibr_gwas/slurm-logs/mark-stderr-%j.txt
#SBATCH -J gemma
#SBATCH --mail-type=END
#SBATCH --mail-user=maktaylor@ucdavis.edu

echo "Starting at `date`"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors."
echo "Current working directory is `pwd`"
echo "SLURM_NODELIST = $SLURM_NODELIST"
echo "SLURM_NODE_ALIASES = $SLURM_NODE_ALIASES"
echo "SLURM_NNODES = $SLURM_NNODES"
echo "SLURM_TASKS_PER_NODE = $SLURM_TASKS_PER_NODE"
echo "SLURM_NTASKS = $SLURM_NTASKS"
echo "SLURM_JOB_ID = $SLURM_JOB_ID"

cd /home/maktaylo/fibr_gwas

readarray -t traitnames < /home/maktaylo/fibr_gwas/phenotypes/phenotypes_list.txt

traitname=${traitnames[$SLURM_ARRAY_TASK_ID]}

mkdir $traitname

path=$traitname
cd $traitname
fibrpath="/home/maktaylo/fibr_gwas/"
traitdash="/"
traitpath=$fibrpath$traitname
traitout=$traitpath$traitdash$traitname

#name the alternative phenotype file
phenodirec="/home/maktaylo/fibr_gwas/phenotypes/" #name directory holding the alternative phenotype file
phenotxtext="_alternate_pheno.txt" #name the alternative phenotype file extension
phenotxt=$phenodirec$traitname$phenotxtext #name the alternative phenotype file name

#make binary PLINK files to input into GEMMA
cp /home/maktaylo/fibr_gwas/1001genomes/k2029.* $traitpath #copy the ped and map files from 1001genomes directory to the phenotype directory since it will otherwise be deleted
cd $traitpath
  /home/maktaylo/plink/plink-1.07-x86_64/plink --file k2029 --noweb --map3 --pheno $phenotxt --make-bed --out $traitout

#this will tell us the allele frequencies in my sample
#/home/maktaylo/plink/plink-1.07-x86_64/plink --file $traitname --noweb --map3 --freq

###########################################################################################
###########################################################################################
#gemma calls 
output="_output"
###########################################################################################
#run Bayesian sparse linear mixed model  (BSLMM)
cd $traitpath #transfer to directory where the PLINK output is!
out=$traitname$output
/home/maktaylo/GEMMA2/gemma -bfile $traitname -bslmm 1 -maf 0.01 -o $out

################################################################################################
#use GEMMA to make relatedness matrix for use in LMM & trait prediction
/home/maktaylo/GEMMA2/gemma -bfile $traitname -gk 1 -o $traitname
kmatext=".cXX.txt" #this is the file extension for the relatedness matrix created by gemma
fibrpath="/home/maktaylo/fibr_gwas/"
kmat=$fibrpath$path

cd output #must find the relatedness matrix in the output directory that GEMMA creates
cp *.cXX.txt $kmat #copy relatedness matrix to directory where PLINK files are

################################################################################################
#run EMMA (LMM in GEMMA)
cd ..
kfile=$traitname$kmatext
lmmout="_output_lmm"
lmout=$traitname$lmmout
/home/maktaylo/GEMMA2/gemma -bfile $traitname -k $kfile -lmm 4 -maf 0.01 -o $lmout

#the no population structure correction
################################################################################################
#run EMMA (LINEAR MODEL not LINEAR MIXED MODEL in GEMMA) WITHOUT KINSHIP CORRECTION
lmmout="_output_lmm_noK"
lmo_noK_out=$traitname$lmmout
/home/maktaylo/GEMMA2/gemma -bfile $traitname -lm 4 -maf 0.01 -o $lmo_noK_out

################################################################################################
#predict traits from BSLMM
pred=$traitname"_predicted"
epmfile=$traitpath"/output/"$traitname"_output.param.txt"
emufile=$traitpath"/output/"$traitname"_output.log.txt"
ebvfile=$traitpath"/output/"$traitname"_output.bv.txt"
kfile=$traitpath"/output/"$traitname".cXX.txt"

/home/maktaylo/GEMMA2/gemma -bfile $traitname -epm $epmfile -emu $emufile -ebv $ebvfile -k $kfile -predict 1 -o $pred

################################################################################################
#clean out giant phenotype/genotype files that were made files that have been made
#rm *ped
#rm *map
#rm *bim
#rm *bed
################################################################################################



