#!/bin/bash

# usage: ./insert_template.sh template.pdb output_directory/
# creates a new pdb_seqres.txt fasta file and a mmcif_files/ folder
# containing a single .cif PDB file. These can be passed to alphafold
# though the --pdb_seqres_database_path and --template_mmcif_dir flags respectively

export PATH=$PATH:/gpfswork/rech/uzu/commun/lib/bin
module load gcc/9.3.0

input=$1
inbase=$(basename $input .pdb)
outdir=$2
outbase=9999
output=$(dirname $input)/$(basename $input .pdb).cif

#cif-tools package needs to be installed
echo "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1 " > $input.2
grep -v CRYST1 $input | grep -v "REMARK" | grep -v TER | grep "\S" >> $input.2
mv $input.2 $input

pdb2cif $input $output

head -3 $output > $output.2
echo "_pdbx_audit_revision_history.revision_date" >> $output.2
echo "2100-01-01" >> $output.2
sed -n '4,$p' $output >> $output.2

sed -i 's/NON-POLYMER/peptide/g' $output.2
mkdir -p $outdir/mmcif_files
mv $output.2 $outdir/mmcif_files/$outbase.cif
pdb2fasta $input | sed "s/$inbase:\([A-Z]\)\s/$outbase\_\1 mol:protein length:/g" >> $outdir/pdb_seqres.txt
