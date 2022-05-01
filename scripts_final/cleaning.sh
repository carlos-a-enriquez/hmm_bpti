#!/bin/bash
pdb_e_fold=$1
rcsb_pdb=$2
pfam=$3

echo "Setting directories" >/dev/stderr
shared=/media/sf_shared_linux/lab_bioinformatics/LB1_project;proj=~/projects/hmm_LB1/

grep -E "^ *[0-9]" $1| awk '{if ($4>=3 && $5<=1) print $(NF-2), $NF}' $proj/proteins/pdb-e_fold_3tgi.txt >$proj/proteins/pdb-e_fold_3tgi_filtered.txt
lines=$(wc $proj/proteins/pdb-e_fold_3tgi_filtered.txt)
echo $lines "chains obtained from pdb-e fold"

$proj/scripts/adv_pdb_parsing.py .$2 >$proj/proteins/rcsb_pdb_cleaned.txt
lines=$(wc $proj/proteins/rcsb_pdb_cleaned.txt)
echo $lines "chains obtained from advanced pdb search"


$proj/scripts/pfam_parse.py $3 >$proj/proteins/pfam_pf00014.clean.txt
lines=$(wc $proj/proteins/pfam_pf00014.clean.txt)
echo $lines "chains obtained from pfam"

echo "Now, lets find the match of all 3 databases" >/dev/stderr

cut -b 8-13 pdb-e_fold_3tgi_filtered.txt |sort >pdb-e_final_list.txt #gets the prot-id:chain and sorts (unique: no repetitions, not needed in my case because I include the chains)
#by checking

#pdb-e fold vs pdb
comm -12 <(awk {'print toupper($0)'} pdb-e_final_list.txt) <(sort rcsb_pdb_cleaned.txt) >comm_pdb_list.txt

#comm_pdb vs pfam
comm -12 comm_pdb_list.txt <(sort -u pfam_pf00014.clean.txt) >three_db_match.txt

count=$(wc three_db_match.txt)
echo $count "chains were obtained from the 3 databases"


: '
blast_file=$1
scripts=$HOME/Desktop/lab_bioinformatics/scripts
tmpfile=`mktemp` #or "$(mktemp)"
#s= ls $scripts
#echo $s
echo "Cleaning blast file in" $tmpfile >/dev/stderr
grep -v "^#" $blast_file | cut -f 1,2,11 > $tmpfile
echo "Processing blast file" >/dev/stderr
$scripts/exercise7.py $tmpfile >$blast_file".bevalue

#$blast_file".tmp"
rm $tmpfile'
