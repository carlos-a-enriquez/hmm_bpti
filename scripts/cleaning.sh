#Setting directories
shared=/media/sf_shared_linux/lab_bioinformatics/LB1_project;proj=~/projects/hmm_LB1/

#Getting the files from pdb-e fold
mv pdb-e_fold_3tgi.txt ~/projects/hmm_LB1/

#How do we extract the useful rows?
less -S pdb-e_fold_3tgi.txt


"""
PDBe Fold v2.59. (src3) 14 Apr 2014 result file.

                  RESULT SUMMARY

  ##      Q-score      P-score     Z-score   RMSD    Nalgn Nsse Ngaps Seq-%     Nmd Nres-Q Nsse->
   1           1       16.18       11.93     0.000    56    4    0           1    0     56      >
   2      0.9986       13.58        10.9     0.112    56    4    0           1    0     56      >
   3      0.9984       13.58        10.9     0.121    56    4    0           1    0     56      >
   4      0.9984       13.58        10.9     0.122    56    4    0           1    0     56      >
   5      0.9975       13.13       10.71     0.151    56    4    0           1    0     56      >
   6      0.9959       10.09       9.336     0.192    56    4    0           1    0     56      >
   7      0.9923       12.25       10.33     0.264    56    4    0      0.8929    0     56      >


"""

grep -E "^ *[0-9]" pdb-e_fold_3tgi.txt |less -S #this gets the useful lines

#A possible command that could be used to obtain the protein ids that match a certain z score and rmsd
grep -E "^ *[0-9]" $proj/proteins/pdb-e_fold_3tgi.txt| awk '{if ($4>=3 && $5<=1) print $(NF-2), $NF}' $proj/proteins/pdb-e_fold_3tgi.txt >$proj/proteins/pdb-e_fold_3tgi_filtered.txt

#filters
#z-score >=3
#rmsd <=1

#final file obtained


#Advanced PDB search
#Processing csv file

###Extracting from pdb advanced search
mv rcsb_pdb_custom_report_20220412080009.csv $proj/proteins


#pandas package used to extract useful data from .csv file
#check jupyter notebook for details
./adv_pdb_parsing.py ../proteins/rcsb_pdb_custom_report_20220412080009.csv |wc #299 rows as expected
./adv_pdb_parsing.py ../proteins/rcsb_pdb_custom_report_20220412080009.csv >../proteins/rcsb_pdb_cleaned.txt





##### Now extracting data from pfam

'''
Note:
The structure data can be copied and pasted directly from the website.
Then, the cleaning procedure involves filling NA rows, exporting as a csv,
selecting the protein-id and chain columns and creating the new concatenated
column.

Omid has saved us some time by obtaining the clean csv. However, we still need
to add the concatanated column. Then, we must print that column into a text file,
which will be used downstream to compare all 3 datasets.
'''

pip install openpyxl
./pfam_parse.py ../proteins/pfam_raw.xlsx >../proteins/pfam_pf00014.clean.txt



#### Now, lets find the match of all 3 databases

cut -b 8-13 pdb-e_fold_3tgi_filtered.txt |sort >pdb-e_final_list.txt #gets the prot-id:chain and sorts (unique: no repetitions, not needed because I include the chains)
#by checking

#pdb-e fold vs pdb
comm -12 <(awk {'print toupper($0)'} pdb-e_final_list.txt) <(sort rcsb_pdb_cleaned.txt) >comm_pdb_list.txt

#comm_pdb vs pfam
comm -12 comm_pdb_list.txt <(sort -u pfam_pf00014.clean.txt) >three_db_match.txt

###Results: 144 id:chain pairs


echo "a, b" |tr -d " "|tr "," "\"
