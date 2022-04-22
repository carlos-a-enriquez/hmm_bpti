#Setting directories
shared=/media/sf_shared_linux/lab_bioinformatics/LB1_project;proj=~/projects/hmm_LB1/;db=~/projects/databases

#Running hmmbuild
hmmbuild bpti.hmm struct_alignment.seq


#Benchmark Downloads

cp 'uniprot-database_(type_pfam+pf00014)+reviewed_yes.fasta.gz' $proj/proteins/uniprot-pf00014.fasta.gz
gunzip uniprot-pf00014.fasta.gz

cp 'uniprot-NOT+database_(type_pfam+pf00014)+reviewed_yes.fasta.gz' $proj/proteins/uniprot-NOT-pf00014.fasta.gz
mv uniprot-NOT-pf00014.fasta.gz $db
gunzip uniprot-NOT-pf00014.fasta.gz


#Cross-validation

grep ">" uniprot-clean-pf00014 |cut -d "|" -f 2 >$proj/proteins/list-clean-pf00014.txt #obtaining identifiers
grep ">" uniprot-NOT-pf00014.fasta |cut -d "|" -f 2 >$proj/proteins/list-not-pf00014.txt


#shuffling
sort -R list-clean-pf00014.txt >shuffle-clean-pf00014.txt

head -n 168 shuffle-clean-pf00014.txt >positive-1.txt
tail -n 168 shuffle-clean-pf00014.txt >positive-2.txt


python -c "print(566629//2)"
'283314'

sort -R list-not-pf00014.txt >shuffle-not-pf00014.txt
head -n 283314 shuffle-not-pf00014.txt >negative-1.txt
tail -n 283315 shuffle-not-pf00014.txt >negative-2.txt
