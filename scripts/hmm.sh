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



### Splitting the fasta files

./cross-validation.py ../proteins/positive-1.txt ../../databases/uniprot-clean-pf00014 |less

./cross-validation.py ../proteins/positive-1.txt ../../databases/uniprot-clean-pf00014 >positive-1.fasta
./cross-validation.py ../proteins/positive-2.txt ../../databases/uniprot-clean-pf00014 >positive-2.fasta

./cross-validation.py ../proteins/negative-1.txt ../../databases/uniprot-NOT-pf00014.fasta >negative-1.fasta
./cross-validation.py ../proteins/negative-2.txt ../../databases/uniprot-NOT-pf00014.fasta >negative-2.fasta


### hmmsearch

hmmsearch --noali --tblout positive-1.match bpti.hmm positive-1.fasta
hmmsearch --noali --tblout positive-2.match bpti.hmm positive-2.fasta

mmsearch --noali --tblout negative-1.match bpti.hmm negative-1.fasta
hmmsearch --noali --tblout negative-2.match bpti.hmm negative-2.fasta

#parsing pid with e-values
grep -v "#" positive-1.match |awk '{print $1,$8}'|less
grep -v "#" positive-1.match |awk '{print $1,$8/168}'|less #e-value independent of the number of values

#oh oh, it seems one pid from the positive group dissapeared
grep -v "#" positive-1.match |wc
'167    3173   24382' #it appears that a +bpti protein had a high e-value and was automatically discarded, which would be a false negative.

#How do we find this missing protein?
comm -13 <(sort <(grep -v "#" positive-1.match |awk '{print $1}')) <(sort positive-1.txt)
'D3GGZ8'

#note: dataset .fasta files were moved to $db to allow for github push
