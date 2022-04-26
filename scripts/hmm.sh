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


#############Continuing

#Turning heuristic filters off
hmmsearch --max --noali --tblout positive-1.tmp bpti.hmm $db/positive-1.fasta

grep -v "#" positive-1.tmp |wc
'168    3192   24528' #Problem solved


#Redoing the hmmsearch with Automatic normalization
hmmsearch -Z 1 -domZ 1 --max --noali --tblout positive-1.tmp bpti.hmm $db/positive-1.fasta
grep -v "#" positive-1.tmp |awk '{print $1,$8}'|less
mv positive-1.tmp positive-1.match

hmmsearch -Z 1 --domZ 1 --max --noali --tblout positive-2.match bpti.hmm $db/positive-2.fasta
hmmsearch -Z 1 --domZ 1 --max --noali --tblout negative-1.match bpti.hmm $db/negative-1.fasta
hmmsearch -Z 1 --domZ 1 --max --noali --tblout negative-2.match bpti.hmm $db/negative-2.fasta


#Classifying our IDs as actual positives or actual negatives
grep -v "#" negative-1.match |awk '{print $1,$8,0}' >negative1.class
grep -v "#" negative-2.match |awk '{print $1,$8,0}' >negative2.class
grep -v "#" positive-1.match |awk '{print $1,$8,1}' >positive1.class
grep -v "#" positive-2.match |awk '{print $1,$8,1}' >positive2.class


#Finding the missing negatives (those eliminated by hmmsearch) and adding them back
comm -23 <(sort negative-1.txt) <(cut -d " " -f 1 negative1.class|sort) |wc
"154899"

comm -23 <(sort negative-1.txt) <(cut -d " " -f 1 negative1.class|sort) |awk '{print $1,10,0}' >>negative1.class
comm -23 <(sort negative-2.txt) <(cut -d " " -f 1 negative2.class|sort) |awk '{print $1,10,0}' >>negative2.class


wc negative1.class negative-1.txt


#Creating the set lists
cat positive1.class negative1.class >set-1.class
cat positive2.class negative2.class >set-2.class


#Optimization table
../scripts/accuracy.py set-1.class 0.001

for i in 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10; do ../scripts/accuracy.py set-1.class $i; done >opt-table-1.txt
for i in 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10; do ../scripts/accuracy.py set-2.class $i; done >opt-table-2.txt


../scripts/accuracy.py <(cat set-1.class set-2.class) 1e-10
