#Setting directories
shared=/media/sf_shared_linux/lab_bioinformatics/LB1_project;proj=~/projects/hmm_LB1/;db=~/projects/databases

#Printing the length of the alignment
grep -v -m 1 ">" struct_alignment.seq |awk '{print length}'
'63'

#writing a awk command to shorten the alignment
grep -v ">" struct_alignment.seq |awk '{print substr($0,4,length($0)-4)}'|awk '{print length}'
'59'

#writing a new command that does not eliminate the headers
awk '{if (substr($0,0,1)==">") {print $0} else {print substr($0,4,length($0)-4)}}' struct_alignment.seq |less

#redirecting into a new alignment file
awk '{if (substr($0,0,1)==">") {print $0} else {print substr($0,4,length($0)-4)}}' struct_alignment.seq >struct_alignment_corrected.seq

#Checking that 3 characters were removed from each sequence
grep -v ">" struct_alignment_corrected.seq |awk '{print length}'|less
'59'

mv struct_alignment.seq struct_alignment.old

: '
All good. the alignment was shortened. Now we proceed with the rest of the
pipeline without any changes (unless explicitly stated down below).

The shortening will now eliminate the first 3 characters instead of the first 2.
2 sequences without gaps is not enough for the model
'

#oops, the corresponding fasta files are already in $db
cp $proj/proteins/*ve-1.txt ../proteins_2/
cp $proj/proteins/*ve-2.txt ../proteins_2/
rm $proj/proteins_2/*ve*

set_pos1=$db/positive-1.fasta;set_pos2=$db/positive-2.fasta;set_neg1=$db/negative-1.fasta;set_neg2=$db/negative-2.fasta

hmmsearch --noali --tblout positive-1.match bpti.hmm $set_pos1
hmmsearch --noali --tblout positive-2.match bpti.hmm $set_pos2
hmmsearch --noali --tblout negative-1.match bpti.hmm $set_neg1
hmmsearch --noali --tblout negative-2.match bpti.hmm $set_neg2

#always the same discarded protein
comm -13 <(sort <(grep -v "#" positive-1.match |awk '{print $1}')) <(sort $proj/proteins/positive-1.txt)
'D3GGZ8'


#Redoing the hmmsearch with Automatic normalization
hmmsearch -Z 1 --domZ 1 --max --noali --tblout positive-1.match bpti.hmm $db/positive-1.fasta
hmmsearch -Z 1 --domZ 1 --max --noali --tblout positive-2.match bpti.hmm $db/positive-2.fasta
hmmsearch -Z 1 --domZ 1 --max --noali --tblout negative-1.match bpti.hmm $db/negative-1.fasta
hmmsearch -Z 1 --domZ 1 --max --noali --tblout negative-2.match bpti.hmm $db/negative-2.fasta


#Classifying our IDs as actual positives or actual negatives
grep -v "#" negative-1.match |awk '{print $1,$8,0}' >negative1.class
grep -v "#" negative-2.match |awk '{print $1,$8,0}' >negative2.class
grep -v "#" positive-1.match |awk '{print $1,$8,1}' >positive1.class
grep -v "#" positive-2.match |awk '{print $1,$8,1}' >positive2.class

#Finding the missing negatives (those eliminated by hmmsearch) and adding them back
comm -23 <(sort $proj/proteins/negative-1.txt) <(cut -d " " -f 1 negative1.class|sort) |wc
'154899'

comm -23 <(sort $proj/proteins/negative-1.txt) <(cut -d " " -f 1 negative1.class|sort) |awk '{print $1,10,0}' >>negative1.class
comm -23 <(sort $proj/proteins/negative-2.txt) <(cut -d " " -f 1 negative2.class|sort) |awk '{print $1,10,0}' >>negative2.class


#the analysis from the creation of positive-1.match was redone because positive1.class has 167 lines instead of 168, as expected
wc negative1.class negative2.class
: '
 283314  849942 3536357 negative1.class
 283315  849945 3535992 negative2.class
 566629 1699887 7072349 total'

 wc positive1.class positive2.class
: '
  168  504 2812 positive1.class
  168  504 2834 positive2.class
  336 1008 5646 total'

#all good

#creating sets
cat positive1.class negative1.class >set-1.class
cat positive2.class negative2.class >set-2.class

#Finding the best threshold for each set
for i in 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10 1e-11 1e-12; do ../scripts/accuracy.py set-1.class $i; done >opt-table-1.txt
for i in 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10 1e-11 1e-12; do ../scripts/accuracy.py set-2.class $i; done >opt-table-2.txt

#concatenated sets
for i in 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10 1e-11 1e-12; do ../scripts/accuracy.py <(cat set-1.class set-2.class) $i; done >opt-table-concatenated.txt

#finding the best thresholds
sort -gk6 opt-table-1.txt|less
sort -gk6 opt-table-2.txt|less
sort -gk6 opt-table-concatenated.txt|less


: ' The best thresholds are as follows:

1e-08 for set 1
1e-08 for set 2

1e-09 for concatenated


'


#Obtaining final alignment (that considered by the HMM model logo, in which only 57 characters are present)
grep -v -m 1 ">" struct_alignment_corrected.seq |awk '{print length}'
grep -v ">" struct_alignment_corrected.seq |awk '{print substr($0,1,length($0)-2)}'|awk '{print length}'
awk '{if (substr($0,0,1)==">") {print $0} else {print substr($0,1,length($0)-2)}}' struct_alignment_corrected.seq >struct_alignment_hmm_corrected.seq #eliminates last 2 chars in each sequence
grep -v ">" struct_alignment_hmm_corrected.seq |awk '{print length}'|less #checking the alignment length is correct
'57 chars, all good'

#simplifying header
awk -F ":" '{if (substr($0,0,1)==">") {print ">"$2":"substr($3,1,1)} else {print $0}}' struct_alignment_hmm_corrected.seq >struct_alignment_hmm_corrected_header.seq #leaves only pid and chain in header
grep -v ">" struct_alignment_hmm_corrected_header.seq |awk '{print length}'|less #checking the alignment length is the same
'57 chars, all good'

#adding consensus sequence to struct_alignment_hmm_corrected_header.seq
hmmemit -c bpti.hmm >> struct_alignment_hmm_corrected_header.seq
tmpfile=$(mktemp)
awk '{if (substr($0,0,1)==">") {print $0} else {print substr($0,1,57)}}' struct_alignment_hmm_corrected_header.seq > $tmpfile #ensuring that the length is mantained.
mv $tmpfile struct_alignment_hmm_corrected_header.seq

: 'In addition, the struct_alignment_hmm_corrected_header.seq file
was manually modified so the header of the consensus sequence is just
"consensus"'




### New procedure: aligning false positive and false negative sequences to model

#########
#Set 1 cross-validation threshold: 1e-08
echo '#Set 1 cross-validation threshold: 1e-08' >set_1_false_positives.txt
echo '#listing false positives:' >>set_1_false_positives.txt
echo '#Set 1 cross-validation threshold: 1e-08' >set_1_false_negatives.txt
echo '#listing false negatives:' >>set_1_false_negatives.txt

#Finding false positives
#awk '{if($NF==0 && $2<1e-08) {print $0}}' <(sort -gk2 set-1.class) |less
awk '{if($NF==0 && $2<1e-08) {print $0}}' <(sort -gk2 set-1.class) >>set_1_false_positives.txt

#finding false negatives
#awk '{if($NF==1 && $2>=1e-08) {print $0}}' <(sort -gk2 set-1.class) |less
awk '{if($NF==1 && $2>=1e-08) {print $0}}' <(sort -gk2 set-1.class) >>set_1_false_negatives.txt


##########
#Set 2 cross-validation threshold: 1e-08
echo '#Set 2 cross-validation threshold: 1e-08' >set_2_false_positives.txt
echo '#listing false positives:' >>set_2_false_positives.txt
echo '#Set 2 cross-validation threshold: 1e-08' >set_2_false_negatives.txt
echo '#listing false negatives:' >>set_2_false_negatives.txt

#Finding false positives
#awk '{if($NF==0 && $2<1e-08) {print $0}}' <(sort -gk2 set-2.class) |less
awk '{if($NF==0 && $2<1e-08) {print $0}}' <(sort -gk2 set-2.class) >>set_2_false_positives.txt

#finding false negatives
#awk '{if($NF==1 && $2>=1e-08) {print $0}}' <(sort -gk2 set-2.class) |less
awk '{if($NF==1 && $2>=1e-08) {print $0}}' <(sort -gk2 set-2.class) >>set_2_false_negatives.txt

###########
#Obtaining the ids to query
grep -v "#" set_1_false_negatives.txt |cut -d " " -f1|less

#Set 1 false positives
tmpfile=$(mktemp)
#../scripts/cross-validation.py <(grep -v "#" set_1_false_positives.txt |cut -d " " -f1) $db/uniprot-NOT-pf00014.fasta |less
../scripts/cross-validation.py <(grep -v "#" set_1_false_positives.txt |cut -d " " -f1) $db/uniprot-NOT-pf00014.fasta >$tmpfile
hmmsearch -Z 1 --domZ 1 --max bpti.hmm $tmpfile >bpti_set_1_false_positives.txt
rm -v $tmpfile

#Set 1 false negatives
tmpfile=$(mktemp)
../scripts/cross-validation.py <(grep -v "#" set_1_false_negatives.txt |cut -d " " -f1) $db/uniprot-clean-pf00014 >$tmpfile
hmmsearch -Z 1 --domZ 1 --max bpti.hmm $tmpfile >bpti_set_1_false_negatives.txt
rm -v $tmpfile

#Set 2 false positives
tmpfile=$(mktemp)
../scripts/cross-validation.py <(grep -v "#" set_2_false_positives.txt |cut -d " " -f1) $db/uniprot-NOT-pf00014.fasta >$tmpfile
hmmsearch -Z 1 --domZ 1 --max bpti.hmm $tmpfile >bpti_set_2_false_positives.txt
rm -v $tmpfile


#Set2 false negatives
tmpfile=$(mktemp)
../scripts/cross-validation.py <(grep -v "#" set_2_false_negatives.txt |cut -d " " -f1) $db/uniprot-clean-pf00014 >$tmpfile
hmmsearch -Z 1 --domZ 1 --max bpti.hmm $tmpfile >bpti_set_2_false_negatives.txt
rm -v $tmpfile


: 'Note:
A glaring error has been found when comparing the e-value of set 1 false negatives
when redoing the hmmsearch.

It appears that positive-1.match was done incorrectly because it has the
following metadata:

"# Option settings: hmmsearch --tblout positive-1.tmp --noali --max bpti.hmm /home/carloskez/projects/databases/positive-1.fasta"

This tells me that the "-Z 1 --domZ 1 --max" options were not included and, as
such, e-values are not normalized.


'


### Repeating accuracy analysis to obtain extended statistics

for i in 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10 1e-11 1e-12; do ../scripts/accuracy.py set-1.class $i; done >ext_opt-table-1.txt
for i in 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10 1e-11 1e-12; do ../scripts/accuracy.py set-2.class $i; done >ext_opt-table-2.txt

#concatenated sets
for i in 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10 1e-11 1e-12; do ../scripts/accuracy.py <(cat set-1.class set-2.class) $i; done >ext_opt-table-concatenated.txt

#finding the best thresholds
sort -gk6 ext_opt-table-1.txt|less
sort -gk6 ext_opt-table-2.txt|less
sort -gk6 ext_opt-table-concatenated.txt|less
