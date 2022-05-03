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
hmmsearch -Z 1 -domZ 1 --max --noali --tblout positive-1.match bpti.hmm $db/positive-1.fasta
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
sort -gk6 opt-table-2.txt|less


: ' The best thresholds are as follows:

1e-09 for set 1
1e-08 for set 2

1e-09 for concatenated


'
