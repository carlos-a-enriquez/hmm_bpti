#Setting directories
shared=/media/sf_shared_linux/lab_bioinformatics/LB1_project;proj=~/projects/hmm_LB1/;db=~/projects/databases

#Printing the length of the alignment
grep -v -m 1 ">" struct_alignment.seq |awk '{print length}'
'63'

#writing a awk command to shorten the alignment
grep -v ">" struct_alignment.seq |awk '{print substr($0,3,length($0)-3)}'|awk '{print length}'

#writing a new command that does not eliminate the headers
awk '{if (substr($0,0,1)==">") {print $0} else {print substr($0,3,length($0)-3)}}' struct_alignment.seq |less

#redirecting into a new alignment file
awk '{if (substr($0,0,1)==">") {print $0} else {print substr($0,3,length($0)-3)}}' struct_alignment.seq >struct_alignment_corrected.seq

#Checking that 3 characters were removed from each sequence
grep -v ">" struct_alignment_corrected.seq |awk '{print length}'|less
'60'

mv struct_alignment.seq struct_alignment.old

: '
All good. the alignment was shortened. Now we proceed with the rest of the
pipeline without any changes (unless explicitly stated down below).'
