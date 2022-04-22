#Setting directories
shared=/media/sf_shared_linux/lab_bioinformatics/LB1_project;proj=~/projects/hmm_LB1/;db=~/projects/databases

### preparing input for e-fold alignment
grep ">" rep-cluster95_v2.fasta |cut -d " " -f 1 |tr -d ">"|tr "_" ":" >sel_chains.txt

#This is extracted from the fasta file
'''
1brb:I
1fak:I
1t8l:B
1zr0:B
3byb:A
3wny:A
4dtg:K
4isl:B
4u32:X
5nx1:C
6gfi:C
6q61:A
6yhy:A

'''

#3tgi:I was manually added
cp sel_chains.txt $shared


#recovering pdb-e fold fasta alignment
cp $shared/alignment/fasta.seq ./struct_alignment.seq



#fixing the alignment to get it all in one line
awk '{if (substr($0,1,1)==">") {print "\n"$0} else {printf "%s",$0}}' ali3d.txt  #not necessary in my case

#NOTE: spaces manually eliminated between aligned sequences
