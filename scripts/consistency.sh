#Setting directories
shared=/media/sf_shared_linux/lab_bioinformatics/LB1_project;proj=~/projects/hmm_LB1/;db=~/projects/databases


#checking the consistency of the bult model with the training alignment
hmmsearch --noali bpti.hmm rep-cluster95_v2.fasta >bpti_consistency.txt
less bpti_consistency.txt
