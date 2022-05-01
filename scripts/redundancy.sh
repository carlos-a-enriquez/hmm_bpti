#Setting directories
shared=/media/sf_shared_linux/lab_bioinformatics/LB1_project;proj=~/projects/hmm_LB1/;db=~/projects/databases

wget https://ftp.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt
mv three_db_match.txt comm_chains.dat


#fixing comm_chains.dat
for i in `cat comm_chains.dat`; do echo $i; done|less

sed -i'.old' 's/:/_/' comm_chains.dat   #replaces ":" with "_"

awk -F '_' '{new_var=tolower($1)"_"$2; print new_var}' comm_chains.dat >comm_chains_lower.dat #protein-id turns lower case


#checking the sequence data pdb_seqres.txt
less pdb_seqres.txt

'''
>101m_A mol:protein length:154  MYOGLOBIN
MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRVKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGNFGADAQGAMNKALELFRKDIAAKYKELGYQG
>102l_A mol:protein length:165  T4 LYSOZYME
MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAAKSELDKAIGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAALINMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRVITTFRTGTWDAYKNL
>102m_A mol:protein length:154  MYOGLOBIN
MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKAGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGNFGADAQGAMNKALELFRKDIAAKYKELGYQG

'''
#as you can see, comm_chains_lower.dat now matches the format that we need

#extracting all sequences matching with our list of protein-id:chain
for i in `cat comm_chains_lower.dat`; do grep -A 1 ">"$i $db/pdb_seqres.txt ; done >ss.fasta &


### online redundancy analysis
#continue analysis at http://weizhong-lab.ucsd.edu/cdhit-web-server/cgi-bin/index.cgi?cmd=cd-hit
cp ss.fasta $shared


#results done
cp *.fas.1 $proj/proteins/rep-cluster95.fasta

#checking the final sequences
less rep-cluster95.fasta
#oops, one sequence appears too long

#checking long sequence
grep -A 1 "^>4bnr" rep-cluster95.fasta |grep -v "^>4bnr" |wc
'1       1     101'                                             #in fact, it is too long. I will delete it.

grep ">" rep-cluster95.fasta |awk '{if(substr($3,8,length($3))+0<100) {print $0}}'|wc  #This is an automated way to check for sequences that are too long (checks the 'length' field in header)

#sequence 4bnr deleted... rep-cluster95_v2.fasta

#Final set: 26 lines, 13 sequences
