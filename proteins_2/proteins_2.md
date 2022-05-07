# Proteins_2: New alignment test

As identified by the consistency test, there are some alignment positions
that seem to add little information. Indeed, they result in low posterior
probabilities when running hmmsearch on the very training sequences.

Therefore, a script will be created to correct this alignment.

The script is alignment_shortener.sh and it contains all the commands that were repeated.

Then, the entire HMM routine will be repeated on the resulting modified
alignment.


***Note:***
Note:
A glaring error has been found when comparing the e-value of set 1 false negatives
when redoing the hmmsearch.

It appears that positive-1.match was done incorrectly because it has the
following metadata:

"# Option settings: hmmsearch --tblout positive-1.tmp --noali --max bpti.hmm /home/carloskez/projects/databases/positive-1.fasta"

This tells me that the "-Z 1 --domZ 1 --max" options were not included and, as
such, e-values are not normalized.

To correct the error, after commit "451a6bca748c36727b0eb73cb2c65c0c8a9871bc",
I will rerun the following lines:

```bash
#Redoing the hmmsearch with Automatic normalization
hmmsearch -Z 1 --domZ 1 --max --noali --tblout positive-1.match bpti.hmm $db/positive-1.fasta
hmmsearch -Z 1 --domZ 1 --max --noali --tblout positive-2.match bpti.hmm $db/positive-2.fasta
hmmsearch -Z 1 --domZ 1 --max --noali --tblout negative-1.match bpti.hmm $db/negative-1.fasta
hmmsearch -Z 1 --domZ 1 --max --noali --tblout negative-2.match bpti.hmm $db/negative-2.fasta

```
