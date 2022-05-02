# Proteins_2: New alignment test

As identified by the consistency test, there are some alignment positions
that seem to add little information. Indeed, they result in low posterior
probabilities when running hmmsearch on the very training sequences. 

Therefore, a script will be created to correct this alignment.

Then, the entire HMM routine will be repeated on the resulting modified
alignment. 
