#!/usr/bin/env python
import sys


def parse_fasta(fasta):
    d={}
    with open(fasta) as f:
        for line in f:
            if line.find('>')==0:  #The header starts with the ">"
                pid=line.split('|')[1] #obtains the identifier
            else:
                d[pid]=d.get(pid,'')+line.rstrip() #If this is not the header, then it is the sequence and most be added to the value of the previously found pid
    return d


if __name__ == '__main__':
    ids=sys.argv[1]
    fasta=sys.argv[2]
    d = parse_fasta(fasta)
    lid=open(ids).read().rstrip().split('\n') #will create a string out of the pid list, then it will split it into a list by line
    for pid in lid:  #check the identifiers. Is it in my dictionary? If yes, return the sequence in fasta format. This gives us the sequences from the list of identifiers
        if pid in d:
            print (">"+pid)
            print(d[pid])
