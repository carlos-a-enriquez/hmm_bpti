#! /usr/bin/env python

import sys
import pandas

def clean_csv(filename, a=2, b=3):
    rcsb = pd.read_csv(filename)
    rcsb2 = rcsb2.iloc[:,[a,b]]



if __name__ == "__main__":
    filename = sys.argv[1]
