#! /usr/bin/env python

import sys
import pandas as pd

def clean_csv(filename, a=2, b=3):
    rcsb = pd.read_csv(filename)
    rcsb2 = rcsb.iloc[:,[a,b]]
    lst_col = list(rcsb2.columns.values)
    rcsb2 = rcsb2.assign(**{lst_col[0]:rcsb2[lst_col[0]].str.split(',')}).explode(lst_col[0], ignore_index=False) #fixing comma separated chain names
    rcsb2['chain_identifiers'] = rcsb2[lst_col[1]]+':'+rcsb2[lst_col[0]].str.strip()
    return rcsb2


def chain_printer(df):
    for index,row in df.iterrows():
        print(row['chain_identifiers'])



if __name__ == "__main__":
    filename = sys.argv[1]
    rcsb2 = clean_csv(filename)
    chain_printer(rcsb2)
