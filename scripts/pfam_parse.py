#! /usr/bin/env python


import sys
import pandas as pd


def pfam_clean(filename, a=2, b=3):
    df = pd.read_excel(filename).fillna(method='ffill', axis=0)
    df = df.iloc[:,a:b+1]
    lst_col = list(df.columns.values)
    df['id:chain'] = df[lst_col[0]].str.strip()+":"+df[lst_col[1]].str.strip()
    return df

def chain_printer(df):
    for index,row in df.iterrows():
        print(row['id:chain'])






if __name__ == "__main__":
    filename = sys.argv[1]
    df = pfam_clean(filename)
    chain_printer(df)
