#! /usr/bin/env python

import sys
import numpy as np


def get_cm(filename, th=1, eval_pos=-2, bin_class=-1):
    cm = np.zeros((2,2))
    with open(filename) as fh:
        for line in fh:
            line = line.rstrip()
            if not line: continue #ignore empty lines
            cols = line.split()
            pc = int(cols[bin_class])  #bin_class is the column position where the ACTUAL binary classification (+/- bpti) was assigned
            e = float(cols[eval_pos])  #eval_pos is the column position where the e-value is stored
            if pc == 0:
                if e >= th:
                    cm[0,0] += 1
                else:
                    cm[1,0] += 1
            elif pc == 1:
                if e >= th:
                    cm[0,1] += 1
                else:
                    cm[1,1] += 1
    return cm


def get_accuracy(cm):
    tp, tn, fp, fn = cm[1,1], cm[0,0], cm[1,0], cm[0,1]
    acc = (tp+tn)/(tp+tn+fp+fn)
    return acc


def get_mcc(cm):
    n = cm[0,0]*cm[1,1]-cm[0,1]*cm[1,0]
    d = np.sqrt((cm[0,0]+cm[0,1])*(cm[0,0]+cm[1,0])*(cm[1,1]+cm[0,1])*(cm[1,1]+cm[1,0]))
    return n/d





if __name__ == "__main__":
    filename = sys.argv[1]
    th = float(sys.argv[2])
    cm = get_cm(filename, th)
    acc = get_accuracy(cm)
    #print(cm, np.sum(cm), acc)
    print('TH:', th, 'Q2:', get_accuracy(cm), 'MCC:', get_mcc(cm), 'TN:', cm[0,0], 'FN:', cm[0,1], 'FP:', cm[1,0], 'TP:', cm[1,1])




"""Creating the confusion matrix
Pred(row)/Actual(col)
 0    1
0 TN FN
1 FP TP

Prediction definition
0 if e_val >= th
1 if e_val < th
"""
