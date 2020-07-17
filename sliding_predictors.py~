#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
from argparse import RawTextHelpFormatter
import os
import numpy as np

def argument_parser():

    parser = argparse.ArgumentParser(description=__doc__, prog='sliding_predictors.py', formatter_class=RawTextHelpFormatter)
    parser.add_argument("-s", "--seq_in", required=True, dest="seq_in",
                            help="File with sequence in FASTA format.")
    parser.add_argument("-i", "--increment", required=False, dest="increment", default=20, type=int,
                            help="Input alignment file. Multiple sequence alignment in one-line FASTA format.")
    parser.add_argument("-w", "--window_size", required=False, dest="window", default=200, type=int,
                            help="Input alignment file. Multiple sequence alignment in one-line FASTA format.")

    args = parser.parse_args() 

    seq_in = args.seq_in
    increment = args.increment
    window = args.window
    
    return seq_in, increment, window

def sliding_prediction():

    seqlen = len(seq)
    num_steps = int(((seqlen-window)/increment)+1)
    print(">sequence")
    print(seq)
    #["RNAfold", "ProbKnot", "Fold", "CentroidFold", "CONTRAfold", "IPknot"]
    
#    predictor = "IPknot"
    all_ss = ""

    for i in range(increment, window, increment):
        print(">first")

        after = '_'*(seqlen-i)

        segment = seq[:i]
        
        create_tempfile(segment)
        
        ss = prediction_method(predictor, 'seq.temp', segment)
        
        segment_ss = ss+after
        all_ss += segment_ss+"\n"
        print(ss+after)
        print(segment+after)
    
    for i in range(0, num_steps):
        print(">middle")

        before = "_"*(i*increment)
        segment = seq[i*increment:window+i*increment]
        after = '_'*(seqlen-(i*increment+len(segment)))

        create_tempfile(segment)
        ss = prediction_method(predictor, 'seq.temp', segment)
        
        segment_ss = before+ss+after
        all_ss += segment_ss+"\n"
        print(before+ss+after)
        print(before+segment+after)

    step = (i+1)*increment

    for i in range(step, seqlen, increment):
        print(">last")

        segment = seq[i:]    
        before = '_'*(seqlen-len(segment)) 

        create_tempfile(segment)
        ss = prediction_method(predictor, 'seq.temp', segment)
        
        segment_ss = before+ss
        all_ss += segment_ss+"\n"
        print(before+segment)
        print(before+ss)
    
    cons_temp = open('cons_temp.txt', 'w')
    cons_temp.write(all_ss)
    cons_temp.close()

    cmd = "find_SS_consensus cons_temp.txt 0.5"
    os.system(cmd)
    cmd = "mv single_consensus.txt "+seq_in.split(".")[0]+"_"+predictor+"_"+str(window)+"_"+str(increment)+".ss"   
    os.system(cmd)

def create_tempfile(segment):
    
    temp = open('seq.temp','w')
    temp.write(">temp\n")
    temp.write(segment)
    temp.close()


def read_files():

    with open(seq_in, 'r') as f:
        for line in f:
            if line[0] !=">":
                seq = line.strip().replace("T","U").replace("t","u").upper()  # sequence with SHAPE casette, read from the sequen
    return seq


def prediction_method(predictor, sequence, segment):

    sequencelen = len(segment)
    print(sequencelen)
    if predictor == "RNAfold":
        #RNAfold
        cmd = "RNAfold -p -d2 --noLP < %s" % (sequence)
        ss = os.popen(cmd).read().splitlines()[2].split(' ')[0]
    elif predictor == "ProbKnot":
        # RNAstructure ProbKnot
        cmd = "ProbKnot --sequence " + sequence + " ProbKnot.ct"
        os.system(cmd)
        cmd = "ct2dot ProbKnot.ct 1 ProbKnot.dot"
        os.system(cmd)
        f = open("ProbKnot.dot", 'r')
        ss = f.read().splitlines()[2]
        cmd = 'rm ProbKnot.ct ProbKnot.dot'
        os.system(cmd)
        f.close()
    elif predictor == "Fold":
        # RNAstructure Fold
        cmd = "Fold " + sequence + " Fold.ct"
        os.system(cmd)
        cmd = "ct2dot Fold.ct 1 Fold.dot"
        os.system(cmd)
        f = open("Fold.dot", 'r')
        ss = f.read().splitlines()[2]
        cmd = 'rm Fold.ct Fold.dot'
        os.system(cmd)
        f.close()
    elif predictor == "CentroidFold":
        # CentroidFold
        cmd = "centroid_fold " + sequence
        try:
            ss = os.popen(cmd).read().splitlines()[2].split(' ')[0]
        except:
            ss = '.'*sequencelen
    elif predictor == "CONTRAfold":
        # CONTRAfold
        cmd = "contrafold predict " + sequence
        ss = os.popen(cmd).read().splitlines()[2]
    elif predictor == "IPknot":
        # ipknot
        cmd = "ipknot " + sequence
        try:
            ss = os.popen(cmd).read().splitlines()[3]
        except:
            ss = '.'*(sequencelen)

    return ss

if __name__ == '__main__':

    predictors = ["RNAfold", "ProbKnot", "Fold", "CentroidFold", "CONTRAfold", "IPknot"]
    seq_in, increment, window = argument_parser()
    seq = read_files()
    for i in range(0,len(predictors)):
        predictor = predictors[i]
        print(predictor)
        sliding_prediction()
 
    os.system("rm seq.temp *.ps cons_temp.txt")
    