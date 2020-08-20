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
    
    all_ss = ""
    all_ss2=">sequence\n"+seq+"\n"


    for i in range(increment, window, increment):
        print(">1-"+str(i))
        after = '-'*(seqlen-i)

        segment = seq[:i]
        
        create_tempfile(segment)
        
        ss = prediction_method(predictor, 'seq.temp', segment)
        
        segment_ss = ss+after
        all_ss += segment_ss+"\n"
        all_ss2 += ">1-"+str(i)+"\n"+segment_ss+"\n"
        print(ss+after)
        #print(segment+after)
    
    for i in range(0, num_steps):
        
        before = "-"*(i*increment)
        
        segment = seq[i*increment:window+i*increment]
        after = '-'*(seqlen-(i*increment+len(segment)))
        start = len(before)+1
        end = start+window-1 
        print(">"+str(start)+"-"+str(end))
        create_tempfile(segment)
        ss = prediction_method(predictor, 'seq.temp', segment)
        
        segment_ss = before+ss+after
        all_ss += segment_ss+"\n"
        all_ss2 += ">"+str(start)+"-"+str(end)+"\n"+segment_ss+"\n"
        print(before+ss+after)

    step = (i+1)*increment

    for i in range(step, seqlen, increment):

        segment = seq[i:]    
        before = '-'*(seqlen-len(segment)) 
        
        start = len(before)+1
        print(">"+str(start)+"-"+str(seqlen))
        create_tempfile(segment)
        ss = prediction_method(predictor, 'seq.temp', segment)
        
        segment_ss = before+ss
        all_ss += segment_ss+"\n"
        all_ss2 += ">"+str(start)+"-"+str(seqlen)+"\n"+segment_ss+"\n"
        print(before+ss)
    
    cons_temp = open('cons_temp.txt', 'w')
    cons_temp.write(all_ss)
    cons_temp.close()

    cmd = "find_SS_consensus cons_temp.txt 0.5"
    os.system(cmd)
    cmd = "mv single_consensus.txt "+seq_in.split(".")[0]+"_"+predictor+"_"+str(window)+"_"+str(increment)+"_consensus.ss"   
    os.system(cmd)
    
    out = open(seq_in.split(".")[0]+"_"+predictor+"_"+str(window)+"_"+str(increment)+"_all.ss", 'w')
    out.write(all_ss2)
    out.close()

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


def spotrna_ncbp(sequence, ss):
#    print("clean spot-RNA")
    
    pair_list = get_pair_list(ss)
    ss_clean = remove_nc(sequence, pair_list, ss)
    
    return ss_clean
    
    

def get_pair_list(ss):
#    print("get pair list")
    
    db_list = [['(',')'],['[',']'],['<','>'],['{','}'],['A','a'],['B','b'],['C','c'],['D','d'],['E','e']]
    allowed_characters = '()[]<>{}AaBbCcDdEe.'

    stack_list =[]
    pairs_list =[]


    # stack-pop for all versions of brackets form the db_list    
    for i in range(0, len(db_list)):
        for c, s in enumerate(ss):
            if s == db_list[i][0]:
                stack_list.append(c)
            elif s == db_list[i][1]:
                if len(stack_list) == 0:
                    sys.exit("There is no opening bracket for nt position "+str(c+1)+'-'+ss[c])
                elif s == db_list[i][1]:
                    pairs_list.append([stack_list.pop(), c])
        if len(stack_list) > 0:
            err = stack_list.pop()
            sys.exit("There is no closing bracket for nt position "+str(err)+'-'+ss[err])

    pairs_list_clean = [x for x in pairs_list if x != []]

    return pairs_list

    

def remove_nc(sequence, pair_list, ss_nc):
#    print("removing non canonical base pairs")
    nt_dict = {"A":"U","U":"A","C":"G","G":"C"}

    pair_list_to_remove = []
    for i in range(0, len(pair_list)):
        nt1 = seq[pair_list[i][0]]
        nt2 = seq[pair_list[i][1]]
        if nt_dict[nt1] != nt2:
            pair_list_to_remove.append(pair_list[i])
    
    l = []

    ss_clean_l = list(ss_nc)
    
    if len(pair_list_to_remove) != 0:
        for i in range(0, len(pair_list_to_remove)):
            ss_clean_l[pair_list_to_remove[i][0]] = "." 
            ss_clean_l[pair_list_to_remove[i][1]] = "."

    ss_clean = ''.join(ss_clean_l)
    
    return ss_clean


    
def prediction_method(predictor, sequence, segment):

    sequencelen = len(segment)
    print(segment, "segment")
    #print(sequencelen)
    
    if predictor == "LinearFold_C":
        cmd = "cat " + sequence + " | linearfold"
        ss = os.popen(cmd)
        #quit()
        ss = os.popen(cmd).read().splitlines()[2].split(' ')[0]
    elif predictor == "LinearFold_V":
        cmd = "cat " + sequence + " | linearfold -V"
        ss = os.popen(cmd).read().splitlines()[2].split(' ')[0]
    elif predictor == "RNAfold":
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
    elif predictor == "ContextFold":
        cmd = "cp "+sequence+" "+PATH_TO_CONTEXTFOLD+". ; cd "+PATH_TO_CONTEXTFOLD+"; java -cp bin contextFold.app.Predict in:" +sequence
        os.system(cmd)
        cmd = "cp "+PATH_TO_CONTEXTFOLD+sequence+".pred ."        
        os.system(cmd)
        ss = os.popen("cat "+sequence+".pred").read().splitlines()[2]
    # for SPOT-RNA i need to add function removing non-cnaonical pairings
    elif predictor == "SPOT-RNA":
        cmd = "SPOT-RNA.py --inputs "+sequence+" --outputs './'"    
        os.system(cmd)
        cmd = "ct2dot temp.ct 1 temp.dot"
        os.system(cmd)
        ss_nc = os.popen("cat temp.dot").read().splitlines()[2]
        ss = spotrna_ncbp(segment, ss_nc)
        
    return ss

if __name__ == '__main__':

    PATH_TO_CONTEXTFOLD = "/home/fryzjer/Apps/ContextFold_1_00/"
    predictors = ['SPOT-RNA']
#    predictors = ["LinearFold_C","LinearFold_V","RNAfold", "ProbKnot", "Fold", "CentroidFold", "CONTRAfold", "IPknot"]
    seq_in, increment, window = argument_parser()
    seq = read_files()
    for i in range(0,len(predictors)):
        predictor = predictors[i]
        print(predictor)
        sliding_prediction()
 
    os.system("rm  *.ps cons_temp.txt seq.temp.pred")
    