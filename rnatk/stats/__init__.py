#!/usr/bin/env python
#TODO: add doctest or unittest to each function

"""
rnatk.stats contains several python functions useful for \
statistical analysis of RNA secondary structure
"""

__all__ = ['randomXijDict', 'numRandomList', 'kernelXij', 'searchXijc',
           'fig_plot','Mij', 'equLen', 'equLen_file', 'similarity',
           'entropy', 'MIxy', 'MIxyz']
__author__ = "Ivan Castro"
__copyright__ = "Copyright 2015, rnatk Project"
__license__ = "GPL"
__version__ = "0.2.0"
__email__ = "sivanc7@gmail.com"

import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import pylab
from scipy import stats
from math import log

def equLen(*args):
    """
    Check for equal length among several strings (parameters).
    If any string has different length with respect to another one, the system exits.
    Otherwise the system continues.
    Modules required:
    - sys
    Usage: <string-1> ... <string-n>
    """
    if all(len(args[0]) == len(_arg) for _arg in args[1:]):
        pass
    else:
        raise ValueError('All sequences muste be the same length')
        
def equLen_file(inFile):
    """
    Evaluate whether one sequence from the file is different to the others. The function requieres
    a list with elements Bio.SeqRecord.SeqRecord.
    Module required:
    - equLen
    Usage: <inFile>
    """
    for item_A in range(len(inFile)):
        for item_B in range(len(inFile)):
            if item_B-item_A >= 1:
                equLen(str(inFile[item_A].seq), str(inFile[item_B].seq))

def similarity(str1, str2):
    """
    Returns the similarity index between two string. Both string must have equal lenght.
    Module required:
    - equLen (from rnatk.stats)
    Usage: <string 1> <string 2>
    """
    equLen(str1, str2)
    matches = 0.0
    for index in range(len(str1)):
        if str1[index] == str2[index]:
            matches += 1
    similarity = matches/len(str1)
    return similarity

def Mij(pos1, pos2, logBase=2, mode='basepair'):
    """
    Calculate the mutual information between two RNA positions
    in an alignment. Each position must be converted into a string
    of RNA. The mode defines if the function must compute only those
    sequences obsered in the basepairs (mode=basepairs) or if the function
    must compute the MI no matters the basepairs (mode=free).
    Modules required:
    - sys
    - math
    - equLen
    Usage: <sequence1> <sequence2> <logBase> <mode (default=basepair)>
    """
    equLen(pos1,pos2)
    MI = 0
    if mode == 'basepair':
        allow = ['CG','GC','AU','UA','GU','UG']
    if mode == 'free':
        allow = ['AA', 'AC', 'AG', 'AU',
                 'CA', 'CC', 'CG', 'CU',
                 'GA', 'GC', 'GG', 'GU',
                 'UA', 'UC', 'UG', 'UU']
    for a in allow:
        Fxy, Fx, Fy = 0, 0, 0
        for b in range(len(pos1)):
            if pos1[b]+pos2[b] == a:
                Fxy += 1
        if Fxy > 0:
            Fx = pos1.count(a[0])/float(len(pos2))
            Fy = pos2.count(a[1])/float(len(pos2))
            if (Fxy/float(len(pos1)))*math.log(((Fxy/float(len(pos1)))/(Fx*Fy)),int(logBase)) <= 0:
                MI += 0
            else:
                MI += (Fxy/float(len(pos1)))*math.log(((Fxy/float(len(pos1)))/(Fx*Fy)),int(logBase))
        else:
            MI += 0
    return MI

def entropy(logbase, *seqs):
    """
    This function find the entropy for one or several sequences (joint entropy)
    Equation used here:
    H(X1, ..., Xn) = -sum{x1}...sum{xn} P(x1, ..., xn)log[P(x1, ..., xn)]
    Module required:
    - math
    Usage: <logbase> <seq1> ... <seqn>
    """
    equLen(*seqs)
    joint_list = []
    for index_A in range(len(seqs[0])):
        paired = ''
        for index_B in range(len(seqs)):
            paired += seqs[index_B][index_A]
        joint_list.append(paired)
    set_pair = set(joint_list)
    Hjoint = 0
    for elem in set_pair:
        prob = joint_list.count(elem)/float(len(joint_list))
        Hjoint -= prob*math.log(prob,logbase)
    return Hjoint

def MIxy(seq1, seq2, logbase):
    """
    Calculate the mutual information between two DNA/RNA sequences
    using the entropy and joint entropy. Both sequence must have equal
    lenght.
    Equation used here:
    MIxy = H(X)+H(Y)-H(X,Y)
    Usage: <seq1> <seq2> <logbase>
    """
    equLen(seq1, seq2)
    if (entropy(logbase, seq1) == 0) or (entropy(logbase, seq2) == 0):
        return 0
    else:
        return entropy(logbase, seq1)+entropy(logbase, seq2)-entropy(logbase, seq1, seq2)

def MIxyz(seq1, seq2, seq3, logbase=2):
    """
    Calculate the mutual information among three DNA/RNA sequences
    using the entropy and joint entropy. All the sequences must have equal
    lenght.
    Equation used here:
    MIxyz = H(X)+H(Y)+H(Z)-[H(X,Y)+H(Y,Z)+H(X,Z)]+H(X,Y,Z)
    Usage: <seq1> <seq2> <seq3> <logbase>
    """
    equLen(seq1, seq2, seq3)
    if (entropy(logbase, seq1) == 0) or (entropy(logbase, seq2) == 0) or (entropy(logbase, seq3) == 0):
        return 0
    else:
        return entropy(logbase, seq1)+entropy(logbase, seq2)+entropy(logbase, seq3)-\
             (entropy(logbase, seq1, seq2)+entropy(logbase, seq2, seq3)+entropy(logbase, seq1, seq3))+\
             entropy(logbase, seq1, seq2, seq3)

def numRandomList(times, lower, upper):
    """
    This fuction creates a list with a determinate amount
    of random float numbers that varying from a lower bound to
    an upper bound.
    Modules required:
    - numpy (as np)
    Usage: <times> <lower> <upper>
    """
    list = []
    for a in range(int(times)):
        list.append(np.random.uniform(float(lower), float(upper)))
    return list

def randomXijDict(lenSeq, upper, bioSense=True):
    """
    This function builds a dictionary of (lenSeq**2) keys. Each key
    is associated with a random number that varies from 0 to an upper
    bound.
    Modules required:
    - numpy (as np)
    Usage: <lenSeq> <upper> 
    """
    dict = {}
    for a in range(int(lenSeq)):
        for b in range(int(lenSeq)):
            if bioSense:
                if b-a>3:
                    dict[a,b] = np.random.uniform(0, float(upper))
            else:
                if b-a>=1:
                    dict[a,b] = np.random.uniform(0, float(upper))
    return dict

def kernelXij(XijList):
    """
    Creates a Kernel Density Function from a list of values.
    Modules required:
    - stats (from scipy)
    - numpy (as np)
    Usage: <list>
    """
    KDE = stats.gaussian_kde(np.array(XijList, dtype=np.float))
    return KDE

def searchXijc(KDE, critical, lower, upper):
    """
    Given a KDE, this function search the integral from <lower> to
    <upper>, where <lower> is fixed, that results in <critical>.
    Module required:
    - numpy (as np)
    - stats (from scipy)
    Usage: <KDE> <critical value> <lower limit> <upper limit>
    """
    stop = True
    search = lower
    while stop:
        if KDE.integrate_box_1d(search+0.00001,upper) > critical:
            search += 0.0001
        else: stop = False
    return search, KDE.integrate_box_1d(search,upper)

def fig_plot(array, lower, upper, fileName, xlabel=''):
    """
    Save a png file for a KDE.
    Module required:
    - numpy (as np)
    - stats (from scipy)
    - matplotlib.plot (as plt)
    - pylab
    Usage: <array> <lower limit> <upper limit> <outFile> <type of covariance (string)>
    """
    x1 = np.array(array, dtype=np.float)
    kde1 = stats.gaussian_kde(x1)
    kde2 = stats.gaussian_kde(x1, bw_method = 'silverman')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x1, np.zeros(x1.shape), 'b+', ms=20)
    x_eval = np.linspace(lower, upper, num=200)
    ax.plot(x_eval, kde1(x_eval), 'k-')
    plt.ylabel('Density')
    plt.xlabel(xlabel)
    pylab.savefig(fileName+'.png')
    plt.clf()
