#!/usr/bin/env python
#TODO: filter position with H(X)=0 as option (add the change to changeLog file)
#TODO: add doctest or unittest to each function
#TODO: make the function generalized
"""
Performs an analysis of Multivariate Mutual Information (MMI) (three variables) to the inputted fasta file.
"""

__all__ = ['mult_mut_info']
__author__ = "Ivan Castro"
__copyright__ = "Copyright 2015, rnatk Project"
__license__ = "GPL"
__version__ = "0.2.0"
__email__ = "sivanc7@gmail.com"

from os.path import expanduser
import rnatk
import rnatk.stats import MIxyz
from rnatk.seq import randomSeq
from rnatk.stats import (kernelXij, searchXijc)

def mult_mut_info(inFile, logbase=2, repeats, siglevel=0.05, bioSense=False):
	if (type(inFile) == str) and (str(type(inFile[0])) != "<class 'Bio.SeqRecord.SeqRecord'>"):
		InFile = rnatk.fasta.openFasta(inFile)
	matrix = rnatk.fasta.transpose(inFile)
	dictResults = {}
	for indexA in range(len(matrix)):
		for indexB in range(len(matrix)):
			for indexC in range(len(matrix)):
				if bioSense:
					if (indexB-indexA>3) and (indexC-indexB>3):
						value = MIxyz(matrix[indexA], matrix[indexB], matrix[indexC], logbase)
						dictResults[indexA+1, indexB+1, indexC+1] = value
				else:
					if indexA<indexB<indexC:
						value = MIxyz(matrix[indexA], matrix[indexB], matrix[indexC], logbase)
						dictResults[indexA+1, indexB+1, indexC+1] = value
	MIxyzRandom = []
	for k in range(repeats):
		seq3Tpl = (randomSeq(len(InFile), 'RNA'), randomSeq(len(InFile), 'RNA'), randomSeq(len(InFile), 'RNA')) # make a 3-tuple of random sequences
		MIxyzRandom.append(MIxyz(seq3Tpl[0], seq3Tpl[1], seq3Tpl[2], logbase)) # compute and store the MIxyz value of the three random sequences
	MIxyzCrit = searchXijc(kernelXij(MIxyzRandom), siglevel, -0.5, 2.0)[0] # here is defined the MIxyz critical value
	return dictResults
