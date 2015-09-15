#!/usr/bin/env python
#TODO: add doctest or unittest to each function

"""
rnatk.pdb contains several python fucntions useful to DNA/RNA sequences \
manipulation
"""

__all__ = []
__author__ = "Ivan Castro"
__copyright__ = "Copyright 2015, rnatk Project"
__license__ = "GPL"
__version__ = "0.2.0"
__email__ = "sivanc7@gmail.com"

import math
import glob
import copy
import subprocess
import matplotlib.pyplot as plt
import moderna
import rnatk

def distance3D_2points(p1, p2):
    """
    Return the distance between two points in a 3D space.
    Module required:
    - math
    Usage: <point1> <point2>
    """
    return math.sqrt(sum((p1 - p2)**2 for p1, p2 in zip(p1, p2)))

def rmsd_calc(distances):
    """
    Return the root mean square deviation given several distances
    Module required:
    - math
    Usage: <distances>
    """
    return math.sqrt((1.0/len(distances))*sum(dist**2 for dist in distances))

def get_pdb_from_folder(pathFolder):
    """
    Return a list with all the path+names of the files contained under the \
    given folder
    Module required:
    - glob
    Usage: <pathFolder>
    """
    return glob.glob(pathFolder)

def get_atom_info_from_pdb(pdb_file):
    """
    Return a dictionary that contains info about each ATOM in the pdb file.
    The index is [<num atom>, <num residue> <chain>].
    Usage: <pdb_file>
    """
    infoDict = {}
    if type(pdb_file) == str:
        pdbFile = open(pdb_file, 'r').readlines()
    else:
        pdbFile = pdb_file
    for elem in pdbFile:
        if elem[0:6].strip(' ') == 'ATOM':
            infoDict[(int(elem[6:11].strip(' ')), int(elem[22:26].strip(' ')), elem[21:22].strip(' '))] = {
                'ATOM':elem[0:6].strip(' '),
                'ATOM number':int(elem[6:11].strip(' ')),
                'type ATOM':elem[12:16].strip(' '),
                'alternate location':elem[16:17].strip(' '),
                'residue':elem[17:20].strip(' '),
                'chain':elem[21:22].strip(' '),
                'residue number':int(elem[22:26].strip(' ')),
                'x position':float(elem[30:38].strip(' ')),
                'y position':float(elem[38:46].strip(' ')),
                'z position':float(elem[46:54].strip(' ')),
                'occupancy':float(elem[54:60].strip(' ')),
                'temperature':float(elem[60: 66].strip(' ')),
                'element':elem[76:78].strip(' '),
                'charge':elem[78:80].strip(' ')
                }
    return infoDict

def get_chain_info(pdb_file, chain='A'):
    """
    This function returns a dictionary that contains info about each ATOM
    of the specified chain from a pdb file. The index is
    [<num atom>, <num residue> <chain>]
    Usage: <pdb_file> <chain='A'>
    """
    pdb_atoms = get_atom_info_from_pdb(pdb_file)
    pdb_atoms_chain = []
    towrite = ''
    result = {}
    for elem in pdb_atoms:
        if elem[2] == chain:
            pdb_atoms_chain.append(elem)
            pdb_atoms_chain.sort()
    else:
        for elem in pdb_atoms_chain:
            result[elem] = pdb_atoms[elem]
        return result

def rmsd_pdb(atom_ref, number_residue, chain, pdb_files):
    """
    This function return the RMSD for one equivalent atom among two or several PDB files.
    The pdb_files parameter must be a tuple o list containing the path to the
    PDBs files or a dictionary with pdb information.
    All the residues must be aligned.
    Usage: <atom_ref> <number_residue> <chain> <pdb_files>
    """
    distances_list = []
    if type(pdb_files[0]) != dict:
        for index in range(len(pdb_files)):
            pdb_files[index] = get_chain_info(pdb_files[index], chain)
    for index in range(len(pdb_files)):
        for subindex in range(len(pdb_files)):
            if subindex-index >= 1:
                for key in pdb_files[index]:
                    for subkey in pdb_files[subindex]:
                        if (subkey[1] == number_residue) and (subkey[2] == chain) and (key[1] == number_residue) and (key[2] == chain):
                            if (pdb_files[index][key]['type ATOM'] == atom_ref) and (pdb_files[subindex][subkey]['type ATOM'] == atom_ref):
                                point1 = (pdb_files[index][key]['x position'], pdb_files[index][key]['y position'], pdb_files[index][key]['z position'])
                                point2 = (pdb_files[subindex][subkey]['x position'], pdb_files[subindex][subkey]['y position'], pdb_files[subindex][subkey]['z position'])
                                distance = distance3D_2points(point1, point2)
                                distances_list.append(distance)
    if distances_list:
        return rmsd_calc(distances_list)

def get_clean_sequence_rna(pdbFile):
    """
    With this function, the residues beloging to the pdb that are represented
    by a dot (.) are changed by the undefinided nucleotide 'N'. Finally, the
    fuction returns that sequence
    Module required:
    - moderna
    """
    clean_sequence = ''
    pdb = moderna.load_model(pdbFile)
    sequence = str(pdb.get_sequence())
    for residue in sequence:
        if residue in ['A', 'G', 'U', 'C']:
            clean_sequence += residue
        if residue == '.':
            clean_sequence += 'N'
    return clean_sequence

def filter_atoms_by_ref(pdb_files, atom_ref="C1'"):
    """
    The aim of this function is to filter the atoms by its type.
    The function returns a list that can be used to write a new pdb file
    Module required:
    - copy
    """
    pdb_files_copy = copy.copy(pdb_files)
    if type(pdb_files_copy) == str:
        pdb_files_tem = []
        pdb_files_tem.append(pdb_files_copy)
        pdb_files_copy = pdb_files_tem
    for index in range(len(pdb_files_copy)):
        filter_lines = []
        try:
            pdb_files_copy[index] = open(pdb_files_copy[index], 'r').readlines()
        except:
            pass
        for line in pdb_files_copy[index]:
            if line[12:16].strip(' ') == atom_ref:
                filter_lines.append(line)
        pdb_files_copy[index] = filter_lines
    return pdb_files_copy

def plot_rmsd_pdbFiles(atom_ref, from_residue, until_residue, chain, pdb_files):
    """
    Plot the RMSD value for the specified residues from its atom of reference.
    It is recommended to use first the function filter_atoms_by_ref in order to save
    time. If the pdb residues are not aligned, it is recommended to use the
    function align_residues_pdb which perform the filter_atoms_by_ref function
    and align the residues. The object returned by each function can be used
    directly here in the parameter pdb_files.
    Module required:
    - matplotlib.pyplot as plt
    """
    rmsd_list = []
    for elem in range(from_residue, until_residue+1):
        print elem
        rmsd_result = rmsd_pdb(atom_ref, elem, chain, pdb_files)
        rmsd_list.append(rmsd_result)
    print rmsd_list
    plt.ylabel('RMSD value')
    plt.xlabel('positions')
    plt.grid(True)
    plt.plot(list(range(from_residue, until_residue+1)), rmsd_list, 'k.-')
    plt.show()

def get_seq_from_pdb(pdb_file):
    """
    This function returns the nucleotide sequence from a pdb file
    """
    seq = ''
    pdb_file = filter_atoms_by_ref(pdb_file, "C1'")
    for index in range(len(pdb_file[0])):
        if int(pdb_file[0][index][22:26].strip(' ')) == index+1:
            if pdb_file[0][index][17:20].strip(' ') == 'UNK':
                seq += 'N'
            else:
                seq += pdb_file[0][index][17:20].strip(' ')
    return seq

def make_fasta_seqs_from_pdbs(pdb_files, outFile):
    """
    This functions build a fasta file from a set (list) of pdb_files.
    The paramater pdb_files must be a list with the path+name of each
    single pdb file. Within each path+name, the function will look for a
    GenBank code access which will be use it as name into each sequences in the fasta
    file.
    Module required:
    - rnatk
    """
    seqList = []
    nameList = []
    for pdb in pdb_files:
        seq = get_seq_from_pdb(pdb)
        seqList.append(seq)
        try:
            name = re.findall(r'[A-Z]{1,3}[0-9]{3,9}.[0-9]{1,9}.[0-9]{1}?',pdb)[0]
        except:
            try:
                name = re.findall(r'[A-Z]{1,3}[0-9]{3,9}', pdb)[0]
            except:
                print 'Warning: Impossible to found a GenBank code access for %s' %(pdb)
                name = pdb
        nameList.append(name)
    rnatk.fasta.makeFasta(outFile, nameList, seqList)

def align_residues_pdb(pdb_files, aligner='clustalo', atom_ref="C1'"):
    """
    To use this funcion, the pdb_files must be a list with the path+name of each
    single pdb file.
    This functions returns an aligned representation (from all the pdb files) for
    each single pdb. This returned object can be used to plot its RMSD (with the
    plot_rmsd_pdbFiles function). To rebuild each pdb file with the aligned
    number of residue, use the function build_aligned_pdb.
    Module required:
    - subprocess
    - copy
    - rnatk
    """
    make_fasta_seqs_from_pdbs(pdb_files, expanduser('~')+'/tmp_pdbFasta')
    subprocess.check_call([aligner, '-i', expanduser('~')+'/tmp_pdbFasta', '-o', expanduser('~')+'/tmp_pdbFasta_align'])
    aligned = rnatk.fasta.openFasta(expanduser('~')+'/tmp_pdbFasta_align')
    nameLstRef, seqLstRef = rnatk.fasta.getNameSeq(aligned)
    for seq_index in range(len(seqLstRef)):
        seqLstRef[seq_index] = list(seqLstRef[seq_index])
    subprocess.check_call(['rm', expanduser('~')+'/tmp_pdbFasta', expanduser('~')+'/tmp_pdbFasta_align'])
    pdb_files_org = copy.copy(pdb_files)
    pdb_files = filter_atoms_by_ref(pdb_files, atom_ref)
    for index in range(len(pdb_files)):
        residue_index = 0
        times = 0
        while residue_index < len(pdb_files[index]):
            if (seqLstRef[index][residue_index-times] == 'N') and (pdb_files[index][residue_index-times] == 'UNK'):
                pass
            elif seqLstRef[index][residue_index-times] == pdb_files[index][residue_index-times]:
                pass
            elif (seqLstRef[index][residue_index-times] != pdb_files[index][residue_index-times][17:20].strip(' ')) and not (seqLstRef[index][residue_index-times] == 'N') and not (pdb_files[index][residue_index-times] == 'UNK'):
                for new_index in range(residue_index-times, len(pdb_files[index])):
                    refer = str(int(pdb_files[index][new_index][22:26].strip(' '))+1)
                    while len(refer) != 4:
                        refer = ' '+refer
                    to_replace = pdb_files[index][new_index][0:22]+refer+pdb_files[index][new_index][26:80]+'\n'
                    pdb_files[index][new_index] = to_replace
                times += 1
                residue_index -= 1
                seqLstRef[index].remove('-')
            residue_index += 1
    return pdb_files

def RNA_file_tertiary_interactions(path_to_rnaview, pdb_file):
    """
    This is a hidden function.
    The RNAview executable should be in the /bin folder.
    The pdb_file parameter must be a string indicating the path and name
    of the pdb file.
    """
    if type(pdb_file) == str:
        process = subprocess.Popen([path_to_rnaview, pdb_file])
        process.wait()

def RNA_folder_tertiary_interactions(path_to_rnaview, path_to_folder):
    """
    This is a hidden function.
    The RNAview executable should be in the /bin folder.
    The list must contain string items indicating the path and name
    of the pdb files.
    """
    list_pdb_files = get_pdb_from_folder(path_to_folder)
    for elem in list_pdb_files:
        RNA_file_tertiary_interactions(path_to_rnaview, elem)

def get_list_num_res(pdb, atom_ref="C1'", filter_atoms=True):
    """
    This function returns a list containing all the residue number of the pdb file
    """ 
    list_res_num = []
    if filter_atoms:
        pdb_ref = filter_atoms_by_ref(pdb, atom_ref)
    else:
        pdb_ref = pdb
    for index, item in enumerate(pdb_ref):
        list_res_num.append(int(pdb_ref[index][22:26].strip(' ')))
    print list_res_num
    return list_res_num
