import pymesh 
import numpy as np
from scipy.spatial import cKDTree
# This script evaluates which pdbs in a list involve a receptor: helix binder. 
# Useful tools : It invokes dssp

from IPython.core.debugger import set_trace
import os
datadir = '/home/gainza/lpdi_fs/masif_design_paper/masif/data/masif_site'
pdbdir = os.path.join(datadir,'data_preparation/01-benchmark_pdbs/')
listfile = os.path.join(datadir,'lists/pdbbind_dimers_w_affinities.txt')

# Parse the pdbbind file. 

pdbbind_lines = open(listfile)
affinity = {}
all_metrics = {'mM':1000000, 'uM':1000, 'nM':1, 'pM':0.001, 'fM':0.000001}
for line in pdbbind_lines.readlines():
    ppi_pair_id = line.split(' ')[0]
    val = line.split(' ')[1].rstrip()
    val = val.replace('~', '=')
    val = val.replace('<', '=')
    val = val.split('=')[1]
    metric = ''.join([val[-2], val[-1]])
    assert(metric in all_metrics)
    val = float(val.replace(metric, ''))
    val = val*all_metrics[metric]
    affinity[ppi_pair_id] = val

# Find helical peptides that bind to the target. 
# Parse the pdb files and find the residues in the interface. 
from Bio.PDB import * 
from Bio.PDB.DSSP import DSSP
from IPython.core.debugger import set_trace

names = []
affinities = []
parser = PDBParser()
#for ppi_pair_id in affinity: 
for ppi_pair_id in affinity: 
    #if '5szk' not in ppi_pair_id: 
    #    continue
    pdbid = ppi_pair_id.split('_')[0]
    chain1 = ppi_pair_id.split('_')[1]
    chain2 = ppi_pair_id.split('_')[2]
    chainCombin = [(chain1, chain2), (chain2, chain1)]
    #chainCombin = [(chain1, chain2)]
    for chain_pair in chainCombin:

        chainA = chain_pair[0]
        chainB = chain_pair[1]

        try:
            struct_nameA = os.path.join(pdbdir,pdbid+'_'+chainA+'.pdb')
            structA = parser.get_structure(struct_nameA, struct_nameA)
        except: 
            continue

        try:
            struct_nameB = os.path.join(pdbdir,pdbid+'_'+chainB+'.pdb')
            structB = parser.get_structure(struct_nameB, struct_nameB)
        except: 
            #print('Error with {}'.format(struct_nameB))
            continue

        modelA = structA[0]
        modelB = structB[0]

        dsspA = DSSP(modelA, struct_nameA)

        # Get interfaceA: 
        atomsA = Selection.unfold_entities(structA, 'A')
        atomsB = Selection.unfold_entities(structB, 'A')
        resA = Selection.unfold_entities(structA, 'R')

        coordsA = [x.get_coord() for x in atomsA]
        coordsB = [x.get_coord() for x in atomsB]

        ckd = cKDTree(coordsB)
        distsA_to_B, r = ckd.query(coordsA)
        
        # Get the residues in the interface
        interfaceA = np.where(distsA_to_B < 5.0)[0]
        residInterface = [atomsA[x].get_parent().get_id()[1] for x in interfaceA]
        residInterface = sorted(residInterface)

        # Get the secondary structure of these residues. 
        helix_sum = []
#        prevIntRes = -1000
        curHelixSum = 0.0
        for ix, elem in enumerate(dsspA):
            resid = resA[ix].get_id()[1]
            if resid in residInterface: 
                if elem[2] == 'H' or  \
                    elem[2] == 'G' or \
                    elem[2] == 'I' or \
                    elem[2] == 'T': 
                    #if resid - prevIntRes < 6: 
                    curHelixSum += residInterface.count(resid)
                else: 
                    helix_sum.append(curHelixSum)
                    curHelixSum = 0.0
            if elem[2] == '-': 
                helix_sum.append(curHelixSum)
                curHelixSum = 0.0
                    #prevIntRes = resid
        helix_sum.append(curHelixSum)
                
        if len(helix_sum) > 0:
            frac = np.max([x/len(residInterface) for x in helix_sum])
            if frac > 0.75: 
                print('fraction: {} receptor:{} helix:{} affinity:{}'.format(frac, pdbid+'_'+chainB, pdbid+'_'+chainA, affinity[ppi_pair_id]))
        

