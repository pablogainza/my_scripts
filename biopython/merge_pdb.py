from Bio.PDB import *
import sys
import os
# Read two PDBs with one chain each, set the first one to have chain = 'A' and the second one chain = 'B'
# Output a file with the two. 
parser  = PDBParser()
struct1_fn = sys.argv[1]
struct2_fn = sys.argv[2]
# 
struct1_bn = os.path.basename(struct1_fn)
struct2_bn = os.path.basename(struct2_fn)
site = struct2_fn.split('/')[-3]

outfilename = "merged_out/{}_{}_{}.pdb".format(struct1_bn.split('.')[0], site, struct2_bn.split('.')[0] )

struct1 = parser.get_structure(struct1_fn, struct1_fn)
struct2 = parser.get_structure(struct2_fn, struct2_fn)

# Select residues to extract and build new structure
model1 = Selection.unfold_entities(struct1, "M")[0]
model2 = Selection.unfold_entities(struct2, "M")[0]

structBuild = StructureBuilder.StructureBuilder()
structBuild.init_structure("output")
structBuild.init_seg(" ")
structBuild.init_model(0)
structBuild.init_chain('A')
structBuild.init_chain('B')
outputStruct = structBuild.get_structure()

for chain1 in model1: 
    for residue in chain1: 
        outputStruct[0]['A'].add(residue)

for chain2 in model2:
    for residue in chain2:
        outputStruct[0]['B'].add(residue)
        
# Output the selected residues
pdbio = PDBIO()
pdbio.set_structure(outputStruct)
pdbio.save(outfilename)
print(outfilename)
