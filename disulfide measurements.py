#!/usr/bin/env python
# coding: utf-8

# This is a test to get the S-S distance and torsion angle for residues 215 and 348. 
# Bad style in programming in that there are no error traps, if for example res1 would have a number 
# pointing to a non-Cys residue which has no SG atom

# In[30]:


from Bio import PDB
from Bio.PDB.vectors import calc_dihedral
from numpy import pi

parser = PDB.PDBParser()

pdb1 ='./AF-K0RA12-F1-model_v4.pdb' 
structure = parser.get_structure("KORA12", pdb1) 
model = structure[0] 
chain = model['A'] 
res1 = chain[215] 
res2 = chain[348]
atom1 = res1['CB'] 
atom2 = res1['SG'] 
atom3 = res2['SG'] 
atom4 = res2['CB'] 

distance = atom2-atom3 # vector subtraction
torsion = calc_dihedral(atom1.get_vector(), 
                        atom2.get_vector(), 
                        atom3.get_vector(), 
                        atom4.get_vector())
#for some reason the torsion calculation did not work without .get_vector()

torsion = torsion * 180 / pi  # in degrees

#print(distance)
#print(torsion)
print("Distance: %.2f Å" % distance) #pdb coordinates are in Ångstrom units
print("Torsion angle: %.2f degrees" % torsion)


# In[ ]:




