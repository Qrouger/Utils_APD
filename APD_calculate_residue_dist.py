from Bio import PDB
import numpy as np
from tabulate import tabulate
import pandas as pd
import csv

print("PDB file:")
pdb_file=input()
print("Protein 1")
prot_1=input()
print("Protein 2")
prot_2=input()
def calculate_residue_distances(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    result = [[prot_1,prot_2,"Distance Ä Cb"]]
    i = 0
    for model in structure:
       list_chain = model.get_list()
       while i<len(list_chain) :
          for chain1 in list_chain : #B and C
             i += 1
             for residue1 in chain1 : #number of residue (len of the chain)
                for atom1 in residue1 : #type of the atom
                   for index in range(i,len(list_chain)) :
                      chain2 = list_chain[index]
                      for residue2 in chain2 :
                         for atom2 in residue2 :
                            if atom1.get_id()=="CB" and atom2.get_id()=="CB" :
                               distance = atom1 - atom2
                               if distance <= 5 :
                                  result.append([residue1.get_resname()+str(residue1.get_id()[1]),residue2.get_resname()+str(residue2.get_id()[1]),str(distance)])
                               else :
                                  pass
    #files = os.listdir(directory)
    #for file in files :
    #   if file == result_*
    fileout = prot_1+"_and_"+prot_2+"_res_int.csv"
    np_table = np.array(result)
    with open(fileout, "w", newline="") as file :
        mywriter = csv.writer(file, delimiter=",")
        mywriter.writerows(np_table)
    print("Write table")

if __name__ =="__main__" :

	calculate_residue_distances(pdb_file)
