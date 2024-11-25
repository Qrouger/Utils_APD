from Bio import PDB
import numpy as np
from tabulate import tabulate
import pandas as pd
import csv
import copy

print("It's a cryo Model ? y/n")
pdb_mod=input()
print("PDB file:")
pdb_file=input()

def calculate_residue_distances(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    dict_int = dict()
    dict_information = dict()
    atom_possible_contact = ["O","OH","NH2","NH1","OG","NE2","ND2","NZ","NE","N","OE1","OE2","OD2","OG1"]

    for model in structure:
       list_chain = model.get_list()
       if pdb_mod == "y" :
          cryo = 2
       else :
          cryo = 1
       print(str(len(list_chain)//cryo)+" chains detected")
       for i in range(0,len(list_chain)//cryo) :
          chain1 = list_chain[i] #B and C
          for residue1 in chain1 : #number of residue (len of the chain)
             for atom1 in residue1 : #type of the atom
                for index in range(i+1,len(list_chain)//cryo) :
                   chain2 = list_chain[index]
                   for residue2 in chain2 :
                      for atom2 in residue2 :
                         #if atom1.get_id()=="CB" and atom2.get_id()=="CB" :
                            distance = atom1 - atom2
                            if distance <= 3.64 :
                               if chain1.get_id()+chain2.get_id() in dict_int.keys() : #to make different table for different interaction  
                                  dict_int[chain1.get_id()+chain2.get_id()].append([chain1.get_id()+":"+residue1.get_resname()+"  "+str(residue1.get_id()[1]),chain2.get_id()+":"+residue2.get_resname()+"   "+str(residue2.get_id()[1]),str(distance)])
                                  dict_information[chain1.get_id()+chain2.get_id()].append([residue1.get_resname()+str(residue1.get_id()[1]),atom1.get_id(),residue2.get_resname()+str(residue2.get_id()[1]),atom2.get_id(),str(distance)])
                               else :
                                  dict_int[chain1.get_id()+chain2.get_id()] = [["Chain "+chain1.get_id(),"Chain "+chain2.get_id(),"Distance Ä"]]
                                  dict_information[chain1.get_id()+chain2.get_id()] = [[" "," ","Chain "+chain1.get_id(),"Chain "+chain2.get_id(),"Distance Ä"]]
                                  dict_int[chain1.get_id()+chain2.get_id()].append([chain1.get_id()+":"+residue1.get_resname()+"   "+str(residue1.get_id()[1]),chain2.get_id()+":"+residue2.get_resname()+"   "+str(residue2.get_id()[1]),str(distance)])
                                  dict_information[chain1.get_id()+chain2.get_id()].append([residue1.get_resname()+str(residue1.get_id()[1]),atom1.get_id(),residue2.get_resname()+str(residue2.get_id()[1]),atom2.get_id(),str(distance)])
                            else :
                               pass

    int_list = list()
    save_dict = copy.deepcopy(dict_int)
    for chains in dict_int.keys() :
      for line in range(len(save_dict[chains])) :
         if dict_information[chains][line][1] not in atom_possible_contact or dict_information[chains][line][3] not in atom_possible_contact and dict_information[chains][line][0] != " " : #sort in function of atom id
             for index in range(len(dict_int[chains])) :
                  if dict_int[chains][index-1] == save_dict[chains][line] :
                     dict_int[chains].pop(index-1)
    save_dict2 = copy.deepcopy(dict_int)
    for chains2 in save_dict2.keys() : #sort in function of the distance
      for line2 in range(len(save_dict2[chains2])) :
         int_list.append(save_dict2[chains2][line2])
    for chains2 in save_dict2.keys() : #sort in function of the distance
      for line2 in range(len(save_dict2[chains2])) :
         for interaction in range(len(int_list)) :
            if int_list[interaction][0] == save_dict2[chains2][line2][0] and int_list[interaction][1] == save_dict2[chains2][line2][1] and float(save_dict2[chains2][line2][2]) > float(int_list[interaction][2]) :
               for index in range(len(dict_int[chains2])) :
                  if dict_int[chains2][index-1] == save_dict2[chains2][line2] :
                     dict_int[chains2].pop(index-1)
            elif int_list[interaction][0] == save_dict2[chains2][line2][0] and int_list[interaction][1] == save_dict2[chains2][line2][1] and float(save_dict2[chains2][line2][2]) < float(int_list[interaction][2]) :
               for index in range(1,len(dict_int[chains2])) :
                  if dict_int[chains2][index-1] == int_list[interaction] :
                     dict_int[chains2].pop(index-1)
            else :
                pass
      fileout = chains2+"_res_int.csv"
      np_table = np.array(dict_int[chains2])
      with open(fileout, "w", newline="") as file :
         mywriter = csv.writer(file, delimiter=",")
         mywriter.writerows(np_table)
      print("Write table")

if __name__ =="__main__" :

	calculate_residue_distances(pdb_file)
	
### Last release
parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('protein', path_int + "/ranked_0.pdb")
        dict_int = dict()
        int_already_know = dict()
        proteins = path_int.split('/')[2].split('_and_')
        for model in structure:
            list_chain = model.get_list()
            for index1 in range(0,len(list_chain)) :
                chain1 = list_chain[index1] #B and C
                for residue1 in chain1 : #number of residue (len of the chain)
                    for atom1 in residue1 : #type of the atom
                        if atom1.get_id() == "CA" :
                            for index2 in range(index1+1,len(list_chain)) :
                                chain2 = list_chain[index2]
                                for residue2 in chain2 :
                                    for atom2 in residue2 :
                                        if atom2.get_id() == "CA" :
                                            distance = atom1 - atom2
                                            if distance <= 6 : #filtered on pLDDT and distance, be stringent to avoid false residue interaction (or maybe use PAE ?)
                                                res_int = chain1.get_id()+":"+residue1.get_resname()+" "+str(residue1.get_id()[1])," "+chain2.get_id()+":"+residue2.get_resname()+" "+str(residue2.get_id()[1])
                                                print(res_int)
                                                if chain1.get_id()+chain2.get_id() in dict_int.keys() : #to make different table for different interaction
                                                    if res_int in int_already_know.keys() and int_already_know[res_int] > str(distance) :
                                                        dict_int[chain1.get_id()+chain2.get_id()].remove([res_int[0],res_int[1]," "+str(int_already_know[res_int])])
                                                        dict_int[chain1.get_id()+chain2.get_id()].append([res_int[0],res_int[1]," "+str(distance)])
                                                        int_already_know[res_int] = str(distance)
                                                    elif res_int in int_already_know.keys() and int_already_know[res_int] < str(distance) : #skip double interaction with differents atoms
                                                        pass
                                                    else :
                                                        dict_int[chain1.get_id()+chain2.get_id()].append([res_int[0],res_int[1]," "+str(distance)])
                                                        int_already_know[res_int] = str(distance)
                                                else :
                                                    dict_int[chain1.get_id()+chain2.get_id()] = [[proteins[0]," "+proteins[1]," Distance Ä"]]
                                                    dict_int[chain1.get_id()+chain2.get_id()].append([res_int[0],res_int[1]," "+str(distance)])
                                                    int_already_know[res_int] = str(distance)
                                            else :
                                                pass
        for chains in dict_int.keys() :
            file.define_interface(dict_int[chains],proteins)
            fileout = chains+"_res_int.csv"
            np_table = np.array(dict_int[chains])
            with open(f"{path_int}/"+fileout, "w", newline="") as file :
                mywriter = csv.writer(file, delimiter=",")
                mywriter.writerows(np_table)
            #print("Write table")
