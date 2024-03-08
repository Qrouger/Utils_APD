import copy

all_lines = str()
print("Enter a fasta file")
file = input()
dict_prot = dict()
save_prot =  None
ban_list =  list()
list_int = []
lenght = 0
y_n = "n"
i1 = 0
i2 = 0
with open(file, 'rb') as handle :
    for lines in handle :
        line = lines.decode('UTF-8')
        new_lines = line.strip().split()
        if new_lines != [] :
           if new_lines[0][0] == ">" :
              dict_prot[save_prot] = lenght
              save_prot = new_lines[0][1:len(new_lines[0])]
              lenght = 0
           else :
             lenght += len(new_lines[0])
    del(dict_prot[None])
    for proteins in dict_prot.keys() :
        i1 += 1
        for proteins2 in dict_prot.keys() :
            i2 += 1
            if proteins != proteins2 and i1 < i2 :
                adition = dict_prot[proteins]+ dict_prot[proteins2]
                if adition <= 2000 :
                   list_int.append([proteins,proteins2])
                else :
                   ban_list.append([proteins,proteins2])
    for interactions in list_int :
       write_l = str(interactions[0])+";"+str(interactions[1])+"\n"
       all_lines = all_lines + write_l
    with open("file_3090_ti_int.txt", "w") as fh :
       fh.write(all_lines)
print("list banned",sorted(ban_list))

