import os
print("Fasta file")
fasta_file = input()
name = fasta_file.split(".")
name = name[0]

def create_multimer_file(fasta_file) :
        lenght = 0
        with open(fasta_file, "r") as handle :
            for line in handle :
               line = line.strip().split()
	       if len(line) > 0 :
                  if line[0][0] != ">" :
                     lenght += len(line)
        num_multimer = 2000//lenght
        if num_multimer > 15 :
           num_multimer = 15
        filename = "multimer_"+name+".txt"
        with open(filename,"w") as fh :
           for nbr_line in range(int(num_multimer)) :
              sentence = name+","+str(nbr_line+1)+"\n"
              fh.write(sentence)

if __name__ == "__main__" :
	create_multimer_file(fasta_file)


