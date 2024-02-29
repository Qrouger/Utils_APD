print("Name of protein")
name = input()
print ("Numbers of homo-oligomers")
n = input()
def create_multimer_file(name,n) :
	filename = "multimer_"+name+".txt"
	with open(filename,"w") as fh :
		for nbr_line in range(int(n)) :
			sentence = name+","+str(nbr_line+1)+"\n"
			fh.write(sentence)

if __name__ == "__main__" :
	create_multimer_file(name,n)


