import itertools

# https://stackoverflow.com/questions/36429507/python-combinations-without-repetitions
t = ["A","B","C","D","E","F","G","H","I","J","K","L"]
c = list(itertools.combinations(t, 3))

print(c)
print(len(c))

pdb = "1G63"
file_lst = []

print(type(c[1]))

for comb in c:
    file_lst.append(pdb+"_"+"".join(comb)+".pdb")

print(file_lst)
print(len(file_lst))