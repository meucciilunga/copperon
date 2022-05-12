from cmath import inf
from numpy import float64


path = "/home/katanga/Coding/copperon/analysis/output/genomes.txt"
output = "/home/katanga/Coding/copperon/analysis/output/tax_ids.log"

def average(list):
    return sum(list) / len(list)

with open(path, 'r') as f:

    CopY_associated_density = []
    Non_CopY_associated_density = []
    CopY_species = set()
    all_species = set()

    for line in f.readlines():
        entry = line.rstrip('\n').split('\t')

        operator_density = float64(entry[11])
        if entry[6] == '0':
            continue

        if '1313' in entry[6]

        if 'Y' in entry[-1]:
            CopY_associated_density.append(operator_density)
            CopY_species.add(entry[1])

        if not 'Y' in entry[-1]:
            Non_CopY_associated_density.append(operator_density)
            if average(Non_CopY_associated_density) == inf:
                print(entry)
            
        

with open(output, 'w') as f:
    for item in CopY_species:
        out = "{}\n".format(item)
        f.write(out)

print(average(Non_CopY_associated_density))
print(average(CopY_associated_density))
print(len(CopY_species))