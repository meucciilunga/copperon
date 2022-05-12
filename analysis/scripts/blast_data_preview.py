path = "/home/katanga/Coding/copperon/blast/results/CopY"


with open(path, 'r') as f:

    CopY_associated = set()

    for line in f.readlines():
        entry = line.rstrip('\n').split('\t')

        print(entry[0])

        CopY_associated.add(entry[0])



    print(len(CopY_associated))

        
