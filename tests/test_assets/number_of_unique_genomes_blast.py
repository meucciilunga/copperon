CopA_blast_results_path = '/home/katanga/Coding/copperon/tests/test_assets/CopA_blast_result.txt'
CopY_blast_results_path = '/home/katanga/Coding/copperon/tests/test_assets/CopY_blast_result.txt'

with open(CopA_blast_results_path, 'r') as CopA_results:
    accessions = []

    for line in CopA_results.readlines():
        accessions.append(line.split('\t')[0])

    msg = "Number of Unique CopA Blast Result Genomes: {}".format(len(set(accessions)))
    print(msg)

with open(CopY_blast_results_path, 'r') as CopY_results:
    accessions = []

    for line in CopY_results.readlines():
        accessions.append(line.split('\t')[0])

    msg = "Number of Unique CopY Blast Result Genomes: {}".format(len(set(accessions)))
    print(msg)