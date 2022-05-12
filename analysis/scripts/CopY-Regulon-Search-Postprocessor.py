tax_cat = "/home/katanga/Coding/copperon/database/data_assets/categories.dmp"
tax_node = "/home/katanga/Coding/copperon/database/data_assets/nodes.dmp"
tax_name = "/home/katanga/Coding/copperon/database/data_assets/names.dmp"
tax_cat_file = open(tax_cat, 'r')
tax_node_file = open(tax_node, 'r')
tax_name_file = open(tax_name, 'r')

# Process data from taxonomic names file and map
# all known tax_ids to a name; more info here: 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3245000/
def scientific_names_table_builder(name_file):
    scientifc_names_table = dict()

    for line in name_file.readlines():
        entry = line[:-3].split('\t|\t') 
        # remove last 2 chars (ruin parsing pattern) then 
        # split entries into an array

        entry_tax_id = int(entry[0])
        name = entry[1]
        name_type = entry[3].lower()

        # Only map a taxid to its scientific name
        if name_type == 'scientific name':
            scientifc_names_table[entry_tax_id] = name

    return(scientifc_names_table)

# Process data from taxonomic node file and create
# a function (in the form of a hash table) that can
# map any given node to its first ancestor node
def lineage_map_builder(node_file):
    lineage_table = dict()
    rank_table = dict()

    for line in node_file.readlines():
        entry = line[:-3].split('\t|\t')
        entry_tax_id = int(entry[0])
        parent_tax_id = int(entry[1])
        rank = entry[2]

        # Pair rank and current tax_id
        rank_table[entry_tax_id] = rank

        # Pair child with parent
        lineage_table[entry_tax_id] = parent_tax_id
    
    return(lineage_table, rank_table)


# Walk up a branch to the root of lineage table
# (due to recursive structure, can't include
# the original tax_id in the node list it returns)
def branch_walk(tax_id, lineage_table):
    parent_id = lineage_table.get(tax_id)

    if parent_id != 1:
        node_vec = branch_walk(parent_id, lineage_table)
        node_vec.append(parent_id)
        return(node_vec)
    else:
        return([parent_id]) # nucleation point

# Given a tax_id, builds its root branch via 
# recursive calls to the lineage table
def build_branch_nodes(tax_id, lineage_table):
    output = branch_walk(tax_id, lineage_table)
    output.append(tax_id)

    return(output)


# Given a tax_id, scientific names table, and a lineage 
# table, find the tax_id's branch (nodes) and build a its 
# full name (superkingdom always starts at third node)
def taxon_name_builder(tax_id, scientific_names, lineage_table):
    nodes = build_branch_nodes(tax_id, lineage_table)
    full_name = []

    for item in nodes:
        sub_name = scientific_names[item]
        full_name.append(sub_name)

    name = ' / '
    print(name.join(full_name[2:]))

# Process data from taxonomic category file 
# (tax_cat file lists only species or lower!)
def build_species_finder_table(tax_cat_file):
    tax_cat_entries = dict()

    for line in tax_cat_file.readlines():
        entry = line.split('\t')
        entry_tax_id = entry[2]
        species_level_id = entry[1]
        
        # will allow us to see how many unique **species** are in the prokaryote genome database
        # by creating a map that links strains/subspecies to the major species they fall under
        tax_cat_entries[entry_tax_id] = species_level_id
    
    return(tax_cat_entries)

names = scientific_names_table_builder(tax_name_file)
lineage = lineage_map_builder(tax_node_file)

# Process a list of tax_ids
path = "/home/katanga/Coding/copperon/analysis/output/tax_ids.log"
with open(path, 'r') as f:
    for line in f.readlines():
        tax_id = int(line.rstrip())
        taxon_name_builder(tax_id, names, lineage[0])