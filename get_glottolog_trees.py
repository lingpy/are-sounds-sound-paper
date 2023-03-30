#%%
from ete3 import Tree
import re
import pandas as pd
import numpy as np
import requests
import os
from tqdm import tqdm



#%%



url = 'https://cdstar.eva.mpg.de//bitstreams/EAEA0-56EE-48A3-6B18-0/tree_glottolog_newick.txt'
glotName = 'tree_glottolog_newick.txt'

r = requests.get(url)

with open(glotName, 'wb') as f:
   f.write(r.content)



#%%

with open(glotName) as f:
    raw = f.readlines()

os.remove(glotName)

trees = []

for i, ln in enumerate(raw):
    ln = ln.strip()
    ln = re.sub(r"\'[A-Z][^[]*\[", "[", ln)
    ln = re.sub(r"\][^']*\'", "]", ln)
    ln = re.sub(r"\[|\]", "", ln)
    ln = ln.replace(":1", "")
    trees.append(Tree(ln, format=1))

#%%


glot = Tree()
for t in trees:
    glot.add_child(t)


nonLeaves = [nd.name for nd in glot.traverse()
             if nd.name != '' and not nd.is_leaf()]

for nm in tqdm(nonLeaves):
    nd = glot & nm
    nd.name = ''
    nd.add_child(name=nm)

taxa = glot.get_leaf_names()

#%%

dbs = os.listdir("languages")

#%%
tree_pth = "glottolog_trees"
isExist = os.path.exists(tree_pth)

if not isExist:
    os.mkdir(tree_pth)
#%%

for db in dbs:
    d = pd.read_table(os.path.join('languages',db))
    glot_db = glot.copy()
    glot_db.prune([l for l in d.Glottocode if l in taxa])
    for l in glot_db.get_leaves():
        doculects = d.Language[d.Glottocode == l.name]
        for dc in doculects:
            l.add_child(name=dc)
    glot_db.write(
        outfile=os.path.join(tree_pth, db.split(".")[0]+"_glottolog.tre"), 
        format=9
        )
