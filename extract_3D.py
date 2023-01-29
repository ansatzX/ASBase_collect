import h5py
import os
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem


def conf_gen(smi, conf_number, smi_charge):
    m = Chem.MolFromSmiles(smi)
    m = Chem.AddHs(m)
    print(smi)
    params = AllChem.ETKDGv3()
    cids = AllChem.EmbedMultipleConfs(m, numConfs=conf_number, params=params)
    res = AllChem.MMFFOptimizeMoleculeConfs(m,
                                            numThreads=40,
                                            mmffVariant="MMFF94")
    energy_low = min(res)
    label = res.index(energy_low)

    conf_mol = Chem.MolToXYZBlock(m, confId=label)
    del cids
    return conf_mol


db = h5py.File("testdb.hdf5", "r")

IDs = list(db.keys())

for id in IDs:
    group = db[id]
    pdb_block = ""
    xyz_block = ""

    if "XTB optimized structure" in group.keys():
        pdb_block = group["XTB optimized structure"][0].decode()
    else:
        if "MMFF94 optimized structure" in group.keys():
            pdb_block = group["MMFF94 optimized structure"][0].decode()
        else:
            smi = group["SMILES"][0].decode()
            xyz_block = conf_gen(smi, 500, 0)

    if xyz_block != "":
        with open("xyzs/" + id + ".xyz", "w") as f:
            f.write(xyz_block)
    else:
        with open("pdbs/" + id + ".pdb", "w") as f:
            f.write(pdb_block)
    del group
    
db.close()
