import json
import os
from pymongo import MongoClient


def get_molecules_collection(db_dir):
    db_path = os.path.join(db_dir, 'tasks_db.json')
    with open(db_path) as f:
        db_creds = json.load(f)
    conn = MongoClient(db_creds['host'], db_creds['port'])
    db = conn[db_creds['database']]
    if db_creds['admin_user']:
        db.authenticate(db_creds['admin_user'], db_creds['admin_password'])
    ipea_coll = db['molecules']
    return ipea_coll

def transform_molecule_doc(mol1):
    mol2 = dict()
    mol2["inchi_root"] = mol1["inchi_root"]
    mol2["can"] = mol1["can"]
    mol2["smiles"] = mol1["smiles"]
    mol2["elements"] = mol1["elements"]
    mol2["nelements"] = mol1["nelements"]
    mol2["user_tags"] = mol1["user_tags"]
    mol2["run_tags"] = mol1["run_tags"]
    mol2["reduced_cell_formula_abc"] = mol1["reduced_cell_formula_abc"]
    mol2["implicit_solvent"] = mol1["implicit_solvent"]
    mol2["pretty_formula"] = mol1["pretty_formula"]
    mol2["pointgroup"] = mol1["pointgroup"]

    mol2["task_id"] = mol1["neutral"]["task_id"]
    mol2["task_id_deprecated"] = mol1["neutral"]["task_id_deprecated"]
    mol2["snlgroup_id_final"] = mol1["neutral"]["snlgroup_id_final"]
    mol2["charge"] = mol1["neutral"]["charge"]
    mol2["spin_multiplicity"] = mol1["neutral"]["spin_multiplicity"]
    mol2["snl_final"] = mol1["neutral"]["snl_final"]
    mol2["molecule"] = mol1["neutral"]["molecule"]
    mol2["xyz"] = mol1["neutral"]["xyz"]
    mol2["inchi"] = mol1["neutral"]["inchi"]

    mol2["IE"] = mol1["IP"]["sol"]
    mol2["EA"] = mol1["EA"]["sol"]

    mol2['solvation_energy'] = mol1["solvation_energy"]
    mol2["svg"] = mol1["svg"]
    return mol2

def copy_collections():
    source_cred_dir = os.environ['DB_LOC_SRC']
    dest_cred_dir = os.environ['DB_LOC_DEST']
    coll_src = get_molecules_collection(source_cred_dir)
    coll_dest = get_molecules_collection(dest_cred_dir)
    molecules = coll_src.find()
    for mol_db in molecules:
        mol_web = transform_molecule_doc(mol_db)
        coll_dest.insert(mol_web)
