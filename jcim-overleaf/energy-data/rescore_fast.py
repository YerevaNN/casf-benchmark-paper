"""Equivalent batched reference loading for rescore.py.
Run: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python docs/energy_analysis/rescore_fast.py
The scientific protocol, identity checks, and cached calculation inputs are unchanged.
"""
from pathlib import Path
import hashlib,json
import numpy as np
import zarr
import rescore
from casf_benchmark.chembl3d import loader
from casf_benchmark.chembl3d.identity import expected_identity_from_mapping
from casf_benchmark.cli import analyze_conformer_sets as analysis

def batched_references(row,mol_id,topology_root,zarr_root):
    group=str(row['chembl3d_group']).zfill(3)
    expected,record_index=expected_identity_from_mapping(row)
    template=loader.load_topology_mol(group,str(row['chembl3d_mol_id']),topology_root,
        expected_smiles=expected,sdf_record_index=record_index)
    assert template is not None
    base=zarr_root/group
    ids=zarr.open_array(str(base/'mol_id'),mode='r')
    indices=loader.find_mol_id_indices(ids,str(row['chembl3d_mol_id']))
    if not indices:return [],None
    # Each coordinate chunk contains 500,000 conformers. Batch reads avoid
    # repeatedly decompressing the same chunk for every candidate stereoisomer.
    coords=zarr.open_array(str(base/'coord'),mode='r').oindex[indices,:,:]
    numbers=zarr.open_array(str(base/'numbers'),mode='r').oindex[indices,:]
    expected_numbers=loader._topology_atomic_numbers(template)
    mols=[]
    for index,xyz,atoms in zip(indices,coords,numbers):
        loader._validate_atomic_numbers(atoms,expected_numbers,index,str(row['chembl3d_mol_id']))
        mol=loader._set_coords_from_row(template,xyz)
        if loader.mol_matches_expected_smiles(mol,expected,from_3d=True):mols.append(mol)
    return analysis.sample_reference_chembl_mols(mols,mol_id)

if __name__=='__main__':
    analysis.get_reference_chembl_mols=batched_references
    rescore.main()
    path=Path(__file__).parent/'provenance.json'
    data=json.loads(path.read_text())
    data['reference_loading_optimization']={
        'source':'rescore_fast.py',
        'sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'description':'Batched Zarr reads; same template, candidate order, atomic-number and stereoisomer checks. Existing validated cache entries reused without changing the energy protocol.'}
    path.write_text(json.dumps(data,indent=2)+'\n')
