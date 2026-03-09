from chembl_webresource_client.new_client import new_client
from chembl_webresource_client.settings import Settings
import pandas as pd
from tqdm import tqdm

Settings.Instance().TIMEOUT = 60
Settings.Instance().CACHING = True
Settings.Instance().CONCURRENT_SIZE = 50

def download_mechanisms(human_only=True, limit=None):
    molecule  = new_client.molecule
    mechanism = new_client.mechanism
    target    = new_client.target

    print("Fetching mechanisms...")
    mechs = list(mechanism.filter(target_organism='Homo sapiens') if human_only else mechanism.filter())
    if limit:
        mechs = mechs[:limit]
    print(f"Found {len(mechs)} mechanisms")

    compound_ids = list({m['molecule_chembl_id'] for m in mechs if m.get('molecule_chembl_id')})
    target_ids   = list({m['target_chembl_id']   for m in mechs if m.get('target_chembl_id')})

    print("Fetching compound data...")
    compound_cache = {}
    for i in tqdm(range(0, len(compound_ids), 100)):
        for c in molecule.filter(molecule_chembl_id__in=compound_ids[i:i+100]):
            compound_cache[c['molecule_chembl_id']] = c

    print("Fetching target data...")
    target_cache = {}
    for i in tqdm(range(0, len(target_ids), 100)):
        for t in target.filter(target_chembl_id__in=target_ids[i:i+100]):
            uniprot_ids = [
                xref['xref_id']
                for comp in t.get('target_components', [])
                for xref in comp.get('target_component_xrefs', [])
                if xref.get('xref_src_db') == 'UniProt'
            ]
            gene_names = [
                comp.get('gene_name', '')
                for comp in t.get('target_components', [])
                if comp.get('gene_name')
            ]
            target_cache[t['target_chembl_id']] = {
                'target_name':        t.get('pref_name', ''),
                'target_type':        t.get('target_type', ''),
                'target_organism':    t.get('organism', ''),
                'target_subunit':     t.get('target_components', [{}])[0].get('component_synonym', '') if t.get('target_components') else '',
                'uniprot_accessions': ';'.join(uniprot_ids),
                'gene_names':         ';'.join(gene_names)
            }

    print("Building records...")
    rows = []
    for m in mechs:
        cid  = m.get('molecule_chembl_id')
        tid  = m.get('target_chembl_id')
        c    = compound_cache.get(cid, {})
        t    = target_cache.get(tid, {})
        refs = m.get('mechanism_refs', [])
        ref  = refs[0] if refs else {}

        rows.append({
            'molecule_chembl_id':   cid,
            'compound_name':        c.get('pref_name', ''),
            'max_phase':            c.get('max_phase', ''),
            'molecule_type':        c.get('molecule_type', ''),
            'parent_drug':          c.get('parent_molecule_chembl_id', ''),
            'first_approval':       c.get('first_approval', ''),
            'indication_class':     c.get('indication_class', ''),
            'mechanism_id':         m.get('mec_id', ''),
            'mechanism_of_action':  m.get('mechanism_of_action', ''),
            'action_type':          m.get('action_type', ''),
            'direct_interaction':   m.get('direct_interaction', ''),
            'molecular_mechanism':  m.get('molecular_mechanism', ''),
            'disease_efficacy':     m.get('disease_efficacy', ''),
            'selectivity_comment':  m.get('selectivity_comment', ''),
            'binding_site_comment': m.get('binding_site_comment', ''),
            'target_chembl_id':     tid,
            'target_name':          t.get('target_name', ''),
            'target_type':          t.get('target_type', ''),
            'target_organism':      t.get('target_organism', m.get('target_organism', '')),
            'target_subunit':       t.get('target_subunit', ''),
            'uniprot_accessions':   t.get('uniprot_accessions', ''),
            'gene_names':           t.get('gene_names', ''),
            'ref_id':               ref.get('ref_id', ''),
            'ref_type':             ref.get('ref_type', ''),
            'ref_url':              ref.get('ref_url', '')
        })

    df = pd.DataFrame(rows).drop_duplicates(
        subset=['molecule_chembl_id', 'target_chembl_id', 'mechanism_of_action']
    ).sort_values('mechanism_of_action')

    df.to_csv('input.txt', index=False, sep='\t')
    print(f"\nSaved input.txt — {len(df)} records")
    return df


# Run
df = download_mechanisms()

# Download
from google.colab import files
files.download('input.txt')
