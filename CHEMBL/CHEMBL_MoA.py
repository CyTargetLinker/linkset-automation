from chembl_webresource_client.new_client import new_client
from chembl_webresource_client.settings import Settings
import pandas as pd
from datetime import datetime

Settings.Instance().TIMEOUT = 60
Settings.Instance().CACHING = True
Settings.Instance().CONCURRENT_SIZE = 20


def get_all_mechanisms_count():
    mechanism = new_client.mechanism
    try:
        all_mechanisms_list = list(mechanism.filter())
        print(f"Total mechanisms in ChEMBL: {len(all_mechanisms_list)}")
        return len(all_mechanisms_list)
    except Exception as e:
        print(f"Error getting mechanisms: {e}")
        return 0


def download_all_mechanism_data(limit_compounds=None, human_only=False, tqdm_class=None):

    print(f"Start time: {datetime.now()}")
    molecule  = new_client.molecule
    mechanism = new_client.mechanism
    target    = new_client.target

    # ── 1. Get mechanisms ──────────────────────────────────────────────
    print("1. Getting mechanisms...")
    try:
        if human_only:
            all_mechanisms_list = list(mechanism.filter(target_organism='Homo sapiens'))
        else:
            all_mechanisms_list = list(mechanism.filter())
        print(f"Found {len(all_mechanisms_list)} mechanisms")
    except Exception as e:
        print(f"Error: {e}")
        return None

    if not all_mechanisms_list:
        return None

    # ── 2. Unique compound IDs ─────────────────────────────────────────
    print("2. Extracting unique compounds...")
    compound_ids = list(set([
        m.get('molecule_chembl_id') for m in all_mechanisms_list
        if m.get('molecule_chembl_id')
    ]))
    print(f"Found {len(compound_ids)} unique compounds")

    if limit_compounds:
        compound_ids = compound_ids[:limit_compounds]
        print(f"Limited to {len(compound_ids)} compounds")

    # ── 3. Fetch compound data ─────────────────────────────────────────
    # NOTE: CAS number is NOT available in ChEMBL molecule API.
    # cross_references only contains: DailyMed, PubChem SIDs, Wikipedia, TG-GATEs.
    # We collect PubChem SID as the best available external compound identifier.
    print("3. Fetching compound data...")
    compound_cache   = {}
    target_cache     = {}
    failed_compounds = 0

    with tqdm_class(total=len(compound_ids), desc="Fetching compounds") as pbar:
        for i in range(0, len(compound_ids), 50):
            for comp_id in compound_ids[i:i + 50]:
                try:
                    comp_data = molecule.get(comp_id)
                    if comp_data:
                        compound_cache[comp_id] = comp_data
                    else:
                        failed_compounds += 1
                except Exception:
                    failed_compounds += 1
                pbar.update(1)

    print(f"Fetched {len(compound_cache)} compounds ({failed_compounds} failed)")

    # ── 4. Process mechanisms ──────────────────────────────────────────
    print("4. Processing mechanism records...")

    mechanisms_by_compound = {}
    for mech in all_mechanisms_list:
        comp_id = mech.get('molecule_chembl_id')
        if comp_id and comp_id in compound_cache:
            mechanisms_by_compound.setdefault(comp_id, []).append(mech)

    all_data = []

    with tqdm_class(total=len(mechanisms_by_compound), desc="Processing mechanisms") as pbar:
        for compound_id, compound_mechanisms in mechanisms_by_compound.items():
            comp_data = compound_cache.get(compound_id, {})

            # ── PubChem SID (best available external ID in ChEMBL) ──
            # CAS is NOT in ChEMBL. PubChem CID can be used to look up CAS externally.
            pubchem_sid = ''
            for xref in comp_data.get('cross_references', []) or []:
                if xref.get('xref_src') == 'PubChem':
                    pubchem_sid = xref.get('xref_id', '')
                    break  # take first PubChem SID

            # ── InChIKey (stable universal structure identifier) ──
            inchi_key = ''
            mol_struct = comp_data.get('molecule_structures') or {}
            if mol_struct:
                inchi_key = mol_struct.get('standard_inchi_key', '') or ''

            for mechanism_data in compound_mechanisms:
                try:
                    target_id   = mechanism_data.get('target_chembl_id')
                    target_info = {}

                    if target_id:
                        if target_id in target_cache:
                            target_info = target_cache[target_id]
                        else:
                            try:
                                target_data = target.get(target_id)
                                if target_data:
                                    uniprot_ids  = []
                                    gene_symbols = []  # from target_component_synonyms GENE_SYMBOL
                                    hgnc_ids     = []  # from target_component_xrefs HGNC

                                    for component in target_data.get('target_components', []) or []:

                                        # ── Gene symbol: confirmed key is target_component_synonyms ──
                                        for syn in component.get('target_component_synonyms', []) or []:
                                            if syn.get('syn_type') == 'GENE_SYMBOL':
                                                gene_symbols.append(syn.get('component_synonym', ''))
                                                break  # one gene symbol per component

                                        # ── xrefs: UniProt and HGNC ──
                                        for xref in component.get('target_component_xrefs', []) or []:
                                            src   = xref.get('xref_src_db', '')
                                            xid   = xref.get('xref_id', '')
                                            xname = xref.get('xref_name', '')
                                            if src == 'UniProt':
                                                uniprot_ids.append(xid)
                                            elif src == 'HGNC':
                                                hgnc_ids.append(xid)
                                                # xref_name is the gene symbol e.g. "CA1"
                                                # use as fallback if target_component_synonyms was empty
                                                if xname and xname not in gene_symbols:
                                                    gene_symbols.append(xname)

                                    target_info = {
                                        'target_chembl_id'  : target_id,
                                        'target_name'       : target_data.get('pref_name', ''),
                                        'target_type'       : target_data.get('target_type', ''),
                                        'target_organism'   : target_data.get('organism', ''),
                                        'uniprot_accessions': ';'.join(filter(None, uniprot_ids)),
                                        'gene_symbol'       : ';'.join(filter(None, gene_symbols)),
                                        'hgnc_id'           : ';'.join(filter(None, hgnc_ids)),
                                    }
                                    target_cache[target_id] = target_info

                            except Exception:
                                target_info = {
                                    'target_chembl_id'  : target_id,
                                    'target_name'       : '',
                                    'target_type'       : '',
                                    'target_organism'   : mechanism_data.get('target_organism', ''),
                                    'uniprot_accessions': '',
                                    'gene_symbol'       : '',
                                    'hgnc_id'           : '',
                                }
                                target_cache[target_id] = target_info

                    mechanism_refs = mechanism_data.get('mechanism_refs', []) or []

                    record = {
                        # ── Compound ──────────────────────────────────────────────
                        'molecule_chembl_id' : compound_id,
                        'compound_name'      : comp_data.get('pref_name', ''),
                        # CAS not available in ChEMBL; use InChIKey or PubChem SID
                        # to look up CAS via PubChem API if needed
                        'inchi_key'          : inchi_key,
                        'pubchem_sid'        : pubchem_sid,
                        'max_phase'          : comp_data.get('max_phase', ''),
                        'molecule_type'      : comp_data.get('molecule_type', ''),
                        'parent_drug'        : comp_data.get('parent_molecule_chembl_id', ''),
                        'first_approval'     : comp_data.get('first_approval', ''),
                        'indication_class'   : comp_data.get('indication_class', ''),
                        # ── Mechanism ─────────────────────────────────────────────
                        'mechanism_id'       : mechanism_data.get('mec_id', ''),
                        'mechanism_of_action': mechanism_data.get('mechanism_of_action', ''),
                        'action_type'        : mechanism_data.get('action_type', ''),
                        'direct_interaction' : mechanism_data.get('direct_interaction', ''),
                        'molecular_mechanism': mechanism_data.get('molecular_mechanism', ''),
                        'disease_efficacy'   : mechanism_data.get('disease_efficacy', ''),
                        'selectivity_comment': mechanism_data.get('selectivity_comment', ''),
                        'binding_site_comment': mechanism_data.get('binding_site_comment', ''),
                        # ── Target ────────────────────────────────────────────────
                        'target_chembl_id'   : target_info.get('target_chembl_id', ''),
                        'target_name'        : target_info.get('target_name', ''),
                        'target_type'        : target_info.get('target_type', ''),
                        'target_organism'    : target_info.get('target_organism', mechanism_data.get('target_organism', '')),
                        'uniprot_accessions' : target_info.get('uniprot_accessions', ''),
                        'gene_symbol'        : target_info.get('gene_symbol', ''),
                        'hgnc_id'            : target_info.get('hgnc_id', ''),
                        # ── References ────────────────────────────────────────────
                        'ref_id'             : mechanism_refs[0].get('ref_id', '')   if mechanism_refs else '',
                        'ref_type'           : mechanism_refs[0].get('ref_type', '') if mechanism_refs else '',
                        'ref_url'            : mechanism_refs[0].get('ref_url', '')  if mechanism_refs else '',
                    }
                    all_data.append(record)

                except Exception:
                    continue

            pbar.update(1)

    # ── 5. Save ────────────────────────────────────────────────────────
    print("5. Saving results...")
    if not all_data:
        print("No data collected!")
        return None

    df = pd.DataFrame(all_data)
    df_unique = df.drop_duplicates(subset=['mechanism_id'])

    filename = 'input.txt'
    df_unique.to_csv(filename, index=False, sep='\t')

    print(f"\nResults Summary:")
    print(f"File              : {filename}")
    print(f"Mechanism records : {len(df_unique)}")
    print(f"Unique compounds  : {df_unique['molecule_chembl_id'].nunique()}")
    print(f"Unique targets    : {df_unique['target_chembl_id'].nunique()}")

    sample_cols = ['molecule_chembl_id', 'compound_name', 'inchi_key',
                   'pubchem_sid', 'gene_symbol', 'hgnc_id', 'uniprot_accessions']
    print("\nSample (first 3 rows):")
    print(df_unique[[c for c in sample_cols if c in df_unique.columns]].head(3).to_string(index=False))

    return df_unique


# ── Entry point ────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("ChEMBL Mechanism Data Download")
    print("=" * 60)

    try:
        from tqdm import tqdm
        tqdm_class = tqdm
    except ImportError:
        class SimpleTqdm:
            def __init__(self, total=None, desc=""):
                self.total = total; self.desc = desc; self.current = 0
                print(f"{desc}...")
            def update(self, n=1):
                self.current += n
                if self.total and self.current % max(1, self.total // 10) == 0:
                    print(f"  {self.desc}: {self.current/self.total*100:.1f}%")
            def __enter__(self): return self
            def __exit__(self, *args): print(f"  {self.desc}: Complete")
        tqdm_class = SimpleTqdm

    total = get_all_mechanisms_count()
    if not total:
        print("Could not get count. Exiting.")
        exit()

    data = download_all_mechanism_data(human_only=True, tqdm_class=tqdm_class)
    if data is not None:
        print(f"\nDone! End time: {datetime.now()}")
    else:
        print("Download failed.")
