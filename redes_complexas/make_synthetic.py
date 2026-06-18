#!/usr/bin/env python3
"""fetch_rcsb_metadata.py — consulta a API de dados do RCSB (GraphQL) e gera:
  - data/6B1T_entities.json : metadados das entidades/assembly
  - data/6B1T_annotation.csv : node,label (CADEIA:RESNAME:RESSEQ -> nome da proteína)
A anotação resultante pode alimentar diretamente:
    python run_pipeline.py --pdb 6B1T --annotation data/6B1T_annotation.csv

Requer acesso a data.rcsb.org. Rode na sua máquina (não no sandbox restrito)."""
import json, sys, os, urllib.request, warnings
warnings.filterwarnings("ignore")

URL = "https://data.rcsb.org/graphql"
QUERY = """
query structure($id: String!) {
  entry(entry_id: $id) {
    struct { title }
    exptl { method }
    em_3d_reconstruction { resolution }
    rcsb_entry_info { deposited_atom_count deposited_polymer_monomer_count
                      polymer_entity_count_protein molecular_weight }
    assemblies {
      rcsb_struct_symmetry { kind type symbol oligomeric_state stoichiometry }
    }
    polymer_entities {
      rcsb_polymer_entity { pdbx_description }
      rcsb_polymer_entity_container_identifiers { entity_id auth_asym_ids uniprot_ids }
      entity_poly { rcsb_sample_sequence_length }
      rcsb_entity_source_organism { scientific_name }
    }
  }
}"""

def fetch(pdb_id="6B1T"):
    body = json.dumps({"query": QUERY, "variables": {"id": pdb_id.upper()}}).encode()
    req = urllib.request.Request(URL, data=body,
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)["data"]["entry"]

def build_annotation(structure_path, entry):
    """Mapeia cada resíduo (CADEIA:RESNAME:RESSEQ) ao nome da proteína da entidade."""
    import protein_net as pn
    from Bio.PDB.Polypeptide import is_aa
    chain2name = {}
    for ent in entry["polymer_entities"]:
        name = ent["rcsb_polymer_entity"]["pdbx_description"]
        for ch in ent["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"]:
            chain2name[ch] = name
    st = pn.load_structure(structure_path, "6B1T")
    rows = ["node,label"]
    for r in list(st)[0].get_residues():
        if not is_aa(r, standard=False):
            continue
        ch = r.get_parent().id
        key = f"{ch}:{r.resname}:{r.id[1]}{r.id[2].strip()}"
        rows.append(f"{key},{chain2name.get(ch,'unknown')}")
    return "\n".join(rows)

if __name__ == "__main__":
    pdb_id = sys.argv[1] if len(sys.argv) > 1 else "6B1T"
    os.makedirs("data", exist_ok=True)
    entry = fetch(pdb_id)
    json.dump(entry, open(f"data/{pdb_id}_entities.json", "w"), indent=2)
    print("Título:", entry["struct"]["title"])
    print("Método:", entry["exptl"][0]["method"],
          "| resolução:", entry.get("em_3d_reconstruction", [{}])[0].get("resolution"))
    for a in entry.get("assemblies", []):
        for s in (a.get("rcsb_struct_symmetry") or []):
            print(f"Simetria: {s['type']} ({s['symbol']}), "
                  f"estado {s['oligomeric_state']}, estequiometria {s['stoichiometry']}")
    print("\nEntidades (proteínas):")
    for ent in entry["polymer_entities"]:
        cid = ent["rcsb_polymer_entity_container_identifiers"]
        print(f"  entidade {cid['entity_id']}: "
              f"{ent['rcsb_polymer_entity']['pdbx_description']} "
              f"| cadeias {cid['auth_asym_ids']} | UniProt {cid.get('uniprot_ids')}")
    # gera CSV de anotação se o .cif estiver disponível
    path = f"data/{pdb_id}.cif"
    if os.path.exists(path):
        open(f"data/{pdb_id}_annotation.csv", "w").write(build_annotation(path, entry))
        print(f"\n-> data/{pdb_id}_annotation.csv (use com run_pipeline.py --annotation)")
