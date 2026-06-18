# Query GraphQL da API de dados do RCSB (data.rcsb.org/graphql)
# Variável: { "id": "6B1T" }
query structure($id: String!) {
  entry(entry_id: $id) {
    rcsb_id
    struct { title }
    exptl { method }
    em_3d_reconstruction { resolution }
    symmetry { space_group_name_H_M }
    rcsb_entry_info {
      deposited_atom_count
      deposited_polymer_monomer_count
      polymer_entity_count_protein
      molecular_weight
    }
    assemblies {
      rcsb_assembly_container_identifiers { assembly_id }
      pdbx_struct_assembly { rcsb_details rcsb_candidate_assembly }
      rcsb_struct_symmetry { kind type symbol oligomeric_state stoichiometry }
      rcsb_assembly_info { modeled_polymer_monomer_count polymer_entity_count_protein }
    }
    polymer_entities {
      rcsb_polymer_entity { pdbx_description }
      rcsb_polymer_entity_container_identifiers { entity_id auth_asym_ids uniprot_ids }
      polymer_entity_instances {
        rcsb_polymer_entity_instance_container_identifiers { auth_asym_id }
      }
      entity_poly { rcsb_sample_sequence_length }
      rcsb_entity_source_organism { scientific_name }
    }
  }
}
