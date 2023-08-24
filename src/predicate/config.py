from omegaconf import OmegaConf

default_yaml = """
settings:
  run_hde: False
  variant_name: ''

mutated_files:
  fasta_sequences: ''
  gisaid_csv_table: ''

reference_files:
  fasta_sequences: ''
  protein_sequences: ''

metabolic_network:
  model: ''
  objective: ''
  metabolites_ids: 
    #energy requirements
    - {'name': 'ADP', 'metabolite': 'adp_c'}
    - {'name': 'H_atom', 'metabolite': 'h_c'}
    - {'name': 'H2O', 'metabolite': 'h2o_c'}
    - {'name': 'PI', 'metabolite': 'pi_c'}
    - {'name': 'PPI', 'metabolite': 'ppi_c'}
    #nucleotides: left-hand terms 
    - {'name': 'ATP', 'metabolite': 'atp_c'}
    - {'name': 'CTP', 'metabolite': 'ctp_c'}
    - {'name': 'GTP', 'metabolite': 'gtp_c'}
    - {'name': 'UTP', 'metabolite': 'utp_c'}
    #amino acids: left-hand terms
    - {'name': 'A', 'metabolite': 'ala__L_c'}
    - {'name': 'R', 'metabolite': 'arg__L_c'}
    - {'name': 'N', 'metabolite': 'asn__L_c'}
    - {'name': 'D', 'metabolite': 'asp__L_c'}
    - {'name': 'C', 'metabolite': 'cys__L_c'}
    - {'name': 'E', 'metabolite': 'glu__L_c'}
    - {'name': 'Q', 'metabolite': 'gln__L_c'}
    - {'name': 'G', 'metabolite': 'gly_c'}
    - {'name': 'H', 'metabolite': 'his__L_c'}
    - {'name': 'I', 'metabolite': 'ile__L_c'}
    - {'name': 'L', 'metabolite': 'leu__L_c'}
    - {'name': 'K', 'metabolite': 'lys__L_c'}
    - {'name': 'M', 'metabolite': 'met__L_c'}
    - {'name': 'F', 'metabolite': 'phe__L_c'}
    - {'name': 'P', 'metabolite': 'pro__L_c'}
    - {'name': 'S', 'metabolite': 'ser__L_c'}
    - {'name': 'T', 'metabolite': 'thr__L_c'}
    - {'name': 'W', 'metabolite': 'trp__L_c'}
    - {'name': 'Y', 'metabolite': 'tyr__L_c'}
    - {'name': 'V', 'metabolite': 'val__L_c'}
    
  additional_metabolites:
    - {}


structural_proteins:
  - {}

  
other_copy_numbers:
  Cnp: 1
  Cg: 1

other_parameters:
  kATP: 4
  kPPi: 1

amino_acids_and_weights:
  - { amino_acid: 'A', molecular_weight: 89.1 }
  - { amino_acid: 'R', molecular_weight: 174.2 }
  - { amino_acid: 'N', molecular_weight: 132.1 }
  - { amino_acid: 'D', molecular_weight: 133.1 }
  - { amino_acid: 'C', molecular_weight: 121.2 }
  - { amino_acid: 'E', molecular_weight: 147.1 }
  - { amino_acid: 'Q', molecular_weight: 146.2 }
  - { amino_acid: 'G', molecular_weight: 75.1 }
  - { amino_acid: 'H', molecular_weight: 155.2 }
  - { amino_acid: 'I', molecular_weight: 131.2 }
  - { amino_acid: 'L', molecular_weight: 131.2 }
  - { amino_acid: 'K', molecular_weight: 146.2 }
  - { amino_acid: 'M', molecular_weight: 149.2 }
  - { amino_acid: 'F', molecular_weight: 165.2 }
  - { amino_acid: 'P', molecular_weight: 115.1 }
  - { amino_acid: 'S', molecular_weight: 105.1 }
  - { amino_acid: 'T', molecular_weight: 119.1 }
  - { amino_acid: 'W', molecular_weight: 204.2 }
  - { amino_acid: 'Y', molecular_weight: 181.2 }
  - { amino_acid: 'V', molecular_weight: 117.1 }

nucleotides_per_virus_mol: {
  A: 135.13,
  U: 112.09,
  G: 151.13,
  C: 111.1
}
"""


def default_config():
    return OmegaConf.create(default_yaml)