# Predicate

- [**1. Overview**](#1-overview)
- [**2. Installation**](#2-installation)
- [**3. Quickstart**](#3-quickstart)
  - [**3.1 Virus biomass calculation**](#31-virus-biomass-calculation)
  - [**3.2 Virus biomass calculation with copy numbers**](#32-virus-biomass-calculation-with-copy-numbers)
  - [**3.3 Host-derived enforcement**](#33-host-derived-enforcement)
- [**4. Advanced settings**](#4-advanced-settings)
  - [**4.1 Adapting metabolites IDs to your model**](#41-adapting-metabolites-ids-to-your-model)
  - [**4.2 Adding additional metabolites to the virus biomass**](#42-adding-additional-metabolites-to-the-virus-biomass)
  - [**4.3 Changing the grams of nucleotide per mole of virus**](#43-changing-the-grams-of-nucleotide-per-mole-of-virus)
  - [**4.4 Number of APT and PPI molecules required**](#44-number-of-apt-and-ppi-molecules-required)
  - [**4.5 Changing amino acids molecular weights**](#45-changing-amino-acids-molecular-weights)
  - [**4.6 Mutations**](#46-mutations)
- [**5. License**](#5-license)

## 1. Overview

Predicate is a command line version application of the PREDICtor of Antiviral TargEts (PREDICATE) application introduced in: [journal.pcbi.1010903](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010903). 
The origianl code can be found at: [pymCADRE repository](https://github.com/draeger-lab/pymCADRE). 
The main purpose of this application is to:
- Calculate a virus biomass function and integrate it into an SBML metabolic network.
- Target selection through reaction-knockout and host-derived enforcement.

## 2. Installation

Predicate depends on [pymCADRE](https://github.com/draeger-lab/pymCADRE),
which requires a Python >= 3.8.5 version.

Predicate can be installed via **pip** package manager:

```bash
pip install git+https://github.com/alexOarga/predicate
```

## 3. Quickstart

### 3.1 Virus biomass calculation

:link: :arrow_right: [This example can be found [here](https://github.com/alexOarga/predicate/tree/main/examples/1_vbof)]

To calculate a virus biomass function, predicate requires a metabolic model in SBML format, a virus genome in FASTA format and the protein sequences in FASTA format.
To run predicate you will need to create a `config.yml` file. The file should look as follows:

```yaml
settings:
  variant_name: variant_name  # This is the name that will be used to name the output files.

reference_files:
  fasta_sequences: '<Path to the fasta file with the virus genome>'
  protein_sequences: '<Path to the fasta file with the protein sequences>'

metabolic_network:
  model: '<Path to the SBML file with the metabolic model>'
  objective: '<Id of the current biomass reaction of the model>' # The growth produced by this reaction will be compared with the growth produced by the virus biomass function.
```

To run predicate, run the following command:

```bash
predicate config.yml
```

This will generate a new SBML metabolic network with the virus biomass function integrated into the model. The id of this reactions will be `VBOF`.

### 3.2 Virus biomass calculation with copy numbers

:link: :arrow_right: [This example can be found [here](https://github.com/alexOarga/predicate/tree/main/examples/2_vbof_copy_numbers)]

The previous examples assumes that each proteins has a copy number of 1 which is generally not the case. To calculate the virus biomass function with copy numbers defined for each protein,
first, create a `config.yml` file as follows:

```yaml
settings:
  variant_name: variant_name 

reference_files:
  fasta_sequences: '<Path to the fasta file with the virus genome>'
  protein_sequences: '<Path to the fasta file with the protein sequences>'

metabolic_network:
  model: '<Path to the SBML file with the metabolic model>'
  objective: '<Id of the current biomass reaction of the model>' 

structural_proteins:
  - { id: 'YOUR_PROTEIN_1_ID', copy_number: <Protein 1 copy number> }  
  - { id: 'YOUR_PROTEIN_2_ID', copy_number: <Protein 2 copy number> }   

other_copy_numbers:
  Cnp: 1  # Copy number of non-structural proteins
  Cg: 1   # Copy number of the genome
```

In this file, you can specify the copy number of each protein in the `structural_proteins` section.
For each structural protein, you need to provide the id of the protein and the copy number. :warning: The id provided for each structural protein must match the id of the protein in the protein sequences file.
Besides, strutural proteins, you can also specify the copy number of non-structural proteins and the copy number of the genome in the `other_copy_numbers` section.
By default, all proteins that do not have a copy number defined in the `structural_proteins` section will be assumed to be non-structural proteins, and hence, will have the copy number defined in the `other_copy_numbers.Cnp` section.

Once you have created the `config.yml` file, you can run predicate as follows:

```bash
predicate config.yml
```

As before, this will generate a new SBML metabolic network with the virus biomass function integrated into the model. The id of this reactions will be `VBOF`.

### 3.3 Host-derived enforcement

:link: :arrow_right: [This example can be found [here](https://github.com/alexOarga/predicate/tree/main/examples/3_hde)]

To run host-derived enforcement on the generated virus biomass model, simply create a config file as before and add the following line: `run_hde: True` in the `settings` section:

```yaml
settings:
  variant_name: variant_name 
  run_hde: True

...
```

## 4. Advanced settings

:link: :arrow_right: [An example containing all the following sections can be found [here](https://github.com/alexOarga/predicate/tree/main/examples/4_full_config)]

### 4.1 Adapting metabolites IDs to your model

Predicate uses by default metabolites identifiers given by the [BiGG database](http://bigg.ucsd.edu/). If your SBML model uses different identifiers, you need to manually specify the IDs of all aminoacids, nucleotides (ATP, CTP, GTP, UTP) and
other energy requirements (ADP, H, H2O, PI, PPI). This can be done by adding the `metabolites_ids` section in the `metabolic_network` section of the `config.yml` file:

```yaml
metabolic_network:
  model: <your path>
  objective: <your biomass objective reaction>
  metabolites_ids: 
     #energy requirements
    - {'name': 'ADP', 'metabolite': 'adp_c'}
    - {'name': 'H_atom', 'metabolite': 'h_c'}
    - {'name': 'H2O', 'metabolite': 'h2o_c'}
    - {'name': 'PI', 'metabolite': 'pi_c'}
    - {'name': 'PPI', 'metabolite': 'ppi_c'}
    #nucleotides
    - {'name': 'ATP', 'metabolite': 'atp_c'}
    - {'name': 'CTP', 'metabolite': 'ctp_c'}
    - {'name': 'GTP', 'metabolite': 'gtp_c'}
    - {'name': 'UTP', 'metabolite': 'utp_c'}
    #amino acids
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
```

For example, if ADP is represented in your model as `my_adp` you can change the following line: `- {'name': 'ADP', 'metabolite': 'adp_c'}` to `- {'name': 'ADP', 'metabolite': 'my_adp'}`.

### 4.2 Adding additional metabolites to the virus biomass

If you want to add additional metabolites to the virus biomass function, you can do so by adding the `additional_metabolites` section in the `metabolic_network` section of the `config.yml` file:

```yaml
metabolic_network:
  model: <your path>
  objective: <your biomass objective reaction>
  additional_metabolites:
    - {'metabolite': 'pchol_hs_c', 'stoichiometry': 0.038400}
    - {'metabolite': 'pe_hs_c', 'stoichiometry': 0.014566}
    - {'metabolite': 'pail_hs_c', 'stoichiometry': 0.006621}
    - {'metabolite': 'ps_hs_c', 'stoichiometry': 0.001986}  
    - {'metabolite': 'chsterol_c', 'stoichiometry': 0.000012} 
    - {'metabolite': 'sphmyln_hs_c', 'stoichiometry': 0.001986}
```

For example, the above example adds 6 more metabolites to the virus biomass function. :warning: A positive stoichiometry means that the metabolite is produced by the virus biomass function, while a negative stoichiometry means that the metabolite is consumed by the virus biomass function.

### 4.3 Changing the grams of nucleotide per mole of virus

To change the number of grams of nucleotide per mole of virus, you can add the `nucleotides_per_virus_mol` section of the `config.yml` file:

```yaml
nucleotides_per_virus_mol: {
  A: 135.13,
  U: 112.09,
  G: 151.13,
  C: 111.1
}
```

In the above example, the number of grams of nucleotide A per virus mole is 135.13, the number of grams of nucleotide U per virus mole is 112.09, the number of grams of nucleotide G per virus mole is 151.13 and the number of grams of nucleotide C per virus mole is 111.1.

### 4.4 Number of APT and PPI molecules required

To change the number of ATP molecules required for the polymerization of a nucleotide, or to change the number of molecules of PPI required to bond 2 nucleotides, you can add the `other_parameters` section of the `config.yml` file:

```yaml
other_parameters:
  kATP: 4
  kPPi: 1
```

In the above example, 4 ATP molecules are required for the polymerization of a nucleotide and 1 PPI molecule is required to bond 2 nucleotides.

### 4.5 Changing amino acids molecular weights

Although this might not be frequent, you can change the molecular weights of the amino acids by adding the `amino_acids_and_weights` section of the `config.yml` file:

```yaml
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
```

In the above example, the molecular weight of amino acid A is 89.1, the molecular weight of amino acid R is 174.2, etc.

### 4.6 Mutations

Mutation visualization and analysis as in the original Predicate code is not
implemented yet.

## 5. License

This code was directly taked from the [pymCADRE repository](https://github.com/draeger-lab/pymCADRE) and shares the same license. Please see the [LICENSE](LICENSE) file for details.
