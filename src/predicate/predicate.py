import os, sys
import cobra
from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis import single_reaction_deletion
from cobra.flux_analysis import flux_variability_analysis
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator
import matplotlib.patches as mpatches
import seaborn as sns
import warnings

from .utils import read_model


def read_sequences(path):
    all_sequences = {}
    # read fasta file
    fasta_sequences = SeqIO.parse(open(path),'fasta')
    for fasta in fasta_sequences:
        # store sequence id and sequence
        name, sequence = fasta.id, str(fasta.seq)
        # add info to dictionary
        all_sequences[name] = sequence
    #print(all_sequences)
    return all_sequences


def analyze_nucleotide_counts(all_sequences, csv_table):
    # if multiple input sequences 
    if csv_table != "":
        raise NotImplementedError("GISAID mutation table is not yet supported")
        
        # extract mutations info
        data_df = pd.read_table(csv_table,sep=';')
        print("data_df: ", data_df)
        
        # create df to store nucleotides count
        nucleotide_counts_df = pd.DataFrame(columns=['GISAID Accession ID', 'A','T','G','C'])

        for name, seq in all_sequences.items():

            new_row = {'GISAID Accession ID':data_df.loc[data_df['strain']==name]['gisaid_epi_isl'].values[0], 'A':seq.count("A"), 'T':seq.count("T"), 'G':seq.count("G"),'C':seq.count("C")}
            nucleotide_counts_df = nucleotide_counts_df.append(new_row, ignore_index=True)

        nucleotide_counts_df = nucleotide_counts_df.sort_values(by=["A","T","G","C"])

        ################################################################################################
        # Additional (optional) Calculations: calculate log2FC between mutated and not-mutated sequence
        ################################################################################################
        ref_seq = ref_rec[0].seq

        ref_seq_ncl_counts = [ref_seq.count("A"),ref_seq.count("T"),ref_seq.count("G"),ref_seq.count("C")]
        print("ref_seq_ncl_counts:", ref_seq_ncl_counts)

        log2FC_A,log2FC_T,log2FC_G,log2FC_C = [],[],[],[]
        for index,row in nucleotide_counts_df.iterrows():

            log2FC_A.append(np.log2(row['A']/ref_seq_ncl_counts[0]))
            log2FC_T.append(np.log2(row['T']/ref_seq_ncl_counts[1]))
            log2FC_G.append(np.log2(row['G']/ref_seq_ncl_counts[2]))
            log2FC_C.append(np.log2(row['C']/ref_seq_ncl_counts[3]))

        nucleotide_counts_df['log2FC_A'] = log2FC_A
        nucleotide_counts_df['log2FC_T'] = log2FC_T
        nucleotide_counts_df['log2FC_G'] = log2FC_G
        nucleotide_counts_df['log2FC_C'] = log2FC_C
        
        print("nucleotide_counts_df: ", nucleotide_counts_df)
        
            
        ##########################
        # Plot nucleotide counts
        ##########################
        mpl.rcParams['axes.spines.right'] = False

        plt.rcParams['font.size'] = '16'

        nucleotide_counts_df[['GISAID Accession ID','log2FC_A','log2FC_T','log2FC_G','log2FC_C']].plot(x='GISAID Accession ID',
                kind='bar',
                stacked=False,
                figsize=(19,10),
                color=["#d7191c", "#fdae61", "#abdda4", "#2b83ba"],
                width=0.7)

        plt.axhline(0.0, color='k', linewidth=1.2)
        plt.legend(loc="best")
        plt.xticks(rotation=89)
        plt.ylabel("Counts")

        plt.tight_layout()

        plt.savefig(path+"all_sequences_nucleotides_counts_with_log2FC.pdf")
        plt.savefig(path+"all_sequences_nucleotides_counts_with_log2FC.png")

        plt.show()
        
        ########################
        # Analyse Mutations
        #######################
        mutation_data = data_df['AA Mutations']
        print("mutation_data: ", mutation_data)
        
        # Important check: do all mutations involve letters related to amino acids or involve also "unknown" letters
        # allowed: all amino acid letters + del, stop, ins, dupl words + numbers
        allowed_chars = amino_acids+['s','t','o','p','d','e','l','i','n','s','u']+['0','1','2','3','4','5','6','7','8','9']

        c = 0
        for d in mutation_data:

            for m in d.split(","):

                for char in m.strip().split(" ")[1]:

                    if char not in allowed_chars:

                        # print character and index of occured sequence
                        print('weird char detected: ', char, c)
            c+=1
        
        
        mutations_df = pd.DataFrame(columns=['Virus Name', 'GISAID Accession ID', 'Mutations'])

        mutations_df['Virus Name'] = data_df['strain']
        mutations_df['GISAID Accession ID'] = data_df['gisaid_epi_isl']
        mutations_df['Mutations'] = data_df['AA Mutations']
        mutations_df.to_csv(path+"mutations_info_"+variant_name+"_mean_VBOF.csv")
        print('mutation_df: ', mutations_df)


        # store types of mutations observed in the 20 downloaded sequences
        mut_types = []

        for data in mutation_data:

            for mut in data.split(","):

                if mut.strip() not in mut_types:

                    mut_types.append(mut.strip())

        print('\nmut_types: ', mut_types)

        # count how oft a mutation appears
        mut_frequencies = {}
        for mut in mut_types:

            count = 0

            for data in mutation_data:

                if mut in data:

                    count +=1

            mut_frequencies[mut] = count

        mut_frequencies = dict(sorted(mut_frequencies.items(), key=lambda item: item[1]))
        print('\nmut_frequencies: ', mut_frequencies)
        
            
        ##############################
        # Plot Mutations Information
        ##############################
        plt.figure(figsize=(19,8)).gca().yaxis.set_major_locator(MaxNLocator(integer=True))
        plt.rcParams['font.size'] = '19'

        names = list(mut_frequencies.keys())
        values = list(mut_frequencies.values())

        barplot = plt.bar(range(len(mut_frequencies)), values, tick_label=names, color="#d7191c")

        # highlight different mutations
        for key,value in mut_frequencies.items():

            if 'del' in key:
                barplot[list(mut_frequencies).index(key)].set_color('#3288bd')

            if 'stop' in key:
                barplot[list(mut_frequencies).index(key)].set_color('#99d594')

            if 'ins' in key:
                barplot[list(mut_frequencies).index(key)].set_color('#fc8d59')

            if 'dupl' in key:
                barplot[list(mut_frequencies).index(key)].set_color('#984ea3')


        plt.xlim([-0.9,len(mut_frequencies)])
        plt.ylim([0,21])

        plt.xticks(rotation=89)

        plt.ylabel("Counts")
        plt.xlabel("Mutation Type")

        blue_patch = mpatches.Patch(color='#3288bd', label='Deletions')
        green_patch = mpatches.Patch(color='#99d594', label='Involves stop codon')
        plt.legend(handles=[blue_patch, green_patch])

        plt.tight_layout()

        plt.savefig(path+variant_name+"_sequences_mutation_counts.pdf")
        plt.savefig(path+variant_name+"_sequences_mutation_counts.png")

        plt.show()
        
        ###########################################
        # Plot Mutations Information per category
        ###########################################
        # store mutations per protein
        spike_prot_muts, n_prot_muts, m_prot_muts, e_prot_muts,nsp_prot_muts, ns_prot_muts = {},{},{},{},{},{}

        for key, value in mut_frequencies.items():
            if "Spike" in key:
                spike_prot_muts[key.split("Spike")[1].strip()] = value

            elif "N " in key:
                n_prot_muts[key.split("N")[1].strip()] = value

            elif "M " in key:
                m_prot_muts[key.split("M")[1].strip()] = value

            elif "E " in key:
                e_prot_muts[key.split("E")[1].strip()] = value

            elif "NSP" in key:
                nsp_prot_muts[key] = value

            else:
                ns_prot_muts[key] = value

        print("spike_prot_muts", spike_prot_muts)
        print("n_prot_muts", n_prot_muts)
        print("m_prot_muts", m_prot_muts)
        print("e_prot_muts", e_prot_muts)
        print("nsp_prot_muts", nsp_prot_muts)
        print("ns_prot_muts", ns_prot_muts)
        print("\n Total length", len(mut_frequencies))
        print("\n Total length", len(spike_prot_muts)+len(n_prot_muts)+len(m_prot_muts)+ len(e_prot_muts)+ len(nsp_prot_muts)+ len(ns_prot_muts))
        
        
        plt.rcParams['font.size'] = 18
        plt.rcParams['xtick.labelsize'] = 18

        colors=["#d53e4f", "#fc8d59", "#fee08b", "#e6f598", "#99d594", "#3288bd"]

        fig, axes = plt.subplots(2, 3, figsize=(30, 20), sharey=True)

        plt.subplots_adjust(left=0.125,
                            bottom=0.1, 
                            right=0.9, 
                            top=0.9, 
                            wspace=0.2, 
                            hspace=0.35)

        sns.barplot(ax=axes[0,0], x=list(spike_prot_muts.keys()), y=list(spike_prot_muts.values()),color=colors[0])
        axes[0,0].set_title("Spike Protein", fontweight="bold")
        axes[0,0].set_xticklabels(list(spike_prot_muts.keys()), rotation=90)

        sns.barplot(ax=axes[0,1], x=list(n_prot_muts.keys()), y=list(n_prot_muts.values()), color=colors[1])
        axes[0,1].set_title("N Protein", fontweight="bold")
        axes[0,1].set_xticklabels(list(n_prot_muts.keys()), rotation=90)

        sns.barplot(ax=axes[0,2], x=list(m_prot_muts.keys()), y=list(m_prot_muts.values()), color=colors[2])
        axes[0,2].set_title("Membrane Protein", fontweight="bold")
        axes[0,2].set_xticklabels(list(m_prot_muts.keys()), rotation=90)

        if len(e_prot_muts) != 0:
            # ! omit if no mutation on the E protein found (e.g., in gamma SARS-CoV-2 variant)
            sns.barplot(ax=axes[1,0], x=list(e_prot_muts.keys()), y=list(e_prot_muts.values()), color=colors[3])
            axes[1,0].set_title("Envelope Protein", fontweight="bold")
            axes[1,0].set_xticklabels(list(e_prot_muts.keys()), rotation=90)

        sns.barplot(ax=axes[1,1], x=list(nsp_prot_muts.keys()), y=list(nsp_prot_muts.values()), color=colors[4])
        axes[1,1].set_title("NSP Proteins",fontweight="bold")
        axes[1,1].set_xticklabels(list(nsp_prot_muts.keys()), rotation=90)

        sns.barplot(ax=axes[1,2], x=list(ns_prot_muts.keys()), y=list(ns_prot_muts.values()), color=colors[5])
        axes[1,2].set_title("NS Ptoteins", fontweight="bold")
        axes[1,2].set_xticklabels(list(ns_prot_muts.keys()), rotation=90)

        plt.tight_layout()

        plt.savefig(path+variant_name+"_sequences_mutation_counts_alltogether.pdf")
        plt.savefig(path+variant_name+"_sequences_mutation_counts_alltogether.png")

        plt.show()

    # if only one input sequence, simply count nucleotides
    else:
        
        nuc_seq = list(all_sequences.values())[0]
        print("Count nucleotide A: ", nuc_seq.count("A"))
        print("Count nucleotide C: ", nuc_seq.count("C"))
        print("Count nucleotide G: ", nuc_seq.count("G"))
        print("Count nucleotide T: ", nuc_seq.count("T"))
        print("Total length: ",len(nuc_seq))  
        return nuc_seq
        
        
def host_derived_enforcement(model, vbof_name, bof_name, growth_rate_vbof, growth_rate_bof, aa_weights):
    
    """ 
    Host-derived Enforcement
    
    Inputs:
        - model: SBML model 
        - vbof_name: string of vbof ID in given model
        - bof_name: string of bof ID in given model
        - growth_rate_vbof: virus growth rate 
        - growth_rate_bof: host maintenance rate
        
    Outputs:
        - data: list of lists with all predicted targets and the respective % reamining BOF
        
    """
    
    model.objective = bof_name
    print('BOF optimized:', model.optimize().objective_value)
    fva_mac = flux_variability_analysis(model,fraction_of_optimum=0.98)
    fva_mac = fva_mac.reset_index()
    fva_mac = fva_mac.rename(columns={'index': 'Reaction_ID', 'minimum': 'H_(F-)', 'maximum': 'H_(F+)'})

    model.objective = vbof_name
    print('VBOF optimized:', model.optimize().objective_value)
    fva_vbof = flux_variability_analysis(model,fraction_of_optimum=0.75)
    fva_vbof = fva_vbof.reset_index()
    fva_vbof = fva_vbof.rename(columns={'index': 'Reaction_ID', 'minimum': 'V_(F-)', 'maximum': 'V_(F+)'})

    fva = fva_mac.merge(fva_vbof, left_on='Reaction_ID', right_on='Reaction_ID')
    fba_vbof_fluxes = model.optimize().fluxes

    hde = {}
    zero_fba_vbof_flux = []

    for i in range(len(fva)):

        if fva.loc[i]['H_(F+)'] == fva.loc[i]['V_(F+)'] and fva.loc[i]['H_(F-)'] == fva.loc[i]['V_(F-)']:
            continue

        elif fva.loc[i]['H_(F+)'] > fva.loc[i]['V_(F+)'] and fva.loc[i]['H_(F-)'] >= fva.loc[i]['V_(F-)']:
            new_lower_bound = fva.loc[i]['H_(F+)'] - ((fva.loc[i]['H_(F+)']-fva.loc[i]['V_(F+)'])/2)
            new_upper_bound = fva.loc[i]['H_(F+)']

        elif fva.loc[i]['H_(F-)'] < fva.loc[i]['V_(F-)'] and fva.loc[i]['H_(F+)'] <= fva.loc[i]['V_(F+)']:
            new_lower_bound = fva.loc[i]['H_(F-)']
            new_upper_bound = fva.loc[i]['H_(F-)'] - ((fva.loc[i]['H_(F-)']-fva.loc[i]['V_(F-)'])/2)

        elif fva.loc[i]['H_(F-)'] > fva.loc[i]['V_(F-)'] and fva.loc[i]['H_(F+)'] < fva.loc[i]['V_(F+)']:
            new_lower_bound = fva.loc[i]['H_(F-)']
            new_upper_bound = fva.loc[i]['H_(F+)']

        elif fva.loc[i]['H_(F-)'] >= fva.loc[i]['V_(F-)'] and fva.loc[i]['H_(F+)'] <= fva.loc[i]['V_(F+)']:
            new_upper_bound = fva.loc[i]['H_(F+)'] - ((fva.loc[i]['H_(F+)']-fva.loc[i]['V_(F+)'])/2)
            new_lower_bound = fva.loc[i]['H_(F-)'] - ((fva.loc[i]['H_(F-)']-fva.loc[i]['V_(F-)'])/2)   

        else:
            print('else: ', fva.loc[i]['Reaction_ID'])


        with model: 
            # store old (initial) reaction bounds
            old_lb = model.reactions.get_by_id(fva.loc[i]['Reaction_ID']).lower_bound
            old_ub = model.reactions.get_by_id(fva.loc[i]['Reaction_ID']).upper_bound

            direction_before = model.reactions.get_by_id(fva.loc[i]['Reaction_ID']).reaction.split(" ")

            # ensure that lower bound will always be < than upper bound
            if new_lower_bound > new_upper_bound:
                model.reactions.get_by_id(fva.loc[i]['Reaction_ID']).lower_bound = round(new_upper_bound,6)
                model.reactions.get_by_id(fva.loc[i]['Reaction_ID']).upper_bound = round(new_lower_bound,6)
            else:
                model.reactions.get_by_id(fva.loc[i]['Reaction_ID']).lower_bound = round(new_lower_bound,6)
                model.reactions.get_by_id(fva.loc[i]['Reaction_ID']).upper_bound = round(new_upper_bound,6)

            direction_after = model.reactions.get_by_id(fva.loc[i]['Reaction_ID']).reaction.split(" ")

            # optimize BOF 
            model.objective = bof_name
            host = model.optimize().objective_value
            # optimize VBOF 
            model.objective = vbof_name
            virus = model.optimize().objective_value
                    
            # report reactions that reduce the virus growth rate to below 50% of its initial growth rate
            if virus < 0.5*growth_rate_vbof:
                
                print(direction_before,direction_after)


                if fba_vbof_fluxes[fva.loc[i]['Reaction_ID']] == 0.0:

                    # ignore those, since they don't influence the virus's growth
                    print('Zero FBA flux when VBOF optimized: ', fva.loc[i]['Reaction_ID'])
                    zero_fba_vbof_flux.append(fva.loc[i]['Reaction_ID'])

                else:

                    if  new_lower_bound >= 0 and new_upper_bound > 0:
                        hde[fva.loc[i]['Reaction_ID']] = [host, virus]
                        print('old bounds:', old_lb, old_ub)
                        print('Forward:', fva.loc[i]['Reaction_ID'], new_lower_bound, new_upper_bound)

                    elif new_lower_bound < 0 and new_upper_bound <= 0:
                        hde[fva.loc[i]['Reaction_ID']] = [host, virus]
                        print('old bounds:', old_lb, old_ub)
                        print('Reverse:', fva.loc[i]['Reaction_ID'], new_lower_bound, new_upper_bound)

                    elif new_lower_bound < 0 and new_upper_bound > 0:
                        hde[fva.loc[i]['Reaction_ID']] = [host, virus]
                        print('old bounds:', old_lb, old_ub)
                        print('Both directions:', fva.loc[i]['Reaction_ID'], new_lower_bound, new_upper_bound)

                    else:
                        print('old bounds:', old_lb, old_ub)
                        print('Unknown:', fva.loc[i]['Reaction_ID'], new_lower_bound, new_upper_bound)
                    
                print("------------------------------------------")        

    
    # remove BOF and VBOF if appear in targets
    if vbof_name in hde:
        hde.pop(vbof_name)
    if bof_name in hde:
        hde.pop(bof_name)
        

    print('Total Host Derived Enformcement (HDE) targets', len(hde))
    print() 
    
    hde_vbof = {}
    for key in hde.keys():
        hde_vbof[key] = (hde[key][1]/growth_rate_vbof)*100
    
    hde_vbof = {k: v for k, v in sorted(hde_vbof.items(), key=lambda item: item[1])}
    
    data = [[key, hde_vbof[key]] for key in hde_vbof.keys()]
    print('HDE Targets sorted:', data)
    print()
    print('HDE - top first target:', data[0])
    
    return data


def create_vbof_and_compute_hde(
    variant_name,    
    all_sequences, 
    amino_acids,
    aa_weights,
    Cg,
    Cnp,
    kATP,
    kPPi,
    grams_perm_mv, 
    records, 
    prot_copy_numbers, 
    model,
    metabolites, 
    additional_metabolites, 
    objective,
    csv_table,
    run_hde
):

    count_iterations = 0

    # store stoichiometric coefficients of all VBOF compounds of all input sequences
    stoichiometric_coeffs = {} 

    # store results of KO and HDE experiments of all input sequences
    ko_results_all_sequences, hde_results_all_sequences = {}, {} 

    # store all counts from all input sequences
    amino_acids_counts = {}
    for aa in amino_acids:
        amino_acids_counts[aa] = []

    # for all input sequences
    for ids, seq in all_sequences.items():
        
        print(f"Seq.Nr. {count_iterations}. Seq.ID. {ids}")
        nucl_aa = [] # store counts of nucleotides and amino acids in the current sequence
        count_dict = {}
        
        #********************************************
        # Nucleotides Investment
        #********************************************
        
        # count nucleotides in all sequences
        A, U, G, C = 0, 0, 0, 0

        # Get nucleotide count from .fasta file
        A = seq.count('A')
        U = seq.count('T')
        G = seq.count('G')
        C = seq.count('C')
        print("    Nucleotide count A,U,G,C: ",A,U,G,C)

        # compute moles of each nulceotide in a mole of virus particle
        AUtot = Cg * 2 * (A + U)
        GCtot = Cg * 2 * (G + C)
        print('    Moles of nucleotide in a model virus particle - AU: ', AUtot)
        print('    Moles of nucleotide in a model virus particle - CG: ', GCtot)

        # convert moles of nucleotides into grams of nucleotide per mole of virus
        GA = AUtot * grams_perm_mv['A']
        GU = AUtot * grams_perm_mv['U']
        GG = GCtot * grams_perm_mv['G']
        GC = GCtot * grams_perm_mv['C']

        # compute the mass of all genome components
        Gi = GA + GU + GG + GC
        print('    Mass of all genome compontents: ', Gi)    
        
        # Count Amino Acids in proteins
        for letter in amino_acids:
            num = 0
            for i in range(len(records)):
                prot_id = records[i].id
                if prot_id in prot_copy_numbers:
                    num += records[i].seq.count(letter) * prot_copy_numbers[prot_id]
                else:
                    num += records[i].seq.count(letter) * Cnp
            
            # for current amino acid
            count_dict[letter] = [num]

        # **********************************************************
        # Introduce Mutations in reference protein sequence 
        # **********************************************************   
        # if table with mutations is given
        if csv_table != '':
            raise NotImplementedError("GISAID mutation table is not yet supported")
            
            for mut in list(data_df[data_df['strain']==ids]['AA Mutations'])[0].split(","):

                if 'del' in mut or 'stop' in mut:

                    to_delete = mut.strip().split(" ")[1][0]

                    # consider different weights of structural proteins when deleting/inserting amino acid
                    if 'Spike' in mut:
                        count_dict[to_delete][0] = count_dict[to_delete][0] - Cs
                    elif "E " in mut: #e
                        count_dict[to_delete][0] = count_dict[to_delete][0] - Ce
                    elif "M " in mut: #m
                        count_dict[to_delete][0] = count_dict[to_delete][0] - Cm
                    elif "N " in mut: #n
                        count_dict[to_delete][0] = count_dict[to_delete][0] - Cn
                    else:
                        count_dict[to_delete][0] = count_dict[to_delete][0] - 1   

                elif 'ins' in mut or 'dupl' in mut:

                    to_delete = mut.strip().split(" ")[1][0]

                    # consider differnt weights of structural proteins when deleting/inserting amino acid
                    if 'Spike' in mut:
                        count_dict[to_delete][0] = count_dict[to_delete][0] + Cs
                    elif "E " in mut: #e
                        count_dict[to_delete][0] = count_dict[to_delete][0] + Ce
                    elif "M " in mut: #m
                        count_dict[to_delete][0] = count_dict[to_delete][0] + Cm
                    elif "N " in mut: #n
                        count_dict[to_delete][0] = count_dict[to_delete][0] + Cn
                    else:
                        count_dict[to_delete][0] = count_dict[to_delete][0] + 1

                else:

                    to_delete = "".join(filter(str.isupper, mut.strip().split(" ")[1]))[0]
                    to_add = "".join(filter(str.isupper, mut.strip().split(" ")[1]))[1]

                    # consider differnt weights of structural proteins when deleting/inserting amino acid
                    if 'Spike' in mut:
                        count_dict[to_delete][0] = count_dict[to_delete][0] - Cs
                        count_dict[to_add][0] = count_dict[to_add][0] + Cs
                    elif "E " in mut: #e
                        count_dict[to_delete][0] = count_dict[to_delete][0] - Ce
                        count_dict[to_add][0] = count_dict[to_add][0] + Ce
                    elif "M " in mut: #m
                        count_dict[to_delete][0] = count_dict[to_delete][0] - Cm
                        count_dict[to_add][0] = count_dict[to_add][0] + Cm

                    elif "N " in mut: #n
                        count_dict[to_delete][0] = count_dict[to_delete][0] - Cn
                        count_dict[to_add][0] = count_dict[to_add][0] + Cn
                    else:
                        count_dict[to_delete][0] = count_dict[to_delete][0] - 1
                        count_dict[to_add][0] = count_dict[to_add][0] + 1
        
        # store counts of amino acids
        for keys, values in count_dict.items():
            amino_acids_counts[keys].append(values[0])
        
        # ********************************************
        # Amino Acid Investment
        # ********************************************
        
        count_aa_weight = {}
        for letter in amino_acids:
            count_aa_weight[letter] = count_dict[letter][0] * aa_weights[letter]
        
        # Compute mass of all proteome components
        Gj = 0    
        for letter in amino_acids:
            Gj += count_aa_weight[letter]

        print('    Mass of all proteome sequences: ', Gj)
        
        #**********************
        # Total Viral Mass
        #**********************
        Mv = Gi + Gj
        print('    Total Viral mass (Genome mass + Proteome mass):', Mv)


        SA = 1000 * (AUtot / Mv) #stoichiometric coefficient of Adenine
        SG = 1000 * (GCtot / Mv) #stoichiometric coefficient of Guanine
        print("    Stoichiometric coefficient of Adenine: ", SA)
        print("    Stoichiometric coefficient of Guanine: ", SG)

        for letter in amino_acids:
            letter_count = count_dict[letter][0]
            count_dict[letter].append(1000 * (letter_count / Mv)) #entry 4 = SjX

        #**********************
        # ATP Requirement
        #**********************

        Xj = 0
        for letter in amino_acids:
            Xj += count_dict[letter][0]

        ATOT = (Xj * kATP) - kATP
        SATP = 1000 * (ATOT / Mv)
        print('    Stoichiometric coefficient of ATP: ', SATP)

        #*****************************
        # Pyrophosphate Liberation
        #*****************************
        PG = ((A + U + G + C) * kPPi) - kPPi
        PR = ((A + U + G + C) * kPPi) - kPPi

        PTOT = Cg * (PG + PR)

        SPPi = 1000 * (PTOT / Mv)
        print('    Stoichiometric coefficient of Pyrophosphate (PPi) Liberation: ', SPPi)
                
        #***************************************************************
        # VBOF Reconstruction + Incorporation in host-virus cell 
        #***************************************************************
        with model: 
            
            # add new VBOF
            reaction = Reaction('VBOF')
            reaction.name = 'Viral biomass objective function'
            reaction.lower_bound = 0 
            reaction.upper_bound = 1000 

            reaction.add_metabolites({
                #energy requirements
                model.metabolites.get_by_id(metabolites['ADP']): SATP,  
                model.metabolites.get_by_id(metabolites['H_atom']): SATP,  
                model.metabolites.get_by_id(metabolites['H2O']): -SATP,  
                model.metabolites.get_by_id(metabolites['PI']): SATP, 
                model.metabolites.get_by_id(metabolites['PPI']): SPPi, 

                #nucleotides: left-hand terms 
                model.metabolites.get_by_id(metabolites['ATP']): -(SATP + SA), 
                model.metabolites.get_by_id(metabolites['CTP']): -SG,
                model.metabolites.get_by_id(metabolites['GTP']): -SG,
                model.metabolites.get_by_id(metabolites['UTP']): -SA,

                #amino acids: left-hand terms
                model.metabolites.get_by_id(metabolites['A']): -count_dict['A'][1],
                model.metabolites.get_by_id(metabolites['R']): -count_dict['R'][1], 
                model.metabolites.get_by_id(metabolites['N']): -count_dict['N'][1],
                model.metabolites.get_by_id(metabolites['D']): -count_dict['D'][1],
                model.metabolites.get_by_id(metabolites['C']): -count_dict['C'][1],
                model.metabolites.get_by_id(metabolites['E']): -count_dict['E'][1],
                model.metabolites.get_by_id(metabolites['Q']): -count_dict['Q'][1],
                model.metabolites.get_by_id(metabolites['G']): -count_dict['G'][1],
                model.metabolites.get_by_id(metabolites['H']): -count_dict['H'][1],
                model.metabolites.get_by_id(metabolites['I']): -count_dict['I'][1],
                model.metabolites.get_by_id(metabolites['L']): -count_dict['L'][1],
                model.metabolites.get_by_id(metabolites['K']): -count_dict['K'][1],
                model.metabolites.get_by_id(metabolites['M']): -count_dict['M'][1],
                model.metabolites.get_by_id(metabolites['F']): -count_dict['F'][1],
                model.metabolites.get_by_id(metabolites['P']): -count_dict['P'][1],
                model.metabolites.get_by_id(metabolites['S']): -count_dict['S'][1],
                model.metabolites.get_by_id(metabolites['T']): -count_dict['T'][1],
                model.metabolites.get_by_id(metabolites['W']): -count_dict['W'][1],
                model.metabolites.get_by_id(metabolites['Y']): -count_dict['Y'][1],
                model.metabolites.get_by_id(metabolites['V']): -count_dict['V'][1], 

            })

            # add extra metabolites
            for mid, stoi in additional_metabolites.items():
                reaction.add_metabolites({
                    model.metabolites.get_by_id(mid): stoi
                })


            model.add_reactions([reaction])

            s_matrix = cobra.util.array.create_stoichiometric_matrix(model=model, array_type='DataFrame')
            # store stoichiometric coeffs as absolute values
            stoichiometric_coeffs[ids] = list(map(abs,list(s_matrix['VBOF'].loc[s_matrix['VBOF'] != 0])))
            # store compounds participating in VBOF
            vbof_compounds = list(s_matrix['VBOF'].loc[s_matrix['VBOF'] != 0].index)

            print("    VBOF compounds: ",vbof_compounds)

            #######################################
            # Save model
            #######################################
            from cobra.io import load_json_model, save_json_model, load_matlab_model, save_matlab_model, read_sbml_model, write_sbml_model

            name = f"Seq_{count_iterations}_Id_{ids}_Var_{variant_name}_VBOF.xml"
            file = os.path.join(os.getcwd(), name)
            write_sbml_model(model, file)
            print("    Model saved as: ", file)
            
            #***************************************
            # Targets Prediction Experiments
            #***************************************

            # (1) Single reaction knock-out experiments
            print("    \n * Single Reaction Knock-outs * \n")
            
            # optimize host
            model.objective = objective
            growth_rate_bec = model.optimize().objective_value
            react_ko_mac = single_reaction_deletion(model)
            react_ko_mac = react_ko_mac.rename(columns={'growth': 'growth_host', 'status': 'status_host'})
            # change indices from numbers to reaction IDs
            new_indices = []
            for i in react_ko_mac.index:
                new_indices.append(list(i)[0])
            react_ko_mac.set_index([pd.Index(new_indices)], inplace=True)

            # optimize virus
            model.objective = 'VBOF'
            growth_rate_vbof = model.optimize().objective_value
            react_ko_vbof = single_reaction_deletion(model)
            react_ko_vbof = react_ko_vbof.rename(columns={'growth': 'growth_virus', 'status': 'status_virus'})
            # change indices from numbers to reaction IDs
            new_indices2 = []
            
            for i in react_ko_vbof.index:
                new_indices2.append(list(i)[0])
            react_ko_vbof.set_index([pd.Index(new_indices2)], inplace=True)

            # merge results from BOF and VBOF knock-out
            growth_virus = react_ko_vbof["growth_virus"]
            merged = react_ko_mac.join(growth_virus)

            # for each reaction: calculate what percentage of host growth remains after single reaction knock-out
            percent_bec = []
            for i in merged.index:
                percent_bec.append(merged.loc[i]['growth_host']/growth_rate_bec*100)
            merged.insert(loc = 0, column='%growth_host', value=percent_bec)

            percent_vbof = []
            for i in merged.index:
                percent_vbof.append(merged.loc[i]['growth_virus']/growth_rate_vbof*100)     
            merged.insert(loc = 0, column='%growth_virus', value=percent_vbof)

            # keep only interesting -- > host 99% and virus 0%
            print('    KO Targets: ', list(merged.loc[(merged['%growth_host'] > 99) & (merged['%growth_virus'] == 0.0)].index))
            ko_results_all_sequences[ids] = list(merged.loc[(merged['%growth_host'] > 99) & (merged['%growth_virus'] == 0.0)].index)

            #######################################
            # Host-derived enforcement
            #######################################
            if run_hde:
                hde_results_all_sequences[ids] = host_derived_enforcement(model, 'VBOF', objective, growth_rate_vbof, growth_rate_bec, aa_weights)

            
            
        print('\n')
        count_iterations += 1


def validate_metabolites_ids_map(model, metabolites_ids):
    for original, mid in metabolites_ids.items():
        try:
            m = model.metabolites.get_by_id(mid)
        except KeyError:
            raise KeyError(f"""
                The metabolite ID {mid} was given for biochemical '{original}' in the configuration,
                however, this metabolite ID does not exist in the model. Please ensure that
                there is a metabolite with this ID in the model. You can also change the
                metabolite ID assigned to the biochemical '{original}' in the .yml configuration file, by setting
                this field:

                    metabolic_network:
                        ...
                        metabolites_ids:
                            - ...
                            - {{'name': '{original}', 'metabolite': '<your metabolite ID>' }}
                            - ...
                        ...
            """)


def validate_additional_metabolites(model, additional_metabolites):
    for mid, _ in additional_metabolites.items():
        try:
            m = model.metabolites.get_by_id(mid)
        except KeyError:
            raise KeyError(f"""
                Metabolite with ID {mid} defined in the 'additional_metabolites' field of the
                configuration file does not exist in the model. Please ensure that there is a
                metabolite with this ID in the model.
            """)


def validate_proteins(records, structural_proteins):
    pro_seq_ids = set([record.id for record in records])
    for entry in structural_proteins:
        prot_id = entry['id']
        prot_cn = entry['copy_number']
        if prot_id not in pro_seq_ids:
            warnings.warn(f"""
                Protein with ID {prot_id} defined in the 'structural_proteins' field of the
                configuration file does not exist in the protein sequences file. The copy number
                that will be used for this protein is 1, instead of the copy number defined in the
                configuration file: {prot_cn}.
            """)


def validate_vbof(model):
    try:
        model.reactions.get_by_id('VBOF')
        warnings.warn(f"""
            The model already contains a reaction with ID 'VBOF'. This reaction will be
            overwritten by a new reaction with this ID.
        """)
        model.remove_reactions(['VBOF'])
        return model
    except KeyError:
        # reaction does not exist, no warning here
        return model


def validate_objective(model, objective):
    try:
        model.reactions.get_by_id(objective)
    except KeyError:
        raise KeyError(f"""
            Objective reaction with ID '{objective}', defined in the configuration file does not exist in the model.
            Please ensure that there is a reaction with this ID in the model.
        """)


def main(data):

    path = os.getcwd()
    
    variant_name = data['settings']['variant_name']
    run_hde = bool(data['settings']['run_hde'])

    # read input fasta files
    all_sequences = read_sequences(data['reference_files']['fasta_sequences'])
    records = list(SeqIO.parse(data['reference_files']['protein_sequences'], "fasta"))
    csv_table = data['mutated_files']['gisaid_csv_table']

    # validate proteins
    structural_proteins = data['structural_proteins']
    validate_proteins(records, structural_proteins)

    # analyze nucleotide counts
    nuc_seq = analyze_nucleotide_counts(all_sequences, csv_table)

    # read model
    model = read_model(data['metabolic_network']['model'])
    metabolites_ids = data['metabolic_network']['metabolites_ids']
    metabolites_ids = dict(zip([x['name'] for x in metabolites_ids], [x['metabolite'] for x in metabolites_ids]))
    additional_metabolites = data['metabolic_network']['additional_metabolites']
    additional_metabolites = dict(zip([x['metabolite'] for x in additional_metabolites], [x['stoichiometry'] for x in additional_metabolites]))
    validate_metabolites_ids_map(model, metabolites_ids)
    validate_additional_metabolites(model, additional_metabolites)
    model = validate_vbof(model)
    validate_objective(model, data['metabolic_network']['objective'])

    # setup variables
    prot_copy_numbers = dict(zip([x['id'] for x in structural_proteins], [x['copy_number'] for x in structural_proteins]))

    amino_acids_and_weights = data['amino_acids_and_weights']
    aa_weights = dict(zip([x['amino_acid'] for x in amino_acids_and_weights], [x['molecular_weight'] for x in amino_acids_and_weights]))
    
    amino_acids = [x['amino_acid'] for x in amino_acids_and_weights]

    create_vbof_and_compute_hde(
        variant_name,
        all_sequences, 
        amino_acids,
        aa_weights,
        data['other_copy_numbers']['Cg'],
        data['other_copy_numbers']['Cnp'],
        data['other_parameters']['kATP'],
        data['other_parameters']['kPPi'],
        data['nucleotides_per_virus_mol'], 
        records, 
        prot_copy_numbers, 
        model,
        metabolites_ids, 
        additional_metabolites, 
        data['metabolic_network']['objective'],
        csv_table,
        run_hde
    )


