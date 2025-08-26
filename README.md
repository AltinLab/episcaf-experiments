

The purpose of this project is to extend the PepSeq assay to folded peptides (folded PepSeq).

To investigate if this is possible, we will work in 3 stages:

## 1. Sourced antibody

Using PDB entries containing human L+H chains and a bound antigen, define the conformational epitope,
remove all other residues from the PDB entry, and use RFdiffusion to generate new scaffolding for the epitope region.

No conformational epitope definition is currently proposed. While epitope residues can be identified, grouping them into distinct epitopes is not solved (AFAIK)

RFdiffusion is run on either

- The labelled subset of residues from the PDB entry (labelled by ABDB) and the antigen chain
- The labelled CDR subset of residues from the PDB entry and the antigen chain (in this case, verify there are no interactions outside of the CDR loops and antigen chain)

AF3 is run on only the antigen scaffold (no antibody) to better assess the filtering pipeline with a realistic use-case
(where the antibody is not known)

After AF3 quality assessment, the goal is to generate 30 epitopes x 30 different scaffolds for each to run
in an assay.

Also, 104AA is the variable region. There will then be 17 fixed flanking residues (2 N-terminal and 15 C-terminal):
GA-[104mer]-GDSLSWLLRLLNGNA

## 2. Censored antibody

Re-run the above experiment but without the antigen in the structure fed to RFdiffusion, the purpose here is 
to determine if we can still generate a bindable scaffold even without feeding the antibody chain to RFdiffusion.

## 3. Antigen-only

Given only a pathogenic protein, determine the likely epitope regions, scaffold them, and then measure reactivity.


## Installation

Modify your `.condarc` as specified here: https://github.com/nrbennet/dl_binder_design?tab=readme-ov-file#conda-environment-