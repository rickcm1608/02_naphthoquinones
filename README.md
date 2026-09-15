# HSA - Naphthoquinones Docking and Molecular Dynamics

Computational pipeline for the study of naphthoquinone derivatives (D01, D02, D03, RWF) as potential binders of Human Serum Albumin (HSA). This repository contains all analysis scripts used in the associated publication.

## Workflow Overview

```
01_prepare_receptor          -> Clean PDB, fix missing atoms/H, convert to PDBQT
02_binding_site_prediction   -> P2Rank binding pocket prediction (pocket placed manually from this output)
03_ligand_optimization       -> ANI-2x geometry optimization (Neural Network Potential)
04_docking                   -> GPU-accelerated docking with UniDock (Vina scoring)
05_build_topology            -> AMBER19SB + GAFF2 + TIP3P topology (OpenMM)
06_interaction_fingerprint   -> Protein-ligand interaction fingerprint (ProLIF)
07_equilibration             -> NVT + NPT equilibration (0.5 ns + 0.5 ns, OpenMM)
08_production_md             -> NPT production MD -- 11 x 5 ns = 55 ns (CUDA)
09_concatenate_trajectories  -> PBC unwrap, complex centering, Ca alignment
10_gnina_rescoring           -> CNN-based rescoring of MD frames with GNINA
11_rmsd_ligand               -> Ligand RMSD (protein Ca fit, reference = frame 0)
12_rmsd_ligand_nofit         -> Ligand RMSD (no local fit, reference = frame 100)
13_rmsd_docking_pose         -> Ligand RMSD vs initial docking pose
14_rmsd_protein_backbone     -> Protein Ca RMSD over production trajectory
15_hbonds                    -> Protein-ligand hydrogen bond counts per frame
```

The pipeline was run independently for each of the four ligands (D01, D02, D03, RWF).

## Software & Force Fields

| Tool | Purpose |
|------|---------|
| OpenMM | MD engine |
| AMBER ff19SB | Protein force field |
| GAFF2 | Ligand force field |
| TIP3P | Water model |
| UniDock | GPU docking |
| P2Rank | Binding site prediction |
| GNINA | CNN rescoring |
| ANI-2x | Ligand geometry optimization |
| MDAnalysis | Trajectory analysis |
| ProLIF | Interaction fingerprints |
| pytraj | Protein backbone RMSD |

## Data Files

All data files are organized by ligand under `data/<LIGAND>/`.

| File | Description |
|------|-------------|
| `gnina_affinity.csv` | Minimized affinity per MD frame (kcal/mol) |
| `gnina_cnn_score.csv` | GNINA CNN score per frame |
| `gnina_cnn_affinity.csv` | GNINA CNN affinity per frame |
| `rmsd_ligand.csv/png` | Ligand RMSD (protein Ca fit, ref = frame 0) |
| `rmsd_ligand_nofit.csv/png` | Ligand RMSD (no local fit, ref = frame 100) |
| `rmsd_docking_pose.csv/png` | Ligand RMSD vs initial docking pose |
| `rmsd_protein_backbone.csv/png` | Protein Ca RMSD over trajectory |
| `hbonds_counts.csv` | H-bond counts per frame (PROT->UNK and UNK->PROT) |
| `ligand_2d_preview.png` | 2D structure of the ligand (from SMILES) |

## Boltz-2 Predicted Poses

Boltz-2 (protein-ligand co-folding, via Rowan Scientific) was used to predict binding poses and estimate binding affinity for three naphthoquinone derivatives (D01, D02, D03) at site IB of HSA. Three replicate runs were performed per compound (9 structures total).

| Compound | Replicate | pTM | ipTM | Affinity Prob. | Predicted pIC50 | Predicted IC50 (M) |
|----------|-----------|-----|------|----------------|-----------------|---------------------|
| D01 | 1 | 0.959 | 0.954 | 0.567 | 4.73 | 1.87e-5 |
| D01 | 2 | 0.950 | 0.964 | 0.542 | 4.58 | 2.64e-5 |
| D01 | 3 | 0.943 | 0.954 | 0.566 | 4.73 | 1.87e-5 |
| D02 | 1 | 0.951 | 0.946 | 0.643 | 5.33 | 4.63e-6 |
| D02 | 2 | 0.954 | 0.941 | 0.630 | 5.27 | 5.41e-6 |
| D02 | 3 | 0.937 | 0.949 | 0.603 | 5.14 | 7.18e-6 |
| D03 | 1 | 0.943 | 0.933 | 0.671 | 5.70 | 2.01e-6 |
| D03 | 2 | 0.949 | 0.923 | 0.682 | 5.76 | 1.73e-6 |
| D03 | 3 | 0.945 | 0.909 | 0.685 | 5.79 | 1.61e-6 |

D03 shows the highest predicted affinity (mean pIC50 5.75), consistent with docking and MD results. Boltz-2 predictions are available for D01, D02, and D03 only (not RWF).

| File | Description |
|------|-------------|
| `boltz2/D01-1_boltz2.pdb` | Boltz-2 pose -- D01, replicate 1 |
| `boltz2/D01-2_boltz2.pdb` | Boltz-2 pose -- D01, replicate 2 |
| `boltz2/D01-3_boltz2.pdb` | Boltz-2 pose -- D01, replicate 3 |
| `boltz2/D02-1_boltz2.pdb` | Boltz-2 pose -- D02, replicate 1 |
| `boltz2/D02-2_boltz2.pdb` | Boltz-2 pose -- D02, replicate 2 |
| `boltz2/D02-3_boltz2.pdb` | Boltz-2 pose -- D02, replicate 3 |
| `boltz2/D03-1_boltz2.pdb` | Boltz-2 pose -- D03, replicate 1 |
| `boltz2/D03-2_boltz2.pdb` | Boltz-2 pose -- D03, replicate 2 |
| `boltz2/D03-3_boltz2.pdb` | Boltz-2 pose -- D03, replicate 3 |
| `boltz2/boltz2_results.csv` | Full Boltz-2 scoring results |

## Key Results

**GNINA rescoring** (mean minimized affinity over MD frames, kcal/mol):

| Ligand | Mean affinity | Best frame |
|--------|--------------|------------|
| D01 | -5.82 | -7.42 |
| D02 | -8.04 | -9.43 |
| D03 | -8.78 | -10.13 |
| RWF | -8.38 | -9.95 |

D03 and RWF show the strongest binding affinity; D01 is markedly weaker.

**Boltz-2 predicted pIC50** (mean over 3 replicates):

| Ligand | Mean pIC50 |
|--------|-----------|
| D01 | 4.67 |
| D02 | 5.25 |
| D03 | 5.75 |

**Ligand stability** (RMSD over MD trajectory, A):

| Ligand | Mean RMSD | Max RMSD |
|--------|----------|---------|
| D01 | 0.15 | 0.27 |
| D02 | 0.64 | 1.11 |
| D03 | 0.62 | 1.18 |
| RWF | 1.15 | 1.53 |

D01 is the most stable in the binding pocket; RWF shows the most conformational flexibility.

**Protein backbone RMSD** (Ca, mean over trajectory): D01 2.77 A, D02 2.40 A, D03 2.93 A, RWF 2.51 A -- all consistent with a stable protein fold.

**Hydrogen bonds** (mean total per frame): D01 1.74, D02 1.60, D03 2.32, RWF 1.91. D03 forms the most persistent H-bond contacts with HSA.

## Notes

- Large trajectory files (`.dcd`, `.rst`, `.chk`, `.prmtop`, `.inpcrd`) are excluded from this repository. Contact the authors for trajectory data.
- Binary paths (P2Rank, UniDock, GNINA) are hardcoded for the HPC environment used. Update these before running locally.

## Citation

> *To be added upon publication.*
