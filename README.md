# HSA – Naphthoquinones Docking and Molecular Dynamics

Computational pipeline for the study of naphthoquinone derivatives (D01, D02, D03, RWF) as potential binders of Human Serum Albumin (HSA). This repository contains all analysis scripts used in the associated publication.

## Workflow Overview

```
01_prepare_receptor          → Clean PDB, fix missing atoms/H, convert to PDBQT
02_binding_site_prediction   → P2Rank binding pocket prediction
03_residue_selection         → Extract pocket residues as individual PDB files
04_ligand_optimization       → ANI-2x geometry optimization (Neural Network Potential)
05_docking                   → GPU-accelerated docking with UniDock (Vina scoring)
06_build_topology            → AMBER19SB + GAFF2.11 + TIP3P topology (OpenMM)
07_interaction_fingerprint   → Protein-ligand interaction fingerprint (ProLIF)
08_equilibration             → NVT + NPT equilibration (0.5 ns + 0.5 ns, OpenMM)
09_production_md             → NPT production MD — 11 × 5 ns = 55 ns (CUDA)
10_concatenate_trajectories  → PBC unwrap, complex centering, Cα alignment
11_gnina_rescoring           → CNN-based rescoring of MD frames with GNINA
12_rmsd_ligand               → Ligand RMSD (protein Cα fit, reference = frame 0)
13_rmsd_ligand_nofit         → Ligand RMSD (no local fit, reference = frame 100)
14_rmsd_docking_pose         → Ligand RMSD vs initial docking pose
15_rmsd_protein_backbone     → Protein Cα RMSD over production trajectory
16_hbonds                    → Protein–ligand hydrogen bond counts per frame
```

The pipeline was run independently for each of the four ligands (D01, D02, D03, RWF).

## Software & Force Fields

| Tool | Purpose |
|------|---------|
| OpenMM | MD engine |
| AMBER ff19SB | Protein force field |
| GAFF2 (v2.11) | Ligand force field |
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
| `rmsd_ligand.csv/png` | Ligand RMSD (protein Cα fit, ref = frame 0) |
| `rmsd_ligand_nofit.csv/png` | Ligand RMSD (no local fit, ref = frame 100) |
| `rmsd_docking_pose.csv/png` | Ligand RMSD vs initial docking pose |
| `rmsd_protein_backbone.csv/png` | Protein Cα RMSD over trajectory |
| `hbonds_counts.csv` | H-bond counts per frame (PROT→UNK and UNK→PROT) |
| `ligand_2d_preview.png` | 2D structure of the ligand (from SMILES) |

## Boltz-2 Predicted Poses

Boltz-2 (deep learning structure prediction, via Rowan Scientific) was used to predict binding poses and estimate binding affinity (IC50 and Kd) for three naphthoquinone derivatives (D01, D02, D03) at site IB of HSA. The predicted affinities were compared against the docking and MD results as an independent validation. Three poses were generated per compound (total: 9 structures).

| File | Description |
|------|-------------|
| `boltz2/D01-1_boltz2.pdb` | Boltz-2 pose 1 — D01 |
| `boltz2/D01-2_boltz2.pdb` | Boltz-2 pose 2 — D01 |
| `boltz2/D01-3_boltz2.pdb` | Boltz-2 pose 3 — D01 |
| `boltz2/D02-1_boltz2.pdb` | Boltz-2 pose 1 — D02 |
| `boltz2/D02-2_boltz2.pdb` | Boltz-2 pose 2 — D02 |
| `boltz2/D02-3_boltz2.pdb` | Boltz-2 pose 3 — D02 |
| `boltz2/D03-1_boltz2.pdb` | Boltz-2 pose 1 — D03 |
| `boltz2/D03-2_boltz2.pdb` | Boltz-2 pose 2 — D03 |
| `boltz2/D03-3_boltz2.pdb` | Boltz-2 pose 3 — D03 |
| `boltz2/boltz2_results.csv` | Boltz-2 scoring results summary |

## Key Results

**GNINA rescoring** (mean minimized affinity over MD frames, kcal/mol):

| Ligand | Mean affinity | Best frame |
|--------|--------------|------------|
| D01 | −5.82 | −7.42 |
| D02 | −8.04 | −9.43 |
| D03 | −8.78 | −10.13 |
| RWF | −8.38 | −9.95 |

D03 and RWF show the strongest 