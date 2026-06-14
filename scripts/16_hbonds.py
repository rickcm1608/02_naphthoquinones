import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
import pandas as pd

# ===============================
# LOAD SYSTEM AND TRAJECTORY
# ===============================
u = mda.Universe("system_amber.prmtop", "prot_lig_prod_1-11_whole.dcd")

# ===============================
# ATOM SELECTIONS
# ===============================
protein_sel = "protein"
ligand_sel  = "resname UNK"

# Both protein and ligand can act as donors and acceptors
donors    = f"({protein_sel}) or ({ligand_sel})"
acceptors = f"({protein_sel}) or ({ligand_sel})"

# ===============================
# HYDROGEN BOND ANALYSIS
# ===============================
hbond = HydrogenBondAnalysis(
    universe=u,
    donors_sel=donors,
    acceptors_sel=acceptors,
    d_a_cutoff=3.2,
    d_h_a_angle_cutoff=135
)

print("Analyzing protein ↔ ligand hydrogen bonds...")
hbond.run()
print("Done.")

# ===============================
# BUILD DATAFRAME
# ===============================
hb = hbond.results.hbonds

df = pd.DataFrame(
    hb,
    columns=["frame", "donor", "H", "acceptor", "dist", "angle"]
)

# ===============================
# IDENTIFY RESIDUES BY INDEX
# ===============================
atoms = u.atoms

def resname(i):
    return atoms[int(i)].resname

df["donor_res"]    = df["donor"].apply(resname)
df["acceptor_res"] = df["acceptor"].apply(resname)

# ===============================
# CLASSIFY DIRECTION: PROT → UNK vs UNK → PROT
# ===============================
df_prot_to_unk = df[
    (df["donor_res"] != "UNK") &
    (df["acceptor_res"] == "UNK")
]

df_unk_to_prot = df[
    (df["donor_res"] == "UNK") &
    (df["acceptor_res"] != "UNK")
]

# ===============================
# COUNT H-BONDS PER FRAME
# ===============================
frames = sorted(df["frame"].unique())

counts = pd.DataFrame({"frame": frames})

counts["PROT_to_UNK"] = (
    df_prot_to_unk.groupby("frame")
    .size()
    .reindex(frames, fill_value=0)
)

counts["UNK_to_PROT"] = (
    df_unk_to_prot.groupby("frame")
    .size()
    .reindex(frames, fill_value=0)
)

counts = counts.astype(int)

# ===============================
# EXPORT
# ===============================
counts.to_csv("hbonds_counts.csv", index=False)

print("Output file: hbonds_counts.csv")
print("Total PROT → UNK:", df_prot_to_unk.shape[0])
print("Total UNK → PROT:", df_unk_to_prot.shape[0])
