# ==========================================
# Compute RMSD of the ligand (in Å)
# Reference: frame 0 of the no-water trajectory (internal fit to protein Cα)
# ==========================================
import MDAnalysis as mda
import MDAnalysis.analysis.rms
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# === INPUT FILES ===
pdb_ref = "prot_lig_prod_nw.pdb"         # reference structure
traj = "prot_lig_prod_1-5_nw.dcd"        # trajectory
workDir = os.getcwd()
Output_name = "rmsd_ligand"              # output file base name

# === LOAD SYSTEM ===
u = mda.Universe(pdb_ref, traj)
ref = mda.Universe(pdb_ref, pdb_ref)
ref.trajectory[0]  # set reference frame

# === COMPUTE RMSD ===
R = MDAnalysis.analysis.rms.RMSD(
    u, ref,
    select="resname UNK and not (name H*)",
    groupselections=["resname UNK and not (name H*)"]
)
R.run()

# === EXTRACT DATA ===
rmsd = R.rmsd.T    # transpose for easy plotting
# 500 frames for 50 ns → 0.1 ns/frame
time_array = np.arange(0, 50.0, 0.1)

# === PLOT RMSD ===
plt.figure(figsize=(6,4))
plt.plot(time_array, rmsd[2], alpha=0.6, color='blue', linewidth=1.0)
plt.xlim(0, 50)
plt.xlabel("Time (ns)", fontsize=14, fontweight='bold')
plt.ylabel("RMSD (Å)", fontsize=14, fontweight='bold')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.savefig(os.path.join(workDir, Output_name + ".png"), dpi=600, bbox_inches='tight')
plt.close()

# === SAVE RAW DATA ===
raw_data = pd.DataFrame({
    "Time (ns)": time_array,
    "RMSD (Å)": rmsd[2]
})
raw_data.to_csv(os.path.join(workDir, Output_name + ".csv"), index=False)

print("✅ RMSD computed successfully!")
print(f"Saved: {Output_name}.png and {Output_name}.csv in {workDir}")
