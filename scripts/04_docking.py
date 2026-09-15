import os, sys
import locale
import py3Dmol

# ============================
# 0. GENERAL CONFIG
# ============================

workDir = os.getcwd()
receptor = "receptor.pdb"
receptor_pdbqt = "receptor.pdbqt"
ligand = os.path.join(workDir, "ligand_min.pdbqt")   # optimized ligand

# Docking box center and size (from P2Rank pocket)
centerX = 32.79
centerY = 13.53
centerZ = 9.61
sizeX = 15
sizeY = 15
sizeZ = 15

# ============================
# 1. UNIDOCK PARAMETERS
# ============================

seed = "0"
scoring = "vina"
search_mode = "detail"
exhaustiveness = 500
max_step = 40
num_modes = 10

locale.getpreferredencoding = lambda do_setlocale=True: "UTF-8"

# ============================
# 2. BUILD AND RUN UNIDOCK COMMAND
# ============================
# Path to UniDock executable (override with UNIDOCK_BIN env var if not on PATH)
unidock_exec = os.environ.get("UNIDOCK_BIN", "unidock")

unidock_cmd = (
    f"{unidock_exec}"
    f" --receptor {receptor_pdbqt} "
    f" --gpu_batch {ligand} "
    f" --search_mode {search_mode} "
    f" --scoring {scoring} "
    f" --center_x {centerX} --center_y {centerY} --center_z {centerZ} "
    f" --size_x {sizeX} --size_y {sizeY} --size_z {sizeZ} "
    f" --num_modes {num_modes} "
    f" --dir {workDir} "
    f" --seed {seed} "
    f" --max_step {max_step} "
    f" --exhaustiveness {exhaustiveness}"
)

# Save as shell script
with open("unidock.sh", "w") as f:
    f.write(unidock_cmd + "\n")

os.system("chmod +x unidock.sh")
print("Running UniDock...")
os.system("./unidock.sh > unidock.txt 2>&1")

# ============================
# 3. CONVERT OUTPUT
# ============================

ligand_out = os.path.join(workDir, "ligand_min_out.pdbqt")
ligand_out_sdf = os.path.join(workDir, "ligand_out.sdf")

os.system(f"obabel -i pdbqt {ligand_out} -o sdf -O {ligand_out_sdf} -xh")

# ============================
# 4. VISUALIZATION (optional)
# ============================

v = py3Dmol.view(js='https://3dmol.org/build/3Dmol.js')
v.addModel(open(receptor).read())
v.setStyle({'cartoon': {}, 'stick': {'colorscheme':'white','radius':.1}})
v.addModelsAsFrames(open(ligand_out_sdf,'rt').read())
v.setStyle({'model':1},{'stick':{'colorscheme':'greenCarbon'}})
v.zoomTo({'model':1})
v.show()
