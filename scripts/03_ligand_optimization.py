# =============================
# IMPORTS
# =============================
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from openbabel import pybel
import os
import torch
import torchani
from ase import Atoms, io
from ase.optimize import BFGS
from torchani.units import HARTREE_TO_KCALMOL
from rdkit.Geometry import Point3D

workDir = os.getcwd()
Type = "smiles"   # smiles, pdb or mol
smiles_or_filename = "C1C=CC2C(C(=CC(=O)C=2C=1)O)=O"

# =============================
# PART 1 — GENERATE ligand.mol AND ligand.pdb
# =============================

print("\n=== Generating initial ligand structure ===")

if Type == "smiles":
    Smiles = smiles_or_filename
    smiles_fig = Chem.MolFromSmiles(Smiles)
    hmol = Chem.AddHs(smiles_fig)
    AllChem.EmbedMolecule(hmol, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(hmol)
    AllChem.MolToMolFile(hmol, os.path.join(workDir, "ligand.mol"))
    AllChem.MolToPDBFile(hmol, os.path.join(workDir, "ligand.pdb"))

elif Type == "pdb":
    pdb_name = os.path.join(workDir, smiles_or_filename)
    mol = [m for m in pybel.readfile(filename=pdb_name, format='pdb')][0]
    out = pybel.Outputfile(filename="mol_tmp.mol", format="mol", overwrite=True)
    out.write(mol)
    out.close()

    mol_rd = Chem.MolFromMolFile("mol_tmp.mol")
    Smiles = Chem.MolToSmiles(mol_rd)

    smiles_fig = Chem.MolFromSmiles(Smiles)
    hmol = Chem.AddHs(smiles_fig)
    AllChem.EmbedMolecule(hmol)
    AllChem.MMFFOptimizeMolecule(hmol)
    AllChem.MolToMolFile(hmol, os.path.join(workDir, "ligand.mol"))
    AllChem.MolToPDBFile(hmol, os.path.join(workDir, "ligand.pdb"))

else:  # MOL input
    mol_name = os.path.join(workDir, smiles_or_filename)
    mol_rd = Chem.MolFromMolFile(mol_name)
    Smiles = Chem.MolToSmiles(mol_rd)

    smiles_fig = Chem.MolFromSmiles(Smiles)
    hmol = Chem.AddHs(smiles_fig)
    AllChem.EmbedMolecule(hmol)
    AllChem.MMFFOptimizeMolecule(hmol)
    AllChem.MolToMolFile(hmol, os.path.join(workDir, "ligand.mol"))
    AllChem.MolToPDBFile(hmol, os.path.join(workDir, "ligand.pdb"))

print("Detected SMILES:", Smiles)

# Save 2D structure image
Draw.MolToFile(smiles_fig, os.path.join(workDir, "ligand_2d_preview.png"))
img = mpimg.imread(os.path.join(workDir, "ligand_2d_preview.png"))
plt.imshow(img)
plt.axis("off")
plt.show()

# =============================
# PART 2 — GEOMETRY OPTIMIZATION WITH TORCHANI (ANI-2x)
# =============================

print("\n=== GEOMETRY OPTIMIZATION — ANI-2x ===")

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = torchani.models.ANI2x(periodic_table_index=True).to(device)
calculator = torchani.models.ANI2x().ase()

def mol2arr(mol, device=device):
    """Convert an RDKit Mol to TorchANI tensors."""
    pos = mol.GetConformer().GetPositions().tolist()
    atomnums = [a.GetAtomicNum() for a in mol.GetAtoms()]
    coordinates = torch.tensor([pos], requires_grad=True, device=device)
    species = torch.tensor([atomnums], device=device)
    return coordinates, species

# Load the ligand.mol file
mol_rd2 = Chem.MolFromMolFile(os.path.join(workDir, "ligand.mol"), removeHs=False)

# Load into ASE
atoms = io.read(os.path.join(workDir, "ligand.mol"))
atoms.center(vacuum=3.0)
atoms.set_calculator(calculator)

# Run geometry optimization
opt = BFGS(atoms)
opt.run(fmax=0.0001)

io.write(os.path.join(workDir, "ligand_min.xyz"), atoms)

# Update RDKit coordinates from the optimized XYZ
xyz = []
with open(os.path.join(workDir, "ligand_min.xyz")) as f:
    lines = f.readlines()[2:]
    for line in lines:
        parts = line.split()
        xyz.append([float(parts[1]), float(parts[2]), float(parts[3])])

conf = mol_rd2.GetConformer()
for i in range(mol_rd2.GetNumAtoms()):
    x, y, z = xyz[i]
    conf.SetAtomPosition(i, Point3D(x, y, z))

# Save final formats
AllChem.MolToMolFile(mol_rd2, os.path.join(workDir, "ligand_min.mol"))
AllChem.MolToPDBFile(mol_rd2, os.path.join(workDir, "ligand_min.pdb"))

# Convert to PDBQT
os.system(
    "obabel -i mol ligand_min.mol -o pdbqt -O ligand_min.pdbqt -xh --partialcharge"
)

# Report final energy
coords, species = mol2arr(mol_rd2, device)
energy = model((species, coords)).energies

print("Energy (hartree):", energy.item())
print("Energy (kcal/mol):", energy.item() * HARTREE_TO_KCALMOL)

print("\nGenerated files:")
print(" - ligand.mol")
print(" - ligand.pdb")
print(" - ligand_min.xyz")
print(" - ligand_min.mol")
print(" - ligand_min.pdb")
print(" - ligand_min.pdbqt")
print(" - ligand_2d_preview.png")
