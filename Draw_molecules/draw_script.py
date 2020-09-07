import sys
import svgwrite
import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole

IPythonConsole.ipython_useSVG = True

import matplotlib.pyplot as plt


def draw_2_molecules(smiles1, smiles2, legend1, legend2, output_name):
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    plot = Draw.MolsToGridImage(
        [mol1, mol2],
        # mols_list,
        molsPerRow=2,
        subImgSize=(500, 300),
        legends=[legend1, legend2],
        useSVG=False,
    )
    plot.save(output_name + ".jpg")
    plot.show()
