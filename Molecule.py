from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

IPythonConsole.ipython_useSVG=True  #< set this to False if you want PNGs instead of SVGs
def mol_with_atom_index(mol):
    m1 = Chem.Mol(mol)
    for atom in m1.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    return m1

def GetCharges(m):
 m2 = Chem.Mol(m)
 AllChem.ComputeGasteigerCharges(m2)
 #m2 = Chem.Mol(m)
 for at in m2.GetAtoms():
    lbl = '%.2f'%(at.GetDoubleProp("_GasteigerCharge"))
    at.SetProp('atomNote',lbl)
 return m2

def DrawToFileWithCharges(smiles, fname='Pic.png'):
    m = GetCharges(Chem.MolFromSmiles(smiles))
    Draw.MolToFile(m, fname)


#mol = Chem.MolFromSmiles('C1=CC=CC=C1')
#Draw.MolToFile(mol, 'benzene.png')
from enum import Enum
from io import BytesIO
import base64

class Molecule():
    class DrawType(Enum):
        JUST = 0
        INDEXES = 1
        CHARGES = 2
    def GetB64(m, h=400, y=400):
        img = Draw.MolToImage(m, size=(h, y))
        buffer = BytesIO()
        img.save(buffer, format="PNG")
        b64_str = base64.b64encode(buffer.getvalue()).decode("utf-8")
        return b64_str
    def __init__(self, smiles):
        self.smiles = smiles
        self.m = Chem.MolFromSmiles(smiles)
        self.m_with_indexes = mol_with_atom_index(self.m)
        self.m_with_charges = GetCharges(self.m)

    def Get(self, Type: DrawType):
        if Type == self.DrawType.JUST:
            return self.m
        elif Type == self.DrawType.INDEXES:
            return self.m_with_indexes
        elif Type == self.DrawType.CHARGES:
            return self.m_with_charges
        else :
            raise NotImplementedError("No implemented")
    def Draw_Molecule(self, fname = 'PIC.png', Type: DrawType = DrawType.JUST):
        Draw.MolToFile(self.Get(Type), fname)
