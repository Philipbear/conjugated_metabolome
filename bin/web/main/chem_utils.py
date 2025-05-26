from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from molmass import Formula
import re
import pubchempy as pcp
import urllib.parse
import requests


def inchikey_to_common_name(inchikey):
    try:
        compound = pcp.get_compounds(inchikey, 'inchikey')
        if compound:
            if compound[0].synonyms:
                name_lst = compound[0].synonyms
                return name_lst
    except:
        return None
    return None


def smiles_to_formula_inchikey(smiles):
    try:
        # Convert SMILES to a molecule object
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, None
        # Get the molecular formula
        formula = rdMolDescriptors.CalcMolFormula(mol)
        # Get the InChIKey
        inchikey = Chem.MolToInchiKey(mol)
        return formula, inchikey
    except:
        return None, None


def calc_monoisotopic_mass(formula):
    """
    Calculate the exact mass for a given formula string
    """
    try:
        f = Formula(formula)
        return f.monoisotopic_mass
    except:
        return None


def neutralize_formula(formula):
    """
    deal with the charge in the formula
    such as C5H5N+ -> C5H4N
    """
    if not formula:
        return formula

    # Split the formula into the chemical part and the charge part
    match = re.match(r'([A-Za-z0-9]+)([-+]?\d*)', formula)
    if not match:
        return None # Invalid formula

    chemical, charge = match.groups()

    # If there's no charge, return the original formula
    if not charge:
        return chemical

    # Convert charge to integer
    if charge == '-':
        charge = -1
    elif charge == '+':
        charge = 1
    else:
        charge = int(charge) if charge else 0

    # Parse the chemical formula
    elements = re.findall(r'([A-Z][a-z]?)(\d*)', chemical)

    # Find H and its count
    h_index = next((i for i, (elem, _) in enumerate(elements) if elem == 'H'), None)
    if h_index is not None:
        h_count = int(elements[h_index][1]) if elements[h_index][1] else 1
    else:
        h_count = 0

    # Calculate new H count
    new_h_count = h_count - charge

    # Update or insert H in the elements list
    if new_h_count > 0:
        if h_index is not None:
            elements[h_index] = ('H', str(new_h_count) if new_h_count > 1 else '')
        else:
            elements.insert(1, ('H', str(new_h_count) if new_h_count > 1 else ''))
    elif h_index is not None:
        elements.pop(h_index)

    # Reconstruct the formula
    new_formula = ''.join(elem + count for elem, count in elements)

    return new_formula



def get_structure_image_pubchem(smiles):
    """
    Generate and display a structure image from a SMILES string
    using the PubChem structure image API.
    """
    # URL encode the SMILES string
    encoded_smiles = urllib.parse.quote(smiles)
    image_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{encoded_smiles}/PNG"
    
    return image_url

    
def get_structure_image_gnps2(smiles):
    """
    Generate and display a structure image from a SMILES string
    using the GNPS2 structure image API.
    
    Parameters:
    -----------
    smiles : str
        SMILES string of the chemical structure
    """
    # URL encode the SMILES string
    encoded_smiles = urllib.parse.quote(smiles)
    image_url = f"https://structure.gnps2.org/structureimg?smiles={encoded_smiles}"
    
    return image_url


def get_compound_description_pubchem(smiles):
    """
    Get compound description text from PubChem using SMILES string.
    Returns a simple string with all descriptions concatenated.
    """
    try:
        # URL encode the SMILES string
        encoded_smiles = urllib.parse.quote(smiles)
        
        # PubChem REST API endpoint for compound description
        description_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{encoded_smiles}/description/JSON"
        
        response = requests.get(description_url)
        response.raise_for_status()
        
        data = response.json()
        
        if 'InformationList' not in data or not data['InformationList']['Information']:
            return None
        
        # Collect all descriptions
        descriptions = []
        for info in data['InformationList']['Information']:
            if 'Description' in info:
                descriptions.append(info['Description'])
        
        # Return all descriptions joined together
        return ' '.join(descriptions) if descriptions else None
        
    except:
        return None
    
    
if __name__ == "__main__":
    print(calc_monoisotopic_mass("C5H5N"))