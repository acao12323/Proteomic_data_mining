# Importamos las librerías a usar
import pandas as pd
import requests

# Creamos el df_1tup
df_1tup = pd.DataFrame()

# Construye la URL sobre la cual se va a llevar a cabo la consulta
url = f'https://rest.uniprot.org/uniprotkb/search?query=1tup'

# Realiza una solicitud GET para obtener los datos de la proteína
response_1tup = requests.get(url)

# bajamos los datos
if response_1tup.status_code == 200:
    data_1tup = response_1tup.json()
    
    # Los agregamos al df_1tup
    for result in data_1tup["results"]:
        # Dado que result[3] no tiene synonyms crea un error
        # Por lo tanto definimos la consulta de synonyms
        synonyms = result.get("genes", [{}])[0].get("synonyms", [])
        
        new_row = {
            'Uniprot_id' : result["uniProtkbId"],
            'Fecha_publicacion' : result["entryAudit"]["firstPublicDate"],
            'Fecha_modificacion' : result["entryAudit"]["lastAnnotationUpdateDate"],
            'Revisado' : ["True" if "Swiss-Prot" in result["entryType"] else "False"],
            'Nombre_del_gen' : result["genes"][0]["geneName"]["value"],
            'Sinónimos': [synonyms if synonyms else None], # Usamos list comprehension para lidiar con el error
            'Organismo' : result["organism"]["scientificName"],
            # Agregamos la columna Protein_name como el nombre completo de la proteína que se pide en el enuncido
            'Protein_name' : result["proteinDescription"]["recommendedName"]["fullName"]["value"],
            'PDB_ids' : result["primaryAccession"],
            }
        df_1tup = pd.concat([df_1tup, pd.DataFrame(new_row, index=[0])], ignore_index=True)
else:
    print(f"Error tipo {response_1tup.status_code}")

df_1tup


import rdkit
from rdkit import Chem

smiles = "O=C(O)C1=CC=CC=C1C(=O)O"
mol = Chem.MolFromSmiles(smiles)
sdf_string = Chem.MolToMolFile(mol, kekuleSmiles=True)

print(sdf_string)  # Print the SDF string to the console

with open("output.sdf", "w") as f:
    f.write(sdf_string)  # Save the SDF string to a file
