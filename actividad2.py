#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 21:49:25 2024

@author: acao
"""
# 1 Obtención de información sobre proteínas y compuestos químicos
"""
A. Crea una función que acepte como entrada una lista de ID del Protein data_1tup Bank y
descargue de esta base de datos los archivos formato .cif de cada uno de estos ID.
Incluye una forma de controlar cuando un archivo ha podido ser descargado correctamente y
cuando no ha sido posible sin que ello afecte a la ejecución del programa. 
"""

def descargar_cifs(lista_ids):
    # Importamos la librería
    import requests
    
    # URL base para las consultas
    base_url = "https://www.rcsb.org/pdb/files/"
        
    # Creamos el bucle para iterar por id
    for id in lista_ids:
        response = requests.get(f"{base_url}{id}.cif")
        # Comprobamos la solicitud
        if response.status_code == 200:
            #Abrimos el .cif y lo leemos por chunks
            with open(f"{id}.cif", "wb") as archivo:
                for chunk in response.iter_content():
                    archivo.write(chunk) # Lo guardamos en un archivo
            # Imprimimos mensaje de confirmación
            print(f"Archivo {id}.cif descargado correctamente")
        
        else:
            # En caso que no se descargue correctamente
            print(f"Error tipo ({response.status_code}) al descargar {id}.cif")

# Ejemplo de uso
lista_ids_ejemplo = ['1tup', '2xyz', '3def', '4ogq', '5jkl', '6mno', '7pqr', '8stu', '9vwx', '10yza']
descargar_cifs(lista_ids_ejemplo)




"""B. Toma el ID 1tup y determina su identificador en UniProt a través de su API.
Una vez obtengas el identificador, consulta la información en UniProt,
investiga el objeto creado por la llamada y extrae la siguiente información:
    
    la fecha de publicación de la entrada,
    la fecha de la última modificación,
    si está manualmente revisada (Swiss-Prot) o si no ha sido revisada (Trembl),
    el nombre del gen y sus sinónimos,
    el organismo al cual pertenece,
    el nombre completo de la proteína,
    la secuencia de aminoácidos de una sola letra y los pdb ids asociados a dicha entrada de UniProt.

Guarda esta información en un data_1tupFrame de Pandas que contenga las siguientes columnas: 
['Uniprot_id', 'Fecha_publicacion', 'Fecha_modificacion', 'Revisado', 'Nombre_del_gen',
'Sinónimos', 'Organismo', 'PDB_ids'] (2,5)"""

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
    
    # Iteramos por cada resultado en results
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




""""C. Investiga en esta entrada de UniProt de la 1tup si hay información sobre algún cofactor.
Indica dónde se identifica al cofactor y utiliza la API de PubChem para extraer información sobre:
    
    su identificador de compuesto en Pubchem (cid),
    su peso molecular exacto,
    su inchi,
    inchikey y
    sobre sus nombres según la iupac.
    
Guarda esta información en un data_1tupFrame que contenga las siguientes columnas:
['Compuesto', 'Pubchem_id', 'Peso_molecular', 'Inchi', 'Inchikey', 'Iupac_name'] (1,5)."""

import pandas as pd
import requests

# El cofactor es:
cofactor = data_1tup.get("results")[5].get("comments")[1].get("cofactors")[0].get("name")
chebi_id = data_1tup.get("results")[5].get("comments")[1].get("cofactors")[0].get("cofactorCrossReference").get("id")

# Al buscarlo en PubChem encontramos que su CID es 32051
url_cofactor = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/32051/JSON"

# Creamos el df para el cofactor
df_cofactor = pd.DataFrame()

# Bajamos la información
response_cofactor = requests.get(url_cofactor)

if response_cofactor.status_code == 200:
    data_cofactor = response_cofactor.json()
    
    # Recolectamos la información solicitada
    new_row = {
        "Compuesto" : "CHEBI:29105 - zinc(2+)",

        "Pubchem_id" : data_cofactor.get("Record").get("RecordNumber"),
        
        "Peso_molecular" : data_cofactor.get("Record").get("Section")[2].get("Section")[0]
                                           .get("Section")[0].get("Information")[0].get("Value")
                                           .get("StringWithMarkup")[0].get("String"),
        
        "Inchi" : data_cofactor.get("Record").get("Section")[1].get("Section")[1].get("Section")[1]
            .get("Information")[0].get("Value").get("StringWithMarkup")[0].get("String"),
        
        "Inchikey" : data_cofactor.get("Record").get("Section")[1].get("Section")[1].get("Section")[2]
            .get("Information")[0].get("Value").get("StringWithMarkup")[0].get("String"),
            
        "Iupac_name" : data_cofactor.get("Record").get("Section")[1].get("Section")[1].get("Section")[0]
            .get("Information")[0].get("Value").get("StringWithMarkup")[0].get("String")
        }
    
    # Agregamos al df
    df_cofactor = pd.concat([df_cofactor, pd.DataFrame(new_row, index=[0])], ignore_index=True)

else:
    print(f"Error tipo {response_cofactor.status_code}")



"""Biopython y rdkit son librerías muy útiles para manipular datos de origen biológico y químico, respectivamente.
Utilízalas para completar estas tareas sobre el archivo 4ogq.cif: """

# Reutilizamos la función creada en el punto 1
descargar_cifs(["4ogq"])

"""A. Parséalo con MMCIFParser() y guarda en una lista todas las heteromoléculas (sin incluir las aguas). (0,5)"""

from Bio.PDB.MMCIFParser import MMCIFParser
parser = MMCIFParser(QUIET=True)
structure = parser.get_structure('4ogq', "./4ogq.cif")

heteromoleculas = []
for modelo in structure:
    for cadena in modelo:
        for residuo in cadena:
            if residuo.get_resname() != 'HOH': #Todas las heteromoleculas distintas a agua
                heteromoleculas.append(residuo.id)

print(f"Hay {len(heteromoleculas)} heteromoleculas en la estructura (sin incluir las aguas)")



"""B. Parseálo con MMCIF2Dict() en un diccionario, extrae la información sobre la clave '_pdbx_entity_nonpoly' 
y guarda en un DataFrame la información sobre el nombre de cada una de las heteromoléculas,
así como su identificador de tres letras. (0,5)"""

from Bio.PDB.MMCIFParser import MMCIF2Dict
import pandas as pd

dic_4ogq = MMCIF2Dict("./4ogq.cif")

df_4ogp = pd.DataFrame({"entitty_id" : dic_4ogq.get("_pdbx_entity_nonpoly.entity_id"),
                        "entity_name" : dic_4ogq.get("_pdbx_entity_nonpoly.name"),
                        "comp_id" : dic_4ogq.get("_pdbx_entity_nonpoly.comp_id")})




"""C. Toma el nombre de cada heteromolécula y mediante la API de Pubchem consigue su SMILES,
para los casos que sea posible.
Transforma estos SMILES con rdkit en SDF y añade ambas columnas al DataFrame del Apartado 2B. (1)"""

import pubchempy as pcp
from rdkit import Chem    

for index, row in df_4ogp.iterrows():
    nombre = row["entity_name"]

    try:
        # Buscamos en pubchem por el nombre de cada compuesto
        compound = pcp.get_compounds(nombre, "name")

        if compound:
            # Tomamos el primer hallazgo de cada busqueda y buscamos los smiles
            smiles = compound[0].isomeric_smiles
            # Agregamos al df
            df_4ogp.loc[index, "smiles"] = smiles
            
            # Convertimos a sdf
            mol_smiles = Chem.MolFromSmiles(smiles)
            sdf = Chem.MolToMolBlock(mol_smiles)
            # Agregamos al df
            df_4ogp.loc[index, "SDF"] = sdf
            
    except:
        print(f"Error en: {nombre}")


"""D. Genera un archivo .sdf en el cual guardes todas las moléculas SDF.
Añade para cada una de ellas un campo que incluya información sobre su 'Molecular_weight'.
Ayúdate de rdkit para calcular dicho campo y asegúrate que el archivo que has creado
puede utilizarse para crear un objeto mol de rdkit con cada una de las moléculas. (1,5)"""

"Pendiente pero ya no me da tiempo"
