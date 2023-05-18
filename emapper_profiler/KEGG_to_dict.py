#import csv
import json

#cargamos todos los pathwyas con sus KOs de KEGG
KOs_KEGGs_file = open("KOs.keg_cleanned", "r")


## generamos los diccionarios 

def build_KEGG_dict_from_file(KOs_KEGGs):
    KEGG_dict = {} #almacena los Keggs pathways ID y una lista con sus KOs
    kos_dict = {}
    pathway = ""
    kegg_count = 0 # esto es un count de los ko
    for line in KOs_KEGGs:
        if line.startswith("C"):
            if pathway != "":
                KEGG_dict[pathway]= [pathway_description,kegg_count, KOS_list]
                kegg_count = 0

            KOS_list = []
            field = line.rstrip().split("    ")[1].split(" ")
            pathway = field[0]
            pathway_description = " ".join(field[1::])
        else:
            field = line.rstrip().split("      ")[1].split("  ")		
            try:
                ko = field[0]
                symbol = field[1].split("; ")[0]
                description = field[1].split("; ")[1]
                kegg_count+=1
            except:
                continue

            KOS_list.append({"KO":ko, "symbol":symbol, "description":description})	
            kos_dict[ko] = {}
            kos_dict[ko]['symbol'] = symbol
            kos_dict[ko]['description'] = description

    return KEGG_dict, kos_dict



# Dictionary to be saved
KEGG_dict, kos_dict = build_KEGG_dict_from_file(KOs_KEGGs_file)

## JSON ##
with open("KEGG_pathway_dict.txt", "w") as file:
    json.dump(KEGG_dict, file)  # encode dict into JSON

with open("KEGG_kos_dict.txt", "w") as file:
    json.dump(kos_dict, file)  # encode dict into JSON


## CSV ##
# with open("KOs_KEGGs_dict.csv", "w", newline="") as file: # Open a csv file for writing
#     # Create a writer object
#     writer = csv.DictWriter(file fieldnames=KEGG_dict.keys())

#     # Write the header row
#     writer.writeheader()

#     # Write the data rows
#     writer.writerow(KEGG_dict)




