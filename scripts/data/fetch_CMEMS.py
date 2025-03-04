import copernicusmarine
#copernicusmarine.login() # Nécessaire pour se logger sur le serveur de Copernicus Marine
from pprint import pprint

#output_directory = '/Users/remyguillermin/Programmation/Stage/IGE-Stochastic/data/RAW/'
output_directory = '/home/guilremy/IGE-Stochastic/data/RAW/'

dataset = 'METOFFICE-GLO-SST-L4-REP-OBS-SST' # SST
dataset = 'SEALEVEL_GLO_PHY_L4_MY_008_047' # SSH

get_result = copernicusmarine.get(
    dataset_id=dataset, # ID du dataset sur https://data.marine.copernicus.eu/products
    dry_run=False, # Si True, affiche les fichiers à télécharger sans les télécharger
    filter="*2018/*", # Filtre les fichiers à télécharger
    output_directory=output_directory, # Répertoire de sortie
    #create_file_list= output_directory + "files_to_download.txt", # Crée un fichier avec la liste des fichiers à télécharger
    no_directories=True, # Si True, télécharge les fichiers dans le répertoire de sortie sans créer de sous-répertoires
)

pprint(f"List of saved files: {get_result}")