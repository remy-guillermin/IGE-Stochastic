import copernicusmarine
#copernicusmarine.login() # Nécessaire pour se logger sur le serveur de Copernicus Marine
from pprint import pprint

#output_directory = '/Users/remyguillermin/Programmation/Stage/IGE-Stochastic/data/RAW/'
output_directory = '/home/guilremy/IGE-Stochastic/data/RAW/SLA_CMEMS'

#dataset = 'METOFFICE-GLO-SST-L4-REP-OBS-SST' # SST
dataset = 'cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.125deg_P1D' # SLA
#dataset = 'cmems_obs-mob_glo_phy-sss_my_multi_P1D' # SSS

get_result = copernicusmarine.get(
    dataset_id=dataset, # ID du dataset sur https://data.marine.copernicus.eu/products
    dry_run=False, # Si True, affiche les fichiers à télécharger sans les télécharger
    filter="*2019/*", # Filtre les fichiers à télécharger
    output_directory=output_directory, # Répertoire de sortie
    #create_file_list= output_directory + "files_to_download.txt", # Crée un fichier avec la liste des fichiers à télécharger
    no_directories=True, # Si True, télécharge les fichiers dans le répertoire de sortie sans créer de sous-répertoires
)

pprint(f"List of saved files: {get_result}")

copernicusmarine.subset(
  dataset_id="cmems_mod_glo_phy_my_0.083deg_P1D-m",
  variables=["uo", "vo"],
  minimum_longitude=20,
  maximum_longitude=70,
  minimum_latitude=-40,
  maximum_latitude=10,
  start_datetime="2017-01-01T00:00:00",
  end_datetime="2021-06-30T00:00:00",
  minimum_depth=0,
  maximum_depth=1,
  output_filename = "CMEMS_Indian_currents_surface.nc",
  output_directory = "copernicus-data"
)