import os
import re
from datetime import datetime

def rename_files(directory):
    pattern = re.compile(r"(\d{8})(\d{6})-UKMO-L4_GHRSST-SSTfnd-OSTIA-GLOB_REP-.*")
    
    for filename in os.listdir(directory):
        match = pattern.match(filename)
        if match:
            date_str = match.group(1)  # Extract YYYYMMDD
            date_obj = datetime.strptime(date_str, "%Y%m%d")
            new_filename = f"{date_obj.strftime('%Y-%m-%d')}-METOFFICE-GLO-SST-L4-REP-OBS-SST"
            old_path = os.path.join(directory, filename)
            new_path = os.path.join(directory, new_filename)
            os.rename(old_path, new_path)
            print(f"Renamed: {filename} -> {new_filename}")

# directory = "/Users/remyguillermin/Programmation/Stage/IGE-Stochastic/data/RAW" 
directory = '/home/guilremy/IGE-Stochastic/data/RAW/'
rename_files(directory)