import os
import subprocess

job_name = 'gls-nan'
ntasks = 1

job_script=f"""#!/bin/bash
#SBATCH --account=c1601279
#SBATCH --job-name="{job_name}"
#SBATCH --constraint=GENOA
#SBATCH --nodes=1
#SBATCH --threads-per-core=1
#SBATCH --ntasks-per-node={ntasks}
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=01:30:00
#SBATCH --output="/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/scripts/job/{job_name}.out"
 
echo Bonjour du noeud `hostname`

source /lus/home/PERSO/grp_rguillermin/rguillermin/.bashrc
source /lus/home/CT1/c1601279/rguillermin/python_environment/bin/activate

srun --ntasks-per-node={ntasks} --cpus-per-task=1 --threads-per-core=1 python /lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/scripts/ensembles/remove_NaN.py
"""

filename=f"job-{job_name}.sh"
filepath = os.path.join('/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/scripts/job', filename)
with open(filepath, "w") as f:
    f.write(job_script)
    
os.chmod(filepath, 0o755)


try:
    result = subprocess.run(["sbatch", filepath], capture_output=True, text=True, check=True)
    print(f'{filename} submitted.')
    print(result.stdout.strip())
except subprocess.CalledProcessError as e:
    print('Error submitting job:')
    print(e.stderr)