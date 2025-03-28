# REGRIDING

On ADASTRA, one should manually install `ESMFpy` by downloading and compiling the source code from [earthsystemmodeling.org](https://earthsystemmodeling.org/download/)
```bash
cd esmf-8.8.0/
make info
```

this should output
```
makefile:13: *** Environment variable ESMF_DIR needs to be set to the top ESMF directory. Please see the README file for examples of setting ESMF_DIR.  Stop.
```

One should execute
```bash
echo 'export ESMF_DIR=/lus/home/CT1/c1601279/rguillermin/esmf-8.8.0' >> ~/.bashrc
module purge 
module load cray-libsci_acc/24.07.0
module load intel/2022.1.0
make
```

This will build `ESMF`. Those two modules works together but maybe other combination would work too.

One should check the build with `make all_tests`.

When `EMSF` has been built and checked, one can proceeds to the installation of `esmpy`

```bash
echo 'export ESMFMKFILE=/lus/home/CT1/c1601279/rguillermin/esmf-8.8.0/lib/libO/Linux.gfortran.64.mpiuni.default/esmf.mk' >> ~/.bashrc
cd src/addon/esmpy
python3 -m pip install .
source ~/.bashrc
```

And now one can import `esmpy` with 
```python
import esmpy
```

