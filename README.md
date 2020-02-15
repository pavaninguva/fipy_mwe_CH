# To set-up

This is a conda-forge installation of FiPy 3.4 with the numpy package downgraded to 1.16.4. 

To reproduce the conda environment: 
```
conda env create -f environment.yml
```

# To run

To run the script: 
```
python mwe.py --petsc
```

To do a mpirun: 
```
mpirun -np 8 python mwe.py --petsc 
```

# Postprocessing

The output .vtk file needs to be edited since the celltypes are currently listed as `41` and should be set to `9`. 

This can be done on the commandline: 
```
sed -i 's/^41$/9/' foo.vtk
```