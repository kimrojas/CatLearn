#!/bin/bash


# OMP_NUM_THREADS=1 NP=1 mpirun -np 1 python catlearn-example-01.py
# OMP_NUM_THREADS=1 NP=2 mpirun -np 2 python catlearn-example-01.py
# OMP_NUM_THREADS=1 NP=3 mpirun -np 3 python catlearn-example-01.py
# OMP_NUM_THREADS=1 NP=4 mpirun -np 4 python catlearn-example-01.py
# OMP_NUM_THREADS=1 NP=5 mpirun -np 5 python catlearn-example-01.py
# OMP_NUM_THREADS=1 NP=6 mpirun -np 6 python catlearn-example-01.py
# OMP_NUM_THREADS=1 NP=7 mpirun -np 7 python catlearn-example-01.py
# OMP_NUM_THREADS=1 NP=8 mpirun -np 8 python catlearn-example-01.py


OMP_NUM_THREADS=1 python catlearn-example-01.py
OMP_NUM_THREADS=2 python catlearn-example-01.py
OMP_NUM_THREADS=3 python catlearn-example-01.py
OMP_NUM_THREADS=4 python catlearn-example-01.py
OMP_NUM_THREADS=5 python catlearn-example-01.py
OMP_NUM_THREADS=6 python catlearn-example-01.py
OMP_NUM_THREADS=7 python catlearn-example-01.py
OMP_NUM_THREADS=8 python catlearn-example-01.py