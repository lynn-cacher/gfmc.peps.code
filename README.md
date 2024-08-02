#PEPS_GFMC
1, Windows or Linux system for X-86/64.

2, An MPI implementation, such as OpenMPI or MPICH.

3, A C++ compiler.

4, A Python interpreter.

5, PyTorch platform.

6, The main file is located in PEPS_GFMC\special\main.cpp.

7, The Makefile is located in PEPS_GFMC Depending on the C++ compiler being used, you may need to modify the "CC" and "CPPFLAGS" in the Makefile accordingly.

8, The executable program will appear in PEPS_GFMC\console\ with the name "exe".

9, The command to submit the program is in PEPS_GFMC\console\slurm.bscc.t6.txt. Make necessary modifications based on the platform.

10, Specific parameters for the calculations can be modified in PEPS_GFMC\special\pars.cpp and PEPS_GFMC\peps\config.py, such as lattice size and other settings.





