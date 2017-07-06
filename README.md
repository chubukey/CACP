This package consists of two part,
1. LAMMPS related source file in lammps_src folder
2. main cpfe and driver code in main_code folder, together with a folder contains input file for an example run

The installation of CACP can be done by following steps
1. download and install LAMMPS version Feb-1-2014
2. copy the source files from lammps_src into your lammps src/ folder, replace if needed
3. make using Makefile.cacp on marcc or your Makefile dependes on the server enviroment. ( Upto here you get a lammps that can run regular MD simulation but also provides necessary interface for coupled code.)
4. compile lammps as static library using " make Makefile.lib " and "make -f Makefile.lib foo" , foo is ur executable name, succesfuly build gives liblammps_foo.a file.
5. modify the Makefile in the main_code folder, change the path accordingly for your lammps source file path.
6. Before the final make, you need to make sure your compiler is intel openmpi, which can be examined by using "which mpicc", change .bashrc files accordingly.
7. Using command "Make" in the main_code folder, successful build gives executable file "fem_uel".

The usage of CACP code: Copy the executable to your input folder, use "qsub run_CACP.batch" to submit to server.

For details on how the CACP code works, and how to adjust the parameters for your simulation, refer to Introduction_to_CACP.pdf


