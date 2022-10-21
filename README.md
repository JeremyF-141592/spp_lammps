# spp_lammps
Self-propelled particles simulation with auto-alignment in lammps

# How to install

Once you have lammps installed (<https://github.com/lammps/lammps>), simply copy/paste the content of the /src folder into lammps/src.
Be informed that some core files of lammps are customized in order to add the new types of particles. If you have a custom lammps and want to keep your changes, careful merging will be necessary.
**The package Brownian is also required.**. Copy/paste the files from lammps/src/BROWNIAN/ to lammps/src.

Once you copied the file, compile lammps using "make serial" or "make mpi" (parallel executable).


# Quick start

You will find two quickstart scripts in the /examples folder. Simply modify the parameters in these scripts to your liking and begin the simulation by using the compiled lammps executable, ex. :

    cd lammps/src
    make serial
    ./lmp_serial -in ../../spp_lammps/examples/in.quickstart.run

You can then use a visualization software such as Ovito (<https://www.ovito.org/>).
