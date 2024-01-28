This repository implements the equations of self-alignment found in several active matter papers such as Lam et al. 2015.

These equations are useful to simulate vibrated polar objects such as Kilobots, Hexbugs, Pogobots and other damped bristle-bot like robots.

On top of the self-alignment property, a differential formulation of distributed learning is implemented, following the equations described in the Wiki.

# Installation

Clone this repository along with the "stable" branch of the Lammps repository

```
git clone https://github.com/JeremyF-141592/spp_lammps.git
git clone -b stable https://github.com/lammps/lammps.git
```

Customization of Atoms in Lammps have changed while this repository was under development.
It is advised to revert to a previous commit :

```
cd lammps
git reset --hard d618b0ffc05dfd86915c0d148c0a72fba995eba4
```

Copy/paste the content of both the spp_lammps/src and spp_lammps/src/BROWNIAN directories into the src directory of Lammps.

From the parent directory of both repositories :

```
cp ./spp_lammps/src/BROWNIAN/* ./lammps/src/
cp ./spp_lammps/src/* ./lammps/src/
```

From the lammps/src directory, compile the serial or parallel executable :

```
cd lammps/src
make serial
make mpi
```

This yields two executable files, `lmp_serial` and `lmp_mpi`, allowing to execute lammps scripts.

# Quick start

Two example scripts are given in the example folder. All parameters, regarding self-alignment and learning, can be modified through script variable defined on top of the example scripts. The lit region of the experiments is defined with the lammps `region` command, as a cylinder in both examples. Once all parameters are set, the script can be passed to the executable, either in serial or parallel :

```
lmp_serial -in in.myscript.run
mpiexec -n NumberOfProcesses lmp_mpi -in in.myscript.run
```

A result folder should be created to store the resulting data.
In order to execute the N=256 preset example with 8 processes, starting from the parent directory of spp_lammps and lammps, this would yield :

```
mkdir ExamplePreset
cd ExamplePreset
mpiexec -n 8 ../lammps/src/lmp_mpi -in ../spp_lammps/examples/in.phototaxis256.run
```

# Notes
By default, the light intensity received in the dark is set to 0.25 and 0.75 in the light.

The type of controller can be chosen from a set of pre-programmed controllers shown on Controllers_annotated.png in this repository. However, any controller can be programmed following the procedure described in the Wiki.

The name of the dump file is set through lammps scripts. In the examples, the following format for each file is used :
evo2D_${N}_${eta}_${D}_${tau_v}_${tau_n}_${J}_${epsilon}_${alpha}_${alphaq}_${comm_radius}_${controller}_${Seed}.cfg


