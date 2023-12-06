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

Two example scripts are given in the example folder.


