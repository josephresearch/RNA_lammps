[![DOI](https://zenodo.org/badge/735131874.svg)](https://zenodo.org/doi/10.5281/zenodo.10428291)
- [1. Source code for the manuscript](#1-source-code-for-the-manuscript)
  - [1.1. System requirements](#11-system-requirements)
  - [1.2. Installation](#12-installation)
  - [1.3. Implementation](#13-implementation)
  - [1.4. Demo](#14-demo)
    - [1.4.1. Single chain](#141-single-chain)
    - [1.4.2. Multichain Condensation](#142-multichain-condensation)
  - [1.5. Main simulations](#15-main-simulations)
# 1. Source code for the manuscript

> **Dilimulati Aierken and Jerelle A. Joseph, *Accelerated Simulations of RNA Phase Separation: A
> Systematic Study of Non-redundant Tandem Repeats*, BioRxiv (2023). https://doi.org/10.1101/2023.12.23.573204**

We are happy to share our implementation of a coarse-grained RNA model. Please feel free to use it
and cite our paper (Preprint: https://doi.org/10.1101/2023.12.23.573204)

## 1.1. System requirements
Linux systems, C++ compilers, and MPI capability. Tested on [Della
Cluster](https://researchcomputing.princeton.edu/systems/della) at Princeton University with Intel
compiler 2019 and an OptiPlex 7050 Tower with Red Hat Enterprise Linux 8.5.

## 1.2. Installation
For more detailed instructions about how to compile LAMMPS, please refer to
https://docs.lammps.org/Build_cmake.html. Here, we show the compilation of LAMMPS with the
coarse-grained RNA model in a multi-CPU environment.


1. Download the stable version `23Jun2022_update3` of LAMMPS to the target directory.
   ```console
   $ wget https://github.com/lammps/lammps/archive/refs/tags/stable_23Jun2022_update3.tar.gz
   ```
2. Download the coarse-grained RNA force field files to the same directory:
   [pair_base_rna.h](./src_code/pair_base_rna.h) and
   [pair_base_rna.cpp](./src_code/pair_base_rna.cpp) from the folder [src_code](./src_code/).
   
3. Uncompress the downloaded package
   ```console
   $ tar zxf stable_23Jun2022_update3.tar.gz
   ```

4. Include the two files [pair_base_rna.h](./src_code/pair_base_rna.h) and
   [pair_base_rna.cpp](./src_code/pair_base_rna.cpp) into the
   `/lammps-stable_23Jun2022_update3/src/EXTRA-PAIR/` folder from the uncompressed LAMMPS source
   code.
   ```console
   $ mv pair_base_rna.h ./lammps-stable_23Jun2022_update3/src/EXTRA-PAIR/
   $ mv pair_base_rna.cpp ./lammps-stable_23Jun2022_update3/src/EXTRA-PAIR/
   ```

5. Compile LAMMPS using `cmake` with `MPI`, `EXTRA-PAIR`, and `MOLECULE` options. Here is the
   compilation options we have tested on a Linux desktop (main compiler: g++ (GCC) 8.5.0 20210514
   (Red Hat 8.5.0-18) with mpicxx):
   ```console
   $ cd lammps-stable_23Jun2022_update3
   $ mkdir build
   $ cd build
   $ cmake3 -D CMAKE_INSTALL_PREFIX=~/.local \
   > -D LAMMPS_MACHINE=rna \
   > -D ENABLE_TESTING=yes \
   > -D BUILD_MPI=yes \
   > -D BUILD_OMP=yes \
   > -D CMAKE_BUILD_TYPE=Release \
   > -D CMAKE_CXX_COMPILER=g++ \
   > -D  CMAKE_CXX_FLAGS_RELEASE="-O2 -march=native -DNDEBUG" \
   > -D PKG_EXTRA-PAIR=yes \
   > -D PKG_MOLECULE=yes \
   > -D PKG_RIGID=yes ../cmake
   $ make -j 8
   $ make install
   ```
This will create a LAMMPS executable called `lmp_rna` in the directory `~/.local/bin/`. Thus, LAMMPS
can be called by `~/.local/bin/lmp_rna`.

## 1.3. Implementation
In LAMMPS, the potential between base pairs in RNA named as `base/rna` and should be overlayed with
the excluded volume potential. In "atom" types, Adenine is defined as type `1`, Cytosine as type
`2`, Guanine as type `3`, and Uracil as type `4`. Here is how the potentials are implemented in
LAMMPS scripts:
```py
#### Excluded volume WCA-LJ potential and Base-Pair Potential
pair_style hybrid/overlay lj/cut 10.0 base/rna 18.0 # cuttoffs


#  pair_coeff for LJ, specify 2:
#    * energy
#    * sigma
pair_coeff * * lj/cut 2.0 8.9089871814
pair_modify pair lj/cut shift yes    # option to ensure energy is calculated corectly


#  pair_coeff for Base Pair Interactions, specify 9:
#    * energy Ubp: A-1, C-2, G-3, U-4
#    * harmonic K_r
#    * Equilibrium distance r_0
#    * harmonic K_theta
#    * dihedral K_phi
#    * angle theta1
#    * angle theta2
#    * dihedral phi1
#    * dihedral phi2
pair_coeff 1 1 base/rna 0.00000 3.0 13.8 1.5 0.5 1.8326 0.9425 1.8326 1.1345 #AA
pair_coeff 1 2 base/rna 0.00000 3.0 13.8 1.5 0.5 1.8326 0.9425 1.8326 1.1345 #AC
pair_coeff 1 3 base/rna 0.00000 3.0 13.8 1.5 0.5 1.8326 0.9425 1.8326 1.1345 #AG
pair_coeff 1 4 base/rna 3.33333 3.0 13.8 1.5 0.5 1.8326 0.9425 1.8326 1.1345 #AU
pair_coeff 2 2 base/rna 0.00000 3.0 13.8 1.5 0.5 1.8326 0.9425 1.8326 1.1345 #CC
pair_coeff 2 3 base/rna 5.00000 3.0 13.8 1.5 0.5 1.8326 0.9425 1.8326 1.1345 #CG
pair_coeff 2 4 base/rna 0.00000 3.0 13.8 1.5 0.5 1.8326 0.9425 1.8326 1.1345 #CU
pair_coeff 3 3 base/rna 0.00000 3.0 13.8 1.5 0.5 1.8326 0.9425 1.8326 1.1345 #GG
pair_coeff 3 4 base/rna 3.33333 3.0 13.8 1.5 0.5 1.8326 0.9425 1.8326 1.1345 #GU
pair_coeff 4 4 base/rna 0.00000 3.0 13.8 1.5 0.5 1.8326 0.9425 1.8326 1.1345 #UU
```
## 1.4. Demo
In the [Demo](./demo/) folder, we included two types of simulations: single chain and multichain.

### 1.4.1. Single chain
There are two files for the single-chain simulation. One file is
[config_single.dat](./demo/config_single.dat) which includes the initial configuration of a single
RNA chain with sequence `AGGCAGCAGCCAAAAGGCAGCAGCCA`. The other file is
[lmp_single.in](./demo/lmp_single.in) which includes the rest of the setup for simulating the chain
at 293K with the Langevin thermostat. After copying the two files into a test folder run the
following command to simulate using 1 CPU:
```console
$ mpirun -np 1 ~/.local/bin/lmp_rna -in lmp_single.in
```
### 1.4.2. Multichain Condensation
There are two files for the multi-chain simulation. One file is
[config_multi.dat](./demo/config_multi.dat) which includes the initial configuration of a condensed
droplet of (CAG)$_{47} \times 64$. The other file is [lmp_multi.in](./demo/lmp_multi.in) that
includes the rest of the setup for simulating the chain at 293K with the Langevin thermostat. After
copying the two files into a test folder run the following command to simulate using 4 CPUs:
```console
$ mpirun -np 8 ~/.local/bin/lmp_rna -in lmp_multi.in
```
## 1.5. Main simulations
In the [mainSimulations](./mainSimulations/) folder, we included three files for the simulation of
(CAG)$_{40} \times 64$, and all other main simulations in the manuscript follow the same procedure.

1. The file [config_cag40.dat](./mainSimulations/config_cag40.dat) sets up the initial configuration
by placing 64 chains of (CAG)$_{40}$ evenly in a cubic box with a concentration of $200\mu M$.

2. The file [lmp_main.in](./mainSimulations/lmp_main.in) sets up the simulation environment for
   LAMMPS.
  
3. The file [readme.md](./mainSimulations/readme.md) explains the simulation script file in detail.
