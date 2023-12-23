## LLPS simulation of (CAG)$_{40} \times 64$ at $200\mu M$ concentration
Here we explain the LAMMPS input file The file [lmp_main.in](./lmp_main.in) sets up the simulation environment for LAMMPS in detail.
 1. The unit, dimension, boundary condition and the atom_style are set initially:
      ```py
      units real
      dimension 3
      boundary p p p
      atom_style full
      ```
 2. Read the initial setup and set up the temperatures for annealing.
      ```py
      read_data  config_cag40.dat

      variable temperature1 equal 373
      variable temperature2 equal 363
      variable temperature3 equal 353
      variable temperature4 equal 343
      variable temperature5 equal 333
      variable temperature6 equal 323
      variable temperature7 equal 313
      variable temperature8 equal 303
      variable temperature9 equal 293
      variable randomSeed equal 13652
      ```
   3. Set up the force field
      ```py
      #### harmonic bond between bonds
      bond_style harmonic
      #  pair_coeff for harmonic bond, specify 3:
      #    * bond type 
      #    * K/2
      #    * equilibrium distance
      bond_coeff 1 7.5 5.9


      #### harmonic angle potential between bonds 
      angle_style harmonic
      #  pair_coeff for harmonic angle, specify 3:
      #    * angle type 
      #    * K/2
      #    * equilibrium angle in degrees
      angle_coeff 1 5.0 150.0


      #### Excluded volume WCA-LJ potential and Base-Pair Potential
      pair_style hybrid/overlay lj/cut 10.0 base/rna 18.0

      # pair_coeff 1 1 -5.0 3.0 13.8
      # pair_style  hybrid/overlay base/test 18.0 lj/cut 10.0

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
      #    * Equilibrium angle theta1
      #    * Equilibrium angle theta2
      #    * Equilibrium angle phi1
      #    * Equilibrium angle phi2
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



      special_bonds lj 0.0 0.0 1.0 # for 1-2, 1-3, and 1-4
      ```
4.  Assign initial velocities,  set up the neighbor list and set the timestep size.
    ```py
    velocity        all create ${temperature1} ${randomSeed}

    # neighbour list settings
    neighbor  25.0 bin

    # Timestep and computational parameters
    comm_style      tiled
    timestep        10

    neigh_modify    every 10 delay 0
    ```
5. Set up the Langevin thermostat for the first melting temperature.
   ```py
    fix             fxnve   all nve
    fix             fxlange all langevin ${temperature1} ${temperature1} 100000.0 ${randomSeed}
    fix             fxbal  all balance 1000 1.05 rcb
   ```
6. Set up the thermo out put, restart file writing frequency and the trajectory dump file frequency
   ```py
   # Thermo output settings
   thermo          1000
   thermo_style    custom step pe
   thermo_modify   flush yes

   restart      1000000 restart

   dump            1 all custom 10000 40CAG_200um_${temperature1}_${temperature9}.lammpstrj id mol type xu yu zu
   dump_modify     1 sort id
   ```
7. Anneal from 373K to 293K by simulating the system at each temperature for 5ns.
    ```py
    # run1
    run             500000

    # run2
    #-------------------
    fix             fxlange all langevin ${temperature2} ${temperature2} 100000.0 ${randomSeed}

    run             500000

    # run3
    #-------------------
    fix             fxlange all langevin ${temperature3} ${temperature3} 100000.0 ${randomSeed}

    run             500000

    # run4
    #-------------------
    fix             fxlange all langevin ${temperature4} ${temperature4} 100000.0 ${randomSeed}

    run             500000

    # run5
    #-------------------
    fix             fxlange all langevin ${temperature5} ${temperature5} 100000.0 ${randomSeed}

    run             500000

    # run6
    #-------------------
    fix             fxlange all langevin ${temperature6} ${temperature6} 100000.0 ${randomSeed}

    run             500000

    # run7
    #-------------------
    fix             fxlange all langevin ${temperature7} ${temperature7} 100000.0 ${randomSeed}

    run             500000

    # run8
    #-------------------
    fix             fxlange all langevin ${temperature8} ${temperature8} 100000.0 ${randomSeed}

    run             500000
    ```
8. Finally, simulate the system at 293K for 2 microseconds.
    ```py
    # run5
    #-------------------
    fix             fxlange all langevin ${temperature9} ${temperature9} 100000.0 ${randomSeed}

    run             200000000
    ```

To simulated this system with 16 CPUs, run the following command after copying the files to the simulation folder:
```sh
$ mpirun -np 8 ~/.local/bin/lmp_rna -in lmp_main.in
```