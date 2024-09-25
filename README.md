# DaRT + RBE Simulation

Simulate a DaRT seed and the associated DNA damage in a two part simulation.

## Origin
Repo created from two repos:
  - git@github.com:lg9783/DaRT.git/main at commit cf85bac0cef4292bd5818d35b9e3141e5a319558 
  - git@github.com:lg9783/RBE.git/main at commit 687700c5cde0efff0275dbc8966887792f094a65

## Overview
The two parts of the simulation are:

1. [The DaRT seed simulation](./DaRT), simulates the seed and saves the particle distribution in rings at different radii from the seed.
1. The particle distribution provides the input to the [DNA damage simulation](./DNA_damage) which simulates the DNA damage in a specified radius from the seed due to the particle distribution.

More details on how to run each are in the readme files: [DaRT](./DaRT/README.md) and [DNA_damage](./DNA_damage/README.md).

## Example of how to run
This is an example of how to run both parts of the simulation. After both components have been compiled (assumed to be in a folder called 'build'), and the DNA geometry has been created by running (./DNA_damage/geometryFiles/writeContinuousGeometry.py).

### DaRT simulation

Create the dart simulation input file (dart.in)
```
/run/verbose 1
/control/verbose 2
/det/cellPos_min 155 micrometer
/det/cellPos_max 255 micrometer
/det/cellPos_Nrings 3
/run/initialize
/gun/particle ion
/gun/ion 88 224 0 0
/run/printProgress 1
/run/beamOn 10
```

Run the dart simulation, from the [results file](results)
```
../DaRT/build/dart -mac dart.in -out example_simulation  -seed 1
```

This should create two files:
1. example_simulation.root - this contains information about the simulation and a few basic plots
1. example_simulation.bin - this contains the saved particle spectrum (in binary format to save space)

### DNA damage simulation

The DNA damage simulation consists of two steps:
1. Geant4 simulation
2. Calculation of DNA damage clusters from the Geant4 simulation results

Run the DNA_damage Geant4 simulation, again from the [results file](results):
```
../DNA_damage/simulation/build/rbe -mac ../DNA_damage/simulation/rbe_PS.in -decayPS example_simulation.bin -out example_simulation_DNA -sugar ../DNA_damage/simulation/geometryFiles/sugarPos_n2_300nm_H75nm_R11nm_33hist_50deg.bin -histone ../DNA_damage/simulation/geometryFiles/histonePos_n2_300nm_H75nm_R11nm_33hist_50deg.bin -seed 1 
```
This should create a new root file called example_simulation_DNA.root

Run the DNA damage clustering algorithm, for all three rings (specified using --copy) in the DaRT simulation, again from the [results file](results). The example below will included all particles in the decay chain in the DNA damage calculation.

```
conda activate clustering

python ../DNA_damage/Clustering/run.py --filename ../root/example_simulation_DNA.root --output dart_19Jun_2555um_seed3216_all.csv --sugar ../DNA_damage/geometryFiles/sugarPos_n2_300nm_H75nm_R11nm_33hist_50deg.bin --continuous 1 --copy 0 --simulationType=decay --part1_particleSource "[alphaRa224,alphaRn220,alphaPo216,alphaBi212,alphaPo212,e-Rn220,e-Po216,e-Pb212,e-Bi212,e-Tl208,e-Po212,gammaRn220,gammaPo216,gammaPb212,gammaBi212,gammaTl208,gammaPo212]"
 
python ../DNA_damage/Clustering/run.py --filename ../root/example_simulation_DNA.root --output dart_19Jun_2755um_seed3216_all.csv --sugar ../DNA_damage/geometryFiles/sugarPos_n2_300nm_H75nm_R11nm_33hist_50deg.bin --continuous 1 --copy 1 --simulationType=decay --part1_particleSource "[alphaRa224,alphaRn220,alphaPo216,alphaBi212,alphaPo212,e-Rn220,e-Po216,e-Pb212,e-Bi212,e-Tl208,e-Po212,gammaRn220,gammaPo216,gammaPb212,gammaBi212,gammaTl208,gammaPo212]"
 
python ../DNA_damage/Clustering/run.py --filename ../root/example_simulation_DNA.root --output dart_19Jun_2955um_seed3216_all.csv --sugar ../DNA_damage/geometryFiles/sugarPos_n2_300nm_H75nm_R11nm_33hist_50deg.bin --continuous 1 --copy 2 --simulationType=decay --part1_particleSource "[alphaRa224,alphaRn220,alphaPo216,alphaBi212,alphaPo212,e-Rn220,e-Po216,e-Pb212,e-Bi212,e-Tl208,e-Po212,gammaRn220,gammaPo216,gammaPb212,gammaBi212,gammaTl208,gammaPo212]"
```