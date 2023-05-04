# BlobCrystallinOligomer

Modeling tools for simulations of oligomers consisting of very-coarse-grained representations of Alpha-B Crystallin.

## Installation

The program requires the Boost library, specifically the program options module.
Building and installation is done with CMake.
To configure, build and install in `installdir`, run
```
CXX=clang++ cmake -S . -B build -DCMAKE_INSTALL_PREFIX=[installdir] -DBUILD_TESTING=OFF
cmake --build build
cmake --install build
```

## Running simulations

Example input files for running simulations can be found in `scripts/examples/`.
To see all configuration options, run

`blobCrystallinOligomer -h`

To run a simulation, enter

`blobCrystallinOligomer -i [configuration file] > [log file]`

## Viewing configurations

Tcl scripts for VMD are available in `scripts/vmd` for viewing configurations.
To watch a simulation live, run

`vmd -e [scripts dir path]/pipe.tcl -args [vmd scripts directory] [output filebase]`

To watch a simulation that has finished, start VMD and run

```
set libdir [vmd scripts directory]
set filebase [output filebase
source [vmd scripts directory]/view_coors.tcl
```
