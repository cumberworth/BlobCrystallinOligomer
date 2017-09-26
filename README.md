# BlobCrystallinOligomer

Modeling tools for simulations of oligomers consisting of very-coarse-grained representations of Alpha-B Crystallin.

## Installation

A Makefile is provided for compilation and installation.
To install with the default settings run

`make && make install`

To set the installation directory, modify the `TARGETDIR` variable.
To change the optimization level or add debug information (`g`), modify the OPTLEVEL variable.
Boost is required by the program, specifically a version with the program options module.
If you are need to install it locally on a cluster, follow the Makefile instructions in the comments to specify the locations of the headers and library; you will also need to modify the path variable:

`export LD_LIBRARY_PATH=${BOOSTLIBRARYLOCATION}:$LD_LIBRARY_PATH`

### Setting up on Aurora

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

`
set libdir [vmd scripts directory]
set filebase [output filebase
source [vmd scripts directory]/view_coors.tcl
`
