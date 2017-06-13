// main.cpp

#include "BlobCrystallinOligomer/param.h"
#include "BlobCrystallinOligomer/config.h"
#include "BlobCrystallinOligomer/energy.h"
#include "BlobCrystallinOligomer/simulation.h"

int main(int argc, char* argv[]) {
    param::InputParams params {argc, argv};
    config::Config conf {params};
    energy::Energy ene {conf, params};
    simulation::NVTMCSimulation sim {conf, ene, params};
    sim.run();
}
