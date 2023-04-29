// main.cpp

#include <memory>

#include "BlobCrystallinOligomer/param.h"
#include "BlobCrystallinOligomer/config.h"
#include "BlobCrystallinOligomer/energy.h"
#include "BlobCrystallinOligomer/random_gens.h"
#include "BlobCrystallinOligomer/simulation.h"

int main(int argc, char* argv[]) {

    using std::unique_ptr;
    using std::make_unique;

    param::InputParams params {argc, argv};
    auto random_num {make_unique<random_gens::RandomGens>()};
    auto conf {make_unique<config::Config>(params, *random_num)};
    auto ene {make_unique<energy::Energy>(*conf, params)};
    simulation::NVTMCSimulation sim {*conf, *ene, params, *random_num};
    sim.run();
}
