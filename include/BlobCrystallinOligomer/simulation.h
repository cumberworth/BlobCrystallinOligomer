// simulation.h

#ifndef SIMULATION_H
#define SIMULATION_H

#include "BlobCrystallinOligomer/param.h"
#include "BlobCrystallinOligomer/config.h"
#include "BlobCrystallinOligomer/energy.h"

namespace simulation {

    using param::InputParams;
    using config::Config;
    using energy::Energy;

    class MCSimulation {
        public:
            virtual void run() = 0;

    };

    class NVTMCSimulation {
        public:
            NVTMCSimulation(
                    Config& config,
                    Energy& energy,
                    InputParams params);
            void run();
    };
}

#endif // SIMULATION_H
