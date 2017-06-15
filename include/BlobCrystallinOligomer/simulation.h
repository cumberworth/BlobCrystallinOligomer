// simulation.h

#ifndef SIMULATION_H
#define SIMULATION_H

#include "BlobCrystallinOligomer/param.h"
#include "BlobCrystallinOligomer/config.h"
#include "BlobCrystallinOligomer/energy.h"
#include "BlobCrystallinOligomer/random_gens.h"

namespace simulation {

    using param::InputParams;
    using config::Config;
    using energy::Energy;
    using random_gens::RandomGens;

    class MCSimulation {
        public:
            virtual void run() = 0;

        private:
            RandomGens m_random_num;

    };

    class NVTMCSimulation {
        public:
            NVTMCSimulation(
                    Config& config,
                    Energy& energy,
                    InputParams params,
                    RandomGens& random_num);
            void run();
    };
}

#endif // SIMULATION_H
