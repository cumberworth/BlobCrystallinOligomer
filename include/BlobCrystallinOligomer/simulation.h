// simulation.h

#ifndef SIMULATION_H
#define SIMULATION_H

#include "BlobCrystallinOligomer/space.h"
#include "BlobCrystallinOligomer/config.h"
#include "BlobCrystallinOligomer/energy.h"

namespace simulation {

    using params::InputParams;
    using energy::Energy;

    class MCSimulation {
        public:
            virtual void run() = 0;

    };

    template<typename configT>
    class NTVMCSimulation {
        public:
            NTVMCSimulation(
                    configT& config,
                    Energy<configT>& energy,
                    InputParams params);
            void run();
    };
}

#endif // SIMULATION_H
