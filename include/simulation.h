// simulation.h

#ifndef SIMULATION_H
#define SIMULATION_H

namespace simulation {

    class MCSimulation {
        public:
            virtual void run() = 0;

    };

    class NTVMCSimulation {
        public:
            void run();
    };
}

#endif // SIMULATION_H
