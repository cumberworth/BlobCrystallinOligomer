// simulation.h

#ifndef SIMULATION_H
#define SIMULATION_H

#include <memory>
#include <string>
#include <vector>

#include "BlobCrystallinOligomer/param.h"
#include "BlobCrystallinOligomer/config.h"
#include "BlobCrystallinOligomer/energy.h"
#include "BlobCrystallinOligomer/movetype.h"
#include "BlobCrystallinOligomer/ofile.h"
#include "BlobCrystallinOligomer/random_gens.h"
#include "BlobCrystallinOligomer/shared_types.h"

namespace simulation {

    using param::InputParams;
    using config::Config;
    using energy::Energy;
    using movetype::MCMovetype;
    using ofile::VTFOutputFile;
    using ofile::PatchOutputFile;
    using random_gens::RandomGens;
    using shared_types::eneT;
    using shared_types::stepT;
    using shared_types::timeT;
    using std::unique_ptr;
    using std::string;
    using std::vector;

    /** Simulation interface */
    class MCSimulation {
        public:
            virtual ~MCSimulation();
            virtual void run() = 0;
    };

    // Eventually migrate much of this to a more general class if other simulation
    // method classes are to be designed

    /** Canonical ensemble simulation */
    class NVTMCSimulation {
        public:
            NVTMCSimulation(
                    Config& conf,
                    Energy& ene,
                    InputParams params,
                    RandomGens& random_num);
            void run();

        private:
            Config& m_config;
            Energy& m_energy;
            RandomGens& m_random_num;
            eneT m_beta;

            vector<unique_ptr<MCMovetype>> m_movetypes;
            vector<double> m_cum_probs;
            vector<stepT> m_move_attempts;
            vector<stepT> m_move_accepts;

            stepT m_steps;
            timeT m_duration;
            stepT m_logging_freq;
            stepT m_config_output_freq;
            stepT m_op_output_freq;

            VTFOutputFile m_vtf_file;
            PatchOutputFile m_patch_file;
            VTFOutputFile m_pipe_vtf_file;
            PatchOutputFile m_pipe_patch_file;

            void construct_movetypes(InputParams params);
            void setup_output_files(InputParams params);
            int select_movetype();
            void log_move(stepT step, string movetype_label, bool accepted);
            void log_summary();
    };
}

#endif // SIMULATION_H
