// param.h

#ifndef PARAM_H
#define PARAM_H

#include <string>

#include "BlobCrystallinOligomer/shared_types.h"

namespace param {
 
    using shared_types::eneT;
    using shared_types::distT;
    using shared_types::stepT;
    using std::string;

    /** Parse fractions in strings
      *
      * E.g., "3/4"
      */
    class Fraction {
        public:
            Fraction(string unparsed_fraction);
            double to_double();
        private:
            double m_double_fraction;
            double m_numerator;
            double m_denominator;
    };


    /** Input parameters from file */
    class InputParams {
        public:
            InputParams(int argc, char* argv[]);

            // System input 
            string m_config_filename;
            string m_energy_filename;
            eneT m_temp;

            // Simulation options
            stepT m_steps;
            distT m_max_cutoff; // Not very nice to put here

            // Movetypes
            distT m_max_disp_tc;
            distT m_max_disp_rc;
            distT m_max_disp_a;
            double m_translation_met;
            double m_rotation_met;
            double m_translation_vmmc;
            double m_rotation_vmmc;
            double m_ntd_flip;
            
            // Output
            string m_output_filebase;
            stepT m_logging_freq;
            stepT m_config_output_freq;
            stepT m_op_output_freq;

        private:
            string m_rotation_met_raw;
            string m_translation_met_raw;
            string m_rotation_vmmc_raw;
            string m_translation_vmmc_raw;
            string m_ntd_flip_raw;

            /*  Process inputs that require more than default constructor */
            void post_process_inputs();
    };
}

#endif // PARAM_H
