// param.h

#ifndef PARAM_H
#define PARAM_H

#include <string>

#include "BlobCrystallinOligomer/shared_types.h"

namespace param {
 
    using shared_types::eneT;
    using shared_types::distT;
    using std::string;

    class InputParams {
        public:
            InputParams(int argc, char* argv[]);

            // System input 
            string m_config_filename;
            string m_energy_filename;
            eneT m_temp;

            // Movetypes
            distT m_max_disp_tc;
            distT m_max_disp_rc;
            distT m_max_disp_a;
    };
}

#endif // PARAM_H
