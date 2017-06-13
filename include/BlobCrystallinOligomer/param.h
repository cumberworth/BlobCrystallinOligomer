// param.h

#ifndef PARAM_H
#define PARAM_H

#include <string>

namespace param {

    using std::string;

    class InputParams {
        public:
            InputParams(int argc, char* argv[]);

            // System input 
            string m_config_filename;
            string m_energy_filename;
    };
}

#endif // PARAM_H
