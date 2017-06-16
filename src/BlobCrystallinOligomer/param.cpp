// param.cpp

#include <iostream>
#include <fstream>
#include <string>

#include "boost/program_options.hpp"

#include "BlobCrystallinOligomer/param.h"
#include "BlobCrystallinOligomer/shared_types.h"

namespace param {

    namespace po = boost::program_options;

    using shared_types::eneT;
    using shared_types::distT;
    using std::string;
    using std::cout;

    InputParams::InputParams(int argc, char* argv[]) {
        
        // Command line options
        po::options_description cl_options {"Command line options"};
        cl_options.add_options()
            ("parameter_filename,i", po::value<string>())
        ;

        // Parameter file options
        po::options_description inp_options {"System input"};
        inp_options.add_options()
            ("config_filename",
                po::value<string>(&m_config_filename)->default_value(""),
                "File for system topology")
            ("energy_filename",
                po::value<string>(&m_energy_filename)->default_value(""),
                "File for system energy")
            ("temp",
                po::value<eneT>(&m_temp)->default_value(300),
                "Maximum dispacement for translations")
        ;
        po::options_description move_options {"Movetype options"};
        move_options.add_options()
            ("max_disp_tc",
                po::value<distT>(&m_max_disp_tc)->default_value(1),
                "Maximum dispacement for translations")
            ("max_disp_rc",
                po::value<distT>(&m_max_disp_rc)->default_value(1),
                "Maximum displacement for selecting center of rotation")
            ("max_disp_a",
                po::value<distT>(&m_max_disp_a)->default_value(1),
                "Maximum displacement for selecting rotation angle")
        ;

        // Displayed options
        po::options_description displayed_options {"Allowed options"};
        displayed_options.add(cl_options).add(inp_options).add(move_options);

        // Parse command line input
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, displayed_options), vm);
        po::notify(vm);
        if (vm.count("help")) {
            cout << displayed_options;
        }
        if (not vm.count("parameter_filename")) {
            cout << "Input parameter file must be provided.\n";
            std::exit(1);
        }
        string param_filename {vm["parameter_filename"].as<string>()};

        // Parse parameter file input
        std::ifstream param_file {param_filename};
        po::store(po::parse_config_file(param_file, displayed_options), vm);
        po::notify(vm);
    }
}
