// param.cpp

#include <iostream>
#include <fstream>
#include <string>

#include "boost/program_options.hpp"

#include "BlobCrystallinOligomer/param.h"

namespace param {

    namespace po = boost::program_options;

    using std::cout;

    Fraction::Fraction(string unparsed_fraction) {
        string delimiter {"/"};
        auto delim_pos {unparsed_fraction.find(delimiter)};

        // Assume is real number
        if (delim_pos == string::npos) {
            m_numerator = stod(unparsed_fraction);
            m_denominator = 1;
        }
        else {
            auto end_pos {unparsed_fraction.size()};
            m_numerator = stod(unparsed_fraction.substr(0, delim_pos));
            m_denominator = stod(unparsed_fraction.substr(delim_pos + 1, end_pos));
        }
        m_double_fraction = m_numerator / m_denominator;
    }

    double Fraction::to_double() {
        return m_double_fraction;
    }

    InputParams::InputParams(int argc, char* argv[]) {
        
        // Displayed options
        po::options_description displayed_options {"Allowed options"};

        // Command line options
        po::options_description cl_options {"Command line options"};
        cl_options.add_options()
            ("parameter_filename,i",
                po::value<string>(),
                "Input file")
            ("help,h",
             "Display available options")
        ;
        displayed_options.add(cl_options);

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
        displayed_options.add(inp_options);

        po::options_description sim_options {"Simulation options"};
        sim_options.add_options()
            ("steps",
                po::value<stepT>(&m_steps)->default_value(0),
                "Number of steps")
            ("max_cutoff",
                po::value<distT>(&m_max_cutoff)->default_value(0),
                "Maximum cutoff value of any included potential")
        ;
        displayed_options.add(sim_options);

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
            ("translation_met",
                po::value<string>(&m_translation_met_raw)->default_value("0"),
                "Probability of performing a translation Metropolis movetype")
            ("rotation_met",
                po::value<string>(&m_rotation_met_raw)->default_value("0"),
                "Probability of performing a rotation Metropolis movetype")
            ("translation_vmmc",
                po::value<string>(&m_translation_vmmc_raw)->default_value("0"),
                "Probability of performing a translation VMMC")
            ("rotation_vmmc",
                po::value<string>(&m_rotation_vmmc_raw)->default_value("0"),
                "Probability of performing a rotation VMMC")
            ("ntd_flip",
                po::value<string>(&m_ntd_flip_raw)->default_value("0"),
                "Probability of performing a NTD flip")
        ;
        displayed_options.add(move_options);

        po::options_description output_options {"Output options"};
        output_options.add_options()
            ("output_filebase",
                po::value<string>(&m_output_filebase)->default_value(""),
                "Filebase for all output files")
            ("logging_freq",
                po::value<stepT>(&m_logging_freq)->default_value(1),
                "Logging frequency")
            ("config_output_freq",
                po::value<stepT>(&m_config_output_freq)->default_value(0),
                "Configuration output frequency")
            ("op_output_freq",
                po::value<stepT>(&m_op_output_freq)->default_value(0),
                "Order parameters output frequency")
        ;
        displayed_options.add(output_options);

        // Parse command line input
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, displayed_options), vm);
        po::notify(vm);
        if (vm.count("help")) {
            cout << "\n";
            cout << displayed_options;
            cout << "\n";
            exit(1);
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

        post_process_inputs();
    }

    void InputParams::post_process_inputs() {

        // Parse fractions for movetype probabilities
        Fraction rotation_met_fraction {m_rotation_met_raw};
        m_rotation_met = rotation_met_fraction.to_double();
        Fraction translation_met_fraction {m_translation_met_raw};
        m_translation_met = translation_met_fraction.to_double();
        Fraction rotation_vmmc_fraction {m_rotation_vmmc_raw};
        m_rotation_vmmc = rotation_vmmc_fraction.to_double();
        Fraction translation_vmmc_fraction {m_translation_vmmc_raw};
        m_translation_vmmc = translation_vmmc_fraction.to_double();
        Fraction ntd_flip_fraction {m_ntd_flip_raw};
        m_ntd_flip = ntd_flip_fraction.to_double();
    }
}
