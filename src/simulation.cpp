// simulation.cpp

#include <iostream>
#include <memory>

#include "BlobCrystallinOligomer/config.h"
#include "BlobCrystallinOligomer/energy.h"
#include "BlobCrystallinOligomer/movetype.h"
#include "BlobCrystallinOligomer/ofile.h"
#include "BlobCrystallinOligomer/param.h"
#include "BlobCrystallinOligomer/random_gens.h"
#include "BlobCrystallinOligomer/simulation.h"

namespace simulation {

    using config::Config;
    using energy::Energy;
    using movetype::RotationVMMCMovetype;
    using movetype::TranslationVMMCMovetype;
    using movetype::NTDFlipMCMovetype;
    using ofile::OutputConfigsFile;
    using param::InputParams;
    using random_gens::RandomGens;
    using shared_types::stepT;
    using std::cout;
    using std::setw;
    using std::unique_ptr;

    NVTMCSimulation::NVTMCSimulation(Config& conf, Energy& ene,
            InputParams params, RandomGens& random_num):
            m_config {conf}, m_energy {ene}, m_random_num {random_num},
            m_beta {1/params.m_temp}, m_logging_freq {params.m_logging_freq},
            m_config_output_freq {params.m_config_output_freq},
            m_op_output_freq {params.m_op_output_freq},
            m_config_file {params.m_output_filebase + ".confs", conf} {

        construct_movetypes(params);
    }

    void NVTMCSimulation::run() {
        for (stepT step {1}; step != (m_steps + 1); step++) {
            // Do a move
            int movetype_i {select_movetype()};
            MCMovetype& movetype {*m_movetypes[movetype_i]};
            bool accepted {movetype.move()};
            m_move_attempts[movetype_i]++;
            m_move_accepts[movetype_i] += accepted;

            // Log
            if (m_logging_freq and step % m_logging_freq) {
                log_move(step, movetype.label(), accepted);
            }

            // Output configuration and order parameters
            if (m_config_output_freq and step % m_config_output_freq) {
                m_config_file.write_timestep(m_config, step);
            }
            //if (m_op_output_freq and step % m_op_output_freq) {
                // Write op to file
            //}
        }
        log_summary();
    }

    void NVTMCSimulation::construct_movetypes(InputParams params) {
        // This is pretty ugly
        // Also creates the cumalitive probability array
        double cum_prob {0};
        if (params.m_rotation_vmmc) {
            MCMovetype* movetype;
            movetype = new RotationVMMCMovetype {m_config, m_energy,
                    m_random_num, params};
            m_movetypes.emplace_back(movetype);
            cum_prob += params.m_rotation_vmmc;
            m_cum_probs.push_back(cum_prob);
        }
        if (params.m_translation_vmmc) {
            MCMovetype* movetype;
            movetype = new TranslationVMMCMovetype {m_config, m_energy,
                    m_random_num, params};
            m_movetypes.emplace_back(movetype);
            cum_prob += params.m_translation_vmmc;
            m_cum_probs.push_back(cum_prob);
        }
        if (params.m_ntd_flip) {
            MCMovetype* movetype;
            movetype = new NTDFlipMCMovetype {m_config, m_energy,
                    m_random_num, params};
            m_movetypes.emplace_back(movetype);
            cum_prob += params.m_ntd_flip;
            m_cum_probs.push_back(cum_prob);
        }
    }

    int NVTMCSimulation::select_movetype() {
        double prob {m_random_num.uniform_real()};
        size_t i;
        for (i = 0; i != m_cum_probs.size(); i++) {
            if (prob < m_cum_probs[i]) {
                break;
            }
        }

        return i;
    }

    void NVTMCSimulation::log_move(stepT step, string label, bool accepted) {
        cout << "Step: " << step << "\n";
        cout << "Movetype : " << label << "\n";
        cout << "Accepted : " << accepted << "\n";
        cout << "\n";
    }

    void NVTMCSimulation::log_summary() {
        cout << "Run summary" << "\n";
        cout << "Movetype" << setw(10);
        cout << "Attempts" << setw(10);
        cout << "Accepts" << setw(10);
        cout << "Frequency" << "\n";
        for (size_t i {0}; i != m_movetypes.size(); i++) {
            cout << m_movetypes[i]->label() << setw(10);
            cout << m_move_attempts[i] << setw(10);
            cout << m_move_accepts[i] << setw(10);
            cout << static_cast<double>(m_move_attempts[i]) /
                    m_move_accepts[i] << "\n";
        }
    }
}
