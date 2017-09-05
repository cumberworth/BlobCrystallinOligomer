// simulation.cpp

#include <iostream>

#include "BlobCrystallinOligomer/simulation.h"

namespace simulation {

    using movetype::RotationVMMCMovetype;
    using movetype::TranslationVMMCMovetype;
    using movetype::NTDFlipMCMovetype;
    using std::cout;
    using std::setw;

    NVTMCSimulation::NVTMCSimulation(Config& conf, Energy& ene,
            InputParams params, RandomGens& random_num):
            m_config {conf}, m_energy {ene}, m_random_num {random_num},
            m_beta {1/params.m_temp}, m_steps {params.m_steps},
            m_logging_freq {params.m_logging_freq},
            m_config_output_freq {params.m_config_output_freq},
            m_op_output_freq {params.m_op_output_freq},
            m_config_file {params.m_output_filebase + ".vtf", conf} {

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
            if (m_logging_freq and step % m_logging_freq == 0) {
                log_move(step, movetype.label(), accepted);
            }

            // Output configuration and order parameters
            if (m_config_output_freq and step % m_config_output_freq == 0) {
                m_config_file.write_step(m_config, step);
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

        // Prepare vectors that track movetype information
        for (size_t i {0}; i != m_movetypes.size(); i++) {
            m_move_attempts.push_back(0);
            m_move_accepts.push_back(0);
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
