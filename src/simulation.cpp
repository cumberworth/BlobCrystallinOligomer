// simulation.cpp

#include <chrono>
#include <iostream>

#include "BlobCrystallinOligomer/simulation.h"

namespace simulation {

using movetype::MetMCMovetype;
using movetype::VMMCMovetype;
using std::cout;
using std::setw;
using std::chrono::steady_clock;

NVTMCSimulation::NVTMCSimulation(
        Config& conf,
        Energy& ene,
        InputParams params,
        RandomGens& random_num):
        m_config {conf},
        m_energy {ene},
        m_random_num {random_num},
        m_beta {1 / params.m_temp},
        m_steps {params.m_steps},
        m_duration {params.m_duration},
        m_logging_freq {params.m_logging_freq},
        m_config_output_freq {params.m_config_output_freq},
        m_op_output_freq {params.m_op_output_freq},
        m_vtf_file {params.m_output_filebase + ".vtf", conf},
        m_patch_file {params.m_output_filebase + ".patch"},
        m_pipe_vtf_file {params.m_output_filebase + "_pipe.vtf", conf},
        m_pipe_patch_file {params.m_output_filebase + "_pipe.patch"} {

    m_pipe_vtf_file.close();
    m_pipe_patch_file.close();
    construct_movetypes(params);
}

void NVTMCSimulation::run() {
    auto start = steady_clock::now();
    for (stepT step {1}; step != (m_steps + 1); step++) {

        // Do a move
        int movetype_i {select_movetype()};
        MCMovetype& movetype {*m_movetypes[movetype_i]};
        bool accepted {movetype.move()};
        m_move_attempts[movetype_i]++;
        m_move_accepts[movetype_i] += accepted;

        // Check if maximum allowed time reached
        std::chrono::duration<double> dt {(steady_clock::now() - start)};
        if (dt.count() > m_duration) {
            cout << "Maximum time allowed reached\n";
            break;
        }

        // Log
        if (m_logging_freq and step % m_logging_freq == 0) {
            log_move(step, movetype.get_label(), accepted);
        }

        // Output configuration and order parameters
        if (m_config_output_freq and step % m_config_output_freq == 0) {
            m_vtf_file.write_step(m_config, step);
            m_patch_file.write_step(m_config);
            m_pipe_vtf_file.open_write_close(m_config, step);
            m_pipe_patch_file.open_write_step_close(m_config);
        }
        // if (m_op_output_freq and step % m_op_output_freq) {
        //  Write op to file
        //}
    }
    log_summary();
}

void NVTMCSimulation::construct_movetypes(InputParams params) {
    // This is pretty ugly
    // Also creates the cumalitive probability array
    double cum_prob {0};
    if (params.m_translation_met) {
        MCMovetype* movetype;
        string label {"TranslationMetMCMovetype"};
        string movemap_type {"translation"};
        movetype = new MetMCMovetype {
                m_config, m_energy, m_random_num, params, label, movemap_type};
        m_movetypes.emplace_back(movetype);
        cum_prob += params.m_translation_met;
        m_cum_probs.push_back(cum_prob);
    }
    if (params.m_rotation_met) {
        MCMovetype* movetype;
        string label {"RotationMetMCMovetype"};
        string movemap_type {"rotation"};
        movetype = new MetMCMovetype {
                m_config, m_energy, m_random_num, params, label, movemap_type};
        m_movetypes.emplace_back(movetype);
        cum_prob += params.m_rotation_met;
        m_cum_probs.push_back(cum_prob);
    }
    if (params.m_translation_vmmc) {
        MCMovetype* movetype;
        string label {"TranslationVMMCMovetype"};
        string movemap_type {"translation"};
        movetype = new VMMCMovetype {
                m_config, m_energy, m_random_num, params, label, movemap_type};
        m_movetypes.emplace_back(movetype);
        cum_prob += params.m_translation_vmmc;
        m_cum_probs.push_back(cum_prob);
    }
    if (params.m_rotation_vmmc) {
        MCMovetype* movetype;
        string label {"RotationVMMCMovetype"};
        string movemap_type {"rotation"};
        movetype = new VMMCMovetype {
                m_config, m_energy, m_random_num, params, label, movemap_type};
        m_movetypes.emplace_back(movetype);
        cum_prob += params.m_rotation_vmmc;
        m_cum_probs.push_back(cum_prob);
    }
    if (params.m_ntd_flip) {
        MCMovetype* movetype;
        string label {"NTDFlipMCMovetype"};
        string movemap_type {"ntdflip"};
        movetype = new MetMCMovetype {
                m_config, m_energy, m_random_num, params, label, movemap_type};
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
    cout << "Movetype: " << label << "\n";
    cout << "Accepted: " << accepted << "\n";
    cout << "Energy: " << m_energy.calc_total_energy() << "\n";
    cout << "\n";
}

void NVTMCSimulation::log_summary() {
    cout << "Run summary"
         << "\n";
    cout << "Movetype" << setw(10);
    cout << "Attempts" << setw(10);
    cout << "Accepts" << setw(10);
    cout << "Frequency"
         << "\n";
    for (size_t i {0}; i != m_movetypes.size(); i++) {
        cout << m_movetypes[i]->get_label() << setw(10);
        cout << m_move_attempts[i] << setw(10);
        cout << m_move_accepts[i] << setw(10);
        cout << static_cast<double>(m_move_accepts[i]) / m_move_attempts[i]
             << "\n";
    }
}
} // namespace simulation
