// movetype.h

#ifndef MOVETYPE_H
#define MOVETYPE_H

#include <cmath>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "BlobCrystallinOligomer/config.h"
#include "BlobCrystallinOligomer/energy.h"
#include "BlobCrystallinOligomer/hash.h"
#include "BlobCrystallinOligomer/monomer.h"
#include "BlobCrystallinOligomer/param.h"
#include "BlobCrystallinOligomer/random_gens.h"

namespace movetype {

    using config::Config;
    using config::Monomer;
    using config::monomerArrayT;
    using energy::Energy;
    using param::InputParams;
    using random_gens::RandomGens;
    using shared_types::eneT;
    using shared_types::distT;
    using shared_types::rotMatT;
    using shared_types::vecT;
    using std::pair;
    using std::set;
    using std::sqrt;
    using std::string;
    using std::vector;

    /** Return a unit vector with uniform distrition across sphere surface */
    vecT random_unit_vector(RandomGens& random_num);

    /** Return uniform random distance between -max_disp and max_disp
      *
      * Range is [-max_disp, max_disp)
      */
    distT random_displacement(distT max_disp, RandomGens& random_num);

    /** General movetype interface */
    class MCMovetype {
        public:
            MCMovetype(Config& conf, Energy& ene, RandomGens& random_num,
                    InputParams params);

            /** Attempt move and return result (accepted or rejectd) */
            virtual bool move() = 0;
            virtual string label() {return "MCMovetype";}

        protected:
            Config& m_config;
            Energy& m_energy;
            RandomGens& m_random_num;
            eneT m_beta;
    };

    /** Flip the NTD conformation
      *
      * There are two available conformations for the NTD. This movetype will
      * attempts to flip the conformation via reflection in a plane of a
      * randomly selected monomer.
      */
    class NTDFlipMCMovetype: public MCMovetype {
        public:
            using MCMovetype::MCMovetype;
            bool move() {}
            string label() {return "NTDFlipMCMovetype";}
    };

    /** Shared implementation for virtual moves */
    class VMMCMovetype: public MCMovetype {
        public:
            using MCMovetype::MCMovetype;
            bool move();
            string label() {return "VMMCMovetype";}

        private:
            monomerArrayT m_cluster;
            int m_frustrated_links {0};
            vector<int> m_frustrated_mis {};
            set<int> m_interacting_mis;
            set<pair<int, int>> m_pair_mis;

            void virtual generate_movemap(Monomer& seed_monomer) = 0;
            void virtual apply_movemap(Monomer& monomer) = 0;

            void add_interacting_pairs(Monomer& monomer1);
            pair<int, int> select_random_pair();
            double calc_prelink_prob(eneT ene1, eneT ene2);
            bool accept_prelink(double prelink_p);
            bool accept_link(double prelink_for_p, double prelink_rev_p);
            bool accept_move();
            void reset_internal();
    };

    /** Translational virtual movetype */
    class TranslationVMMCMovetype: public VMMCMovetype {
        public:
            TranslationVMMCMovetype(Config& conf, Energy& ene, RandomGens& random_num,
                    InputParams params);
            string label() {return "TranslationVMMCMovetype";}
        private:
            void generate_movemap(Monomer&);
            void apply_movemap(Monomer& monomer);

            distT m_max_disp_tc;
            vecT m_disp_v;
    };

    /** Rotational virtual movetype */
    class RotationVMMCMovetype: public VMMCMovetype {
        public:
            RotationVMMCMovetype(Config& conf, Energy& ene, RandomGens& random_num,
                    InputParams params);
            string label() {return "RotationVMMCMovetype";}
        private:
            void generate_movemap(Monomer& seed_monomer);
            void apply_movemap(Monomer& monomer);

            distT m_max_disp_rc;
            distT m_max_disp_a;
            vecT m_rot_c; // Centre of rotation
            rotMatT m_rot_mat; // Rotation matrix
            
    };
}

#endif // MOVETYPE_H
