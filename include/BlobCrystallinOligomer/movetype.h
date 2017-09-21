// movetype.h

#ifndef MOVETYPE_H
#define MOVETYPE_H

#include <cmath>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "BlobCrystallinOligomer/config.h"
#include "BlobCrystallinOligomer/energy.h"
#include "BlobCrystallinOligomer/hash.h"
#include "BlobCrystallinOligomer/monomer.h"
#include "BlobCrystallinOligomer/param.h"
#include "BlobCrystallinOligomer/shared_types.h"
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
    using std::unique_ptr;
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

    class Movemap {
        public:
            Movemap(RandomGens& random_num);
            virtual ~Movemap() {}

            virtual void generate_movemap(Monomer& monomer) = 0;
            virtual void apply_movemap(Monomer& monomer) = 0;

        protected:
            RandomGens& m_random_num;
    };

    class TranslationMovemap:
            public Movemap {

        public:
            TranslationMovemap(distT max_disp_tc, RandomGens& random_num);

            void generate_movemap(Monomer&);
            void apply_movemap(Monomer& monomer);

        private:
            distT m_max_disp_tc;

            vecT m_disp_v;
    };

    class RotationMovemap:
            public Movemap {

        public:
            RotationMovemap(
                    distT max_disp_rc,
                    distT max_disp_a,
                    RandomGens& random_num);

            void generate_movemap(Monomer& monomer);
            void apply_movemap(Monomer& monomer);

        private:
            distT m_max_disp_rc;
            distT m_max_disp_a;

            vecT m_rot_c; // Centre of rotation
            rotMatT m_rot_mat; // Rotation matrix
    };

    /** General movetype interface */
    class MCMovetype {
        public:
            MCMovetype(Config& conf, Energy& ene, RandomGens& random_num,
                    InputParams params, string label);
            virtual ~MCMovetype() {}

            /** Attempt move and return result (accepted or rejectd) */
            virtual bool move() = 0;

            string get_label();

        protected:
            Config& m_config;
            Energy& m_energy;
            RandomGens& m_random_num;
            eneT m_beta;
            string m_label;
            unique_ptr<Movemap> m_movemap;
    };

    /** Flip the NTD conformation
      *
      * There are two available conformations for the NTD. This movetype will
      * attempts to flip the conformation via reflection in a plane of a
      * randomly selected monomer.
      */
    class NTDFlipMCMovetype:
            public MCMovetype {

        public:
            using MCMovetype::MCMovetype;
            bool move() {}

        private:
    };

    /** Metropolis move */
    class MetMCMovetype:
            public MCMovetype {

        public:
            MetMCMovetype(Config& conf, Energy& ene, RandomGens& random_num,
                    InputParams params, string label, string movemap_type);

            bool move();

        protected:
            bool accept_move(eneT de);
    };

    /** Virtual move */
    class VMMCMovetype:
            public MCMovetype {

        public:
            VMMCMovetype(Config& conf, Energy& ene, RandomGens& random_num,
                    InputParams params, string label, string movemap_type);

            bool move();

        private:

            monomerArrayT m_cluster;
            int m_frustrated_links {0};
            vector<int> m_frustrated_mis {};
            set<pair<int, int>> m_proposed_pairs {};
            set<int> m_interacting_mis;
            set<pair<int, int>> m_pair_mis;
            unique_ptr<Movemap> m_movemap;

            void add_interacting_pairs(Monomer& monomer1);
            pair<int, int> pop_random_pair();
            double calc_prelink_prob(eneT ene1, eneT ene2);
            bool accept_prelink(double prelink_p);
            bool accept_link(double prelink_for_p, double prelink_rev_p);
            bool accept_move();
            void reset_internal();
    };
}

#endif // MOVETYPE_H
