#!/usr/env python3

"""Create coordinates hard sphere system with non-infinite energy.

Output in pdb and json format
"""

import argparse

import numpy as np

from crystallinpy import model


def main():
    args = parse_cl()
    monomers = []
    type = 0
    space = model.CuboidPBC(args.box_len)

    # Randomly insert into box and accept if there are not overlaps
    for i in range(args.num_particles):
        overlap = True
        particle = model.SimpleParticle(i, 'PAR', type)
        monomer = model.SingleParticleMonomer([particle], args.diameter/2, i)
        while overlap:
            monomer.pos = args.box_len*(np.random.rand(3) - 0.5)
            overlap = model.hard_sphere_overlap(monomer, monomers,
                    args.diameter, space)
        
        monomers.append(monomer)

    # Write files
    pdb_file = model.PDBConfigOutputFile(args.output_filebase + '.pdb')
    pdb_file.write(monomers)
    json_file = model.JSONConfigOutputFile(args.output_filebase + '.json')
    json_file.write(monomers, args.box_len)

    return overlap


def parse_cl():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'num_particles',
        type=int,
        help='Number of particles to place in system')
    parser.add_argument(
        'diameter',
        type=float,
        help='Diameter of particles')
    parser.add_argument(
        'box_len',
        type=float,
        help='Length of box to put particles into')
    parser.add_argument(
        'output_filebase',
        type=str,
        help='Output filebase')

    return parser.parse_args()


if __name__ == '__main__':
    main()
