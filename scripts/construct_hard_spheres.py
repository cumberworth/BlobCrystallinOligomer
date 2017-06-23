#!/usr/env python3

"""Create coordinates hard sphere system with non-infinite energy.

Output in pdb and json format
"""

import argparse

import numpy as np

import model


def main():
    args = parse_cl()
    monomers = []
    type = 0
    for i in range(args.num_particles):
        overlap = True
        particle = model.SimpleParticle(i, 'PAR', type)
        monomer = model.SingleParticleMonomer([particle], args.diameter/2, i)
        while overlap:
            monomer.pos = args.box_len*(np.random.rand(3) - 0.5)
            overlap = hard_sphere_overlap(monomer, monomers, args.diameter)
        
        monomers.append(monomer)

    pdb_file = model.PDBConfigOutputFile(args.output_filebase + '.pdb')
    pdb_file.write(monomers)
    json_file = model.JSONConfigOutputFile(args.output_filebase + '.json')
    json_file.write(monomers)


def hard_sphere_overlap(new_hs, hses, diameter):
    """Check for overlap between new hard sphere and system of hard spheres"""
    overlap = False
    for hs in hses:
        if np.linalg.norm(new_hs.pos - hs.pos) < diameter:
            overlap = True
            break

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
        type=int,
        help='Length of box to put particles into')
    parser.add_argument(
        'output_filebase',
        type=str,
        help='Output filebase')

    return parser.parse_args()


if __name__ == '__main__':
    main()
