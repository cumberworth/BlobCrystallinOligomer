#!/usr/env python3

"""Create random configuration of alphaB monomers with no overlaps

Output in pdb and json format
"""

import argparse
import copy
import math
import random

import numpy as np
import transforms3d as trans

import construct_spherical_oligomer as cso
from crystallinpy import model


def main():
    args = parse_cl()
    monomers = []
    type = 0
    space = model.CuboidPBC(args.box_len)

    # Randomly insert into box and accept if there are no overlaps
    radius = args.diameter/2
    acd_ntd_angle = 0.6236561373425675
    blob_angle1 = -0.318827603743
    blob_angle2 = -0.534300804011
    s_monomers = cso.construct_monomers(radius, radius, 2, 2, args.num_monomers,
            acd_ntd_angle, blob_angle1, blob_angle2)
    monomers = []
    for m in s_monomers:
        overlap = True
        while overlap:
            x1r = trans.axangles.axangle2aff([1, 0, 0], 2*math.pi*random.random(), m.center)
            yr = trans.axangles.axangle2aff([0, 1, 0], 2*math.pi*random.random(), m.center)
            x2r = trans.axangles.axangle2aff([1, 0, 0], 2*math.pi*random.random(), m.center)
            m.apply_transformation(x1r)
            m.apply_transformation(yr)
            m.apply_transformation(x2r)
            pos = np.concatenate([args.box_len*(np.random.rand(3) - 0.5), [0]])
            ref_pos = copy.deepcopy(m.particles[0].pos)
            m.particles[0].pos = pos
            for p in m.particles[1:]:
                pdiff = p.pos - ref_pos
                p.pos = space.wrap(pos + pdiff)

            overlap = model.hard_sphere_overlap(m, monomers,
                    args.diameter, space)
        
        monomers.append(m)

    # Write files
    pdb_file = model.PDBConfigOutputFile(args.output_filebase + '.pdb')
    pdb_file.write(monomers)
    json_file = model.JSONConfigOutputFile(args.output_filebase + '.json')
    json_file.write(monomers, args.box_len)

    return overlap


def random_unit_vector():
    """Return a vector on the unit sphere with uniform probability
    
    Taken from Daan's book, which is taken from Allen and Tildesley
    """
    ransq = 2.;
    while ransq >= 1:
        ran1 = 1 - 2*random.random()
        ran2 = 1 - 2*random.random()
        ransq = ran1**2 + ran2**2

    ranh = 2*np.sqrt(1 - ransq)
    x = ran1*ranh
    y = ran2*ranh
    z = 1 - 2*ransq

    return np.array([x, y, z])


def parse_cl():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'num_monomers',
        type=int,
        help='Number of monomers to place in system')
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
