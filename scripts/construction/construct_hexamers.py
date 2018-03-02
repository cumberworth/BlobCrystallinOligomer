#!/usr/env python3

"""Create random configuration of alphaB hexamers with no overlaps

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
    type = 0
    space = model.CuboidPBC(args.box_len)

    # Randomly insert into box and accept if there are no overlaps
    radius = args.diameter/2
    acd_ntd_angle = 0.6236561373425675
    blob_angle1 = -0.318827603743
    blob_angle2 = -0.534300804011
    set_monomers = []
    monomers = cso.construct_monomers(radius, radius, 2, 2, 6*args.num_hexamers,
            acd_ntd_angle, blob_angle1, blob_angle2)
    for i in range(args.num_hexamers):
        h_monomers = [m for m in monomers[i*6:i*6 + 6]]
        dimers = []
        for j in range(0, 5, 2):
            monomer1 = h_monomers[j]
            monomer2 = h_monomers[j + 1]
            dimer = cso.construct_dimer(monomer1, monomer2)
            dimers.append(dimer)

        hexamer = cso.construct_hexamer(dimers[0], dimers[1], dimers[2])

        overlap = [True]*6
        while np.any(overlap):
            #angle1 = 2*math.pi*random.random()
            #angle2 = 2*math.pi*random.random()
            #angle3 = 2*math.pi*random.random()
            delta_pos = args.box_len * (np.random.rand(3) - 0.5)
            delta_pos = np.concatenate([delta_pos, [0]])
            for j, m in enumerate(h_monomers):
                #x1r = trans.axangles.axangle2aff([1, 0, 0], angle1, m.center)
                #yr = trans.axangles.axangle2aff([0, 1, 0], angle2, m.center)
                #x2r = trans.axangles.axangle2aff([1, 0, 0], angle3, m.center)
                #m.apply_transformation(x1r)
                #m.apply_transformation(yr)
                #m.apply_transformation(x2r)
                for p in m.particles:
                    p.pos = space.wrap(p.pos + delta_pos)

                overlap[j] = model.hard_sphere_overlap(m, set_monomers,
                        args.diameter, space)
        
        set_monomers.extend(h_monomers)

    # Write files
    json_file = model.JSONConfigOutputFile(args.output_filebase + '.json')
    json_file.write(set_monomers, args.box_len)


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
        'num_hexamers',
        type=int,
        help='Number of hexamers to place in system')
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
