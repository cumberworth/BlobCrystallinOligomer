#!/usr/env python3

"""Count the number of each type of possible interface in alphaB simulations"""

import argparse

import numpy as np

from crystallinpy import analysis
from crystallinpy.config import Config
from crystallinpy.ifile import VTFInputFile

def main():
    args = parse_args()
    sys = Config(args.sys_filename)
    traj = VTFInputFile(args.vtf_filename)
    primaries, secondaries = analysis.calc_interface_timeseries(sys, traj)
    primaries = np.array(primaries)
    secondaries = np.array(secondaries)
    print((primaries/secondaries).mean())

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            'sys_filename',
            type=str,
            help='System file')
    parser.add_argument(
            'vtf_filename',
            type=str,
            help='VTF file')
    return parser.parse_args()


if __name__ == '__main__':
    main()
