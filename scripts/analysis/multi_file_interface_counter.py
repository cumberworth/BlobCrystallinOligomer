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
    prim_means = []
    sec_means = []
    prim_stds = []
    sec_stds = []
    for param_i in range(args.params):
        prim_reps = []
        sec_reps = []
        for rep in range(args.reps):
            vtf_filename = '{}/{}_rep-{}-{}.vtf'.format(args.vtf_dir, args.filebase, rep, param_i)
            traj = VTFInputFile(vtf_filename)
            prim, sec = analysis.calc_interface_timeseries(sys, traj)
            prim_reps.append(np.mean(prim))
            sec_reps.append(np.mean(sec))

        prim_means.append(np.mean(prim_reps))
        prim_stds.append(np.std(prim_reps))
        sec_means.append(np.mean(sec_reps))
        sec_stds.append(np.std(sec_reps))

    np.savetxt('analysis/{}.prim'.format(args.filebase), np.transpose([prim_means, prim_stds]))
    np.savetxt('analysis/{}.sec'.format(args.filebase), np.transpose([sec_means, sec_stds]))

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            'sys_filename',
            type=str,
            help='System file')
    parser.add_argument(
            'filebase',
            type=str,
            help='File base for VTF files and output')
    parser.add_argument(
            'vtf_dir',
            type=str,
            help='Directory that VTF files are located')
    parser.add_argument(
            'params',
            type=int,
            help='Number of parameter sets')
    parser.add_argument(
            'reps',
            type=int,
            help='Number of reps')
    return parser.parse_args()


if __name__ == '__main__':
    main()
