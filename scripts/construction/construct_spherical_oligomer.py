#!/usr/env python3

"""Create coordinates for coarse-grained alphaB-crystallin 24-mer

Output in pdb and json format
"""

import argparse
import math
import sys
import copy

import numpy as np
import scipy.optimize as opt
import transforms3d as trans
import transforms3d.reflections as refls

from crystallinpy import model


def main():
    args = parse_cl()
    hexamer_edge_length = 65
    tet_arm_length = args.arm_to_edge * hexamer_edge_length

    # Calculate radius of spheres
    acd_radius = hexamer_edge_length / (4*args.num_acd_spheres + 1)
    ntd_radius = acd_radius

    # Calculate ACD-NTD angle
    ntd_length = 2*args.num_ntd_spheres * ntd_radius
#    min_angle = 0
#    max_angle = math.pi - math.acos(-1/3) # This only allows ntds that bulge, unlike in Jehle model
    min_angle = 0.62
    max_angle = 0.63
    constants = (acd_radius, ntd_radius, args.num_acd_spheres,
            args.num_ntd_spheres, tet_arm_length)
    acd_ntd_angle = opt.brentq(construct_oligomer_check_ntd, min_angle,
            max_angle, args=constants)

    # Calculate blob angles
    #guess = [0, 0]
    guess = [-0.319, -0.534]
    constants = (acd_ntd_angle, acd_radius, ntd_radius, args.num_acd_spheres,
            args.num_ntd_spheres, tet_arm_length)
    blob_angle1, blob_angle2 = opt.fmin(construct_oligomer_check_blobs, guess,
            args=constants)

    # Build oligomer and write to file
    monomers = construct_oligomer(acd_radius, ntd_radius,
            args.num_acd_spheres, args.num_ntd_spheres, tet_arm_length,
            acd_ntd_angle, blob_angle1, blob_angle2)
    monomers = modify_structure(monomers)
    pdb_file = model.PDBConfigOutputFile(args.output_filebase + '.pdb')
    pdb_file.write(monomers)
    json_file = model.JSONConfigOutputFile(args.output_filebase + '.json')
    json_file.write(monomers, args.box_len)
    print(acd_ntd_angle, blob_angle1, blob_angle2)

    # Simple tetrahedral with same arm length (for aligning in VMD)
    #tetrahedral = construct_tetrahedral(tet_arm_length)
    #pdb_file = model.PDBConfigOutputFile(args.output_filebase + '_centers.pdb')
    #pdb_file.write(tetrahedral)


def construct_oligomer_check_ntd(acd_ntd_angle, acd_radius, ntd_radius,
            num_acd_spheres, num_ntd_spheres, tet_arm_length):
    """Build construct and output distance between adjacent NTDs

    The distance is dependent on the ACD-NTD angle.
    """

    # Build oligomer
    monomers = construct_oligomer(acd_radius, ntd_radius,
            num_acd_spheres, num_ntd_spheres, tet_arm_length, acd_ntd_angle,
            0, 0)

    # Select two particles that should be touching at the NTD-NTD nexus
    particle1 = monomers[0].ntd_particles[-1]
    particle2 = monomers[23].ntd_particles[-1]
    dist = np.linalg.norm(particle2.pos - particle1.pos)
    res = dist - 2*ntd_radius

    return res


def construct_oligomer_check_blobs(angles, acd_ntd_angle, acd_radius, ntd_radius,
            num_acd_spheres, num_ntd_spheres, tet_arm_length):
    """Build construct and output distance between overlapping blobs 

    The distance is dependent on the ACD-NTD angle.
    """
    blob_angle1, blob_angle2 = angles

    # Build oligomer
    monomers = construct_oligomer(acd_radius, ntd_radius,
            num_acd_spheres, num_ntd_spheres, tet_arm_length, acd_ntd_angle,
            blob_angle1, blob_angle2)

    blob1 = monomers[0].blob_particles[0]
    blob2 = monomers[3].blob_particles[0]
#    blob3 = monomers[7].blob_particles[0]
#    blob4 = monomers[10].blob_particles[0]
#    blob5 = monomers[20].blob_particles[0]
    blob6 = monomers[23].blob_particles[0]
    dist1 = np.linalg.norm(blob1.pos - blob2.pos)
    dist2 = np.linalg.norm(blob1.pos - blob6.pos)

    return dist1 + dist2


def construct_monomers(acd_radius, ntd_radius, num_acd_spheres,
        num_ntd_spheres, num_monomers, acd_ntd_angle, blob_angle1, blob_angle2):
    """Construct monomers given number of particles for ACD and NTD.

    Internaly are some descisions about the kind of particles involved (i.e.,
    how many patches on each particle are where and how they oriented.). So
    probably if the number of particles for each domain changes, then will
    probably want to reconsider those patch choices.
    """
    monomers = []
    for monomer_i in range(num_monomers):
        acd_particles = []
        particle_i = 0
        ptype = 0
        for j in range(num_acd_spheres):

            # Note patch orientation is tied to the way the particles are later
            # oriented
            if j == 0:
                particle = model.OrientedPatchyParticle(particle_i, 'ACD',
                        ptype, patch_norm=np.array([0., -1., 0.]),
                        patch_orient=np.array([1., 0., 0.]))
            else:
                particle = model.PatchyParticle(particle_i, 'ACD',
                        ptype, patch_norm=np.array([-1., 0., 0.]))

            particle_i += 1
            acd_particles.append(particle)
            ptype += 1

        ntd_particles = []
        for j in range(num_ntd_spheres):
            if j == 0:
                particle = model.PatchyParticle(particle_i, 'NTD',
                        ptype, patch_norm=[-1., 0., 0.])
            else:
                particle = model.SimpleParticle(particle_i, 'NTD', ptype)

            particle_i += 1
            ntd_particles.append(particle)
            ptype += 1

        blob_particles = [model.SimpleParticle(particle_i, 'BLB', ptype)]

        monomer = model.AlphaBMonomer(acd_particles, ntd_particles, acd_radius,
                ntd_radius, blob_particles, monomer_i)
        monomer = orient_monomer(monomer, acd_ntd_angle, blob_angle1, blob_angle2)
        monomers.append(monomer)

    return monomers


def orient_monomer(monomer, acd_ntd_angle, blob_angle1, blob_angle2):
    """Orient a monomer to have given ACD-NTD angle"""

    # Place particles along y axis with edge of first sphere touching origin
    p_i = 0
    for particle in monomer.acd_particles:
        particle.pos[1] = 2*p_i*monomer.acd_radius + monomer.acd_radius
        p_i += 1

    # Create rotation transformation to rotate NTDs
    # Do as a negative rotation around the z axis followed by a positive rotation
    # around the axis perpindicular to the direction of the NTD in the xy plane
    z_axis = np.array([0, 0, 1])
    last_acd_pos = monomer.acd_particles[-1].pos
    z_axis_angle = -math.pi/6
    acd_z_rotation = trans.axangles.axangle2aff(z_axis, z_axis_angle,
            last_acd_pos)

    # Rotate last ACD for initial patch placement
    back_z_axis_angle = math.pi/6
    acd_z_rotation = trans.axangles.axangle2aff(z_axis, z_axis_angle,
            last_acd_pos)
    monomer.acd_particles[1].apply_transformation(acd_z_rotation)

    # Make the next axis
    z_rotation = trans.axangles.axangle2aff(z_axis, z_axis_angle)
    ntd_xy_axis = np.dot(z_rotation, np.array([1, 0, 0, 1]))[:3]
    acd_ntdxy_rotation = trans.axangles.axangle2aff(ntd_xy_axis,
            acd_ntd_angle, last_acd_pos)

    ntd_transform = acd_ntdxy_rotation.dot(acd_z_rotation)

    # Continue placing particles along y axis then rotate
    for particle in monomer.ntd_particles:
        pos = particle.pos
        pos[1] += 2*(p_i - 1)*monomer.ntd_radius + 2*monomer.acd_radius + monomer.ntd_radius
        particle.pos = pos
        particle.apply_transformation(ntd_transform)
        p_i += 1

    # Place the blob particle along y axis then rotate
    last_ntd_pos = monomer.ntd_particles[-1].pos
    blob_rotation1 = trans.axangles.axangle2aff([1, 0, 0],
            blob_angle1, last_ntd_pos)
    blob_rotation2 = trans.axangles.axangle2aff([0, 1, 0],
            blob_angle2, last_ntd_pos)
    blob_transform = blob_rotation1.dot(blob_rotation2)
    particle = monomer.blob_particles[0]
    pos = particle.pos
    pos[1] += (2*(p_i - 1)*monomer.ntd_radius + 2*monomer.acd_radius +
            1.5*monomer.ntd_radius)
    particle.pos = pos
    particle.apply_transformation(ntd_transform)
    particle.apply_transformation(blob_transform)

    # Prepare patches
    #z_axis_angle = -math.pi/2
    #y_axis = np.array([0, 1, 0])
    #x_comp = monomer.ntd_particles[0].pos[0]
    #z_comp = monomer.ntd_particles[0].pos[2]
    #xz_proj = np.array([x_comp, 0, z_comp])
    #xz_proj_norm = np.linalg.norm(xz_proj)
    #y_axis_angle = -math.acos(np.dot(z_axis, xz_proj)/xz_proj_norm)
    #for particle in monomer:
    #    pos = particle.pos
    #    z_rotation = trans.axangles.axangle2aff(z_axis, z_axis_angle, pos)
    #    y_rotation = trans.axangles.axangle2aff(y_axis, y_axis_angle, pos)
    #    particle.apply_transformation(z_rotation)
    #    particle.apply_transformation(y_rotation)

    return monomer


def construct_dimer(monomer1, monomer2):
    """Construct dimer.

    Orient two monomers in bound position. Move second relative to first.
    Assumes they are currently in the same position and the ACDs in the xy
    plane.
    """
    z_axis = np.array([0, 0, 1])
    z_axis_angle = math.pi
    z_rotation = trans.axangles.axangle2aff(z_axis, z_axis_angle)

    # Need to rotate around y-axis to get spheres in place
    # Get angle by taking dot product of z axis and zx plane projection of the
    # ntd position, dividing by the zx plane projection norm, and taking the
    # inverse cos of the result
    #y_axis = np.array([0, 1, 0])
    #x_comp = monomer1.ntd_particles[0].pos[0]
    #z_comp = monomer1.ntd_particles[0].pos[2]
    #xz_proj = np.array([x_comp, 0, z_comp])
    #xz_proj_norm = np.linalg.norm(xz_proj)
    #y_axis_angle = 2*math.acos(np.dot(z_axis, xz_proj)/xz_proj_norm)
    #y_rotation = trans.axangles.axangle2aff(y_axis, y_axis_angle)

    #rotation = np.dot(y_rotation, z_rotation)
    #monomer2.apply_transformation(rotation)

    reflection = refls.rfnorm2aff([0, 1, 0])
    monomer2.apply_transformation(reflection)

    return (monomer1, monomer2)


def construct_hexamer(dimer1, dimer2, dimer3):
    """Construct hexamer.

    Orient three dimers in bound position. One corner on y axis.
    Assumes they are currently in the same position and the acds are in the
    y axis in the xy plane.
    """

    # Rotate 2 and 3 pi around z
    z_axis = np.array([0, 0, 1])
    z_axis_angle = math.pi
    z_rotation = trans.axangles.axangle2aff(z_axis, z_axis_angle)
    for dimer in [dimer2, dimer3]:
        for monomer in dimer:
            monomer.apply_transformation(z_rotation)

    # First rotate all dimers to be 30 degrees off y axis
    axis = np.array([0, 0, 1])
    angle30 = math.pi/6

    # Rotate around the corner of the to-be-formed triangle, which is one sphere
    # out from the last ACD domain
    origin1 = dimer1[0].acd_particles[-1].pos.copy()
    origin1[1] += 2*dimer1[0].acd_radius
    am1 = trans.axangles.axangle2aff(axis, angle30, origin1)
    for dimer in [dimer1, dimer2, dimer3]:
        for monomer in dimer:
            monomer.apply_transformation(am1)

    # Rotate second dimer to be in position across y-axis
    angle60 = math.pi/3
    am2 = trans.axangles.axangle2aff(axis, -angle60, origin1)
    for monomer in dimer2:
        monomer.apply_transformation(am2)

    # Rotate third dimer to form bottom of triangle
    origin2 = dimer1[1].acd_particles[-1].pos.copy()
    origin2 += origin2 - dimer1[1].acd_particles[0].pos
    am3 = trans.axangles.axangle2aff(axis, angle60, origin2)
    for monomer in dimer3:
        monomer.apply_transformation(am3)

    # Center triangle on origin
    coor_center = 0
    for dimer in [dimer1, dimer2, dimer3]:
        for monomer in dimer:
            coor_center += monomer.acd_center

    coor_center /= 6
    af = trans.affines.compose(-coor_center[:3], np.eye(3), np.ones(3))
    for dimer in [dimer1, dimer2, dimer3]:
        for monomer in dimer:
            monomer.apply_transformation(af)

    return (dimer1, dimer2, dimer3)


def construct_tetrahedral(arm_length):
    """Construct tetrahedral with atoms at center of hexamers"""

    arm = np.array([0, 0, -arm_length])
    axisx = np.array([1, 0, 0])
    translate_down_y = trans.affines.compose(arm, np.eye(3), np.ones(3))
    angleTet = math.acos(-1/3.)
    rotate_tet_x = trans.axangles.axangle2aff(axisx, angleTet)
    angle60 = math.pi/3
    axisz = np.array([0, 0, 1])
    rotate_60_z = trans.axangles.axangle2aff(axisz, angle60)
    angle120 = 2*math.pi/3
    rotate_120_z = trans.axangles.axangle2aff(axisz, angle120)

    trans_1 = translate_down_y
    trans_2 = rotate_60_z.dot(rotate_tet_x).dot(trans_1)
    trans_3 = rotate_120_z.dot(trans_2)
    trans_4 = rotate_120_z.dot(trans_3)
    trans_mons = [trans_1, trans_2, trans_3, trans_4]
    monomers = []
    for trans_mon, i in zip(trans_mons, range(4)):
        particle = model.SimpleParticle(i, 'TET')
        monomer = model.AlphaBMonomer([particle], [], 1, 1, i)
        monomer.apply_transformation(trans_mon)
        monomers.append(monomer)

    return monomers


def construct_tetrahedral_hexamers(hexamers, arm_length):
    """Construct 24mer (tetraheral hexamers).

    Orient four hexamers in bound position.
    """
    arm = np.array([0, 0, -arm_length])
    translate_down_y = trans.affines.compose(arm, np.eye(3), np.ones(3))
    for hexamer in hexamers:
        for dimer in hexamer:
            for monomer in dimer:
                monomer.apply_transformation(translate_down_y)

    axisx = np.array([1, 0, 0])

    angleTet = math.acos(-1/3.)
    rotate_tet_x = trans.axangles.axangle2aff(axisx, angleTet)

    axisz = np.array([0, 0, 1])
    angle60 = math.pi/3
    rotate_60_z = trans.axangles.axangle2aff(axisz, angle60)
    angle120 = 2*math.pi/3
    rotate_120_z = trans.axangles.axangle2aff(axisz, angle120)

    trans_hex2 = rotate_60_z.dot(rotate_tet_x)
    trans_hex3 = rotate_120_z.dot(trans_hex2)
    trans_hex4 = rotate_120_z.dot(trans_hex3)
    trans_hexes = [trans_hex2, trans_hex3, trans_hex4]
    for hexamer, trans_hex in zip(hexamers[1:], trans_hexes):
        for dimer in hexamer:
            for monomer in dimer:
                monomer.apply_transformation(trans_hex)

    return hexamers


def construct_oligomer(acd_radius, ntd_radius, num_acd_spheres,
        num_ntd_spheres, tet_arm_length, acd_ntd_angle, blob_angle1, blob_angle2):

    monomers = construct_monomers(acd_radius, ntd_radius, num_acd_spheres,
            num_ntd_spheres, 24, acd_ntd_angle, blob_angle1, blob_angle2)
    dimers = []
    for i in range(0, len(monomers), 2):
        monomer1 = monomers[i]
        monomer2 = monomers[i + 1]
        dimer = construct_dimer(monomer1, monomer2)
        dimers.append(dimer)

    hexamers = []
    for i in range(0, len(dimers), 3):
        hexamer = construct_hexamer(dimers[i], dimers[i + 1], dimers[i + 2])
        hexamers.append(hexamer)

    construct_tetrahedral_hexamers(hexamers, tet_arm_length)

    return monomers


def modify_structure(monomers):
    # Add another dimer
    m1 = copy.deepcopy(monomers[0])
    m1.index = 24
    m2 = copy.deepcopy(monomers[1])
    m2.index = 25
    rot_c = m1.blob_particles[0].pos
    rot_axis = [1, 0, 0]
    rot = trans.axangles.axangle2aff(rot_axis, math.pi,
            rot_c)
    m1.apply_transformation(rot)
    m2.apply_transformation(rot)
    monomers.extend([m1, m2])

    return monomers


def parse_cl():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'num_acd_spheres',
        type=int,
        help="Number of spheres for one monomer's ACD")
    parser.add_argument(
        "num_ntd_spheres",
        type=int,
        help="Number of spheres for one monomer's NTD")
    parser.add_argument(
        "arm_to_edge",
        type=float,
        help="Tetrahedral arm to triangle edge ratio")
    parser.add_argument(
        "box_len",
        type=float,
        help="Box length")
    parser.add_argument(
        "output_filebase",
        type=str,
        help="Output filebase")

    return parser.parse_args()


if __name__ == '__main__':
    main()
