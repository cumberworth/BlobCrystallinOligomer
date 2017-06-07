"""Classes for CG modeling of alpha crystallin monomers.

Some parallels to the classes used in the C++ implementation, but should not
assume the structure is the same. These classes are meant to be used in analysis
scripts and file conversions.

The altentaive is wrapping the relevant parts of the C++ version in cython,
which is probably more work.
"""

import json
import string

import numpy as np


class Monomer:
    """Coarse-grain AlphaB crystallin monomer.

    Contains ACD partilces and NTD particles . ACDs protrude from origin with NTD
    extending in positive z.
    """

    def __init__(self, acd_particles, ntd_particles, acd_radius, ntd_radius,
            index):
        self._acd_particles = acd_particles
        self._ntd_particles = ntd_particles
        self._acd_radius = acd_radius
        self._ntd_radius = ntd_radius
        self._index = index
        self._num_particles = len(acd_particles) + len(ntd_particles)

        self._cur_i = -1 # For iterating over contained particles
        self._num_particles = len(ntd_particles) + len(acd_particles)

    # Following two functions are iterator protocol to allow iteration over
    # particles
    def __iter__(self):
        return self

    def __next__(self):
        self._cur_i += 1
        if self._cur_i >= self._num_particles:
            self._cur_i = -1
            raise StopIteration
        else:
            return (self._acd_particles + self._ntd_particles)[self._cur_i]

    @property
    def index(self):
        return self._index

    @property
    def acd_radius(self):
        return self._acd_radius

    @property
    def ntd_radius(self):
        return self._ntd_radius

    @property
    def acd_particles(self):
        return self._acd_particles

    @property
    def ntd_particles(self):
        return self._ntd_particles
    
    @property
    def center(self):
        center = np.zeros(4)
        for particle in self:
            center += particle.pos

        center /= self._num_particles

        return center

    def apply_transformation(self, M):
        """Apply the transformation matrix to all position vectors."""

        for particle in self:
            particle.apply_transformation(M)


class SimpleParticle:
    """Spherical particle.

    Defaults to homogenous coordinates."""

    def __init__(self, index, domain, pos=None):
        self._index = index # Unique particle index
        self._domain = domain # Domain particle is representing (ACD or NTD)
        if pos is None:
            self.pos = np.array([0, 0, 0, 1.]) # Position

        self._label = 'SimpleParticle' # Type label

    @property
    def index(self):
        return self._index

    @property
    def type(self):
        return self._domain

    @property
    def domain(self):
        return self._domain

    @property
    def label(self):
        return self._label

    def apply_transformation(self, M):
        """Apply the transformation matrix to all position vectors."""

        self.pos = np.dot(M, self.pos)


class PatchyParticle(SimpleParticle):
    """Spherical particle with one patch."""

    def __init__(self, *args, patch_norm=None, **kwargs):
        super().__init__(*args, **kwargs)
        if patch_norm is None: # Unit vector normal to patch
            self._patch_norm = np.array([0, 0, 0, 1.])
        else:
            self._patch_norm = patch_norm

        self._label = 'PatchyParticle' # Type label

    @property
    def patch_norm(self):
        return self._patch_norm

    def apply_transformation(self, M):
        """Apply the transformation matrix to all position vectors."""

        super().apply_transformation(M)
        self._patch_norm = np.dot(M, self._patch_norm)


class OrientedPatchyParticle(PatchyParticle):
    """Spherical particle with one orientationally specific patch."""

    def __init__(self, *args, patch_orient=None, **kwargs):
        super().__init__(*args, **kwargs)
        if patch_orient is None: # Unit vector perpindicular to patch norm
            self._patch_orient = np.array([0, 0, 0, 1.])
        else:
            self._patch_orient = patch_orient

        self._label = 'OrientedPatchyParticle' # Type label

    @property
    def patch_orient(self):
        return self._patch_orient

    def apply_transformation(self, M):
        """Apply the transformation matrix to all position vectors."""

        super().apply_transformation(M)
        self._patch_orient = np.dot(M, self._patch_orient)


class PDBConfigOutputFile:
    """Configuration (positions only) in PDB format"""

    # Put spacers after field
    _atom_template = (
            'ATOM  '
            '{serial:>5} '
            ' {name:<3}' # The online specs explain how the space is used
            ' ' # Never use altLoc
            '{resName:<3} '
            '{chainID:1}'
            '{resSeq:>4}'
            '    ' # Never use iCode
            '{x:>8.3f}'
            '{y:>8.3f}'
            '{z:>8.3f}'
            '{occupancy:>6.2f}' +
            #'{tempFactor:>6.2f}'
            ' '*6 +
            ' '*12 +
            '{element:}'
            '  ' # Never use chage
    )

    _ter_template = (
            'TER   '
            '{serial:>5}' +
            ' '*6 +
            '{resName:<3} '
            '{chainID:1}'
            '{resSeq:>4}'
            ' '
    )

    def __init__(self, filename):
        self._filename = filename

    def write(self, monomers):
        """Write simple pdb file of coarse-grained alphaB crystallin model.

        Each monomer is treated as a residue, and the NTD and ACD as different
        atom types.
        """

        output_file = open(self._filename, 'w')
        for monomer in monomers:
            for particle in monomer:
                if particle.domain == 'ACD':
                    radius = monomer.acd_radius
                elif particle.domain == 'NTD':
                    radius = monomer.ntd_radius
                else:
                    radius = 1

                atom_fields = {
                        'serial': particle.index,
                        'name': particle.domain,
                        'resName': 'ABC',
                        'chainID': string.ascii_uppercase[monomer.index],
                        'resSeq': monomer.index,
                        'x': particle.pos[0],
                        'y': particle.pos[1],
                        'z': particle.pos[2],
                        'occupancy': radius,
                        #'tempFactor': '',
                        'element': 'CG'
                }
                line = self._atom_template.format(**atom_fields)
                output_file.write(line + '\n')

        ter_fields = {
                'serial': particle.index + 1,
                'resName': 'ABC',
                'chainID': string.ascii_uppercase[monomer.index],
                'resSeq': monomer.index,
        }
        line = self._ter_template.format(**ter_fields)
        output_file.write(line + '\n')


class JSONConfigOutputFile:
    """Configuration in PDB format"""

    def __init__(self, filename):
        self._filename = filename

    def write(self, monomers):
        """Write json format of coarse-grained alphaB crystallin model."""
        cgmonomer_json = {}
        cgmonomer_json['cgmonomer'] = {'config': []}
        for monomer in monomers:
            monomer_json = {}
            monomer_json['index'] = monomer.index
            monomer_json['particles'] = []
            for particle in monomer:
                particle_json = {}
                particle_json['index'] = particle.index
                particle_json['domain'] = particle.domain
                particle_json['type'] = particle.label
                particle_json['pos'] = particle.pos.tolist()[:3]
                if 'Patchy' in particle.label:
                    particle_json['patch_norm'] = particle.patch_norm.tolist()[:3]

                if 'Oriented' in particle.label:
                    particle_json['patch_orient'] = particle.patch_orient.tolist()[:3]

                monomer_json['particles'].append(particle_json)

            cgmonomer_json['config'].append(monomer_json)

        json.dump(cgmonomer_json, open(self._filename, 'w'), indent=4,
                    separators=(',', ': '))
