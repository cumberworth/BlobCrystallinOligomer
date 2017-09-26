"""Classes for CG modeling of alpha crystallin monomers.

Some parallels to the classes used in the C++ implementation, but should not
assume the structure is the same. These classes are meant to be used in analysis
scripts and file conversions.

The alternative is wrapping the relevant parts of the C++ version in cython,
which is probably more work.
"""

import json
import string

import numpy as np
import transforms3d.affines as affines
import transforms3d.zooms as zooms


def hard_sphere_overlap(new_monomer, monomers, diameter, space):
    """Check for overlap between new hard sphere and system of hard spheres"""
    overlap = False
    new_particles = new_monomer.particles
    for old_monomer in monomers:
        for i in range(len(new_particles)):
            new_pos = new_particles[i].pos
            old_particles = old_monomer.particles
            for j in range(i, len(old_particles)):
                old_pos = old_particles[j].pos
                if space.calc_dist(new_pos, old_pos) < diameter:
                    overlap = True
                    return overlap

    return overlap


class CuboidPBC:
    """Cuboid periodic boundary conditions"""

    def __init__(self, box_len):
        self.r = box_len / 2

    def calc_diff(self, pos1, pos2):
        diff = np.zeros(3)
        for i in range(3):
            comp_diff = pos1[i] - pos2[i]
            if comp_diff > self.r:
                comp_diff = -2*self.r + comp_diff
            elif comp_diff < -self.r:
                comp_diff = 2*self.r + comp_diff;

            diff[i] = comp_diff;

        return diff;

    def calc_dist(self, pos1, pos2):
        return np.linalg.norm(self.calc_diff(pos1, pos2))

    def wrap(self, pos):
        for i in range(3):
            if pos[i] > self.r:
                pos[i] = -2*self.r + pos[i]
            elif pos[i] < -self.r:
                pos[i] = 2*self.r + pos[i]

        return pos


class Monomer:
    """General monomer class"""

    def __init__(self, particles, index):
        self._cur_i = -1 # For iterating over contained particles
        self._num_particles = len(particles)
        self._particles = particles

    # Following two functions are iterator protocol to allow iteration over
    # particles.
    def __iter__(self):
        return self

    def __next__(self):
        self._cur_i += 1
        if self._cur_i >= self._num_particles:
            self._cur_i = -1
            raise StopIteration
        else:
            return (self._particles)[self._cur_i]

    @property
    def center(self):
        center = np.zeros(4)
        for particle in self:
            center += particle.pos

        center /= self._num_particles

        return center

    @property
    def particles(self):
        return self._particles

    def apply_transformation(self, M):
        """Apply the transformation matrix to all position vectors."""

        for particle in self:
            particle.apply_transformation(M)


class SingleParticleMonomer(Monomer):
    """Single particle monomer."""

    def __init__(self, particles, radius, index):
        self._particles = particles
        self._radius = radius
        self._index = index
        self._num_particles = len(particles)

        self._cur_i = -1 # For iterating over contained particles

    @property
    def index(self):
        return self._index

    @property
    def radius(self):
        return self._radius

    @property
    def pos(self):
        return self._particles[0].pos

    @pos.setter
    def pos(self, pos):
        self._particles[0].pos = pos


class AlphaBMonomer(Monomer):
    """Coarse-grain AlphaB crystallin monomer.

    Contains ACD particles and NTD particles. ACDs protrude from origin with NTD
    extending in positive z.
    """

    def __init__(self, acd_particles, ntd_particles, acd_radius, ntd_radius,
            blob_particles, index):
        self._acd_particles = acd_particles
        self._ntd_particles = ntd_particles
        self._blob_particles = blob_particles
        self._particles = acd_particles + ntd_particles + blob_particles
        super().__init__(self._particles, index)
        self._acd_radius = acd_radius
        self._ntd_radius = ntd_radius
        self._index = index

        self._num_acd_particles = len(acd_particles)
        self._num_ntd_particles = len(ntd_particles)
        self._num_blob_particles = len(blob_particles)
        self._num_particles = (self._num_acd_particles + self._num_ntd_particles +
                self._num_blob_particles)

    @property
    def index(self):
        return self._index

    @index.setter
    def index(self, i):
        self._index = i

    @property
    def radius(self):
        return self._acd_radius

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
    def blob_particles(self):
        return self._blob_particles
    
    @property
    def acd_center(self):
        center = np.zeros(4)
        for particle in self.acd_particles:
            center += particle.pos

        center /= self._num_acd_particles

        return center


class SimpleParticle:
    """Spherical particle.

    Defaults to homogeneous coordinates."""

    def __init__(self, index, domain, type, pos=None):
        self._index = index # Unique particle index
        self._domain = domain # Domain particle is representing (ACD or NTD)
        self._type = type
        if pos is None:
            self.pos = np.array([0, 0, 0, 1.]) # Position

        self._form = 'SimpleParticle' # Type form

    @property
    def index(self):
        return self._index

    @property
    def type(self):
        return self._type

    @property
    def domain(self):
        return self._domain

    @property
    def form(self):
        return self._form

    def apply_transformation(self, M):
        """Apply the transformation matrix to all position vectors."""

        self.pos = np.dot(M, self.pos)


class PatchyParticle(SimpleParticle):
    """Spherical particle with one patch."""

    def __init__(self, *args, patch_norm=None, **kwargs):
        super().__init__(*args, **kwargs)
        if patch_norm is None: # Unit vector normal to patch
            self._patch_norm = np.array([0., 0., 0.])
        else:
            self._patch_norm = patch_norm

        self._form = 'PatchyParticle' # Type form

    @property
    def patch_norm(self):
        return self._patch_norm

    def apply_transformation(self, M):
        """Apply the transformation matrix to all position vectors."""

        super().apply_transformation(M)
        translation, rotation, zoom, shear = affines.decompose(M)
        transform = affines.compose(np.zeros(3), rotation, zoom, shear)[:3, :3]
        self._patch_norm = np.dot(transform, self._patch_norm)


class OrientedPatchyParticle(PatchyParticle):
    """Spherical particle with one orientationally specific patch."""

    def __init__(self, *args, patch_orient=None, **kwargs):
        super().__init__(*args, **kwargs)
        if patch_orient is None: # Unit vector perpindicular to patch norm
            self._patch_orient = np.array([0., 0., 0.])
        else:
            self._patch_orient = patch_orient

        self._form = 'OrientedPatchyParticle' # Type form

    @property
    def patch_orient(self):
        return self._patch_orient

    def apply_transformation(self, M):
        """Apply the transformation matrix to all position vectors."""

        super().apply_transformation(M)
        translation, rotation, zoom, shear = affines.decompose(M)
        transform = affines.compose(np.zeros(3), rotation, zoom, shear)[:3, :3]
        self._patch_orient = np.dot(transform, self._patch_orient)


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
                    radius = monomer.radius

                atom_fields = {
                        'serial': particle.index,
                        'name': particle.domain,
                        'resName': 'ABC',
                        #'chainID': string.ascii_uppercase[monomer.index],
                        'chainID': 'A',
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
                if type(particle) in [PatchyParticle, OrientedPatchyParticle]:
                    atom_fields = {
                            'serial': particle.index,
                            'name': 'PAT',
                            'resName': 'ABC',
                            #'chainID': string.ascii_uppercase[monomer.index],
                            'chainID': 'A',
                            'resSeq': monomer.index,
                            'x': particle.pos[0] + particle.patch_norm[0],
                            'y': particle.pos[1] + particle.patch_norm[1],
                            'z': particle.pos[2] + particle.patch_norm[2],
                            'occupancy': 1,
                            #'tempFactor': '',
                            'element': 'CG'
                    }
                    line = self._atom_template.format(**atom_fields)
                    output_file.write(line + '\n')
                if type(particle) is OrientedPatchyParticle:
                    atom_fields = {
                            'serial': particle.index,
                            'name': 'PAT',
                            'resName': 'ABC',
                            #'chainID': string.ascii_uppercase[monomer.index],
                            'chainID': 'A',
                            'resSeq': monomer.index,
                            'x': particle.pos[0] + particle.patch_orient[0],
                            'y': particle.pos[1] + particle.patch_orient[1],
                            'z': particle.pos[2] + particle.patch_orient[2],
                            'occupancy': 1,
                            #'tempFactor': '',
                            'element': 'CG'
                    }
                    line = self._atom_template.format(**atom_fields)
                    output_file.write(line + '\n')

        ter_fields = {
                'serial': particle.index + 1,
                'resName': 'ABC',
                #'chainID': string.ascii_uppercase[monomer.index],
                'chainID': 'A',
                'resSeq': monomer.index,
        }
        line = self._ter_template.format(**ter_fields)
        output_file.write(line + '\n')


class JSONConfigOutputFile:
    """Configuration in JSON format"""

    def __init__(self, filename):
        self._filename = filename

    def write(self, monomers, box_len):
        """Write json format of coarse-grained alphaB crystallin model."""
        cgmonomer_json = {}
        cgmonomer_json['cgmonomer'] = {'config': []}
        cgmonomer_json['cgmonomer']['radius'] = monomers[0].radius
        cgmonomer_json['cgmonomer']['box_len'] = box_len
        for monomer in monomers:
            monomer_json = {}
            monomer_json['index'] = monomer.index
            monomer_json['particles'] = []
            for particle in monomer:
                particle_json = {}
                particle_json['index'] = particle.index
                particle_json['domain'] = particle.domain
                particle_json['form'] = particle.form
                particle_json['type'] = particle.type
                particle_json['pos'] = particle.pos.tolist()[:3]
                if 'Patchy' in particle.form:
                    particle_json['patch_norm'] = particle.patch_norm.tolist()[:3]

                if 'Oriented' in particle.form:
                    particle_json['patch_orient'] = particle.patch_orient.tolist()[:3]

                monomer_json['particles'].append(particle_json)

            cgmonomer_json['cgmonomer']['config'].append(monomer_json)

        json.dump(cgmonomer_json, open(self._filename, 'w'), indent=4,
                    separators=(',', ': '))
