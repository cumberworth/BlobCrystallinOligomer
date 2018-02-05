cimport numpy as np
import numpy as np

from config cimport Config
from config import Config
from ifile cimport VTFInputFile
from ifile import VTFInputFile
from monomer cimport Monomer
from monomer import Monomer

def calc_interface_timeseries(Config config, VTFInputFile trajfile):
    """Calculate number of primary and secondary interfaces per timestep"""
    cdef double cutoff = 2*config.get_radius() + 1
    cdef double blob_cutoff = 15
    cdef list primary_interfaces = []
    cdef list secondary_interfaces = []
    cdef list ntd_interfaces = []
    cdef int i
    cdef int num_configs = trajfile.num_configs
    for i in range(num_configs - 1):
        positions = trajfile.get_config_positions(i)
        config.update_config_positions(positions)

        # Calculate number of primary interfaces
        dists0 = config.calc_dist_pairs(0)
        interacts0 = dists0 < cutoff
        primary_interfaces.append(np.sum(interacts0))

        # Calculate number of secondary interfaces
        dists1 = config.calc_dist_pairs(1)
        dists2 = config.calc_dist_pairs(2)
        interacts1 = dists1 < cutoff
        interacts2 = dists2 < cutoff
        secondary_count = np.sum(np.logical_and(interacts1, interacts2))
        secondary_interfaces.append(secondary_count)

        # Calculate number of NTD interfaces
        dists4 = config.calc_dist_pairs(4)
        interacts4 = dists4 < blob_cutoff
        ntd_interfaces.append(np.sum(interacts4))

    return primary_interfaces, secondary_interfaces, ntd_interfaces
