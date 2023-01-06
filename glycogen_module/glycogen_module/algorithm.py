
from glycogen_module.core import GlycogenStructure
import random
import math


def gillespie_step(structure: GlycogenStructure, k_gs: float = 1.0, k_gp: float = 1.0, k_gbe: float = 1.0, k_gde: float = 1.0):
    ''' This functions takes concentrations of the enzymes and the structure info of a glycogen granules and
    return what is the next reaction to occurs and which time has been spent. (Following a gillespie algorithm)
    '''
    # propensity assuming mass action kinetics
    h_gs = k_gs*structure.gs*len(structure.find_chains_for_gs())
    h_gp = k_gp*structure.gp*len(structure.find_chains_for_gp())
    h_gbe = k_gbe*structure.gbe*len(structure.find_chains_for_gbe())
    h_gde = k_gde*structure.gde*len(structure.find_chains_for_gde())

    a = h_gs + h_gp + h_gbe + h_gde

    if a == 0:
        return "no reaction can be proceed, all propensities are zero", 0
    r2 = random.uniform(0, a)
    r1 = random.uniform(0, 1)

    d_t = (1/a)*math.log(1/r1)
    if r2 < h_gs:
        return "Act_gs()", d_t
    if r2 >= h_gs and r2 < h_gs + h_gp:
        return "Act_gp()", d_t
    if r2 >= h_gs + h_gp and r2 < h_gs + h_gp + h_gbe:
        return "Act_gbe()", d_t
    if r2 >= h_gs + h_gp + h_gbe and r2 < h_gs + h_gp + h_gbe+h_gde:
        return "Act_gde()", d_t
