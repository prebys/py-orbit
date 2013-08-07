"""
This module is a collimator node class for TEAPOT lattice
"""

import os
import math

# import the auxiliary classes
from orbit.utils import orbitFinalize, NamedObject, ParamsDictObject

# import general accelerator elements and lattice
from orbit.lattice import AccNode, AccActionsContainer, AccNodeBunchTracker

# import teapot drift class
from orbit.teapot import DriftTEAPOT

# import Collimator class
from collimator import Collimator

class TeapotFieldTrackerNode(DriftTEAPOT):
    """ 
    The collimator node class for TEAPOT lattice
    """
    def __init__(self, order, bx, by, ax, ay, ex, epx, l, zi, zf, ds, niters,
                  resid, xrefi, yrefi, eulerai, eulerbi, eulergi, apflag, b):
        """
        Constructor. Creates the FieldTracker TEAPOT element.
        """
        DriftTEAPOT.__init__(self,name)
        self.fieldtracker = FieldTracker(order, bx, by, ax, ay, ex, epx, l, zi, zf, ds, niters,
                  resid, xrefi, yrefi, eulerai, eulerbi, eulergi, apflag, b)
        self.setType("fieldtracker teapot")
        self.setLength(l)

    def track(self, paramsDict):
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        self.FieldTracker.trackbunch(bunch)