import figs

from pyraf import iraf
from iraf import stsdas,dither,slitless,axe

import pyfits

import drizzlepac
from drizzlepac import astrodrizzle

import os

import numpy as np

INDEF = iraf.INDEF
no = iraf.no
yes = iraf.yes

def astrodrizzle_run():
    """
    Performs an astrodrizzle run on a set of input images
    """


