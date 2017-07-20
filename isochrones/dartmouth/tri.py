#!/usr/bin/env python

import sys, os
import pandas as pd
import numpy as np
try:
    import cPickle as pickle
except ImportError:
    import pickle
from scipy.interpolate import LinearNDInterpolator as interpnd

from ..config import ISOCHRONES

from .grid import DartmouthModelGrid

def tri_filename(afe="afep0", y=""):
    afe_suffix = "" if (afe=="afep0") else "_{}".format(afe)
    y_suffix = "" if (y=="") else "_{}".format(y)
    filename=os.path.join(ISOCHRONES,'dartmouth{}{}.tri'.format(afe_suffix,y_suffix))
    return filename

def write_tri(filename=None,afe="afep0",y='',bands=['g']):
    if filename is None:
        filename = tri_filename(afe, y)
    print "Writing afe={} y={} to {}".format(afe, y, filename)
    df = DartmouthModelGrid(bands, afe=afe, y=y).df
    N = len(df)
    pts = np.zeros((N,3))
    pts[:,0] = np.array(df['MMo'])
    pts[:,1] = np.array(df['age'])
    pts[:,2] = np.array(df['feh'])
    mags = np.array(df[bands[0]])

    magfn = interpnd(pts,mags)

    with open(filename,'wb') as f:
        pickle.dump(magfn.tri,f)
