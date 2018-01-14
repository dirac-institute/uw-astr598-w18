"""Data Lab helpers for (local) positional cross-matching."""

from __future__ import print_function

__authors__ = 'Robert Nikutta <nikutta@noao.edu>, NOAO Data Lab Team <datalab@noao.edu>'
__version__ = '20171219'

from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import pylab as plt


def xmatch(ra1,dec1,ra2,dec2,maxdist=None,units='deg',method='astropy',**kwargs):

    """Cross-match two sets of ra & dec coordinates locally (i.e. all coordinates are in RAM).

    The function will search for counterparts of ra1/dec1 coordinates
    in the in ra2/dec2 coordinate set, i.e. one can consider ra2/dec2
    to be the catalog that will be searched.

    Parameters
    ----------
    ra1, dec1: 1-d array-like sequences
        RA and declination of first coordinate set, in units of `units`

    ra2, dec2: 1-d array-like sequences
        RA and declination of second coordinate set, in units of `units`

    maxdist : float or None
        If not `None`, then it is the maximum angular distance (in
        units of `units`) to be considered. All distances greater than
        that will be considered non-matches. If `None`, then all
        ra1/dec1 will have matches in ra2/dec2.

    units : str
        Units of `ra1`, `dec1`, `ra2`, `dec2`. Default: 'deg' (decimal degrees).

    method : str
        Currently only astropy's :func:`match_to_catalog_sky()` method
        is supported, i.e. the default 'astropy'.

    Other Parameters
    ----------------
    nthneighbor : int, optional
        If ``method='astropy'``. Which closest neighbor to search for.
        Typically ``1`` is desired here, as that is correct for
        matching one set of coordinates to another. The next likely
        use case is ``2``, for matching a coordinate catalog against
        *itself* (``1`` is inappropriate because each point will find
        itself as the closest match).

    Returns
    -------
    idx : 1-d array
        Index values of the ra1/dec1 counterparts found in
        ra2/dec2. Thus ra2[idx], dec2[idx] will select from the
        ra2/dec2 catalog the matched counterparts of the ra1/dec1
        coordinate pairs.

        If `maxdist` was not `None` but a number instead, then 'idx'
        only contains the objects matched up to the `maxdist` radius.

    dist2d : 1-d array
        The angular distances of the matches found in the ra2/dec2
        catalog. In units of `units`.

        If `maxdist` was not `None` but a number instead, then
        'dist2d' only contains the objects matched up to the `maxdist`
        radius.

    """

    # turn coordinate sequences into 1d arrays
    ra1 = _arrayify(ra1)
    dec1 = _arrayify(dec1)
    ra2 = _arrayify(ra2)
    dec2 = _arrayify(dec2)

    if method == 'astropy':
        unit = getattr(u,units)
        c1 = SkyCoord(ra=ra1*unit, dec=dec1*unit)
        c2 = SkyCoord(ra=ra2*unit, dec=dec2*unit)
        idx, dist2d, dist3d = c1.match_to_catalog_sky(c2,**kwargs)

        if maxdist is not None:
            sel = (dist2d <= maxdist*unit)
            idx = idx[sel]
            dist2d = dist2d[sel]
            
    elif method == 'q3cpy':  # serialize Adam's code in Python first
        raise NotImplementedError("Method '%s' not yet implemented." % method)
            
    else:
        raise Exception("'%s' is not a valid method." % method)

    return idx, dist2d


def _arrayify(a):

    """Turn seq into a 1-d numpy array."""
    
    try:
        arr = np.array(a)
    except:
        raise
    
    return arr



# Diagnostics

def make_catalog(n,min,max):
    ra = np.random.uniform(min,max,n)
    dec = np.random.uniform(min,max,n)
    return ra,dec

def test():

    from time import time

    sizes = np.logspace(2,7,6,dtype=int)
    catalogs1 = [make_catalog(n,10,30) for n in sizes]
    catalogs2 = [make_catalog(n,20,40) for n in sizes]
    
    for j1,c1 in enumerate(catalogs1):
        ra1, dec1 = c1
        for j2,c2 in enumerate(catalogs2):
            ra2, dec2 = c2
            start = time()
            idx, dist = crossmatch.xmatch_local(ra1,dec1,ra2,dec2,units='arcsec',maxdist=None)
            stopp = time()
            delta = stopp-start
            times[j1,j2] = delta
            print(j1,j2,ra1.size,ra2.size,delta,'\n')


def plot_time_vs_catalogsizes(times,extent=(2,7)):

    plt.clf()
    im = plt.imshow(np.log10(times.T),origin='lower',cmap=matplotlib.cm.jet,interpolation='bicubic')
    tickmarks = ('2','3','4','5','6','7')
    plt.xticks(range(6),tickmarks)
    plt.yticks(range(6),tickmarks)
    plt.xlabel('log10(size catalog 1)')
    plt.ylabel('log10(size catalog 2)')
    cb = plt.colorbar(im)
    plt.contour(np.log10(times.T),(0.,),colors='k',linestyles='dashed',lw=2,extend='neither')
    cb.set_label('log10(time) [seconds]')
