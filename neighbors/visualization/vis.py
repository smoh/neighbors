"""Stuff related to visualization"""
import numpy as np
import palettable
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Ellipse
from matplotlib import pyplot as plt

from isochrones.mist import MIST_Isochrone

__all__ = [
    'phasecolor',
    'get_isodf',
    'plotcmd_mist_isochrones',
    'violinplot',
    'pairplot',
    'plot_cov_ellipse'
]

# phase color dictionary for MIST isochrones
# from palettable.colorbrewer.sequential.Blues_7
phasecolor = {-1:'#084594',
              0:'#2171B5',
              2:'#4292C6',
              3:'#6BAED6',
              4:'#9ECAE1',
              5:'#C6DBEF',
              6:'#EFF3FF'}

labels = {
    'G' : r'$G$',
    'J' : r'$J$',
    'H' : r'$H$',
    'K' : r'$K$',
    'PanSTARRS_g' : r'$g_\mathrm{PS}$',
    'PanSTARRS_r' : r'$r_\mathrm{PS}$',
    'PanSTARRS_i' : r'$i_\mathrm{PS}$',
    'PanSTARRS_z' : r'$z_\mathrm{PS}$',
    'PanSTARRS_y' : r'$y_\mathrm{PS}$',
}

mist = MIST_Isochrone(
    bands=['G', 'J', 'H', 'K',
           'PanSTARRS_g', 'PanSTARRS_r', 'PanSTARRS_i', 'PanSTARRS_z', 'PanSTARRS_y'])

# MIST ages are in 0.05 dex increment, but the values are not exact.
# So if I query `mist.df.loc[0, 8.7]` for feh=0, age=8.7 it will incur
# index not found error because the actual value is 8.7000...1, which is kind of stupid.
def get_isodf(feh, logage):
    """Get sub-dataframe for given feh and logage"""
    iage = np.abs(mist.ages-logage).argmin()
    return mist.df.loc[feh,mist.ages[iage]]


def plotcmd_mist_isochrones(band1, band2, feh, logage, ax=None,
                            parallax=None, **kwargs):
    """
    Plot MIST isochrones on band1-band2, band1 color-magnitude diagram

    feh : float
        [Fe/H] in dex
    logage : float
        log (Age in Gyr)
    ax : matplotlib axes (optional)
        if None, current axes is used
    parallax : float
        in mas
    kwargs : dict
        passed to matplotlib.pyplot.plot; should not contain color
    """
    if ax == None:
        ax = plt.gca()
    iso = get_isodf(feh, logage)
    for phase, points in iso.groupby('phase'):
        mag = points[band1]
        if parallax:
            mag = mag - 5*(np.log10(parallax*1e-3)+1)
        ax.plot(points[band1]-points[band2], mag,
                c=phasecolor[phase], **kwargs)
    ax.set_xlabel('{:s} - {:s}'.format(labels[band1], labels[band2]))
    ax.set_ylabel(labels[band1])
    return ax


def violinplot(*args, **kwargs):
    """
    Same as matplotlib violin but with only medians marked in black

    facecolor : matplotlib color
        color of violin patches
    """
    facecolor = kwargs.pop('facecolor', 'steelblue')
    default = dict(showmedians=True, showextrema=False)
    default.update(kwargs)
    violins = plt.violinplot(*args, **default)
    for pc in violins['bodies']:
        pc.set_facecolor(facecolor)
    violins['cmedians'].set_color('k')
    return violins


def pairplot(x1, y1, x2, y2, ax=None, colors=None, linewidths=1, **kwargs):
    """
    Plot pairs connected by lines

    x1, y1, x2, y2 : 1d arrays
        coordinates of two end points of pairs
    colors : one color or 1d array of colors, dafault: None
        Colors used for lines and markers.
        If None (default), random color from tableau20 palette is assigned.
        Note that hex colors are not acceptable by LineCollection
    linewidths : float, 1d array, default: 1
        Linewidths used for lines.
    kwargs : dict
        additional keyword arguments for `matplotlib.collection.LineCollection`
    """
    tab20 = palettable.tableau.Tableau_20.mpl_colors
    if not (len(x1) == len(y1) == len(x2) == len(y2)):
        raise ValueError("length of x1, y1, x2, y2 should be the same")

    if not colors:
        colors = [ tab20[i%20] for i in range(len(x1))]
    elif len(list(colors)) == 1:
        colors = [colors] * len(x1)
    if len(colors) != len(x1):
        raise ValueError("length of colors array should match length of coordinates")

    try:
        lw = float(linewidths)
        linewidths = [lw] * len(x1)
    except TypeError:
        if len(linewidths) != len(x1):
            raise ValueError("length of linewidths array should match length of coordinates")

    segs = list(map(lambda x1,y1,x2,y2:[(x1,y1),(x2,y2)], x1,y1,x2,y2))
    lines = LineCollection(segs, colors=colors, facecolors=colors, linewidths=linewidths, **kwargs)
    ax = plt.gca() if not ax else ax
    ax.add_collection(lines)
    return lines


def plot_cov_ellipse(cov, pos, nstd=2, ax=None, **kwargs):
    """
    Plots an `nstd` sigma error ellipse based on the specified covariance
    matrix (`cov`). Additional keyword arguments are passed on to the
    ellipse patch artist.
    Parameters
    ----------
        cov : The 2x2 covariance matrix to base the ellipse on
        pos : The location of the center of the ellipse. Expects a 2-element
            sequence of [x0, y0].
        nstd : The radius of the ellipse in numbers of standard deviations.
            Defaults to 2 standard deviations.
        ax : The axis that the ellipse will be plotted on. Defaults to the
            current axis.
        Additional keyword arguments are pass on to the ellipse patch.
    Returns
    -------
        A matplotlib ellipse artist
    """
    def eigsorted(cov):
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:,order]

    if ax is None:
        ax = plt.gca()

    vals, vecs = eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))

    # Width and height are "full" widths, not radius
    width, height = 2 * nstd * np.sqrt(vals)
    ellip = Ellipse(xy=pos, width=width, height=height, angle=theta, **kwargs)

    ax.add_artist(ellip)
    return ellip
