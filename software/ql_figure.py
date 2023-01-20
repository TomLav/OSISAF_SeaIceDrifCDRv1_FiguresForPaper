import os.path
from datetime import datetime, time
from math import sqrt, log, ceil
import argparse
from argparse import RawDescriptionHelpFormatter
import numpy as np
from numpy import ma
import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
import matplotlib.image as image
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cmocean
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pyresample import utils, AreaDefinition
from netCDF4 import Dataset

from grid_info import region_params, valid_regions

col_land =  'Tan'
col_sea = '#04629a'
col_iceconc50 = '#9b9b9b'
col_iceconc100 = 'white'

def add_regionalbox(ax, lat_n, lat_s, lon_w, lon_e, col='r', lw=1):

    # Add regional boxes
    # e.g. Barents Sea latlon box:  72-82N, 10-60E
    # lat_n,lat_s,lon_w,lon_e = 72., 82., 10., 60.
    # Fram Strait latlon box: lat_n,lat_s,lon_w,lon_e = 70., 82., -20., 15.
    # Barents Sea latlon box: lat_n,lat_s,lon_w,lon_e = 72., 82., 10., 60.
    # Svalbard latlon box:    lat_n,lat_s,lon_w,lon_e = 72., 85., 0., 40.
    # Laptev Sea latlon box:  lat_n,lat_s,lon_w,lon_e = 69., 80., 100., 145.
    # Chukchi Sea latlon box: lat_n,lat_s,lon_w,lon_e = 65.5, 80., -180., -156.5
    bb_lon_n = np.linspace(lon_w, lon_e)
    bb_lat_n = lat_n* np.ones_like(bb_lon_n)
    bb_lat_e = np.linspace(lat_n, lat_s)
    bb_lon_e = lon_e* np.ones_like(bb_lat_e)
    bb_lon_s = np.linspace(lon_e, lon_w)
    bb_lat_s = lat_s* np.ones_like(bb_lon_s)
    bb_lat_w = np.linspace(lat_s, lat_n)
    bb_lon_w = lon_w* np.ones_like(bb_lat_w)
    bb_lon = np.concatenate((bb_lon_n, bb_lon_e, bb_lon_s, bb_lon_w))
    bb_lat = np.concatenate((bb_lat_n, bb_lat_e, bb_lat_s, bb_lat_w))

    ax.plot(bb_lon, bb_lat, '-', color=col, lw=lw,
            transform=ccrs.PlateCarree())

    return


def parse_args():

    drifttypechoice = ['simple', 'average', 'anomaly']
    valid_flgfmt = ['final', 'proc']

    p = argparse.ArgumentParser("ql_figure",
                                formatter_class=RawDescriptionHelpFormatter)
    p.add_argument('-o', '--output', required=False, default='.',
                   help="Either the output directory (in which case the "
                   "output filename is automatically determined) or the full "
                   "output filepath. If no output is specified, the file will "
                   "be put in the current working directory")
    p.add_argument('-d', '--driftfile', required=False, default=None,
                   help="Filepath to the drift file to plot")
    p.add_argument('-bg', '--bgvar', required=False, default=None,
                   help="Variable name of the first background to plot")
    p.add_argument('-bgf', '--bgfile', required=False, default=None,
                   help="Filepath to the background file to plot")
    p.add_argument('-bgy', '--bgyrot', required=False, default=None,
                   help='Variable name of the y-component with which to '
                   'rotate background (x-component) to the plot coordinate '
                   'reference system')
    p.add_argument('-bgx', '--bgxrot', required=False, default=None,
                   help='Variable name of the x-component with which to '
                   'rotate background (y-component) to the plot coordinate '
                   'reference system')
    p.add_argument('-bg2', '--bgvar2', required=False, default=None,
                   help="Variable name of the second background to plot")
    p.add_argument('-bgf2', '--bgfile2', required=False, default=None,
                   help="Filepath to the second background file to plot")
    p.add_argument('-dt', '--drifttype', required=False, default='simple',
                   choices=drifttypechoice,
                   help="Type of ice drift to plot (default='simple', "
                   "choices={}".format(drifttypechoice))
    p.add_argument('-df', '--driftflags', action='store_true', default=False,
                   help="Colour drift arrows with ice drift quality flags?")
    p.add_argument('--scale', required=False,
                   help="Scale factor to apply to drift arrows (default is "
                   "set by region)")
    p.add_argument('--skip', required=False,
                   help="Skip factor to apply to drift arrows (number of "
                   "arrows which are plotted out of the possible drift "
                   "vectors (default is set by region)")
    p.add_argument('-r', '--region', required=False, default='arc',
                   choices=valid_regions,
                   help="Region to plot, default is 'arc', choices are {}"
                   "".format(valid_regions))
    p.add_argument('--plt_transparent', action='store_true', default=False,
                   help="Plot on a transparent background")
    p.add_argument('--high_dpi', action='store_true', default=False,
                   help="Use high DPI when plotting (300 rather than 180)")
    p.add_argument('--title', action='store_true', default=False,
                   help="Plot a title")
    p.add_argument('--custom', required=False, default=None,
                   help='Add a custom title')
    p.add_argument('--mylab', required=False, default=None,
                   help="Plot a label in the plot?")
    p.add_argument('--colbar', action='store_true', default=False,
                   help="Plot a colour bar?")
    p.add_argument('--nocoast', action='store_true', default=False,
                   help="Do not plot coastlines")
    p.add_argument('--regbox', action='store_true', default=False,
                   help="Plot a regional box?")
    p.add_argument('--flatbox', required=False, default=None,
                   help="Plot a box in current projection, in terms of "
                   "comma-separated list of lat-lon pairs")
    p.add_argument('--transectline', action='store_true', default=False,
                   help="Plot a transect line?")
    p.add_argument('--latlongrid', action='store_true', default=False,
                   help="Plot a lat/lon grid?")
    p.add_argument('--concband', action='store_true', default=True,
                   help="Plot ice concentration as open water/open ice/close "
                   "ice")
    p.add_argument('--inv_y', action='store_true', default=False,
                   help="Plot -1 times the y-variable (depending on format of "
                   "files)")
    p.add_argument('-ff', '--flgfmt', required=False, default='final',
                   help="Format of status flags for plotting, default is "
                   "'final', choices are {}".format(valid_flgfmt))

    args = p.parse_args()

    return args


def nc_read(ncfile, var, skip=None):

    ncdata = {}

    if isinstance(var, list):
        v = var[0]
    else:
        v = var

    # Reading from the NetCDF file
    with Dataset(ncfile, 'r') as dataset:
        ncdata['grid_mapping'] = dataset.variables[v].__dict__['grid_mapping']
        try:
            ncdata['proj4_string'] = dataset.variables[
                ncdata['grid_mapping']].__dict__['proj4_string']
        except:
            ncdata['proj4_string'] = dataset.variables[
                ncdata['grid_mapping']].__dict__['proj4']
        try:
            ncdata['proj_dict'] = utils._proj4.proj4_str_to_dict(
                ncdata['proj4_string'])
        except:
            ncdata['proj_dict'] = utils.proj4.proj4_str_to_dict(
                ncdata['proj4_string'])
        ncdata['xc'] = dataset['xc'][:]
        ncdata['yc'] = dataset['yc'][:]

        try:
            ncdata['fv'] = dataset.variables[v].__dict__['_FillValue']
        except:
            pass

        if isinstance(var, list):
            varlist = var
        else:
            varlist = [var]
        for item in varlist:
            vardata = dataset[item][:]
            if len(vardata.shape) == 3:
                if skip:
                    vardata = vardata[0, ::skip, ::skip]
                else:
                    vardata = vardata[0, :, :]
            else:
                if skip:
                    vardata = vardata[::skip, ::skip]
                else:
                    vardata = vardata[:, :]

            # This is for reducing flags to the lowest integers
            if var in ['statusflag', 'status_flag', 'flag']:
                vardata = np.asarray(vardata, float) + 100
                uniques = np.unique(vardata)
                uniques = uniques[np.logical_not(np.isnan(uniques))]
                for newval, origval in enumerate(uniques):
                    vardata[vardata == origval] = newval
                ncdata['sf_labs'] = [str(int(u - 100)) for u in uniques]

            # NOTE: Be very careful with the fill value here. Trying to
            # use ncdata['fv'] as the fill_value for a status array such
            # as 'flag' means that flag values of 0 (i.e. nominal) are
            # masked out
            if 'fv' in ncdata.keys() and item not in ['statusflag',
                                                      'status_flag', 'flag']:
                ncdata[item] = ma.array(vardata, fill_value=ncdata['fv'])
            else:
                ncdata[item] = ma.array(vardata)

            if var in ['statusflag', 'status_flag', 'flag']:
                #ncdata[item].mask = nanmask
                ncdata[item].mask = None
            else:
                ncdata[item].mask = ncdata[item].data == ncdata[item].fill_value

        if skip:
            ncdata['lon'] = dataset.variables['lon'][::skip, ::skip]
            ncdata['lat'] = dataset.variables['lat'][::skip, ::skip]
        else:
            ncdata['lon'] = dataset.variables['lon'][:]
            ncdata['lat'] = dataset.variables['lat'][:]

        # Try fetching time info
        try:
            try:
                d0 = datetime.strptime(dataset.start_date,'%Y-%m-%d %H:%M:%S')
                d1 = datetime.strptime(dataset.stop_date,'%Y-%m-%d %H:%M:%S')
            except:
                try:
                    d0 = datetime.strptime(dataset.start_date_and_time,
                                           '%Y-%m-%dT%H:%M:%SZ')
                    d1 = datetime.strptime(dataset.end_date_and_time,
                                           '%Y-%m-%dT%H:%M:%SZ')
                except:
                    d0 = datetime.strptime(dataset.time_coverage_start,
                                           '%Y-%m-%dT%H:%M:%SZ')
                    d1 = datetime.strptime(dataset.time_coverage_end,
                                           '%Y-%m-%dT%H:%M:%SZ')
            d0_00 = datetime.combine(d0.date(),time(0))
            d1_00 = datetime.combine(d1.date(),time(0))
            ncdata['sdate'] = d0_00
            ncdata['edate'] = d1_00
            ncdata['tspan_hours'] = (d1_00 - d0_00).total_seconds() / (60.*60.)
        except:
            pass

        ncdata['time'] = dataset.variables['time'][0]
        ncdata['time_bnds0'] = dataset.variables['time_bnds'][0][0]
        ncdata['time_bnds1'] = dataset.variables['time_bnds'][0][1]

    if ncdata['yc'][-1] < ncdata['yc'][0]:
        print("FLIPPING!!!")
        ncdata['yc'] = np.flip(ncdata['yc'])

    # Grid spacing
    sorted_xc = np.sort(ncdata['xc'])
    sorted_yc = np.sort(ncdata['yc'])
    smallest_xc = sorted_xc[0]
    second_smallest_xc = sorted_xc[1]
    smallest_yc = sorted_yc[0]
    second_smallest_yc = sorted_yc[1]
    ncdata['ax'] = second_smallest_xc - smallest_xc
    ncdata['ay'] = second_smallest_yc - smallest_yc

    # Area definitions and extents
    if abs(float(ncdata['xc'][0])) < 10000.:
        sf = 1000.
    else:
        sf = 1.
    ncdata['area_extent'] = (float(ncdata['xc'][0] * sf),
                             float(ncdata['yc'][-1] * sf),
                             float(ncdata['xc'][-1] * sf),
                             float(ncdata['yc'][0] * sf))
    ncdata['mpl_extent'] = (ncdata['area_extent'][0],
                            ncdata['area_extent'][2],
                            ncdata['area_extent'][3],
                            ncdata['area_extent'][1])
    # Shifting by half a pixel for divergences and convergences, which
    # are calculated between pixels
    ncdata['mpl_extent_divs'] = (ncdata['area_extent'][0]
                                 + (0.5 * sf * ncdata['ax']),
                                 ncdata['area_extent'][2]
                                 + (0.5 * sf * ncdata['ax']),
                                 ncdata['area_extent'][3]
                                 + (0.5 * sf * ncdata['ay']),
                                 ncdata['area_extent'][1]
                                 + (0.5 * sf * ncdata['ay']))
    ncdata['area_def'] = AreaDefinition('data', 'data', 'data',
                                        ncdata['proj_dict'],
                                        ncdata['xc'].shape[0],
                                        ncdata['yc'].shape[0],
                                        ncdata['area_extent'])
    ncdata['data_crs'] = ncdata['area_def'].to_cartopy_crs()

    if ncdata['grid_mapping'] in ['Polar_Stereographic_Grid',
                                  'projection_stere']:
        data_globe = ccrs.Globe(semimajor_axis=ncdata['proj_dict']['a'],
                                semiminor_axis=ncdata['proj_dict']['b'])
        if ncdata['lat'][0, 0] > 0:
            ncdata['data_ccrs'] = ccrs.NorthPolarStereo(central_longitude=-45.0,
                                                        globe=data_globe)
            ncdata['hemi'] = 'nh'
        else:
            ncdata['data_ccrs'] = ccrs.SouthPolarStereo(central_longitude=0.0,
                                                        globe=data_globe)
            ncdata['hemi'] = 'sh'
    elif ncdata['grid_mapping'] in ['LambertAzimuthalEqualArea',
                                    'Lambert_Azimuthal_Equal_Area',
                                    'Lambert_Azimuthal_Grid',
                                    'projection_laea']:
        if ncdata['lat'][0, 0] > 0:
            ncdata['data_ccrs'] = ccrs.LambertAzimuthalEqualArea(
                central_longitude=0, central_latitude=90,
                false_easting=0, false_northing=0)
            ncdata['hemi'] = 'nh'
        else:
            ncdata['data_ccrs'] = ccrs.LambertAzimuthalEqualArea(
                central_longitude=0, central_latitude=-90,
                false_easting=0, false_northing=0)
            ncdata['hemi'] = 'sh'
    else:
        raise ValueError("Unrecognised grid mapping {}".format(
            ncdata['grid_mapping']))


    return ncdata


def rotate_bg_components(xdata, ydata, londata, latdata, datacrs, plotcrs):
    """Rotate/regrid background components from x, y in the data coordinate
    system to x, y in the plot coordinate system.

    Assume that masking is the same in the x and y components."""

    pc = ccrs.PlateCarree()

    xplot = xdata.copy()
    yplot = ydata.copy()
    for i in range(xdata.shape[0]):
        for j in range(xdata.shape[1]):
            dxd = xdata[i][j]
            dyd = ydata[i][j]
            lond = londata[i][j]
            latd = latdata[i][j]
            x0p, y0p = plotcrs.transform_point(lond, latd, src_crs=pc)

            x0d, y0d = datacrs.transform_point(lond, latd, src_crs=pc)
            x1d = x0d + dxd
            y1d = y0d + dyd
            x1p, y1p = plotcrs.transform_point(x1d, y1d, src_crs=datacrs)
            xplot[i][j] = x1p - x0p
            yplot[i][j] = y1p - y0p

    return xplot, yplot


def flag_arrow_col(flag, flgfmt='final'):

    # Use colours from https://sashamaps.net/docs/maps/roman-roads-index/
    fblack = '#000000'
    fmaroon = '#800000'
    forange = '#f58231'
    fnavy = '#000075'
    fblue = '#4363d8'
    flavender = '#dcbeff'
    fgrey = '#a9a9a9'

    fbrown = '#9A6324'
    fteal = '#469990'
    fgreen = '#3cb44b'
    fcyan = '#42d4f4'
    fmagenta = '#f032e6'

    flag_cols = {}
    flag_cols['final'] = {30: fblack,  # nominal quality
                          20: fbrown, # single-sensor, with smaller pattern
                                      # block
                          21: forange, # single-sensor, with neighbours as
                                       # constraint
                          22: fmaroon,   # interpolated
                          23: fcyan, # Gap filling in wind drift
                          24: fteal, # Vector replaced by wind drift
                          25: fnavy, # Blended satellite and wind drift
}
    flag_cols['proc'] = {0: 'black',   # nominal quality
                         16: 'purple',  # interpolated
                         13: 'red',     # single-sensor, with neighbours as
                                        # constraint
                         17: 'green',   # single-sensor, with smaller
                                        #pattern block
}

    default = 'black'
    if flag in flag_cols[flgfmt].keys():
        return flag_cols[flgfmt][flag]
    else:
        print("WARNING: unsupported flag value {}. Default to {}."
              "".format(flag, default))
        return default


def colbar_discrete_greys():

   # Discrete colorbar
    cmap = cm.Greys
    cmaplist = [cmap(i) for i in range(cmap.N)]
    collen = len(cmaplist)
    col0 = cmaplist[0]
    col1 = cmaplist[int(collen / 6) - 1]
    col2 = cmaplist[int(collen / 3) - 1]
    whitenum = 0
    for i in range(collen):
        if i <= whitenum:
            cmaplist[i] = col0 # Leave colourmap at white
        elif (i > whitenum) and (i < (collen - 1)):
            cmaplist[i] = col1
        else:
            cmaplist[i] = col2
    # Create the new map
    discrete_greys = matplotlib.colors.LinearSegmentedColormap.from_list(
        'discrete_greys', cmaplist, cmap.N)
    return discrete_greys


def colbar_sic_discrete():

    # SIC discrete colorbar with bottom 5% as white
    cmap_conc = cm.get_cmap('gray', 24)
    colors_i = np.linspace(0.0, 1.0, 24)
    c24 = cmap_conc(colors_i)
    c20 = c24[4:]
    # Below 5% conc set to white
    c20[0] = colors.colorConverter.to_rgba_array('white')
    # This is only needed for 10% conc
    cmap_conc20 = colors.LinearSegmentedColormap.from_list('listconc20',
                                                           c20, 20)
    return cmap_conc20


def colbar_simask_discrete():

    # Sea ice mask discrete colorbar with close ice as white
    cmap_conc = cm.get_cmap('PuBu_r', 3)
    colors_i = np.linspace(0.0, 1.0, 3)
    c3 = cmap_conc(colors_i)
    # Set high value set to white
    c3[-1] = colors.colorConverter.to_rgba_array('white')
    cmap_mask3 = colors.LinearSegmentedColormap.from_list('listmask3',
                                                           c3, 3)
    return cmap_mask3


def colbar_simask_discrete2():

    # Sea ice mask discrete colorbar with close ice as white, matching
    # new ice conc quicklooks

    cmap_sim = colors.LinearSegmentedColormap.from_list('listsim',
                                [col_sea, col_iceconc50, col_iceconc100], 3)

    return cmap_sim


def colbar_sit_discrete():

    # SIT discrete colorbar
    cmap_sit = cm.get_cmap(cmocean.cm.ice, 24)
    colors_i = np.linspace(0.0, 1.0, 24)
    csit28 = cmap_sit(colors_i)
    csit20 = csit28[8:]
    # Below 5% sit set to white
    cmap_sit20 = colors.LinearSegmentedColormap.from_list('listsit20',
                                                           csit20, 20)
    return cmap_sit20


def colbar_sity_discrete():

    # Sea-ice type discrete colorbar
    cmap_sity= cm.get_cmap(cmocean.cm.ice_r, 4)
    colors_i = np.linspace(0.0, 1.0, 4)
    csity = cmap_sity(colors_i)
    cmap_sity = colors.LinearSegmentedColormap.from_list('listsity4',
                                                           csity, 4)
    return cmap_sity


def colbar_ffsf_discrete(numcol):

    # Final format status flags discrete colorbar
    cmap_ffsf= cm.get_cmap(cm.Spectral, numcol)
    colors_i = np.linspace(0.0, 1.0, numcol)
    cffsf = cmap_ffsf(colors_i)
    cmap_ffsf = colors.LinearSegmentedColormap.from_list('listffsf',
                                                           cffsf, numcol)
    return cmap_ffsf


def colbar_type_status_flags():

    # Sea-ice type flags discrete colorbar
    cmap_sitysf = cm.get_cmap(cmocean.cm.deep, 4)
    colors_i = np.linspace(0.0, 1.0, 4)
    csitysf = cmap_sitysf(colors_i)
    cmap_sitysf = colors.LinearSegmentedColormap.from_list('listsitysf4',
                                                           csitysf, 4)
    return cmap_sitysf


def bg_plot_setup(bgvar=None, bgnc=None):
    '''Set up background quantities for plotting'''

    bg = {}
    bg['alpha'] = 1
    bg['norm'] = None
    bg['bounds'] = None
    bg['colbar_type'] = 'cont'

    # Average drift in background
    if bgvar in ['avdX', 'avdY']:
        bg['lbl']  = 'Average {} component over 24 hours [km]'.format(bgvar)
        bg['lim'] = 30
        bg['max'] = bg['lim']
        bg['min'] = -bg['lim']
        bg['cmap_lvl']  = [bg['min'], 0, bg['max']]
        bg['cmap'] = cm.RdBu
        bg['colbar_label'] = "metres"

    elif bgvar == 'displacement':
        bg['lbl']  = 'Displacement over 24 hours [km]'
        #bg['fld']  = pow((pow(u,2) + pow(v,2)), 0.5) / 1000
        bg['max'] = 30
        bg['min']  = 0
        bg['cmap'] = cmocean.cm.thermal
        bg['cmap_lvl'] = np.arange(bg['min'], bg['max'], 5)
        bg['cmap_fmt'] = '%2.0f'
        bg['colbar_label'] = bg['lbl']

    # Final format status flag array
    elif bgvar in ['statusflag', 'status_flag', 'flag']:
        bg['lbl']  = 'Status flags'
        bg['max'] = np.nanmax(bgnc[bgvar]) + 1
        bg['min'] = 0
        bg['cmap'] = colbar_ffsf_discrete(int(bg['max']))
        bg['cmap_labs'] = bgnc['sf_labs']
        bg['cmap_lvl'] = np.arange(bg['min'], bg['max'])
        bg['cmap_lvl'] = [x + 0.5 for x in bg['cmap_lvl']]
        bg['colbar_label'] = ""
        bg['colbar_type'] = 'discrete'

    elif bgvar in ['uncertainty']:
        bg['lbl']  = 'Uncertainty'
        bg['max'] = 1
        bg['min'] = 0
        bg['cmap'] = cmocean.cm.dense
        bg['cmap_lvl'] = np.arange(bg['min'],bg['max'], 0.1)
        bg['colbar_label'] = "Uncertainty"

    elif bgvar in ['uncert_dX_and_dY']:
        bg['lbl']  = 'Uncertainty'
        bg['max'] = ceil(np.nanmax(bgnc[bgvar]))
        bg['min'] = 0
        bg['cmap'] = cmocean.cm.dense
        bg['cmap_lvl'] = np.arange(bg['min'],bg['max'], 1)
        bg['colbar_label'] = "Uncertainty"

    elif bgvar in ['icetype', 'ice_type']:
        bg['lbl'] = 'Sea-ice type'
        bg['min'] = 0.5
        bg['max'] = 4.5
        bg['colbar_type'] = 'discrete'
        bg['cmap_lvl'] = [1, 2, 3, 4]
        bg['cmap_labs'] = ['Open Water', 'FYI', 'Ambiguous', 'MYI']
        bg['cmap'] = colbar_sity_discrete()
        bg['colbar_label'] = "Sea-ice type"

    elif bgvar in ['ice_edge']:
        bg['lbl'] = 'Sea-ice edge'
        bg['min'] = 0.5
        bg['max'] = 3.5
        bg['cmap'] = colbar_simask_discrete()
        bg['cmap_lvl'] = [1, 2, 3]
        bg['cmap_labs'] = ['Open Water', 'Open Ice', 'Close Ice']
        bg['colbar_label'] = "Sea-ice edge"
        bg['colbar_type'] = 'discrete'

    elif bgvar == 'sit':
        bg['lbl'] = 'Sea-ice thickness [m]'
        bg['min'] = 0
        bg['max'] = 2.5
        bg['cmap'] = colbar_sit_discrete()
        bg['cmap_lvl'] = np.arange(bg['min'], bg['max'], 0.5)
        bg['cmap_fmt'] = '%2.0f'
        bg['colbar_label'] = "Sea-ice thickness (m)"

    elif bgvar == 'sit_anom':
        bg['lbl'] = 'Sea-ice thickness anomaly [m]'
        bg['min'] = -1.2
        bg['max'] = 1.2
        bg['cmap'] = cmocean.cm.balance_r
        bg['cmap_lvl'] = [-1, -0.5, 0, 0.5, 1]
        bg['cmap_fmt'] = '%2.0f'
        bg['norm'] = matplotlib.colors.TwoSlopeNorm(vmin=bg['min'],
                                                         vcenter=0,
                                                         vmax=bg['max'])
        bg['colbar_label'] = "Sea-ice thickness anomaly (m)"

    elif bgvar == 'ice_conc':
        bg['lbl'] = 'Sea-ice concentration [%]'
        bg['min'] = 0
        bg['max'] = 100
        bg['cmap'] = colbar_sic_discrete()
        bg['cmap_lvl']  = np.arange(bg['min'],bg['max'],10)
        bg['cmap_fmt']  = '%2i'
        bg['colbar_label'] = "Sea-ice concentration (%)"

    elif bgvar == 'ice_mask':
        bg['lbl'] = 'Sea-ice concentration'
        bg['min'] = 0.5
        bg['max'] = 3.5
        bg['cmap'] = colbar_simask_discrete2()
        bg['cmap_lvl'] = [1, 2, 3]
        bg['cmap_labs'] = ['Open Water', 'Open Ice', 'Close Ice']
        bg['colbar_label'] = "Sea-ice concentration"
        bg['colbar_type'] = 'discrete'

    elif bgvar == 'conc_anom':
        bg['lbl'] = 'Sea-ice concentration anomaly [%]'
        bg['min'] = -100
        bg['max'] = 100
        bg['cmap'] = cmocean.cm.balance_r
        bg['cmap_lvl'] = np.arange(bg['min'],bg['max'],20)
        bg['cmap_fmt'] = '%2i'
        bg['norm'] = matplotlib.colors.TwoSlopeNorm(vmin=bg['min'],
                                                    vcenter=0,
                                                    vmax=bg['max'])
        bg['colbar_label'] = "Sea-ice concentration anomaly (%)"

    elif bgvar == 'sst':
        bg['lbl'] = 'Sea surface temperature [K]'
        bg['min'] = 270
        bg['max'] = 290
        bg['cmap'] = cmocean.cm.speed
        bg['cmap_lvl'] = np.arange(bg['min'],bg['max'],5)
        bg['cmap_fmt'] = '%2i'
        bg['colbar_label'] = "Sea surface temperature (K)"

    elif bgvar == 'sst_anom':
        bg['lbl'] = 'Sea surface temperature anomaly [K]'
        bg['min'] = -6
        bg['max'] = 6
        bg['cmap'] = cmocean.cm.balance
        bg['cmap_lvl'] = [-4, -2, 0, 2, 4]
        bg['cmap_fmt'] = '%2i'
        bg['norm'] = matplotlib.colors.TwoSlopeNorm(vmin=bg['min'],
                                                    vcenter=0,
                                                    vmax=bg['max'])
        bg['colbar_label'] = "Sea surface temperature anomaly (K)"

    elif bgvar in ['uwind_avg', 'vwind_avg']:
        bg['lbl'] = 'Wind component {}'.format(bgvar)
        bg['min'] = -12
        bg['max'] = 12
        bg['cmap'] = cmocean.cm.balance
        bg['cmap_lvl'] = [-12, -6, 0, 6, 12]
        bg['cmap_fmt'] = '%2i'
        bg['norm'] = matplotlib.colors.TwoSlopeNorm(vmin=bg['min'],
                                                    vcenter=0,
                                                    vmax=bg['max'])
        bg['colbar_label'] = "Wind component {}".format(bgvar)

    elif bgvar in ['driftX', 'driftY', 'dX', 'dY']:
        bg['lbl'] = '{} component over 24 hours [km]'.format(bgvar)
        bg['min'] = -30
        bg['max'] = 30
        bg['cmap'] = cmocean.cm.balance
        bg['cmap_lvl'] = [-30, -15, 0, 15, 30]
        bg['cmap_fmt'] = '%2i'
        bg['norm'] = matplotlib.colors.TwoSlopeNorm(vmin=bg['min'],
                                                    vcenter=0,
                                                    vmax=bg['max'])
        bg['colbar_label'] = "Ice drift component {}".format(bgvar)

    elif bgvar in ['t0', 't1']:
        if bgvar == 't0':
            bg['lbl'] = 'Start time [hours UTC]'
        elif bgvar == 't1':
            bg['lbl'] = 'Stop time [hours UTC]'
        bg['min'] = -12
        bg['max'] = 12
        bg['cmap'] = cmocean.cm.balance
        bg['cmap_lvl'] = [-12, -9, -6, -3, 0, 3, 6, 9, 12]
        bg['cmap_fmt'] = '%2i'
        bg['colbar_label'] = bg['lbl']

    elif bgvar in ['tdiff']:
        bg['lbl'] = 'Time span [hours]'
        bg['min'] = 21
        bg['max'] = 27
        bg['cmap'] = cmocean.cm.balance
        bg['cmap_lvl'] = [21, 22, 23, 24, 25, 26, 27]
        bg['cmap_fmt'] = '%2i'
        bg['colbar_label'] = bg['lbl']

    return bg


def ql_figure(output, driftfile=None, bgvar=None, bgfile=None, bgyrot=None,
              bgxrot=None, bgvar2=None, bgfile2=None, drifttype=None,
              driftflags=False, scale=None, skip=None, region='nh',
              plt_transparent=False, high_dpi=True, title=False, custom=None,
              mylab=None, colbar=False, nocoast=False, regbox=False,
              flatbox=False,
              transectline=False, latlongrid=False, concband=True,
              inv_y=False, flgfmt='final'):

    # Finding the region parameters
    rp = region_params(region)

    # Scale and skip for the drift
    if scale is None:
        scale = rp['scale']
    if skip is None:
        skip = rp['skip']

    # Read the drift file
    if driftfile is not None:
        with Dataset(driftfile, 'r') as temp:
            if drifttype == 'simple':
                if 'dX' in temp.variables:
                    dxvar = 'dX'
                    dyvar = 'dY'
                    sflag = 'status_flag'
                    varlist = [dxvar, dyvar, sflag]
                elif 'driftX' in temp.variables:
                    dxvar = 'driftX'
                    dyvar = 'driftY'
                    sflag = 'flag'
                    varlist = [dxvar, dyvar, sflag]
                # Allow plotting of wind vectors also
                elif 'uwind_avg' in temp.variables:
                    dxvar = 'uwind_avg'
                    dyvar = 'vwind_avg'
                    varlist = [dxvar, dyvar]
            elif drifttype == 'average':
                dxvar = 'avdX'
                dyvar = 'avdY'
                varlist = [dxvar, dyvar]

        driftnc = nc_read(driftfile, varlist, skip=skip)

        if inv_y:
            driftnc[dyvar][:] = driftnc[dyvar][:] * -1.

    # Reading in background files
    if bgvar == 'displacement':
        bgnc = nc_read(bgfile, ['dX', 'dY'])
    else:
        if (bgfile is not None) and (bgvar is not None):
            bgnc = nc_read(bgfile, bgvar)
    if (bgfile2 is not None) and (bgvar2 is not None):
        bgnc2 = nc_read(bgfile2, bgvar2)

    if bgvar == 'displacement':
        bgnc['displacement'] = np.sqrt((bgnc['dX'] * bgnc['dX'])
                                       + (bgnc['dY'] * bgnc['dY']))

    if bgvar == 't0' or bgvar == 't1':
        if bgvar == 't0':
            bgnc[bgvar] = bgnc[bgvar] - bgnc['time_bnds0']
        elif bgvar == 't1':
            bgnc[bgvar] = bgnc[bgvar] - bgnc['time_bnds1']
        bgnc[bgvar] = bgnc[bgvar] / 3600.
    if bgvar2 == 't0' or bgvar2 == 't1':
        if bgvar2 == 't0':
            bgnc2[bgvar2] = bgnc2[bgvar2] - bgnc2['time_bnds0']
        elif bgvar2 == 't1':
            bgnc2[bgvar2] = bgnc2[bgvar2] - bgnc2['time_bnds1']
        bgnc2[bgvar2] = bgnc2[bgvar2] / 3600.

    if bgvar == 't0' and bgvar2 == 't1':
        bgnc['tdiff'] = bgnc2['t1'] - bgnc['t0'] + 24
        bgvar = 'tdiff'
        bgnc2 = None
        bgvar2 = None

    # Defining the hemisphere and grid type
    if rp['lllat'] > 0.:
        hemi = 'nh'
    else:
        hemi = 'sh'
    if region.startswith('ease'):
        gridtype = 'ease'
    elif region.startswith('pol'):
        gridtype = 'polstere'
    else:
        raise ValueError("Unrecognised grid type for region {} - expected to "
                         "start with 'ease' or 'pol'".format(region))
    # Define grid based on region
    if gridtype == 'polstere':
        if hemi == 'nh':
            plot_proj4_params = {'proj': 'stere',
                                 'lat_0': 90.,
                                 'lat_ts' : 70.,
                                 'lon_0': -45.0,
                                 'a': 6378273,
                                 'b': 6356889.44891}
            plot_globe = ccrs.Globe(semimajor_axis=plot_proj4_params['a'],
                                    semiminor_axis=plot_proj4_params['b'])
            plot_crs = ccrs.NorthPolarStereo(
                central_longitude=plot_proj4_params['lon_0'], globe=plot_globe)
        else:
            plot_proj4_params = {'proj': 'stere',
                                 'lat_0': -90.,
                                 'lat_ts' : -70.,
                                 'lon_0': 0.,
                                 'a': 6378273,
                                 'b': 6356889.44891}
            plot_globe = ccrs.Globe(semimajor_axis=plot_proj4_params['a'],
                                    semiminor_axis=plot_proj4_params['b'])
            plot_crs = ccrs.SouthPolarStereo(
                central_longitude=plot_proj4_params['lon_0'], globe=plot_globe)
    elif gridtype == 'ease':
        if hemi == 'nh':
            plot_crs = ccrs.LambertAzimuthalEqualArea(central_longitude=0,
                                                      central_latitude=90,
                                                      false_easting=0,
                                                      false_northing=0)
        else:
            plot_crs = ccrs.LambertAzimuthalEqualArea(central_longitude=0,
                                                      central_latitude=-90,
                                                      false_easting=0,
                                                      false_northing=0)
    else:
        raise ValueError("Unrecognised region {}".format(region))

    pc = ccrs.PlateCarree()

    llx, lly = plot_crs.transform_point(rp['lllon'], rp['lllat'], src_crs=pc)
    urx, ury = plot_crs.transform_point(rp['urlon'], rp['urlat'], src_crs=pc)

    # Set the concentrations as open water/open ice/close ice
    if concband and bgvar == 'ice_conc':
        ic = bgnc['ice_conc']
        if ic.max() > 50:
            ic /= 100.
        bgnc['ice_mask'] = np.ma.empty((ic.shape),dtype=int)
        bgnc['ice_mask'].mask = ic.mask
        bgnc['ice_mask'][ic < 0.4] = 1
        bgnc['ice_mask'][np.logical_and(ic >= 0.4, ic < 0.7)] = 2
        bgnc['ice_mask'][(ic >= 0.7)] = 3
        # Stop the ocean from displaying underneath the land
        bgnc['ice_mask'][ic < 0] = 3
        bgvar = 'ice_mask'

    # If the background quantity is something that needs to rotate to the
    # plot coordinate system, rotate this
    if bgyrot is not None:
        bgncy = nc_read(bgfile, bgyrot)
        bgnc[bgvar], _ = rotate_bg_components(bgnc[bgvar], bgncy[bgyrot],
                                              bgnc['lon'], bgnc['lat'],
                                              bgnc['data_ccrs'], plot_crs)
    elif bgxrot is not None:
        bgncx = nc_read(bgfile, bgxrot)
        _, bgnc[bgvar] = rotate_bg_components(bgncx[bgxrot], bgnc[bgvar],
                                              bgnc['lon'], bgnc['lat'],
                                              bgnc['data_ccrs'], plot_crs)


    # Setting up background quantities for the plot
    if bgvar is not None:
        bgdict = bg_plot_setup(bgvar=bgvar, bgnc=bgnc)
    else:
        bgdict = None
    if bgvar2 is not None:
        bgdict2 = bg_plot_setup(bgvar=bgvar2, bgnc=bgnc2)
    else:
        bgdict2 = None

    # The projection keyword determines how the plot will look
    if colbar:
        fig = plt.figure(figsize=(8, 7))
    else:
        fig = plt.figure(figsize=(8, 6))
    ax = plt.axes(projection=plot_crs)
#    ax.set_extent(map_extent, crs=plot_crs)
    ax.set_extent([llx, urx, lly, ury], crs=plot_crs)

    if bgvar2 is not None:
        bg2 = ax.imshow(bgnc2[bgvar2], extent=bgnc2['mpl_extent'],
                        transform=bgnc2['data_ccrs'],
                        cmap=bgdict2['cmap'], vmax=bgdict2['max'],
                        vmin=bgdict2['min'], norm=bgdict2['norm'],
                        interpolation='none')

    # For icetype, we reverse the MYI and ambig numbers (want colour of
    # ambig ice type between FYI and MYI)
    if bgvar in ['icetype', 'ice_type']:
        ammask = bgnc[bgvar] == 4
        mymask = bgnc[bgvar] == 3
        bgnc[bgvar][ammask] = 3
        bgnc[bgvar][mymask] = 4


    if bgvar is not None:
        bg = ax.imshow(bgnc[bgvar], extent=bgnc['mpl_extent'],
                       transform=bgnc['data_ccrs'],
                       cmap=bgdict['cmap'], vmax=bgdict['max'],
                       vmin=bgdict['min'], norm=bgdict['norm'],
                       interpolation='none')

    if colbar:
        if bgdict['colbar_type'] == 'discrete':
            cb = plt.colorbar(bg, ticks=bgdict['cmap_lvl'],
                              orientation='horizontal', pad=0.05, shrink=0.7
            )
            cb.ax.set_xticklabels(bgdict['cmap_labs'])
            cb.set_label(bgdict['colbar_label'])
        else:
            plt.colorbar(bg, ticks=bgdict['cmap_lvl'],
                         orientation='horizontal', pad=0.05, shrink=0.65
            ).set_label(bgdict['colbar_label'])


    if not nocoast:
        ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land',
            '10m', edgecolor='black', facecolor=col_land, linewidth=0.2))
        if hemi == 'sh':
            # Fixing the weird line from the South Pole
            x0, y0 = plot_crs.transform_point(180., -90., src_crs=pc)
            x1, y1 = plot_crs.transform_point(180., -78.1, src_crs=pc)
            ax.plot([x0, x1], [y0, y1], color=col_land, linewidth=1.5)
            # Ice shelves
            ax.add_feature(cfeature.NaturalEarthFeature('physical',
                'antarctic_ice_shelves_polys', '10m', facecolor=col_land,
                edgecolor='black', linewidth=0.2))

    if driftfile is not None:
        pc = ccrs.PlateCarree()
        if drifttype in ['simple', 'average']:
            # Plotting the drift arrows
            ar_col = 'black'

            # Convert drift from km to m
            driftnc[dxvar] = driftnc[dxvar] * 1000.
            driftnc[dyvar] = driftnc[dyvar] * 1000.

            for i in range(driftnc[dxvar].size):
                try:
                    x0, y0 = plot_crs.transform_point(
                        driftnc['lon'][~driftnc[dxvar].mask][i],
                        driftnc['lat'][~driftnc[dxvar].mask][i],
                        src_crs=pc)

                    adx = driftnc[dxvar][~driftnc[dxvar].mask][i]
                    ady = driftnc[dyvar][~driftnc[dyvar].mask][i]
                    len_arrow = sqrt(adx**2 + ady**2)

                    # Calculate the endpoints and therefore dx, dy
                    # components of the drift arrows in the plot
                    # coordinate system
                    xorig, yorig = driftnc['data_ccrs'].transform_point(
                        driftnc['lon'][~driftnc[dxvar].mask][i],
                        driftnc['lat'][~driftnc[dxvar].mask][i],
                        src_crs=pc)
                    xarr = xorig + adx
                    # The dy component needs to be + for wind and - for ice
                    # drift
                    yarr = yorig - ady
                    x1, y1 = plot_crs.transform_point(xarr, yarr,
                        src_crs=driftnc['data_ccrs'])
                    pdx = scale * (x1 - x0)
                    pdy = scale * (y1 - y0)

                    # Set the colour of the drift arrows if this should be
                    # done with status flags
                    if driftflags:
                        myflag = driftnc[sflag][~driftnc[dxvar].mask][i]
                        ar_col = flag_arrow_col(myflag, flgfmt=flgfmt)

                    # If the arrow is too small, mark a symbol instead
                    if len_arrow * scale < 2000:
                        plt.plot(x0, y0, 's', color=ar_col, markersize=1)
                    else:
                        head_length = 0.3 * scale * len_arrow
                        plt.arrow(x0, y0, pdx, pdy, color=ar_col,
                                  shape='full', head_length=head_length,
                                  head_width=15000,
                                  fill=True, length_includes_head=True,
                                  width=4000)


                except:
                    pass # Outside the range of points

#            # Scale arrow
#            # 17280m is distance travelled in 2 days at 10cm/s
#            len_ar_sym = 17280. * scale
#            head_length_sym = 0.3 * len_ar_sym
#            plt.arrow(1350000, 2100000, len_ar_sym, 0, color=ar_col,
#                      shape='full', head_length=head_length_sym,
#                      head_width=30000,
#                      fill=True, length_includes_head=True,
#                      width=7000)

    if regbox:
        add_regionalbox(ax, 70, 80, 100, 145, col='dimgrey', lw=2)

    if flatbox is not None:
        pairs = [float(x) for x in flatbox.split(',')]
        mlats = []
        mlons = []
        while len(pairs) > 0:
            mlat = pairs.pop(0)
            mlon = pairs.pop(0)
            mlats.append(mlat)
            mlons.append(mlon)
        mxes = []
        myes = []
        for la, lo in zip(mlats, mlons):
            mx, my = plot_crs.transform_point(lo, la, src_crs=pc)
            mxes.append(mx)
            myes.append(my)
        ax.plot(mxes, myes, '-', color='orange', lw=2)

    if transectline:
        x0, y0 = plot_crs.transform_point(127.3255, 74.1267, src_crs=pc)
        x1, y1 = plot_crs.transform_point(112.5572, 84.4160, src_crs=pc)
        ax.plot([x0, x1], [y0, y1], color='dimgrey', linewidth=2)

    if latlongrid:
        llats = [0.1 * x for x in range(3600)]
        for llon in [-80. + 5. * x for x in range(33)]:
            llons = [llon] * 3600
            xyes = [plot_crs.transform_point(lalo[0], lalo[1], src_crs=pc)
                    for lalo in zip(llats, llons)]
            xes = []
            yes = []
            for item in xyes:
                xes.append(item[0])
                yes.append(item[1])
            ax.plot(xes, yes, color='black', linewidth=0.2, linestyle='--')
        llats = [x * 20. for x in range(18)]
        for llat in llats:
            x0, y0 = plot_crs.transform_point(llat, 80., src_crs=pc)
            x1, y1 = plot_crs.transform_point(llat, -80., src_crs=pc)
            ax.plot([x0, x1], [y0, y1], color='black', linewidth=0.2,
                    linestyle='--')

    # Label
    ftype = 'italic'
    labpos = {'nh': (1200000, 1900000),
              'sh': (1800000, -3400000),
              'arc': (1500000, 1800000)}
    if mylab is not None:
        plt.text(labpos[region][0], labpos[region][1], mylab, fontsize=14)

    if title:
        if custom is None:
            custom = "OSI-455 SIDrift 75.0km {} \n {} - {}".format(
                rp['long_name'],
                datetime.strftime(driftnc['sdate'], '%Y-%m-%d 12:00:00'),
            datetime.strftime(driftnc['edate'], '%Y-%m-%d 12:00:00'))
        plt.title(custom, fontsize='small', style=ftype)

    # Print area and copyright
    copyrightstr = '© '+ '(' + str(datetime.today().year).zfill(4) + ') ' + 'EUMETSAT'
    txt1x = 0.8
    if 'labelpos' in rp:
        txt1x = rp['labelpos']
    txt1y = -0.03
    copyfontsize = 6
    if 'copyfontsize' in rp:
        copyfontsize = rp['copyfontsize']
    ax.annotate(copyrightstr, xy=(txt1x,txt1y), xycoords='axes fraction',
                fontsize=copyfontsize, style=ftype)

    # Logo
    dirname = os.path.abspath(os.path.dirname(__file__))
    logo = os.path.join(dirname, 'OSISAF_Name_Colour.png')
    # Position as left, bottom, width, height - could move to grid_info
    logo_left = 0.73
    logo_bottom = 0.75
    logo_width = 0.15
    logo_height = 0.15
    if 'logo_left' in rp:
        logo_left = rp['logo_left']
    if 'logo_bottom' in rp:
        logo_bottom = rp['logo_bottom']
    if 'logo_width' in rp:
        logo_width = rp['logo_width']
    if 'logo_height' in rp:
        logo_height = rp['logo_height']
    logo_xy = [logo_left, logo_bottom, logo_width, logo_height]
    img = image.imread(logo)
    ax2 = fig.add_axes(logo_xy)
    ax2.imshow(img, zorder=1)
    ax2.axis('off')


    if colbar and bgvar2 is not None:
        ins_ax = inset_axes(ax, width='20%', height='5%', loc=1, borderpad=2)
        cb2 = plt.colorbar(bg2, cax=ins_ax, ticks=bgdict2['cmap_lvl'],
                           orientation='horizontal')#.set_label(bgdict2['colbar_label'])
        cb2.set_label(bgdict2['colbar_label'])

    if os.path.isdir(output):
        name = None
        if isinstance(bgvar, str):
            name = bgvar
            if type(drifttype) == str:
                name = name + '_' + drifttype
        elif isinstance(drifttype, str):
            name = drifttype
        else:
            print("WARNING! No output name")
        if type(region) == str:
            name = name + '_' + region
        foutput = os.path.basename(driftfile).replace('.nc',
                                                  '_' + name + '.png')
        output = os.path.join(output,foutput)
        print("Output file at {}".format(output))

    dpi = 180
    if high_dpi:
        dpi = 300

    #plt.tight_layout()

    plt.savefig(output, bbox_inches='tight', transparent=plt_transparent,
                dpi=dpi)
#    plt.show()
    plt.close()
    print("Figure is in {}".format(output))


def main():

    args = parse_args()

    driftfile = args.driftfile
    output = args.output
    bgvar = args.bgvar
    bgfile = args.bgfile
    bgyrot = args.bgyrot
    bgxrot = args.bgxrot
    bgvar2 = args.bgvar2
    bgfile2 = args.bgfile2
    drifttype = args.drifttype
    driftflags = args.driftflags
    scale = None
    if args.scale:
        scale = int(args.scale)
    skip = None
    if args.skip:
        skip = int(args.skip)
    region = args.region
    plt_transparent = args.plt_transparent
    high_dpi = args.high_dpi
    title = args.title
    custom = args.custom
    mylab = args.mylab
    colbar = args.colbar
    nocoast = args.nocoast
    regbox = args.regbox
    flatbox = args.flatbox
    transectline = args.transectline
    latlongrid = args.latlongrid
    concband = args.concband
    inv_y = args.inv_y
    flgfmt = args.flgfmt

    ql_figure(output, driftfile=driftfile, bgvar=bgvar, bgfile=bgfile,
              bgyrot=bgyrot, bgxrot=bgxrot, bgvar2=bgvar2, bgfile2=bgfile2,
              drifttype=drifttype, driftflags=driftflags, scale=scale,
              skip=skip, region=region, plt_transparent=plt_transparent,
              high_dpi=high_dpi, title=title, custom=custom, mylab=mylab,
              colbar=colbar, nocoast=nocoast, regbox=regbox, flatbox=flatbox,
              transectline=transectline, latlongrid=latlongrid,
              concband=concband, inv_y=inv_y, flgfmt=flgfmt)


if __name__ == '__main__':

    main()
