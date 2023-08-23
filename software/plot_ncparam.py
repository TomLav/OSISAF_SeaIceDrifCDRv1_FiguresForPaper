'''Plot the parameter files from the icedrift from winds.

Updated version to use combined parameter files with all months'''
import os
import sys
import math
import argparse
from argparse import RawDescriptionHelpFormatter
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import matplotlib
import matplotlib.style
matplotlib.use('Agg')
matplotlib.style.use('classic')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pyresample import utils, AreaDefinition
import cmocean

from grid_info import region_params

tickfont = 30

month_dict = {1: 'jan',
              2: 'feb',
              3: 'mar',
              4: 'apr',
              5: 'may',
              6: 'jun',
              7: 'jul',
              8: 'aug',
              9: 'sep',
              10: 'oct',
              11: 'nov',
              12: 'dec'}

SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 15
BIG_SIZE = 18
HUGE_SIZE = 30
plt.rc('font', size=BIG_SIZE)          # controls default text sizes
#plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
#plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
#plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
#plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
#plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
#plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

def parse_args():

    p = argparse.ArgumentParser("plot_ncparam",
                                formatter_class=RawDescriptionHelpFormatter)
    p.add_argument('-i', '--infile', required=True,
                   help="Input file with field to plot")
    p.add_argument('-m', '--month', required=True,
                   help="Month to plot as an integer (1 for Jan, 2 for Feb...)")
    p.add_argument('-p', '--param', required=True,
                   help="Parameter to plot")
    p.add_argument('-p2', '--param2', required=False, default=None,
                   help="Second parameter to average or mean square with first "
                   "parameter")
    p.add_argument('-s', '--meansq', action='store_true', default=False,
                   help="Use mean square rather than averaging two params. "
                   "This is used for speed.")
    p.add_argument('-g', '--geocurr', action='store_true', default=False,
                   help="Plot the geostrophic current arrows")
    p.add_argument('--gskip', required=False,
                   help="Skip factor to plot sparser geostrophic current "
                   "arrows")
    p.add_argument('--gscale', required=False, default=1,
                   help="Scale factor for geostropic current arrows")
    p.add_argument('--nobg', required=False, default=False, action='store_true',
                   help="Skip plotting the background variable and just plot "
                   "geostrophic currents")
    p.add_argument('-o', '--outdir', required=True,
                   help="Output plot directory")
    p.add_argument('-c', '--cline', action='store_true', default=False,
                   help="Plot a contour line around the un-gapfilled data")
    p.add_argument('-l', '--label', default=None,
                   help="Label for the colorbar (default is to use -p <param>)")
    args = p.parse_args()

    return args


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def cmap_norm_to_lims(cmap, mapmin, mapmax, vmin, vmax):
    '''Return a colourmap that will be normalised from mapmin to mapmax (so
    as to get zero in the right place) but which is cropped to vmin and vmax'''

    if vmin < mapmin or vmax > mapmax:
        raise ValueError("vmin and vmax must be within the bounds of mapmin "
                         "and mapmax")
    minval = (vmin - mapmin) / (mapmax - mapmin)
    maxval = (vmax - mapmin) / (mapmax - mapmin)
    new_cmap = truncate_colormap(cmap, minval=minval, maxval=maxval)

    return new_cmap


def colbar_flags_discrete():

    # The flags we were interested in are
    # 0 - Good data
    # 1 - Gapfilled
    # 2 - Land
    # 3 - No ice
    # 4 - No drift vectors

    colorlist = ["firebrick", "goldenrod", "grey", "midnightblue", "teal"]
    cmap_flags = colors.ListedColormap(colorlist)

    return cmap_flags


def plot_ncparam(infile, month, param, param2=None, meansq=False,
                 geocurr=False, gskip=None, gscale=1, nobg=False, outdir='.',
                 cline=False, label=None):

    if param.endswith('_gapfill'):
        gapfill = True
    else:
        param = param + '_gapfill'
        gapfill = False

    if param2 is not None:
        if not param2.endswith('_gapfill'):
            param2 = param2 + '_gapfill'

    print("Param: {}, Param2: {}, Gapfill: {}".format(param, param2, gapfill))

    with Dataset(infile, 'r') as dataset:
        try:
            grid_mapping = dataset.variables[param].__dict__['grid_mapping']
        except:
            grid_mapping = dataset.variables['A_real'].__dict__['grid_mapping']

        proj4_string = dataset.variables[grid_mapping].__dict__['proj4_string']
        data_proj_dict = utils.proj4.proj4_str_to_dict(proj4_string)
#        data_proj_dict = utils._proj4.proj4_str_to_dict(proj4_string)
        xc = dataset['xc']
        yc = dataset['yc']
        flip = False
        # Flip if necessary
        if yc[-1] < yc[0]:
            flip = True
            yc = np.flip(yc)
        data_area_extent = (float(xc[0]*1000.),
                            float(yc[-1]*1000.),
                            float(xc[-1]*1000.),
                            float(yc[0]*1000.))
        mpl_data_extent = (data_area_extent[0],
                           data_area_extent[2],
                           data_area_extent[3],
                           data_area_extent[1])
        data_area_def = AreaDefinition('data', 'data', 'data',
                                       data_proj_dict,
                                       xc.shape[0],
                                       yc.shape[0],
                                       data_area_extent)

        if grid_mapping in ['Polar_Stereographic_Grid',
                            'projection_stere']:
            data_globe = ccrs.Globe(semimajor_axis=ncdata['proj_dict']['a'],
                                    semiminor_axis=ncdata['proj_dict']['b'])
            if dataset['lat'][0, 0] > 0:
                data_ccrs = ccrs.NorthPolarStereo(central_longitude=-45.0,
                                                  globe=data_globe)
                region = 'polstere-nh'
                hemi = 'nh'
            else:
                data_ccrs = ccrs.SouthPolarStereo(central_longitude=0.0,
                                                  globe=data_globe)
                region = 'polstere-sh'
                hemi = 'sh'
        elif grid_mapping in ['LambertAzimuthalEqualArea',
                              'Lambert_Azimuthal_Equal_Area',
                              'Lambert_Azimuthal_Grid',
                              'projection_laea']:
            if dataset['lat'][0, 0] > 0:
                data_ccrs = ccrs.LambertAzimuthalEqualArea(central_longitude=0,
                                                           central_latitude=90,
                                                           false_easting=0,
                                                           false_northing=0)
                region = 'ease-nh-wide'
                hemi = 'nh'
            else:
                data_ccrs = ccrs.LambertAzimuthalEqualArea(central_longitude=0,
                                                           central_latitude=-90,
                                                           false_easting=0,
                                                           false_northing=0)
                region = 'ease-sh-wide'
                hemi = 'sh'
        else:
            raise ValueError("Unrecognised grid mapping {}".format(
                grid_mapping))

        data = dataset.variables[param][month - 1, :, :]
        dmask = dataset.variables[param][month - 1, :, :].mask
        if param2 is not None:
            data2 = dataset.variables[param2][month - 1, :, :]
            dmask2 = dataset.variables[param2][month -1, :, :].mask

        # Reading in the geostrophic current info
        gcxvar = 'C_real_gapfill'
        gcyvar = 'C_imag_gapfill'
        if gskip is not None:
            lats = dataset.variables['lat'][::gskip, ::gskip]
            lons = dataset.variables['lon'][::gskip, ::gskip]
            geocx = dataset.variables[gcxvar][month - 1, ::gskip, ::gskip]
            geocxmask = dataset.variables[gcxvar][month - 1, ::gskip, ::gskip].mask
            geocy = dataset.variables[gcyvar][month - 1, ::gskip, ::gskip]
            geocymask = dataset.variables[gcyvar][month - 1, ::gskip, ::gskip].mask
        else:
            lats = dataset.variables['lat'][:]
            lons = dataset.variables['lon'][:]
            geocx = dataset.variables[gcxvar][month - 1, :, :]
            geocxmask = dataset.variables[gcxvar][month - 1, :, :].mask
            geocy = dataset.variables[gcyvar][month - 1, :, :]
            geocymask = dataset.variables[gcyvar][month - 1, :, :].mask
        # If the parameter file should not be gapfilled, then we want
        # to use the flags to reduce the number of geostrophic current
        # vectors
        if not gapfill:
            if gskip is not None:
                flags = dataset.variables['flags'][month - 1, ::gskip, ::gskip]
            else:
                flags = dataset.variables['flags'][month - 1, :, :]
            geocxmask[flags != 0] = 1
            geocymask[flags != 0] = 1

        # Some of the data is FLOAT nans, not numpy nans...
        if dmask.sum() == 0:
            dmask_new = np.zeros_like(data, dtype=np.int8)
            for iy, ix in np.ndindex(data.shape):
                if math.isnan(data[iy, ix]):
                    dmask_new[iy, ix] = 1
            dmask = dmask_new
        if param2 is not None:
            if dmask2.sum() == 0:
                dmask2_new = np.zeros_like(data2, dtype=np.int8)
                for iy, ix in np.ndindex(data2.shape):
                    if math.isnan(data2[iy, ix]):
                        dmask2_new[iy, ix] = 1
                dmask2 = dmask2_new
            dmask = np.logical_and(dmask, dmask2)
        if geocxmask.sum() == 0:
            geocxmask_new = np.zeros_like(data, dtype=np.int8)
            for iy, ix in np.ndindex(data.shape):
                if math.isnan(data[iy, ix]):
                    geocxmask_new[iy, ix] = 1
            geocxmask = geocxmask_new
        if geocymask.sum() == 0:
            geocymask_new = np.zeros_like(data, dtype=np.int8)
            for iy, ix in np.ndindex(data.shape):
                if math.isnan(data[iy, ix]):
                    geocymask_new[iy, ix] = 1
            geocymask = geocymask_new

        # If the second parameter is set, average these
        if param2 is not None:
            if meansq:
                dsq = data * data
                d2sq = data2 * data2
                datasq = dsq + d2sq
                data = np.sqrt(datasq)
                data = np.ma.array(data)
            else:
                data = np.mean(np.array([data, data2]), axis=0)
                data = np.ma.array(data)

        if param != 'flags':
            data[dmask] = np.nan
            data.mask = dmask
            count = data.count()
            thresh = math.ceil(count / 100.)
            data1d = data.flatten()
            data1dclean = [x for x in data1d if x != np.nan]
            datasort = sorted(data1dclean)
            print("-- min(data) = ", np.nanmin(data))
            print("-- max(data) = ", np.nanmax(data))
            print("-- lower threshold data (99%) = ", datasort[thresh - 1])
            print("-- upper threshold data (99%) = ", datasort[-thresh])
            print("-- mean(data) = ", np.nanmean(data))

        if param == 'flags':
            # Make a plotting array with
            # 0 - OK (0)
            # 1 - Gapfill (19)
            # 2 - Land (3)
            # 3 - No ice (4)
            # 4 - No drift (20)
            data[data == 19] = 1
            data[data == 3] = 2
            data[data == 4] = 3
            data[data == 20] = 4

        # Information for the contour line
        ff = dataset.variables['flags'][month - 1, :, :]
        ffmask = dataset.variables['flags'][month - 1, :, :].mask
        ff[ff == 0] = 100
        ff[ff != 100] = 0
        if ffmask:
            ff.mask = False
            ff[ffmask] = 0
        if flip:
            ff = np.flip(ff, axis=0)

        # if we do not want to plot the gapfilled field, mask the gapfilled cells
        if not gapfill:
            fmask = dataset.variables['flags'][month - 1, :, :]
            data[fmask == 19] = np.ma.masked


    # Plotting setup
    colbar_type = 'cont'
    inv_yax = False
    colbar_label = label
    if param in  ['A_real', 'A_imag', 'A_real_gapfill', 'A_imag_gapfill']:
        datacmap = cmocean.cm.haline
        vmin = -0.02
        vmax = 0.05
        data_cmap_lvl  = [-0.02, -0.01, 0.0, 0.01, 0.02, 0.03, 0.04, 0.05]
        norm = colors.TwoSlopeNorm(vmin=-0.05, vcenter=0, vmax=vmax)
    if param in  ['C_real', 'C_imag', 'C_real_gapfill', 'C_imag_gapfill']:
        # Speed rather than x- or y- component
        if param2 is not None:
            datacmap = cmocean.cm.haline
            vmin = 0.0
            vmax = 0.05
            data_cmap_lvl  = [0, 0.025, 0.05, ]
            norm = None
        else:
            datacmap = cmocean.cm.balance
            vmin = -0.06
            vmax = 0.06
            data_cmap_lvl  = [-0.06, -0.03, 0.0, 0.03, 0.06]
            norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
        if label is None:
            colbar_label = r'(c)  $U_{wg}\ [m.s^{-1}]$'
    if param in  ['RMS_Res_real', 'RMS_Res_imag']:
        datacmap = cmocean.cm.haline
        vmin = 0.0
        vmax = 0.1
        data_cmap_lvl  = [0, 0.02, 0.04, 0.06, 0.08, 0.1]
        norm = None
    if param in  ['absA', 'absA_gapfill']:
        # express the tranfer coefficient as %
        data *= 100
        datacmap = cmocean.cm.haline
        vmin = 0
        vmax = 3.
        data_cmap_lvl  = [0.0, 1, 2, 3]
        norm = None
        if label is None:
            colbar_label = r"(a)  $|A|\ [\%]$"
    if param in  ['ThetaA', 'ThetaA_gapfill']:
        if hemi == 'nh':
            vmin = -40.
            vmax = 0.
            data_cmap_lvl  = np.arange(vmin, vmax+0.01,10)
            norm = None
            datacmap = cm.Blues_r
        else:
            vmin = 0.
            vmax = +40.
            data_cmap_lvl  = np.arange(vmin, vmax+0.01,10)
            norm = None
            inv_yax = True
            datacmap = cm.Reds
        if label is None:
            colbar_label = r"(b)  $\theta\ [ ^{\circ}]$"
    if param in  ['lon']:
        datacmap = cmocean.cm.balance
        vmin = -180
        vmax = 180
        data_cmap_lvl  = [-150, -100, -50, 0, 50, 100, 150]
        norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    if param in  ['flags']:
        datacmap = colbar_flags_discrete()
        vmin = 0
        vmax = 5
        data_cmap_lvl  = [0.5, 1.5, 2.5, 3.5, 4.5]
        cmap_labs = ['Param', 'Gapfilled', 'Land', 'No ice', 'No drift']
        norm = None
        colbar_type = 'discrete'

    # Fetching the region parameters from the file and converting
    # from lat/lon
    pc = ccrs.PlateCarree()
    rp = region_params(region)
    if rp is None:
        raise ValueError("No grid information for region {}".format(region))

    llx, lly = data_ccrs.transform_point(rp['lllon'], rp['lllat'],
                                         src_crs=pc)
    urx, ury = data_ccrs.transform_point(rp['urlon'], rp['urlat'],
                                         src_crs=pc)
    plot_data_extent = [llx, urx, lly, ury]

    # Setting up the plot
    colbar = True
    if colbar:
        plt.figure(figsize=(8.5, 6))
    else:
        plt.figure(figsize=(8, 6))
    ax = plt.axes(projection=data_ccrs)
    ax.set_extent(plot_data_extent, crs=data_ccrs)

    # Plotting the data
    if not nobg:
        datacbar = ax.imshow(data, cmap=datacmap, extent=mpl_data_extent,
                             transform=data_ccrs, vmin=vmin, vmax=vmax,
                             norm=norm, interpolation='none')
        #cl = ax.contour(data, colors='k', levels=data_cmap_lvl, extent=mpl_data_extent,
        #        transform=data_ccrs, origin='image')
        #ax.clabel(cl, inline=True, fontsize=10)

    # Plotting the geostrophic currents
    if geocurr:
        pc = ccrs.PlateCarree()
        exscale = 1000000
        scale = exscale * gscale
        for i in range(geocx.size):
            try:
                if i % 2: continue
                # This was plot_crs
                x0, y0 = data_ccrs.transform_point(lons[~geocxmask][i],
                                                   lats[~geocxmask][i],
                                                   src_crs=pc)

                adx = geocx[~geocxmask][i]
                ady = geocy[~geocymask][i]
                len_arrow = math.sqrt(adx**2 + ady**2)
                if len_arrow >= 0.05:
                    continue

                # Calculate the endpoints and therefore dx, dy components of
                # the geostrophic current arrows in the plot coordinate system
                xorig, yorig = data_ccrs.transform_point(lons[~geocxmask][i],
                                                         lats[~geocxmask][i],
                                                         src_crs=pc)
                xarr = xorig + adx
                yarr = yorig + ady
                # This was plot_crs
                x1, y1 = data_ccrs.transform_point(xarr, yarr,
                                                   src_crs=data_ccrs)
                pdx = scale * (x1 - x0)
                pdy = scale * (y1 - y0)

                # If the arrow is too small, mark a symbol instead
                if len_arrow * scale < 2000:
                    plt.plot(x0, y0, 's', color='black', markersize=1)
                else:
                    head_length = 0.3 * scale * len_arrow
                    plt.arrow(x0, y0, pdx, pdy, color='black',
                              shape='full', head_length=head_length,
                              head_width=0.005 * scale,
                              fill=True, length_includes_head=True,
                              width=0.001 * scale)

            except:
                pass # Outside the range of points

    # Adding a contour line
    if cline:
        cs = ax.contour(ff, [95], colors='red', linewidths=2,
                        transform=data_ccrs, extent=mpl_data_extent)

    if not nobg:
        if colbar:
            if colbar_type == 'discrete':
                cb = plt.colorbar(datacbar, ticks=data_cmap_lvl,
                                  orientation='horizontal', pad=0.05, shrink=0.6)
                cb.ax.set_yticklabels(cmap_labs)
                cb.set_label(colbar_label)
                cb.ax.tick_params(labelsize=tickfont)
            else:
                #plt.yticks(fontsize=tickfont)
                plt.colorbar(datacbar, ticks=data_cmap_lvl,
                             orientation='horizontal', pad=0.05, shrink=0.6
                         ).set_label(colbar_label)
                #cb.ax.tick_params(labelsize=tickfont)
                #ticklabs = cb.ax.get_yticklabels()
                #cb.ax.set_yticklabels(ticklabs, fontsize=tickfont)

    if param != 'flags':
        ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land',
                                                    '50m', edgecolor='grey',
                                                    facecolor='grey'))
        if hemi == 'sh':
            ax.add_feature(cfeature.NaturalEarthFeature('physical',
                    'antarctic_ice_shelves_polys', '10m', facecolor='grey',
                                                        edgecolor='grey'))

    if not gapfill:
        param = param.replace('_gapfill','')
        if param2:
            param2 = param2.replace('_gapfill','')

    if param2 is not None:
        param = '{}_{}'.format(param, param2)
    outbname = "{}_{}_{}.png".format(os.path.basename(infile)[:-3],
                                     month_dict[month], param)
    outname = os.path.join(outdir, outbname)
    plt.savefig(outname, bbox_inches='tight', dpi=150)
    plt.close()
    print("Figure is in {}".format(outname))


def main():


    args = parse_args()
    infile = args.infile
    month = int(args.month)
    param = args.param
    param2 = args.param2
    meansq = args.meansq
    geocurr = args.geocurr
    if args.gskip:
        gskip = int(args.gskip)
    else:
        gskip = None
    gscale = int(args.gscale)
    nobg  = args.nobg
    outdir = args.outdir
    cline = args.cline
    label = args.label

    plot_ncparam(infile, month, param, param2=param2, meansq=meansq,
                 geocurr=geocurr, gskip=gskip, gscale=gscale, nobg=nobg,
                 outdir=outdir, cline=cline, label=label)


if __name__ == '__main__':

    main()
