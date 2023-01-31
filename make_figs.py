#!/usr/bin/env python

import os
import shutil
from subprocess import check_call
from datetime import datetime, timedelta

"""
   Main script to prepare the figures for the EGU TC manuscript
"""

# Configuration
tmp_dir = './tmp_figs'
if os.path.exists(tmp_dir):
    shutil.rmtree(tmp_dir)
os.mkdir(tmp_dir)

final_dir = './final_figs'
if os.path.exists(final_dir):
    shutil.rmtree(final_dir)
os.mkdir(final_dir)


# This is the main figure counter, to create fig01, fig02, etc...
fig_cnt = 0

# Maps of buoy coverage (6 panels)
# =======================================================
fig_cnt += 1
# 1. run the notebook/script

# 2. assemble the panels into a publication-ready file
nh_files = ['figs/mrg_wind_pamb_ssmiamsr_{a:}_{p:}_map.png'.format(p=pstr,a='NH') for pstr in ('1991-2000','2001-2010','2011-2020')]
cmd = 'convert {} +append {}/fig_buoycover_nh.png'.format(" ".join(nh_files), tmp_dir)
check_call(cmd, shell=True)
sh_files = ['figs/mrg_wind_pamb_ssmiamsr_{a:}_{p:}_map.png'.format(p=pstr,a='SH') for pstr in ('1991-2000','2001-2010','2011-2020')]
cmd = 'convert {} +append {}/fig_buoycover_sh.png'.format(" ".join(sh_files), tmp_dir)
check_call(cmd, shell=True)
cmd = 'convert {}/fig_buoycover_nh.png {}/fig_buoycover_sh.png -append {}/fig{:02d}.png'.format(tmp_dir, tmp_dir, final_dir, fig_cnt)
check_call(cmd, shell=True)


# Example maps free-drift parameters (1 for NH, 1 for SH)
# =======================================================
fig_cnt += 1
# 1. run the notebook/script
cmd = "python software/plot_ncparam.py -i https://thredds.met.no/thredds/dodsC/metusers/thomasl/SIDrift_CDR_v1pre/auxiliary_files/inv_params_osi455_nh_200301-202012_1day.nc -m 7 -o figs -p absA "
check_call(cmd, shell=True)
cmd = "python software/plot_ncparam.py -i https://thredds.met.no/thredds/dodsC/metusers/thomasl/SIDrift_CDR_v1pre/auxiliary_files/inv_params_osi455_nh_200301-202012_1day.nc -m 7 -o figs -p ThetaA "
check_call(cmd, shell=True)
cmd = "python software/plot_ncparam.py -i https://thredds.met.no/thredds/dodsC/metusers/thomasl/SIDrift_CDR_v1pre/auxiliary_files/inv_params_osi455_nh_200301-202012_1day.nc -m 7 -o figs -p C_real -p2 C_imag -s -g --gskip 2 --gscale 15 "
check_call(cmd, shell=True)

cmd = "python software/plot_ncparam.py -i https://thredds.met.no/thredds/dodsC/metusers/thomasl/SIDrift_CDR_v1pre/auxiliary_files/inv_params_osi455_sh_200208-202007_1day.nc -m 12 -o figs -p absA "
check_call(cmd, shell=True)
cmd = "python software/plot_ncparam.py -i https://thredds.met.no/thredds/dodsC/metusers/thomasl/SIDrift_CDR_v1pre/auxiliary_files/inv_params_osi455_sh_200208-202007_1day.nc -m 12 -o figs -p ThetaA "
check_call(cmd, shell=True)
cmd = "python software/plot_ncparam.py -i https://thredds.met.no/thredds/dodsC/metusers/thomasl/SIDrift_CDR_v1pre/auxiliary_files/inv_params_osi455_sh_200208-202007_1day.nc -m 12 -o figs -p C_real -p2 C_imag -s -g --gskip 2 --gscale 15 "
check_call(cmd, shell=True)

# 2. assemble the panels into a publication-ready file
cmd = "montage -label '(a) |A|' figs/inv_params_osi455_nh_200301-202012_1day_jul_absA.png -label '(b) Turning angle' figs/inv_params_osi455_nh_200301-202012_1day_jul_ThetaA.png -label '(c) Geostrophic current' figs/inv_params_osi455_nh_200301-202012_1day_jul_C_real_C_imag.png -tile 3x1 -geometry '622x586>+5+5' -pointsize 23 figs/freedriftparameter_nh_3panels.png"
check_call(cmd, shell=True)

cmd = "montage -label '(a) |A|' figs/inv_params_osi455_sh_200208-202007_1day_dec_absA.png -label '(b) Turning angle' figs/inv_params_osi455_sh_200208-202007_1day_dec_ThetaA.png -label '(c) Geostrophic current' figs/inv_params_osi455_sh_200208-202007_1day_dec_C_real_C_imag.png -tile 3x1 -geometry '622x559>+5+5' -pointsize 23 figs/freedriftparameter_sh_3panels.png"
check_call(cmd, shell=True)

cmd = 'cp figs/freedriftparameter_nh_3panels.png {}/fig{:02d}.png'.format(final_dir, fig_cnt)
check_call(cmd, shell=True)

fig_cnt += 1
cmd = 'cp figs/freedriftparameter_sh_3panels.png {}/fig{:02d}.png'.format(final_dir, fig_cnt)
check_call(cmd, shell=True)


# Example maps of t0, t1, and tdiff
# =================================
fig_cnt += 1
# 1. run the notebook/script
rdate = {'nh': datetime.strptime('20150226', '%Y%m%d'),
         'sh': datetime.strptime('20150826', '%Y%m%d')}
sens = 'amsr2-gw1'
for hemi in ['nh', 'sh']:
    if sens != 'multi-oi':
        fname = "https://thredds.met.no/thredds/dodsC/osisaf/met.no/reprocessed/ice/drift_455m_files/single_sensor/{}/{:%Y}/{:%m}/ice_drift_{}_ease2-750_cdr-v1p0-{}_24h-{:%Y%m%d}1200.nc".format(sens, rdate[hemi], rdate[hemi], hemi, sens, rdate[hemi])
    else:
        fname="https://thredds.met.no/thredds/dodsC/osisaf/met.no/reprocessed/ice/drift_455m_files/merged/{:%Y}/{:%m}/ice_drift_{}_ease2-750_cdr-v1p0_24h-{:%Y%m%d}1200.nc".format(rdate[hemi], rdate[hemi], hemi, rdate[hemi])
    sdate = rdate[hemi] - timedelta(days=1)
    titlestr="{} / {:%Y-%m-%d} to {:%Y-%m-%d}".format(sens.upper(), sdate, rdate[hemi])
    cmd = "python software/ql_figure.py -r ease-{}-wide -bg t0 -bgf {} -o figs --colbar --title --custom '{}' --inv_y".format(hemi, fname, titlestr)
    check_call(cmd, shell=True)
    cmd = "python software/ql_figure.py -r ease-{}-wide -bg t1 -bgf {} -o figs --colbar --title --custom '{}' --inv_y".format(hemi, fname, titlestr)
    check_call(cmd, shell=True)
    cmd = "python software/ql_figure.py -r ease-{}-wide -bg t0 -bgf {} -bg2 t1 -bgf2 {} -o figs --colbar --title --custom '{}' --inv_y".format(hemi, fname, fname, titlestr)
    check_call(cmd, shell=True)

# 2. assemble the panels into a publication-ready file
nh_files = ['figs/ice_drift_{a:}_ease2-750_cdr-v1p0-amsr2-gw1_24h-{r:%Y%m%d}1200_{t:}_simple_ease-{a:}-wide.png'.format(t=tstr, a='nh', r=rdate['nh']) for tstr in ('t0','t1','tdiff')]
cmd = 'convert {} +append {}/fig_time_nh.png'.format(" ".join(nh_files), tmp_dir)
check_call(cmd, shell=True)
sh_files = ['figs/ice_drift_{a:}_ease2-750_cdr-v1p0-amsr2-gw1_24h-{r:%Y%m%d}1200_{t:}_simple_ease-{a:}-wide.png'.format(t=tstr, a='sh', r=rdate['sh']) for tstr in ('t0','t1','tdiff')]
cmd = 'convert {} +append {}/fig_time_sh.png'.format(" ".join(sh_files), tmp_dir)
check_call(cmd, shell=True)
cmd = 'convert {}/fig_time_nh.png {}/fig_time_sh.png -smush -150 {}/fig{:02d}.png'.format(tmp_dir, tmp_dir, final_dir, fig_cnt)
check_call(cmd, shell=True)

# Example maps of dX and dY
# =================================
fig_cnt += 1
# 1. run the notebook/script
rdate = {'nh': datetime.strptime('20100117', '%Y%m%d'),
         'sh': datetime.strptime('20100814', '%Y%m%d')}
reg = {'nh': 'ease-eur-trim',
       'sh': 'ease-wed'}
sens = 'multi-oi'
for hemi in ['nh', 'sh']:
    if sens != 'multi-oi':
        fname = "https://thredds.met.no/thredds/dodsC/osisaf/met.no/reprocessed/ice/drift_455m_files/single_sensor/{}/{:%Y}/{:%m}/ice_drift_{}_ease2-750_cdr-v1p0-{}_24h-{:%Y%m%d}1200.nc".format(sens, rdate[hemi], rdate[hemi], hemi, sens, rdate[hemi])
    else:
        fname="https://thredds.met.no/thredds/dodsC/osisaf/met.no/reprocessed/ice/drift_455m_files/merged/{:%Y}/{:%m}/ice_drift_{}_ease2-750_cdr-v1p0_24h-{:%Y%m%d}1200.nc".format(rdate[hemi], rdate[hemi], hemi, rdate[hemi])
    sdate = rdate[hemi] - timedelta(days=1)
    titlestr="{} / {:%Y-%m-%d} to {:%Y-%m-%d}".format(sens.upper(), sdate, rdate[hemi])
    cmd = "python software/ql_figure.py -r {} -d {} -bg dX -bgf {} -o figs --colbar --title --scale 10 --skip 1 --custom '{}' --latlongrid --inv_y".format(reg[hemi], fname, fname, titlestr)
    check_call(cmd, shell=True)
    cmd = "python software/ql_figure.py -r {} -d {} -bg dY -bgf {} -o figs --colbar --title --scale 10 --skip 1 --custom '{}' --latlongrid --inv_y".format(reg[hemi], fname, fname, titlestr)
    check_call(cmd, shell=True)

# 2. assemble the panels into a publication-ready file
nh_files = ['figs/ice_drift_{a:}_ease2-750_cdr-v1p0_24h-{r:%Y%m%d}1200_{d:}_simple_{q:}.png'.format(d=dstr, a='nh', r=rdate['nh'], q=reg['nh']) for dstr in ('dX','dY',)]
cmd = 'convert {} +append {}/fig_dxdy_nh.png'.format(" ".join(nh_files), tmp_dir)
check_call(cmd, shell=True)
sh_files = ['figs/ice_drift_{a:}_ease2-750_cdr-v1p0_24h-{r:%Y%m%d}1200_{d:}_simple_{q:}.png'.format(d=dstr, a='sh', r=rdate['sh'], q=reg['sh']) for dstr in ('dX','dY',)]
cmd = 'convert {} +append {}/fig_dxdy_sh.png'.format(" ".join(sh_files), tmp_dir)
check_call(cmd, shell=True)
cmd = 'convert {}/fig_dxdy_nh.png {}/fig_dxdy_sh.png -smush -150 {}/fig{:02d}.png'.format(tmp_dir, tmp_dir, final_dir, fig_cnt)
check_call(cmd, shell=True)

# Validation results (hexdXdY) for ssmi-f11 and amsr2-gw1 (winter)
# =================================
fig_cnt += 1
# 1. run the notebook/script

# 2. assemble the panels into a publication-ready file
nh_files = ['figs/{s:}_{A:}_winter{a:}_hexdXdY.png'.format(s=sid, a='nh', A='NH') for sid in ('ssmi-f11', 'ssmis-f17', 'amsr2_tb37')]
cmd = 'convert {} +append {}/fig_winterhexdXdY_nh.png'.format(' '.join(nh_files), tmp_dir)
check_call(cmd, shell=True)
sh_files = ['figs/{s:}_{A:}_winter{a:}_hexdXdY.png'.format(s=sid, a='sh', A='SH') for sid in ('ssmi-f11', 'ssmis-f17', 'amsr2_tb37')]
cmd = 'convert {} +append {}/fig_winterhexdXdY_sh.png'.format(' '.join(sh_files), tmp_dir)
check_call(cmd, shell=True)
cmd = 'convert {}/fig_winterhexdXdY_nh.png {}/fig_winterhexdXdY_sh.png -append {}/fig{:02d}.png'.format(tmp_dir, tmp_dir, final_dir, fig_cnt)
check_call(cmd, shell=True)

# Validation results (hexdXdY) for wind-driven motion (summer)
# =================================
fig_cnt += 1
# 1. run the notebook/script

# 2. assemble the panels into a publication-ready file
cmd = 'convert figs/wind_pamb_ssmiamsr_NH_summernh_hexdXdY.png figs/wind_pamb_ssmiamsr_SH_summersh_hexdXdY.png -append {}/fig{:02d}.png'.format(final_dir, fig_cnt)
check_call(cmd, shell=True)

# Validation results (hexdXdY) for multi-source motion (year-round)
# =================================
fig_cnt += 1
# 1. run the notebook/script

# 2. assemble the panels into a publication-ready file
cmd = 'convert figs/mrg_wind_pamb_ssmiamsr_NH_hexdXdY.png figs/mrg_wind_pamb_ssmiamsr_SH_hexdXdY.png -append {}/fig{:02d}.png'.format(final_dir, fig_cnt)
check_call(cmd, shell=True)

# Validation results (timeseries) for dx (Arctic)
# =================================
fig_cnt += 1
# 1. run the notebook/script

# 2. assemble the panels into a publication-ready file
cmd = 'convert figs/timeseries_winter_nh_dx_cropped.png figs/timeseries_summer_nh_dx_cropped.png +append {}/fig{:02d}.png'.format(final_dir, fig_cnt)
check_call(cmd, shell=True)
