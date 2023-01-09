#!/usr/bin/env python

import os
import shutil
from subprocess import check_call

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
# 1. run the notebook

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
# 1. run the notebook

# 2. assemble the panels into a publication-ready file
cmd = 'cp figs/freedriftparameter_nh_3panels.png {}/fig{:02d}.png'.format(final_dir, fig_cnt)
check_call(cmd, shell=True)

fig_cnt += 1
cmd = 'cp figs/freedriftparameter_sh_3panels.png {}/fig{:02d}.png'.format(final_dir, fig_cnt)
check_call(cmd, shell=True)


# Example maps of t0, t1, and tdiff
# =================================
fig_cnt += 1
# 1. run the notebook

# 2. assemble the panels into a publication-ready file
nh_files = ['figs/ice_drift_{a:}_ease2-750_amsr2-gw1_201502251200-201502261200_{t:}_simple_{a:}-ease-wide.png'.format(t=tstr,a='nh') for tstr in ('t0','t1','tdiff')]
cmd = 'convert {} +append {}/fig_time_nh.png'.format(" ".join(nh_files), tmp_dir)
check_call(cmd, shell=True)
sh_files = ['figs/ice_drift_{a:}_ease2-750_amsr2-gw1_201502251200-201502261200_{t:}_simple_{a:}-ease-wide.png'.format(t=tstr,a='sh') for tstr in ('t0','t1','tdiff')]
cmd = 'convert {} +append {}/fig_time_sh.png'.format(" ".join(sh_files), tmp_dir)
check_call(cmd, shell=True)
cmd = 'convert {}/fig_time_nh.png {}/fig_time_sh.png -append {}/fig{:02d}.png'.format(tmp_dir, tmp_dir, final_dir, fig_cnt)
check_call(cmd, shell=True)

# Example maps of dX and dY
# =================================
fig_cnt += 1
# 1. run the notebook

# 2. assemble the panels into a publication-ready file
nh_files = ['figs/ice_drift_{a:}_ease2-750_cdr-v1p0_24h-201001171200_{d:}_simple_{a:}-ease-eur.png'.format(d=dstr,a='nh') for dstr in ('dX','dY',)]
cmd = 'convert {} +append {}/fig_dxdy_nh.png'.format(" ".join(nh_files), tmp_dir)
check_call(cmd, shell=True)
sh_files = ['figs/ice_drift_{a:}_ease2-750_cdr-v1p0_24h-201008141200_{d:}_simple_{a:}-ease-wed.png'.format(d=dstr,a='sh') for dstr in ('dX','dY',)]
cmd = 'convert {} +append {}/fig_dxdy_sh.png'.format(" ".join(sh_files), tmp_dir)
check_call(cmd, shell=True)
cmd = 'convert {}/fig_dxdy_nh.png {}/fig_dxdy_sh.png -append {}/fig{:02d}.png'.format(tmp_dir, tmp_dir, final_dir, fig_cnt)
check_call(cmd, shell=True)

# Validation results (hexdXdY) for ssmi-f11 and amsr2-gw1 (winter)
# =================================
fig_cnt += 1
# 1. run the notebook

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
# 1. run the notebook

# 2. assemble the panels into a publication-ready file
cmd = 'convert figs/wind_pamb_ssmiamsr_NH_summernh_hexdXdY.png figs/wind_pamb_ssmiamsr_SH_summersh_hexdXdY.png -append {}/fig{:02d}.png'.format(final_dir, fig_cnt)
check_call(cmd, shell=True)

# Validation results (hexdXdY) for multi-source motion (year-round)
# =================================
fig_cnt += 1
# 1. run the notebook

# 2. assemble the panels into a publication-ready file
cmd = 'convert figs/mrg_wind_pamb_ssmiamsr_NH_hexdXdY.png figs/mrg_wind_pamb_ssmiamsr_SH_hexdXdY.png -append {}/fig{:02d}.png'.format(final_dir, fig_cnt)
check_call(cmd, shell=True)

# Validation results (timeseries) for dx (Arctic) 
# =================================
fig_cnt += 1
# 1. run the notebook

# 2. assemble the panels into a publication-ready file
cmd = 'convert figs/timeseries_winter_nh_dx_cropped.png figs/timeseries_summer_nh_dx_cropped.png +append {}/fig{:02d}.png'.format(final_dir, fig_cnt)
check_call(cmd, shell=True)
