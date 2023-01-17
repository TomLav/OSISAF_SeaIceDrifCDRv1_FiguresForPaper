python plot_ncparam.py -i https://thredds.met.no/thredds/dodsC/metusers/thomasl/SIDrift_CDR_v1pre/auxiliary_files/inv_params_osi455_nh_200301-202012_1day.nc -m 7 -o ../figs -p absA_gapfill -c
python plot_ncparam.py -i https://thredds.met.no/thredds/dodsC/metusers/thomasl/SIDrift_CDR_v1pre/auxiliary_files/inv_params_osi455_nh_200301-202012_1day.nc -m 7 -o ../figs -p ThetaA_gapfill -c
#python plot_ncparam.py -i https://thredds.met.no/thredds/dodsC/metusers/thomasl/SIDrift_CDR_v1pre/auxiliary_files/inv_params_osi455_nh_200301-202012_1day.nc -m 7 -o ../figs -p C_real_gapfill -c
python plot_ncparam.py -i https://thredds.met.no/thredds/dodsC/metusers/thomasl/SIDrift_CDR_v1pre/auxiliary_files/inv_params_osi455_nh_200301-202012_1day.nc -m 7 -o ../figs -p RMS_Res_real -p2 RMS_Res_imag

python plot_ncparam.py -i https://thredds.met.no/thredds/dodsC/metusers/thomasl/SIDrift_CDR_v1pre/auxiliary_files/inv_params_osi455_sh_200208-202007_1day.nc -m 12 -o ../figs -p absA_gapfill -c
python plot_ncparam.py -i https://thredds.met.no/thredds/dodsC/metusers/thomasl/SIDrift_CDR_v1pre/auxiliary_files/inv_params_osi455_sh_200208-202007_1day.nc -m 12 -o ../figs -p ThetaA_gapfill -c
#python plot_ncparam.py -i https://thredds.met.no/thredds/dodsC/metusers/thomasl/SIDrift_CDR_v1pre/auxiliary_files/inv_params_osi455_sh_200208-202007_1day.nc -m 12 -o ../figs -p C_real_gapfill -c
python plot_ncparam.py -i https://thredds.met.no/thredds/dodsC/metusers/thomasl/SIDrift_CDR_v1pre/auxiliary_files/inv_params_osi455_sh_200208-202007_1day.nc -m 12 -o ../figs -p RMS_Res_real -p2 RMS_Res_imag

montage -label '(a) |A| parameter (gapfilled)' ../figs/inv_params_osi455_nh_200301-202012_1day_jul_absA_gapfill.png -label '(b) Turning angle' ../figs/inv_params_osi455_nh_200301-202012_1day_jul_ThetaA_gapfill.png -label '(c) Average RMS error' ../figs/inv_params_osi455_nh_200301-202012_1day_jul_RMS_Res_real_RMS_Res_imag.png -tile 3x1 -geometry '622x546>+5+5' -pointsize 23 ../figs/freedriftparameter_nh_3panels.png

montage -label '(a) |A| parameter (gapfilled)' ../figs/inv_params_osi455_sh_200208-202007_1day_dec_absA_gapfill.png -label '(b) Turning angle' ../figs/inv_params_osi455_sh_200208-202007_1day_dec_ThetaA_gapfill.png -label '(c) Average RMS error' ../figs/inv_params_osi455_sh_200208-202007_1day_dec_RMS_Res_real_RMS_Res_imag.png -tile 3x1 -geometry '622x546>+5+5' -pointsize 23 ../figs/freedriftparameter_sh_3panels.png

# 622 546
