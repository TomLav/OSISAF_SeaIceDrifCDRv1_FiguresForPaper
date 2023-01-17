edate=$1
hemi=$2
sensor=$3
sdate=$(date --date="${edate} -1 day" +%Y%m%d)

ey=$(date --date="$edate" +%Y)
em=$(date --date="$edate" +%m)

sdatestr=$(date --date="${sdate}" +%Y-%m-%d)
edatestr=$(date --date="${edate}" +%Y-%m-%d)

titlestr="${sensor^^} / ${sdatestr} to ${edatestr}"

if [[ $sensor != "multi-oi" ]]
then 
fname="https://thredds.met.no/thredds/dodsC/osisaf/met.no/reprocessed/ice/drift_455m_files/single_sensor/amsr2-gw1/${ey}/${em}/ice_drift_${hemi}_ease2-750_cdr-v1p0-${sensor}_24h-${edate}1200.nc"
else
fname="https://thredds.met.no/thredds/dodsC/osisaf/met.no/reprocessed/ice/drift_455m_files/merged/${ey}/${em}/ice_drift_${hemi}_ease2-750_cdr-v1p0_24h-${edate}1200.nc"
fi

if [[ $hemi == "nh" ]]
then
reg="ease-eur"
else
reg="ease-wed"
fi

python ql_figure.py -r ${reg} -d ${fname} -bg dX -bgf ${fname} -o ../figs --colbar --title --scale 10 --skip 1 --custom "${titlestr}" --latlongrid --inv_y

python ql_figure.py -r ${reg} -d ${fname} -bg dY -bgf ${fname} -o ../figs --colbar --title --scale 10 --skip 1 --custom "${titlestr}" --latlongrid --inv_y
