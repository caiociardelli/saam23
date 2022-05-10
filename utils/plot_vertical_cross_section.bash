#!/usr/bin/env bash

# SphGLLTools

# Author: Caio Ciardelli, University of São Paulo, October 2020

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

#-----------------------------------------------------------------------------------------------

# PLOT_VERTICAL_CROSS_SECTION

# USAGE
#   ./utils/plot_vertical_cross_section.bash PARAMETER CBMIN CBMAX

# EXAMPLE
#   ./utils/plot_vertical_cross_section.bash vpvs

# COMMAND LINE ARGUMENTS
#   PARAMETER              - model parameter to be plotted (cb, vpvs, ti, dcb, dvsv, etc.)
#   CBMIN (optional)       - minimum of the color bar
#   CBMAX (optional)       - maximum of the color bar

# DESCRIPTION
#   Plots the vertical cross-section for the grid files created by GLL2DD or by
#   CREATE_EXTRA_VERTICAL_CROSS_SECTION.

#-----------------------------------------------------------------------------------------------

if [ $1 = "vp" ]; then
  label="@~\141@~ (km s@+-1@+)"
elif [ $1 = "vpv" ]; then
  label="@~\141@~@-v@- (km s@+-1@+)"
elif [ $1 = "vph" ]; then
  label="@~\141@~@-h@- (km s@+-1@+)"
elif [ $1 = "vs" ]; then
  label="@~\142@~ (km s@+-1@+)"
elif [ $1 = "vsv" ]; then
  label="@~\142@~@-v@- (km s@+-1@+)"
elif [ $1 = "vsh" ]; then
  label="@~\142@~@-h@- (km s@+-1@+)"
elif [ $1 = "eta" ]; then
  label="@~\150@~"
elif [ $1 = "rho" ]; then
  label="@~\162@~ (g cm@+-3@+)"
elif [ $1 = "qmu" ]; then
  label="Q@-@~\155@~@-"
elif [ $1 = "cb" ]; then
  label="C@-Bulk@- (km s@+-1@+)"
elif [ $1 = "ti" ]; then
  label="T@-i@- (%)"
elif [ $1 = "vpvs" ]; then
  label="@~\141@~/@~\142@~"
elif [ $1 = "dvp" ]; then
  label="@~\144\141@~ (%)"
elif [ $1 = "dvpv" ]; then
  label="@~\144\141@~@-v@- (%)"
elif [ $1 = "dvph" ]; then
  label="@~\144\141@~@-h@- (%)"
elif [ $1 = "dvs" ]; then
  label="@~\144\142@~ (%)"
elif [ $1 = "dvsv" ]; then
  label="@~\144\142@~@-v@- (%)"
elif [ $1 = "dvsh" ]; then
  label="@~\144\142@~@-v@- (%)"
elif [ $1 = "deta" ]; then
  label="dln@~\150@~ (%)"
elif [ $1 = "drho" ]; then
  label="dln@~\162@~ (%)"
elif [ $1 = "dqmu" ]; then
  label="dQ@-@~\155@~@- (%)"
elif [ $1 = "dcb" ]; then
  label="dlnC@-Bulk@- (%)"
else
  echo "Error: $1 parameter is not available!"
  exit 1
fi

help()
{
  echo " PLOT_VERTICAL_CROSS_SECTION"
  echo ""
  echo " USAGE"
  echo "   ./utils/plot_vertical_cross_section.bash PARAMETER"
  echo ""
  echo " EXAMPLE"
  echo "   ./utils/plot_vertical_cross_section.bash vp"
  echo ""
  echo " COMMAND LINE ARGUMENTS"
  echo "   PARAMETER              - model parameter to be plotted (cb, vpvs, ti, dcb, dvsv, etc.)"
  echo "   CBMIN (optional)       - minimum of the color bar"
  echo "   CBMAX (optional)       - maximum of the color bar"
  echo ""
  echo " DESCRIPTION"
  echo "   Plots the vertical cross-section for the grid files created by GLL2DD or by"
  echo "   CREATE_EXTRA_VERTICAL_CROSS_SECTION."
  echo ""
}

REARTH=6371
R410=5961
R650=5721
R1000=5371
KM2DEGREE=111.2
TOPOSCALE=0.08

filename="$1"_VCS.dat
grdname="$1"_VCS.grd
output="$1"_VCS

nd=$(head -n1 $filename | cut -f4 -d' ')
nr=$(head -n1 $filename | cut -f3 -d' ')

line=$(sed '2q;d' $filename)

lat1=$(echo $line | gawk '{print $5}')
lon1=$(echo $line | gawk '{print $6}')
lat2=$(echo $line | gawk '{print $7}')
lon2=$(echo $line | gawk '{print $8}')

lat0=$(echo "($lat1 + $lat2) / 2" | bc -l)
lon0=$(echo "($lon1 + $lon2) / 2" | bc -l)

info=$(./bin/getinfo $filename)

mean=$(echo $info | gawk '{print $5}')
stdv=$(echo $info | gawk '{print $6}')

if [ "$#" -eq 3 ]; then
  cbmin=$2
  cbmax=$3
else
  cbmin=$(echo "$mean - 3 * $stdv" | bc -l)
  cbmax=$(echo "$mean + 3 * $stdv" | bc -l)
fi

string=$(gmt info $filename)

dstring=$(echo "$string" | cut -f3 | gawk -F '[</>]' '{print $2" "$3}')
rstring=$(echo "$string" | cut -f2 | gawk -F '[</>]' '{print $2" "$3}')
vstring=$(echo "$string" | cut -f4 | gawk -F '[</>]' '{print $2" "$3}')

drange=$(echo $dstring | gawk '{print $2 - $1}')
rrange=$(echo $rstring | gawk '{print $2 - $1}')

dd=$(echo "$drange / ($nd - 1)" | bc -l)
dr=$(echo "$rrange / ($nr - 1)" | bc -l)

dmin=$(echo "$dstring" | gawk '{print $1}')
dmax=$(echo "$dstring" | gawk '{print $2}')
rmin=$(echo "$rstring" | gawk '{print $1}')
rmax=$(echo "$rstring" | gawk '{print $2}')
vmin=$(echo "$vstring" | gawk '{print $1}')
vmax=$(echo "$vstring" | gawk '{print $2}')

ermax=$(echo "1.01 * $rmax" | bc -l)

dmean=$(echo "($dmin + $dmax) / 2" | bc -l)
rmean=$(echo "($rmin + $rmax) / 2" | bc -l)

dlabeld=$(echo "$dmean" | bc -l)
dlabelr=$(echo "0.75 * $rmin" | bc -l)
rlabeld=$(echo "$dmin - 0.25 * $drange" | bc -l)
rlabelr=$(echo "1.2 * $rmean" | bc -l)

step=$(echo "111 * $dmax / 6" | bc -l)

gmt project -C$lon1/$lat1 -E$lon2/$lat2 -G1 -Q > great_circle_points.xyp
gmt project -C$lon1/$lat1 -E$lon2/$lat2 -G$step -Q > great_circle_points_dots.xyp
gmt grdtrack -Gextra/earth_relief_15m.grd great_circle_points.xyp | gawk '{print $3" "$4}' > profile.dat

gawk -F"," '$3 != "NaN" {print $1-360,$2,$3}' extra/sam_slab2_dep_02.23.18.xyz | gmt project -C$lon1/$lat1 -E$lon2/$lat2 -Q -W-0.2/0.2 | gawk '{print $4,-$3}' > extra/slab_temp.gmt

gawk '{print $1/111.1,6371-$2}' extra/slab_temp.gmt > extra/slab.gmt

string=$(gmt info profile.dat | cut -f3 | gawk -F '[</>]' '{print $2" "$3}')

emin=$(echo $string | gawk '{print $1}')
emax=$(echo $string | gawk '{print $2}')
erange=$(echo $string | gawk '{print $2 - $1}')
scale=$(echo "$erange / ($TOPOSCALE * $rrange)" | bc -l)
shift=$(echo "((0.2 * $erange) - $emin) / $scale" | bc -l)
step=$(echo "($dmax - $dmin) / 4000" | bc -l)
Rmax=$(echo "1.015 * ($rmax + $shift)" | bc -l)

gawk -v km2dg=$KM2DEGREE -v rmax=$rmax -v shift=$shift -v scale=$scale\
    '{print $1 / km2dg" "rmax + shift + $2 / scale}' profile.dat > profile_rescaled1.dat
cp profile_rescaled1.dat profile_rescaled2.dat

gawk -v dmin=$dmin -v dmax=$dmax -v step=$step -v rmax=$rmax\
    'BEGIN {for (delta = dmax; delta >= 0; delta -= step) printf "%.3f %.1f\n", delta, rmax}' >> profile_rescaled1.dat
gawk -v dmin=$dmin -v dmax=$dmax -v step=$step -v rmax=$Rmax\
    'BEGIN {for (delta = dmax; delta >= 0; delta -= step) printf "%.3f %.1f\n", delta, rmax}' >> profile_rescaled2.dat

gawk -v min=$vmin -v max=$vmax 'BEGIN {printf "Min = %E Max = %E\n", min, max}'
echo 'Creating figure...'

gmt xyz2grd $filename -G$grdname -R$dmin/$dmax/$rmin/$rmax -I$dd/$dr -: -h4
gmt makecpt -Cextra/c002.cpt -T-8/8 > tomo.cpt

gmt begin $output pdf
  gmt set FONT_LABEL 14p,100
  gmt set FONT_ANNOT 14p,100
  gmt set MAP_FRAME_PEN 0.5p
  gmt set MAP_LABEL_OFFSET 0.1c

  gmt grdimage $grdname -JPa10c/$dmean -R$dmin/$dmax/$rmin/$Rmax -BweS -Bxa10+l"Distance (km)" -By+l"Depth (km)" -Ctomo.cpt

  gmt plot profile_rescaled1.dat -W0.02c -L -Glightgray
  gmt plot profile_rescaled2.dat -W0.02c -L -Gwhite

  gmt plot -W0.015c,0,-- << DISCONTINUITIES
$(for distance in $(seq 0 10000); do echo $distance $R410; done)
$(for distance in $(seq 0 10000); do echo $distance $R650; done)
$(for distance in $(seq 0 10000); do echo $distance $R1000; done)
DISCONTINUITIES

  step=$(echo "$dmax / 6" | bc -l)

  for distance in $(seq 0.0 $step $dmax);
    do
      echo "$distance $Rmax" | sed -r 's/,+/./g' | gmt plot -Sc0.15c -Ggold -W0.02c,black -N;
  done

  echo "$dmin $Rmax" | gmt plot -Sc0.2c -Ggreen -W0.02c,black -N
  echo "$dmax $Rmax" | gmt plot -Sc0.2c -Gmagenta -W0.02c,black -N

  gmt text -F+a$dmean+jRM -N <<< "-0.5 $R410 410"
  gmt text -F+a$dmean+jRM -N <<< "-0.5 $R650 650"
  gmt text -F+a$dmean+jRM -N <<< "-0.5 $R1000 1000"

  gmt colorbar -Ctomo.cpt -Ba3f -DJBC+e -B+l"@:13p:$label@::"

  gmt inset begin -DJTR+w10c -M0.1c
    gmt coast -JG$lon0/$lat0/2.5c -Rg -Bafg -Dc -A10000 -Glightgray -Wthinnest

    gmt plot extra/plate_boundaries.dat -JG$lon0/$lat0/2.5c -: -W0.6p,brown
    gmt plot great_circle_points.xyp -JG$lon0/$lat0/2.5c -W0.5p
    gmt plot great_circle_points_dots.xyp -JG$lon0/$lat0/2.5c -Sc0.075c -Ggold -W0.02c,black

    echo "$lon1 $lat1" | gmt plot -JG$lon0/$lat0/2.5c -Sc0.1c -Ggreen -W0.02c,black -N
    echo "$lon2 $lat2" | gmt plot -JG$lon0/$lat0/2.5c -Sc0.1c -Gmagenta -W0.02c,black -N
  gmt inset end
gmt end

echo 'Figure saved!'
