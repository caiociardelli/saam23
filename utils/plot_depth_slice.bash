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

# PLOT_DEPTH_SLICE

# USAGE
#   ./utils/plot_depth_slice.bash PARAMETER DEPTH CBMIN CBMAX

# EXAMPLE
#   ./utils/plot_depth_slice.bash vsv 100

# COMMAND LINE ARGUMENTS
#   PARAMETER              - model parameter to be plotted (cb, vpvs, ti, dcb, dvsv, etc.)
#   DEPTH                  - depth in which you created the depth slice
#   CBMIN (optional)       - minimum of the color bar
#   CBMAX (optional)       - maximum of the color bar

# DESCRIPTION
#   Plots the depth slices for the grid files created by GLL2LL or by CREATE_EXTRA_DEPTH_SLICE.

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
  echo " PLOT_DEPTH_SLICE"
  echo ""
  echo " USAGE"
  echo "   ./utils/plot_depth_slice.bash PARAMETER DEPTH CBMIN CBMAX"
  echo ""
  echo " EXAMPLE"
  echo "   ./utils/plot_depth_slice.bash vsv 100"
  echo ""
  echo " COMMAND LINE ARGUMENTS"
  echo "   PARAMETER              - model parameter to be plotted (cb, vpvs, ti, dcb, dvsv, etc.)"
  echo "   DEPTH                  - depth in which you created the depth slice"
  echo "   CBMIN (optional)       - minimum of the color bar"
  echo "   CBMAX (optional)       - maximum of the color bar"
  echo " DESCRIPTION"
  echo "   Plots the depth slices for the grid files created by GLL2LL or by CREATE_EXTRA_DEPTH_SLICE."
  echo ""
}

if [ "$#" -ne 2 ] && [ "$#" -ne 4 ]; then
  help
  exit 1
fi

filename="$1"_"$2"_DS.dat
grdname="$1"_"$2"_DS.grd
output="$1"_"$2"_DS

np=$(head -n1 $filename | cut -f7 -d' ')
nt=$(head -n1 $filename | cut -f6 -d' ')

dp=$(echo "122 / ($np - 1)" | bc -l)
dt=$(echo " 92 / ($nt - 1)" | bc -l)

info=$(./bin/getinfo $filename)

vmin=$(echo $info | gawk '{printf "%lf", $3}')
vmax=$(echo $info | gawk '{printf "%lf", $4}')
mean=$(echo $info | gawk '{printf "%lf", $5}')
stdv=$(echo $info | gawk '{printf "%lf", $6}')

cbmin=$(echo "$mean - 3.5 * $stdv" | bc -l)
cbmax=$(echo "$mean + 3.5 * $stdv" | bc -l)

if [ "$#" -eq 4 ]; then
  cbmin=$3
  cbmax=$4
fi

gmt xyz2grd $filename -G$grdname -R-110/12/-67/25 -I$dp/$dt -: -h3
gmt makecpt -Cextra/tomo_improved.cpt -T$cbmin/$cbmax > tomo.cpt

gawk -v min=$vmin -v max=$vmax 'BEGIN {printf "Min = %E Max = %E\n", min, max}'
echo 'Creating figure...'

gawk -F"," -v depth="$2" '{if ($3 + depth < 0.2 && $3 + depth > -0.2) {print $1-360,$2}}' extra/sam_slab2_dep_02.23.18.xyz > extra/slab_"$2".xy

gmt begin $output pdf
  gmt set FONT_LABEL 16p,100
  gmt set MAP_FRAME_PEN 0.5p

#  gmt grdimage $grdname -JN10c -R-88/-10/-60/16 -Bxa30fg30 -Bya30fg30 -BWeSn -Ctomo.cpt
  gmt grdimage $grdname -JQ10c -R-88/-10/-60/16 -Bxa30fg30 -Bya30fg30 -BWeSn -Ctomo.cpt
  gmt coast -W0.2,50 -A30000

  gmt plot extra/plate_boundaries_clean.gmt -: -W0.7p
  gmt plot extra/slab_"$2".xy -W1p,magenta

  gmt colorbar -Ctomo.cpt -Baf -DJBC+e -B+l"$label - $2 km"
gmt end

