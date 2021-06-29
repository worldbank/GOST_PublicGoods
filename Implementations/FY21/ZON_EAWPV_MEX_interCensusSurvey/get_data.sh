#! /bin/bash
# Author: Diego Valle-Jones
# Web: https://www.diegovalle.net
# LICENSE: Apache 2.0
# Purpose: Download shapefiles of manzanas (blocks), agebs (census areas), ejes
# viales (streets), interesting areas and a whole bunch of other stuff from
# the Encuesta Intercensal 2015


# As of now, this script has not been tested on Windows, 
# only on Ubuntu and Macs. The script will create a directory
# called 'shps' where all the shapefiles are located, if something
# goes wrong when dowloading be sure to delete it and try again

#Exit on error
set -e

# Projection compatible with Google Maps
PROJECTION="+proj=longlat +ellps=WGS84 +no_defs +towgs84=0,0,0"
# wget command
CURL="curl"

#The way the INEGI download url starts
URL="http://internet.contenidos.inegi.org.mx/contenidos/productos//prod_serv/contenidos/espanol/bvinegi/productos/geografia/Cinter_2015/"
#File ends in
SUFFIX="_s.zip"

#States as numbers
#I have no idea why the INEGI used these numbers to represent the states
declare -a states_numbers=("702825209025" "702825209032" "702825209049" "702825209056" "702825209087" "702825209094" "702825209063" "702825209070" "702825209100" "702825209117" "702825209124" "702825209131" "702825209148" "702825209155" "702825209162" "702825209179" "702825209186" "702825209193" "702825209209" "702825209216" "702825209223" "702825209230" "702825209247" "702825209254" "702825209261" "702825209278" "702825209285" "702825209292" "702825209308" "702825209315" "702825209322" "702825209339")

#State names in the URL
declare -a state_names=("Aguascalientes" "Baja_California" "Baja_California_Sur" "Campeche" "Coahuila_de_Zaragoza" "Colima" "Chiapas" "Chihuahua" "Distrito_Federal" "Durango" "Guanajuato" "Guerrero" "Hidalgo" "Jalisco" "Mexico" "Michoacan_de_Ocampo" "Morelos" "Nayarit" "Nuevo_Leon" "Oaxaca" "Puebla" "Queretaro" "Quintana_Roo" "San_Luis_Potosi" "Sinaloa" "Sonora" "Tabasco" "Tamaulipas" "Tlaxcala" "Veracruz_de_Ignacio_de_la_Llave" "Yucatan" "Zacatecas" )

# The list of shapefiles of manzanas, agebs, etc
declare -a files=("ar.shp" "ent.shp" "lpr.shp" "m.shp" "fm.shp" "sia.shp" "sip.shp"
                   "a.shp" "e.shp" "l.shp" "mun.shp" "sil.shp" "territorioinsular.shp");

#The INEGI filenames are unreadable
declare -a files_nice=("ageb_rural.shp" "entidad.shp" "localidad_rural_no_amanzanada.shp" "manzana.shp"  "frente_de_manzana.shp" "servicios_area.shp" "servicios_puntual.shp"
                   "ageb_urbana.shp" "eje_vial.shp" "localidad_urbana_y_rural_amanzanada.shp" "municipio.shp" "servicios_linea.shp" "territorio_insular.shp");


# State abbreviations
declare -a states=("ags" "bc" "bcs" "camp" "coah" "col" "chis" "chih"
    "df" "dgo" "gto" "gro" "hgo" "jal" "mex" "mich" "mor" "nay" "nl" "oax"
    "pue" "qro" "qroo" "slp" "sin" "son" "tab" "tamps" "tlax" "ver" "yuc"
    "zac");

declare -a state_num=(`seq -s " " -w 1 32`);

# Use gdal to reproject, and then rename the shapefiles to include
# a user friendly abbreviation instead of a number
# First argument: directory of shapefiles shps/state_abbreviation
# Second argument: the state abbreviation
# Third argument: the shapefiles inside the zip file as an array
# Fourth argument: the state number
# TODO: convert the encoding from windows-1252 to utf-8
function reproject {
  name=$3[@]
  arr=("${!name}")
  len=`expr ${#files[*]} - 1`
  # Has to match he number of files in the array
  for i in `seq 0  $len`
  do
    if [ -f $1/conjunto_de_datos/$4${files[$i]} ];
    then
      echo "Creating... " "$1"/$2_${files_nice[$i]}
      ogr2ogr "$1"/$2_${files_nice[$i]} $1/conjunto_de_datos/$4${files[$i]} -t_srs "$PROJECTION"
    else
      echo "No territorio insular in " "$1"/$2_${files_nice[$i]}
    fi
  done
  rm -rf "$1"/conjunto_de_datos
}


counter=0
for state_number in "${states_numbers[@]}"
do
   FILENAME="$URL${state_names[$counter]}/$state_number$SUFFIX"
   # download the files from the inegi server.
   $CURL $FILENAME -o ${states[$counter]}.zip
   # Extract the shapefiles from the zip file
   mkdir -p shps/${states[$counter]}
   unzip -o -L ${states[$counter]}.zip -d shps/${states[$counter]}

   reproject shps/${states[$counter]} ${states[$counter]} files ${state_num[$counter]}
   counter=`expr $counter + 1`
done
# Delete the downloaded zip files
rm -rf *.zip

# You could use the code below to merge all the states into one giant
# shapefile of Mexico. Change "*localidad_urbana.shp" to
# '*ageb_urbana.shp' or '*eje_vial.shp' or whatever
#
# FILEOUT="municipio.shp"
# TYPE="*municipio.shp"
# for i in $(find . -maxdepth 3  -name $TYPE)
# do
#     if [ -f "$FILEOUT" ]
#     then
#         echo "adding state $i to $FILEOUT"
#         ogr2ogr -f 'ESRI Shapefile' -update -append $FILEOUT $i -nln $(basename -s .shp $FILEOUT)
#     else
#         echo "startin merge..."
#         echo "adding state $i to $FILEOUT"
#         ogr2ogr -f 'ESRI Shapefile' $FILEOUT $i
#     fi
# done