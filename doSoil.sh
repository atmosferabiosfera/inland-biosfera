#!/bin/bash

# Script para trocar os parametros de solo 
# Troca pelos parametros da Emily menos o SHC(saturated hydraulic conductivity)

# Ordem dos parametros: 
# Porosity, FC(field capacity), WP(wilting point), bexp, AEP(air entry potential)
# 1 - Agricultura; 2 - Cerrado; 3 - Floresta; 4 - Pastagem

ip=0
if [ "$1" == "agricultura" ]
then
    ip=1
fi
if [ "$1" == "cerrado" ]
then
    ip=2
fi
if [ "$1" == "floresta" ]
then
    ip=3
fi
if [ "$1" == "pastagem" ]
then
    ip=4
fi

param=(
""
"0.4108    0.1719   0.1177    5.299   0.1196   1.72790E-05"
"0.4813    0.1274   0.0848    4.175   0.0708   4.52397E-05"
"0.4619    0.1789   0.1109    5.103   0.0868   4.01978E-05"
"0.3762    0.1663   0.1052    5.001   0.1796   1.41598E-05"
)

inicio=(
"0.92   0.05   0.03    "
"0.81   0.12   0.07    "
"0.65   0.25   0.10    "
"0.42   0.40   0.18    "
"0.20   0.65   0.15    "
"0.60   0.13   0.27    "
"0.32   0.34   0.34    "
"0.09   0.58   0.33    "
"0.53   0.07   0.40    "
"0.10   0.45   0.45    "
"0.20   0.20   0.60    "
)

final=(
"  0.0495		! Sand		  1" 
"  0.0613		! Loamy Sand	  2"
"  0.1101		! Sandy Loam	  3"
"  0.0889		! Loam		  4"
"  0.1668		! Silty Loam	  5"
"  0.2185		! Sandy Clay Loam	  6"
"  0.2088		! Clay Loam	  7"
"  0.2730		! Silty Clay Loam	  8"
"  0.2390		! Sandy Clay	  9"
"  0.2922		! Silty Clay	 10"
"  0.3163		! Clay		 11"
)

for i in {0..10};
do
    sed -i -r "s/${inicio[i]}.*${final[i]}/${inicio[i]}${param[$ip]}${final[i]}/"  params/soil
done
