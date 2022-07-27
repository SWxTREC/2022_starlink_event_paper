#!/bin/tcsh
# Script to grab the GRACE-FO data

foreach year ( 2022 )
  mkdir -p `printf "%04d" $year`
  foreach mm ( `seq 1 2` )
    foreach dd ( `seq 1 31` )
      wget --no-check-certificate -q -nc "https://swarm-diss.eo.esa.int/?do=download&file=swarm/Multimission/GRACE-FO/DNS/Sat_1/"`printf "%04d/GF_OPER_DNS1ACC_2__%04d%02d%02dT000000_%04d%02d%02dT235959_0001.cdf" $year $year $mm $dd $year $mm $dd` -O `printf "%04d/GF_OPER_DNS1ACC_2__%04d%02d%02dT000000_%04d%02d%02dT235959_0001.cdf" $year $year $mm $dd $year $mm $dd`
    end
  end
end
