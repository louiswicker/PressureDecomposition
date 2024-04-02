#!/bin/bash

st=20  # Starting CM1 output number for looping
en=20  # Ending CM1 output number for looping

for(( i=$st; i<=$en; i++ )); do
  time

  ncks -O -v zh,zf,th,qv,qc,qr,w,th0,qv0,u,u0,v,v0,prs,rho,prs0 -d time,$i cm1out.nc cm1tmp0.nc
  ncwa -O -a time cm1tmp0.nc cm1tmp.nc

# cp /path/to/pdcomp/directory/run/def.pdcomp.input /path/to/pdcomp/directory/pdcomp/run/pdcomp.input

# sed -i "s/cm1out_000001/cm1out_$i/g" /path/to/pdcomp/directory/pdcomp/run/pdcomp.input 
# sed -i "s/pdcomp_000001/pdcomp_$i/g" /path/to/pdcomp/directory/pdcomp/run/pdcomp.input

 ./pdcomp.exe
 
  mv pdcomp.nc pdcomp_$i.nc
  wait
done

time
