#!/bin/sh

#---set parameters---
infile=response_functions.csv
psfile=vtf.ps

sizes=0.125
sizee=0.2
widthl=1
cimzx=255/0/255
crezx=255/0/0
crezy=0/0/255
cimzy=0/255/255

#---MAKE <FILE_ROOT>.ps---
# GMT 5
gmtset FONT_TITLE 14p
gmtset FONT_LABEL 14p
gmtset MAP_TICK_LENGTH 0.1c
gmtset MAP_FRAME_PEN 0.02c
# GMT 4
#gmtset HEADER_FONT_SIZE 14p
#gmtset LABEL_FONT_SIZE 14p
#gmtset ANOT_FONT_SIZE 14p

#---PLOT COHERENCES ---
awk -F',' 'NR > 1 {print $2, $7}' ${infile} | \
psxy -JX9l/2 -R0.003/30000/0/1 -Ba1f3l:"Log(Period,sec)":/a0.2f0.2:"Coh@+2@+ Hz":Wsne \
-X3 -Y18 -W2 -K > ${psfile}

# FOR Re(ZX)-COMPONENT 
# PLOT ERROR BARS
awk -F',' 'NR > 1 {print $2, $3, $8*1.96}' ${infile} | \
psxy -R0.003/30000/-1/1 -JX9l/6 -Ba1f3l/a0.5f0.5:"Real":Wsne \
-Ey${sizee}/1,${crezx} -Y-6.5 -O -K >> $psfile
# PLOT CIRCLES
awk -F',' 'NR > 1 {print $2, $3}' ${infile} | \
psxy -R -JX -Sc$sizes -G${crezx} -O -K  >> $psfile

# FOR Re(ZY)-COMPONENT 
# PLOT ERROR BARS
awk -F',' 'NR > 1 {print $2, $5, $9*1.96}' ${infile} | \
psxy -R -JX -Ey${sizee}/1,${crezy} -O -K  >> $psfile
# PLOT CIRCLES
awk -F',' 'NR > 1 {print $2, $5}' ${infile} | \
psxy -R -JX -Sc$sizes -G${crezy} -O -K  >> $psfile

# FOR Im(ZX)-COMPONENT 
# PLOT ERROR BARS
awk -F',' 'NR > 1 {print $2, $4, $8*1.96}' ${infile} | \
psxy -R0.003/30000/-1/1 -JX9l/6 -Ba1f3l:"Log(Period,sec)":/a0.5f0.5:"Imaginary":WSne \
-Ey${sizee}/1,${cimzx} -Y-6.5 -O -K  >> $psfile
# PLOT CIRCLES
awk -F',' 'NR > 1 {print $2, $4}' ${infile} | \
psxy -R -JX -Sc$sizes -G${cimzx} -O -K  >> $psfile

# FOR Im(ZY)-COMPONENT
# PLOT ERROR BARS
awk -F',' 'NR > 1 {print $2, $6, $9*1.96}' ${infile} | \
psxy -R -JX -Ey${sizee}/1,${cimzy} -O -K  >> $psfile
# PLOT CIRCLES
awk -F',' 'NR > 1 {print $2, $6}' ${infile} | \
psxy -R -JX -Sc$sizes -G${cimzy} -O -K  >> $psfile

echo  0 0 BL 0.0 13p,Helvetica,${crezx} "Re(ZX)" > tmp
echo  5 0 BL 0.0 13p,Helvetica,${cimzx} "Im(ZX)" >> tmp
echo 10 0 BL 0.0 13p,Helvetica,${crezy} "Re(ZY)" >> tmp
echo 15 0 BL 0.0 13p,Helvetica,${cimzy} "Im(ZY)" >> tmp
pstext tmp -R0/12/0/2 -JX4/2 -F+j+a+f -N -O -X1 -Y-2 >> $psfile
