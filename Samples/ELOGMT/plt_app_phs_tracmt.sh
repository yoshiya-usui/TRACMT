#!/bin/sh

#---set parameters---
infile=apparent_resistivity_and_phase.csv
psfile=app_phs.ps

sizes=0.125
sizee=0.2
widthl=1
red=255/0/0
blue=0/0/255
green=0/255/0
orange=255/128/0

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
psxy -JX9l/2 -R0.003/30000/0/1 -Ba1f3l:"Log(Period,sec)":/a0.2f0.2:"Coh@+2@+ Ex":Wsne \
-X3 -Y18 -W2 -K > ${psfile}

awk -F',' 'NR > 1 {print $2, $12}' ${infile} | \
psxy -JX9l/2 -R0.003/30000/0/1 -Ba1f3l:"Log(Period,sec)":/a0.2f0.2:"Coh@+2@+ Ey":Wsne \
-Y-2.5 -W2 -O -K >> ${psfile}

#---apparent resistivity---
# FOR XX-COMPONENT 
# PLOT ERROR BARS
awk -F',' 'NR > 1 {print $2, $3, $13*1.96}' ${infile} | \
psxy -R0.003/30000/0.1/10000 -JX9l/6l -Ba1f3l/a1f3l:"Log(App. Resistivity,@~W@~m)":Wsne \
-Ey${sizee}/1,${orange} -Y-6.5 -O -K  >> $psfile
# PLOT CIRCLES
awk -F',' 'NR > 1 {print $2, $3}' ${infile} | \
psxy -R -JX -Sc$sizes -G${orange} -O -K  >> $psfile

# FOR YY-COMPONENT 
# PLOT ERROR BARS
awk -F',' 'NR > 1 {print $2, $10, $19*1.96}' ${infile} | \
psxy -R -JX -Ey${sizee}/1,${green} -O -K  >> $psfile
# PLOT CIRCLES
awk -F',' 'NR > 1 {print $2, $10}' ${infile} | \
psxy -R -JX -Sc$sizes -G${green} -O -K  >> $psfile

# FOR YX-COMPONENT 
# PLOT ERROR BARS
awk -F',' 'NR > 1 {print $2, $8, $17*1.96}' ${infile} | \
psxy -R -JX -Ey${sizee}/1,${blue} -O -K  >> $psfile
# PLOT CIRCLES
awk -F',' 'NR > 1 {print $2, $8}' ${infile} | \
psxy -R -JX -Sc$sizes -G${blue} -O -K  >> $psfile

# FOR XY-COMPONENT 
# PLOT ERROR BARS
awk -F',' 'NR > 1 {print $2, $5, $15*1.96}' ${infile} | \
psxy -R -JX -Ey${sizee}/1,${red} -O -K  >> $psfile
# PLOT CIRCLES
awk -F',' 'NR > 1 {print $2, $5}' ${infile} | \
psxy -R -JX -Sc$sizes -G${red} -O -K  >> $psfile

#---phase---
# FOR XX-COMPONENT 
# PLOT ERROR BARS
awk -F',' 'NR > 1 {print $2, $4, $14*1.96}' ${infile} | \
psxy -R0.003/30000/-180/180 -JX9l/6 -Ba1f3l:"Log(Period,sec)":/a45f45:"Phase(deg.)":WSne \
-Ey${sizee}/1,${orange} -Y-6.5 -K -O  >> $psfile
# PLOT CIRCLES
awk -F',' 'NR > 1 {print $2, $4}' ${infile} | \
psxy -R -JX -Sc$sizes -G${orange} -O -K  >> $psfile

# FOR YY-COMPONENT 
# PLOT ERROR BARS
awk -F',' 'NR > 1 {print $2, $11, $20*1.96}' ${infile} | \
psxy -R -JX -Ey${sizee}/1,${green} -K -O  >> $psfile
# PLOT CIRCLES
awk -F',' 'NR > 1 {print $2, $11}' ${infile} | \
psxy -R -JX -Sc$sizes -G${green} -O -K  >> $psfile

# FOR YX-COMPONENT 
# PLOT ERROR BARS
awk -F',' 'NR > 1 {print $2, $9, $18*1.96}' ${infile} | \
psxy -R -JX -Ey${sizee}/1,${blue} -K -O  >> $psfile
# PLOT CIRCLES
awk -F',' 'NR > 1 {print $2, $9}' ${infile} | \
psxy -R -JX -Sc$sizes -G${blue} -O -K  >> $psfile

# FOR XY-COMPONENT 
# PLOT ERROR BARS
awk -F',' 'NR > 1 {print $2, $6, $16*1.96}' ${infile} | \
psxy -R -JX -Ey${sizee}/1,${red} -K -O  >> $psfile
# PLOT CIRCLES
awk -F',' 'NR > 1 {print $2, $6}' ${infile} | \
psxy -R -JX -Sc$sizes -G${red} -O -K  >> $psfile

echo 0 0 BL 0.0 13p,Helvetica,${orange} XX > tmp
echo 3 0 BL 0.0 13p,Helvetica,${red}    XY >> tmp
echo 6 0 BL 0.0 13p,Helvetica,${blue}   YX >> tmp
echo 9 0 BL 0.0 13p,Helvetica,${green}  YY >> tmp
pstext tmp -R0/12/0/2 -JX4/2 -F+j+a+f -N -O -X2.5 -Y-2 >> $psfile
