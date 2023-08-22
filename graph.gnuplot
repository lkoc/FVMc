#!/usr/local/bin/gnuplot -persist
# set terminal pngcairo  transparent enhanced font "arial,10" fontscale 1.0 size 600, 400 
# set output 'heatmaps.7.png'
set format cb "%4.1f" 
unset parametric
set view map scale 1
set samples 25, 25
set isosamples 50, 50
set contour base
set xyplane relative 0
set cbtics border in scale 0,0 mirror norotate  autojustify
set title "(2D Heat Map of MV Cables)\nZ is contoured. Independent value is color-mapped" 
set title  offset character 0, 1, 0 font "" textcolor lt -1 norotate
set urange [ 20.00000 : 60.0000 ] noreverse nowriteback
set vrange [ 20.00000 : 60.0000 ] noreverse nowriteback
set xlabel "x" 
set xlabel  offset character 3, 0, 0 font "" textcolor lt -1 norotate
set xrange [ * : * ] noreverse writeback
set x2range [ * : * ] noreverse writeback
set ylabel "y" 
set ylabel  offset character -1, 0, 0 font "" textcolor lt -1 norotate
set yrange [ * : * ] noreverse writeback
set y2range [ * : * ] noreverse writeback
set zlabel "z" 
set zlabel  offset character 2, 0, 0 font "" textcolor lt -1 norotate
set zrange [ * : * ] noreverse writeback
set cbrange [ * : * ] noreverse writeback
set rrange [ * : * ] noreverse writeback
set pm3d implicit at s
set colorbox user
set colorbox vertical origin screen 0.9, 0.2 size screen 0.03, 0.6 front  noinvert noborder

#read data of three columns x,y, temperature of external file

splot "results_graph_ex3.txt" using 1:2:3 with pm3d title "XXX"
# asign column 1: to x axis, column 2: to y axis, column 3: to z axis and column 4: to color

color(x,y) = 10. * (1.1 + sin((x-20.)/5.)*cos((y-20.)/10.))
NO_ANIMATION = 1
## Last datafile plotted: "++"
splot '++' using 1:2:3:(color($1,$2)) with lines nosurface  title "4 data columns x/y/z/color"