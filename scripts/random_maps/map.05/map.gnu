set term png crop  font 'Helvetica, 11'
set out 'map.png'
unset key
unset colorbox

set colors classic
set palette maxcolors 2
set view 0,0

set xlabel 'L_x (nm)' #offset 5,0
set ylabel 'L_y (nm)' offset 2,0

#splot [:] [:] [:] 'map.otp' u 1:2:3 with image
plot [:] [:] [:] 'map.otp' u 1:2:3 with image

set term x11

