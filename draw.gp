if (exist("i")==0 || i<0) i=0

if (i%10==0) print i
set xrange [-0.5e0:2.7e0]
set yrange [-0.5e0:1.0e0]
set zrange [ 0.0e0:3.0e0]
splot "fort.999" u 1:2:3 index i ti "" pointtype 6
pause 0.01
i = i+1
if (i<=200) reread

i = -1
