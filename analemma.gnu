set size ratio 1 
#set terminal png size 1000,900  font ",20" 
set terminal png
#set autoscale fix 
#set xrange [-1.2:1.2]
#set yrange [-1.2:1.2]
set xzeroaxis 
set yzeroaxis 
set output 'analemmar.png' 
plot "analemma.dat" u ($2/1):($1/1) w l  notitle
pause -1
