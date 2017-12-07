set size ratio 1
set terminal png size 1000,800 font '12'
set xzeroaxis
set yzeroaxis
set terminal png
set output 'analemma.png
plot analemma.dat u ($1/1):($21) w l  notitle
set output 'analemma.png
plot analemma.dat u ($1/1):($21) w l  notitle
pause -1