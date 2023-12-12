set term png
set output 'graph.png'
set multiplot
plot 'data.txt' with lines
plot 'oldData.txt' with lines
unset multiplot