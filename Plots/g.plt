set term pdf enhanced color dashed size 4,4 font "Helvetica,20"

reset
set output "g.pdf"

set xrange [-2.25:2.25]
set xlabel "g" offset 0,0.5
set ylabel "GS energy" offset 1
plot "g.dat" u 1:2 title "PT" w l, "" u 1:3 title "diag" w l
