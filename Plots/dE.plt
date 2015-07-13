set term pdf enhanced color dashed size 4,4 font "Helvetica,20"

reset
set output "dE.pdf"

#set xrange [-2.25:2.25]
set xlabel "g" offset 0,0.5
set ylabel "ΔE" offset 1
set key bottom center
plot "dE.dat" u 1:2 title "PT" w l, "" u 1:3 title "diag" w l,"" u 1:4 title "CCD" w l,"" u 1:5 title "GF" w l
