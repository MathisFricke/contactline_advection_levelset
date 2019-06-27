set datafile separator ","
set term pdf
set output "results.pdf"
set key above
set xlabel "time / s"
set ylabel "Contact Angle / deg"

set style line 1 \
    linecolor rgb '#dd181f' \
    linetype 1 linewidth 1 \
    pointtype 2 pointsize 0.5

plot "results_fs3d.csv" using 1:2 title "Reference" with linespoints linestyle 1,\
     "contactAngle.csv" using 1:2 title "LevelSet" with points pointsize 0.5,\
     "contactAngle.csv" using 1:3 title "reference" with points pointsize 0.1,\
