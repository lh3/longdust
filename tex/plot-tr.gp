set t po eps co "Helvetica,20"

set size 1,0.9

set style histogram rowstacked
set xtics rotate by 40 right nomirror font "Helvetica,20"
set boxwidth 0.7 relative
set style data histograms
set style fill solid 1.0 border lt -1
set border 3

set key invert

set out "len.eps"
set yran [0:450]
set key top left
set mytics 5
set ylab "Length (Mb)" off +0.0,0
plot "<grep -v TANTAN-85 eval.txt" u ($2*1e-6):xtic(1) t "CenSat" ls 1, \
	"" u ($3*1e-6) t "longdust" ls 2, \
	"" u ($4*1e-6) t "TRF" ls 3, \
	"" u ($5*1e-6) t "SDUST" ls 4, \
	"" u ($6*1e-6) t "remainder" ls 5
