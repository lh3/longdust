set t po eps co enh "Helvetica,20"

set size 1,0.9

set style histogram rowstacked
set xtics rotate by 40 right nomirror font "Helvetica,20"
set boxwidth 0.7 relative
set style data histograms
set style fill solid 1.0 border lt -1
set border 3

set key invert

set out "len.eps"
set yran [0:460]
set key top left
set rmargin 2
set mytics 5
set ylab "Length (Mb)" off +0.0,0
plot "<grep -v TANTAN-85 eval.txt" u ($2*1e-6):xtic(1) t "CenSat" ls 1, \
	"" u ($3*1e-6) t "longdust" ls 2, \
	"" u ($4*1e-6) t "TRF" ls 3, \
	"" u ($5*1e-6) t "SDUST" ls 4, \
	"" u ($6*1e-6) t "remainder" ls 5

reset

set out "f-func.eps"
set xlab "Sequence length (L)"
set xran [0:5000]
set yran [0:700]
set mxtics 5
set mytics 5
set key bot right
g(x) = 0.6 * x
plot "f-func.txt" u 1:2 t "f(L/4^7) at 50% GC" w l sm cs lw 3, \
	"" u 1:3 t "f(L/4^7) at 40% GC" w l sm cs lw 3, \
	g(x) t "0.6L" w l lw 3
