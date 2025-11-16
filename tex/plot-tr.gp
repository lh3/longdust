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
plot "<grep -v FasTAN eval.txt" u ($2*1e-6):xtic(1) t "CenSat" ls 1, \
	"" u ($3*1e-6) t "longdust" ls 2, \
	"" u ($4*1e-6) t "TRF" ls 3, \
	"" u ($5*1e-6) t "SDUST" ls 4, \
	"" u ($6*1e-6) t "remainder" ls 5

set out "rev.eps"
set yran [0:*]
set ylab "Percent difference between strands"
plot "<grep -v FasTAN rev.txt" u ((1-$2/$3)*100):xtic(1) not ls 1, \
	"" u 0:((1-$2/$3)*100):(sprintf("%.2f",(1-$2/$3)*100)) not w labels offset 0,.7 font ',16'

reset

set size 0.75,1

set out "f-func.eps"
set xlab "Sequence length (L)"
set xran [0:5000]
set yran [0:1300]
set mxtics 5
set mytics 4
set key top left
set key offset -3,0
set ylab "Function value" offset 2,0
set lmargin 7

g(x) = 0.6 * x
plot "f-func.txt" u 1:2 t "f(L) at 50% GC" w l sm cs lw 3, \
	"" u 1:3 t "f(L) at 40% GC" w l sm cs lw 3, \
	"" u 1:4 t "f(L) at 30% GC" w l sm cs lw 3, \
	g(x) t "0.6L" w l lw 3

set out "f-func-diff.eps"
set key offset -3.5,0
set ylab "First derivative" offset 1,0
set yran [0:0.6]
set mytics 5
plot "f-func-diff.txt" u 1:2 t "f'(L) at 50% GC" w l sm cs lw 3, \
	"" u 1:3 t "f'(L) at 40% GC" w l sm cs lw 3, \
	"" u 1:4 t "f'(L) at 30% GC" w l sm cs lw 3
