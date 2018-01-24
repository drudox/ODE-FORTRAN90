#set terminal pdfcairo enhanced dashed font "Helvetica, 14 " size 12cm,12cm
#set output "WF_EWT.pdf"

# 
set size 1.25,1

set title "Confronto metodi risolutivi ODE " font "Helvetica, 20" 

#{/Simbol e} per epsilon --> greca

set xlabel "$x$" font "Helvetica, 20" #asse x
#set xlabel "time [$s$]" font "Helvetica, 20" #asse x
set ylabel "$y(x)$" font "Helvetica, 20" #asse y
#set ylabel "$y(t)$" font "Helvetica, 20" #asse y
 

set style line 1 	dashtype 1  linecolor rgb "red"  	lw 6 #ps 1.5
set style line 10	dashtype 2 linecolor rgb "red"	  	lw 4 #ps 1.5

set style line 2 	dashtype 1 linecolor rgb "blue"  	lw 5 #ps 1.5
set style line 3 	dashtype 1 linecolor rgb "green" 	lw 3 

set style line 4 	dashtype 1 linecolor rgb "blue" 		lw 6 #ps 1.5
set style line 40	dashtype 3 linecolor rgb "blue" 		lw 6 #ps 1.5

set style line 5 	dashtype 1 linecolor rgb "magenta"  lw 5
set style line 6 	dashtype 1 linecolor rgb "orange" 	lw 5
set style line 7 	dashtype 1 linecolor rgb "purple" 	lw 6 #ps 1.5
set style line 70 dashtype 2 linecolor rgb "purple" 	lw 5 #ps 1.5
set style line 8 	dashtype 1 linecolor rgb "orange" 	lw 6

set style line 9 	dashtype 1 linecolor rgb "black"    lw 4
set style line 90	dashtype 2 linecolor rgb "black"    lw 4

set style line 11 dashtype 1 linecolor rgb "blue"  lw 8 
set style line 12 dashtype 1 pointtype 12 linecolor rgb "green"  lw 1 ps 2
set style line 13 dashtype 1  pointtype 6 linecolor rgb "black"  lw 6 ps 1.5
set style line 14	dashtype 2 linecolor rgb "blue" lw 4

#set label 1 at 0.3, -1
#set label 1 "Suction side" 

#set label 2 at 0.3, 0.35
#set label 2 "Pressure side"




set grid
#set yrange [0:1.2] #gli assi sono già invertiti da qui
#set xtic 2e-3
set mxtics 5
set xrange [-5:5]
#set ytic 0.1
set mytics 5

#set key top left
set key top right
#set key bottom right

plot \
	  'bwdEuler004.out' u 1:2 w l ls 4 t "BWD Euler" ,\
	  'fwdEuler004.out' u 1:2 w l ls 40 t "FWD Euler" ,\
	  'realSol_004.out' u 1:2 w p ls 13 t "Soluzione Analitica" ,\
	  'AB2_004.out'     u 1:2 w l ls 7 t "AB 2th",\
	  'AB3_004.out'     u 1:2 w l ls 5 t "AB 3th",\
        'rk4_004.out'     u 1:2 w l ls 3 t "RK 4th"
# 'fe001.out' u 1:2 w l ls 4  t "Foward Euler" ,\
# 'be001.out' u 1:2 w l ls 1  t "Backward Euler"
#        'AM3_003.out' u 1:2 w l ls 1  t "AM 3nd (P-C)",\

