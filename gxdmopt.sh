#!/bin/bash
# Copyright (c) 2013-2025 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
# Kyle R. Bryenton <kyle.bryenton@gmail.com>, Felix Kannemann <felix.kannemann@dal.ca>,
# Erin R. Johnson <erin.johnson@dal.ca>, Ross M. Dickson <ross.dickson@dal.ca>,
# Hartmut Schmider <hs7@post.queensu.ca>, and Axel D. Becke <axel.becke@dal.ca>

## modify this
chf="blyp"
c1=0.5942
c2=1.4555
basis="6-31+g*"
ecp=""
G09="g09"
POSTG="postg"
verbose=""
guesstrick=""
cat > $2.route <<EOF
%mem=2GB
%nprocs=4
#p blyp int(grid=ultrafine)
EOF
#########

if [ -f $basis ] ; then
    basisfile=$basis
    basis="gen"
fi
if [ -n "$ecp" ] && [ -f "$ecp" ] ; then
    ecpfile=$ecp
    ecp="pseudo=read"
else
    ecp=""
fi

# read info from gaussian "output"
read atoms derivs charge spin < $2

# route string
sroute="units=au output=wfx"
if [ $derivs == "2" ] ; then
    sroute="${sroute} freq(noraman) punch=derivatives"
    do1="1"
    do2="1"
elif [ $derivs == "1" ] ; then
    sroute="${sroute} force punch=derivatives"
    do1="1"
fi

# prepare the G09 input file ($2.in) for the single-point "Force" calculation
rm -f $2.in >& /dev/null
if [ -n "$guesstrick" ] ; then
    echo "%chk=${guesstrick%.chk}.chk" >> $2.in
fi
if [ -n "$guesstrick" ] && [ -f "${guesstrick%.chk}.chk" ] ; then
    guess="guess=read"
else
    guess=""
fi
cat $2.route >> $2.in
cat >> $2.in <<EOF
# $sroute $basis $ecp $guess

title

$charge $spin
$(sed -n 2,$(($atoms+1))p < $2 | cut -c 1-72)
EOF

if [ "x$basis" == "xgen" ] ; then
    echo "" >> $2.in
    awk '/^ *$/{next}{print}' $basisfile >> $2.in
fi

if [ -n "$ecp" ] ; then
    echo "" >> $2.in
    awk '/^ *$/{next}{print}' $ecpfile >> $2.in
fi

cat >> $2.in <<EOF

$2.wfx

EOF

# run G09
$G09 < $2.in > $2.out
if [ -n "$do1" ] ; then
    head -n $atoms fort.7 | tr D e > $2.fgauss
fi
if [ -n "$do2" ] ; then
    tail -n+$((${atoms}+1)) fort.7 | tr D e > $2.qgauss
fi
rm -f fort.7 >& /dev/null

# run postG
$POSTG $c1 $c2 $2.wfx $chf > $2.outg
grep 'WARNING -- inconsistent nelec' $2.outg && exit

# energy
e=$(grep 'total energy' $2.outg | awk '{print $NF}')
printf "%20.12e%20.12e%20.12e%20.12e\n" $e 0 0 0 | tr e D > $3

# forces
if [ -n "$do1" ] ; then
    grep -A $(($atoms+1)) 'dispersion forces' $2.outg | awk '{print $2, $3, $4}' | tail -n $atoms > $2.fpostg
    paste $2.fgauss $2.fpostg > $2.forces
    awk '{printf("%20.12e%20.12e%20.12e\n",$1-$4,$2-$5,$3-$6)}' $2.forces | tr e D >> $3
fi

# frequencies
if [ -n "$do2" ] ; then
    printf "%20.12e%20.12e%20.12e\n" 0 0 0 | tr e D >> $3 # polarizability
    printf "%20.12e%20.12e%20.12e\n" 0 0 0 | tr e D >> $3 
    for ((i=1;i<=3*${atoms};i++)) ; do
	printf "%20.12e%20.12e%20.12e\n" 0 0 0 | tr e D >> $3 # dip ders
    done

    grep -A $((3*$atoms*(3*$atoms+1)/2+1)) 'dispersion force constant matrix' $2.outg | \
	tail -n $((3*$atoms*(3*$atoms+1)/2)) | \
	awk '{printf "%s ",$NF}NR%3==0{printf"\n"}' > $2.qpostg
    paste $2.qgauss $2.qpostg > $2.freqs
    awk '{printf("%20.12e%20.12e%20.12e\n",$1+$4,$2+$5,$3+$6)}' $2.freqs | tr e D >> $3
fi

if [ -n "$verbose" ] ; then
    # output for debug
    echo "#XDM# gaussian input" >> $4
    cat $2.in >> $4
    echo "#XDM# gaussian output" >> $4
    cat $2.out >> $4
    echo "#XDM# wfx file" >> $4
    cat $2.wfx >> $4
    echo "#XDM# postg output" >> $4
    cat $2.outg >> $4
else
    echo "#XDM# tail -n 50 logfile" >> $4
    tail -n 50 $2.out >> $4
    echo "#XDM# grep 'Delta-E' logfile" >> $4
    grep 'Delta-E' $2.out >> $4
    echo "#XDM# energies from postg output (hartree)" >> $4
    grep 'energy' $2.outg >> $4
fi

rm -f $2.* >& /dev/null
