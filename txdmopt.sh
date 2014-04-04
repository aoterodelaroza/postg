#!/bin/bash

## modify this
chf="lcwpbe"
c1=0.9453
c2=1.1217
export TeraChem=/home/cisborn/manchego_home/working_tc/production/terachem
export LD_LIBRARY_PATH=/usr/local/cuda-5.0/lib64:$LD_LIBRARY_PATH
TERACHEM="/home/cisborn/manchego_home/working_tc/production/terachem/int"
POSTG="/home/cisborn/manchego_home/xdm/postg/postg"
verbose=""
cat > $2.route <<EOF
# charge, spin, run, coordinates set by the script
basis        6-31g* 
method       wpbe
rc_w         0.4
maxit        100
gpus         4  
EOF
#########

function toxyz {
    declare -a aa=("Bq" "H" "He" "Li" "Be" "B" "C" "N" "O" "F"\
                   "Ne" "Na" "Mg" "Al" "Si" "P"  "S"  "Cl" "Ar"\
                   "K"  "Ca" "Sc" "Ti" "V"  "Cr" "Mn" "Fe" "Co"\
                   "Ni" "Cu" "Zn" "Ga" "Ge" "As" "Se" "Br" "Kr"\
                   "Rb" "Sr" "Y"  "Zr" "Nb" "Mo" "Tc" "Ru" "Rh"\
                   "Pd" "Ag" "Cd" "In" "Sn" "Sb" "Te" "I"  "Xe"\
                   "Cs" "Ba" "La" "Ce" "Pr" "Nd" "Pm" "Sm" "Eu"\
                   "Gd" "Tb" "Dy" "Ho" "Er" "Tm" "Yb" "Lu" "Hf"\
                   "Ta" "W"  "Re" "Os" "Ir" "Pt" "Au" "Hg" "Tl"\
                   "Pb" "Bi" "Po" "At" "Rn" "Fr" "Ra" "Ac" "Th"\
                   "Pa" "U"  "Np" "Pu" "Am" "Cm" "Bk" "Cf" "Es"\
                   "Fm" "Md" "No" "Lr")
    while read line; do
	f1=$(awk '{print $1}' <<< $line ) 
	f2=$(awk '{printf "%.12f %.12f %.12f",$2*.52917720859,$3*.52917720859,$4*.52917720859}' <<< $line ) 
	echo ${aa[$f1]} $f2
    done
}

# read info from gaussian "output"
read atoms derivs charge spin < $2

# route string
run="energy"
if [ $derivs == "2" ] ; then
    do1="1"
    do2="1"
    echo "No second derivatives!!!!" > $4
    exit
elif [ $derivs == "1" ] ; then
    do1="1"
    run="gradient"
fi

# prepare the geometry
cat > $2.xyz <<EOF
$atoms

$(sed -n 2,$(($atoms+1))p < $2 | cut -c 1-72 | toxyz)
EOF

# prepare the terachem input file ($2.inp)
cat $2.route > $2.inp
cat >> $2.inp <<EOF
charge  $charge
spinmult $spin
run $run
coordinates $2.xyz
end
EOF

# run terachem
$TERACHEM $2.inp >& $2.inp.out

if [ -n "$do1" ] ; then
    grep -A $atoms 'Gradients' postg.tcfchk | tail -n $atoms > $2.fgauss
fi

# run postG
$POSTG $c1 $c2 postg.tcfchk $chf > $2.pgout

# energy
e=$(grep 'total energy' $2.pgout | awk '{print $NF}')
printf "%20.12e%20.12e%20.12e%20.12e\n" $e 0 0 0 | tr e D > $3

# forces
if [ -n "$do1" ] ; then
    grep -A $(($atoms+1)) 'dispersion forces' $2.pgout | awk '{print $2, $3, $4}' | tail -n $atoms > $2.fpostg
    paste $2.fgauss $2.fpostg > $2.forces
    awk '{printf("%20.12e%20.12e%20.12e\n",$1-$4,$2-$5,$3-$6)}' $2.forces | tr e D >> $3
fi

if [ -n "$verbose" ] ; then
    # output for debug
    echo "#XDM# terachem input" >> $4
    cat $2.inp >> $4
    echo "#XDM# terachem output" >> $4
    cat $2.inp.out >> $4
    echo "#XDM# tcfchk file" >> $4
    cat postg.tcfchk >> $4
else
    echo "#XDM# tail -n 50 logfile" >> $4
    tail -n 10 $2.inp.out >> $4
    echo "#XDM# energies from postg output (hartree)" >> $4
    grep 'energy' $2.pgout >> $4
fi

rm -f $2.* >& /dev/null
