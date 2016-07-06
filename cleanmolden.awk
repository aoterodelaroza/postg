#! /usr/bin/awk -f

FNR==1{getline line < FILENAME}
{getline line < FILENAME}

/^ *Sym=.*$/{
    if (inocc)
	flushblock()
    inocc = 0
    sym = $0
    next
}
/^ *Ene=.*$/{
    if (inocc)
	flushblock()
    inocc = 0
    ene = $0
    next
}
/^ *Spin=.*$/{
    if (inocc)
	flushblock()
    inocc = 0
    spin = $0
    next
}
/^ *Occup=.*$/{
    if (inocc)
	flushblock()
    inocc = 1
    occl = $0
    occ = $NF+0
    next
}
(inocc == 1) && ($1+0 == $1){
    coef[$1+0] = $2
    ncoef = $1+0
    next
}
inocc == 1{
    inocc = 0
}
{print}

function flushblock(){
    if (occ > 1e-5){
	print sym
	print ene
	print spin
	print occl
	for (i=1;i<=ncoef;i++)
	    print i, coef[i]
    }
}
