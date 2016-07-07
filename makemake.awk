#! /usr/bin/awk -f
# Copyright (c) 2013-2016 Alberto Otero de la Roza
# <aoterodelaroza@gmail.com>, Felix Kannemann
# <felix.kannemann@dal.ca>, Erin R. Johnson <erin.johnson@dal.ca>,
# Ross M. Dickson <ross.dickson@dal.ca>, Hartmut Schmider
# <hs7@post.queensu.ca>, and Axel D. Becke <axel.becke@dal.ca>
#
# postg is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# makemake.awk -- generate object dependencies for the Makefile.
# Assumes files and modules share name and good old common sense.
#
# Use: makemake.awk *.f90

{ 
    if (FILENAME != f[fs]){
        fs++
        f[fs] = FILENAME
        gsub(/\.f90$/,"",f[fs])
    }
}
/^( |\t)*module( |\t)*[^ \t\n]*( |\t)*$/{
    ismodule[fs] = 1
}
/^( |\t)*use( |\t)*[^ \t\n]*/{
    nm = tolower($2)
    idx = index(nm,",")
    if (idx != 0)
       nm = substr(nm,0,idx-1)
    for (i=1;i<=uses[nm];i++){
       if (use[nm,i] == f[fs])
          next
    }
    uses[nm]++
    use[nm,uses[nm]] = f[fs]
}

END{
    for (i=1;i<=fs;i++){
        if (ismodule[i] && uses[tolower(f[i])]){
            str = sprintf(": %s.mod",f[i])
            for (j=1;j<=uses[tolower(f[i])];j++)
                str = sprintf("%s.o %s",use[f[i],j],str)
            print str 
        }
    }
}
