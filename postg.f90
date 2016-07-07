!                          _        
!          _ __   ___  ___| |_ __ _ 
!         | '_ \ / _ \/ __| __/ _` |
!         | |_) | (_) \__ \ || (_| |
!         | .__/ \___/|___/\__\__, |
!         |_|                 |___/ 
!    
! postg is a program that calculates the exchange-hole dipole moment
! (XDM) model dispersion coefficients and energy using gaussian-style
! wfn and wfx files.
!
! Copyright (c) 2013-2016 Alberto Otero de la Roza
! <aoterodelaroza@gmail.com>, Felix Kannemann
! <felix.kannemann@dal.ca>, Erin R. Johnson <erin.johnson@dal.ca>,
! Ross M. Dickson <ross.dickson@dal.ca>, Hartmut Schmider
! <hs7@post.queensu.ca>, and Axel D. Becke <axel.becke@dal.ca>
!
! postg is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
program postg
  use meshmod
  use wfnmod
  use param
  use atomicdata

  implicit none

  type(molecule) :: mol
  type(tmesh) :: mesh
  character*(mline) :: line, wfnfil, hfword
  logical :: ok, usec9
  integer :: i, narg, idx, lp
  real*8 :: qpro, c1br, c2br, egauss, ntotal
  character*11 :: wfnext

  ! init
  chf = 0d0
  call atomic_init()

  ! read command line
  narg = command_argument_count()
  if (narg < 3) goto 999
  call getarg(1,line)
  read (line,*,err=999) c1br
  call getarg(2,line)
  read (line,*,err=999) c2br
  call getarg(3,line)
  wfnfil = adjustl(trim(line))
  inquire(file=wfnfil,exist=ok)
  if (.not.ok) then
     write (iout,'("wfn file not found: ",A)') trim(wfnfil)
     stop 1
  endif

  ! header and convert a2
  WRITE(IOUT,'("* POSTG OUTPUT")')
  WRITE(IOUT,'("a1          ",F12.6)') c1br
  WRITE(IOUT,'("a2(ang)     ",F12.6)') c2br
  c2br=c2br/0.52917720859d0

  if (narg >= 4) then
     call getarg(4,line)
     read (line,*,err=999) hfword
     lp = 1
     if (.not.isreal(chf,hfword,lp)) then
        if (trim(lower(hfword)) == "blyp") then
           chf = chf_blyp
           write(iout,'("a_hf        blyp")')
        elseif (trim(lower(hfword)) == "b3lyp") then
           chf = chf_b3lyp
           write(iout,'("a_hf        b3lyp")')
        elseif (trim(lower(hfword)) == "bhandhlyp" .or. trim(lower(hfword)) == "bhandh"&
           .or. trim(lower(hfword)) == "bhah" .or. trim(lower(hfword)) == "bhahlyp") then
           chf = chf_bhahlyp
           write(iout,'("a_hf        bhandhlyp")')
        elseif (trim(lower(hfword)) == "camb3lyp" .or. trim(lower(hfword)) == "cam-b3lyp" ) then
           chf = chf_camb3lyp
           write(iout,'("a_hf        camb3lyp")')
        elseif (trim(lower(hfword)) == "pbe") then
           chf = chf_pbe
           write(iout,'("a_hf        pbe")')
        elseif (trim(lower(hfword)) == "pbe0") then
           chf = chf_pbe0
           write(iout,'("a_hf        pbe0")')
        elseif (trim(lower(hfword)) == "lcwpbe" .or. trim(lower(hfword)) == "lc-wpbe") then
           chf = chf_lcwpbe
           write(iout,'("a_hf        lcwpbe")')
        elseif (trim(lower(hfword)) == "pw86" .or. trim(lower(hfword)) == "pw86pbe") then
           chf = chf_pw86
           write(iout,'("a_hf        pw86pbe")')
        elseif (trim(lower(hfword)) == "b971" .or. trim(lower(hfword)) == "b97-1") then
           chf = chf_b971
           write(iout,'("a_hf        b971")')
        elseif (trim(lower(hfword)) == "hf") then
           chf = 1.0d0
           write(iout,'("a_hf        hf")')
        else
           call error("postg","unknown functional",2)
        endif
     else
        write(iout,'("a_hf        ",f12.6)') chf
     endif
  else
     write(iout,'("a_hf        ",f12.6)') chf
  endif

  ! optional keywords
  usec9 = .false.
  if (narg >= 5) then
     do i = 5, narg
        call getarg(i,line)
        if (trim(lower(line)) == 'c9') then
           usec9 = .true.
        else
           call error("postg","unknown keyword " // trim(line),2)
        endif
     end do
  end if

  ! read wfn and output some info
  idx = index(wfnfil,'.',.true.)
  if (wfnfil(idx+1:idx+3) == "wfx") then
     wfnext = "wfx"
     mol = readwfx(wfnfil,egauss)
  else if (wfnfil(idx+1:idx+3) == "wfn") then
     wfnext = "wfn"
     mol = readwfn(wfnfil,egauss)
  else if (wfnfil(idx+1:idx+4) == "fchk") then
     wfnext = "fck"
     mol = readfchk(wfnfil,egauss)
  else if (wfnfil(idx+1:idx+6) == "tcfchk") then
     wfnext = "tck"
     mol = readtck(wfnfil,egauss)
  else if (wfnfil(idx+1:idx+6) == "molden") then
     wfnext = "molden"
     mol = readmolden(wfnfil,egauss)
  else
     wfnext = "assumed wfn"
     mol = readwfn(wfnfil,egauss)
  endif
  write (iout,'("file          ",A)') adjustl(trim(mol%name))
  write (iout,'("filetype      ",A)') adjustl(trim(wfnext))
  write (iout,'("natoms        ",I6)') mol%n
  write (iout,'("use ECPs?     ",L5)') mol%useecp
  if (any(mol%z(:) < 0)) call error('postg','unrecognized atomic symbol',2)
  write (iout,'("# n  At           x               y               z         nr    nl")')
  do i = 1, mol%n
     write (iout,'(I4,X,A2,X,3(F15.7,X),I6,1X,I6)') &
        i, ptable(mol%z(i)), mol%x(:,i), z2nr(mol%z(i)), z2nang(mol%z(i))
  enddo
  write (iout,'("#")')
  if (mol%wfntyp == 0) then
     write (iout,'("wfn type      ","closed-shell")') 
  elseif (mol%wfntyp == 1) then
     write (iout,'("wfn type      ","open-shell")') 
  elseif (mol%wfntyp == 2) then
     write (iout,'("wfn type      ","restricted open-shell")') 
  elseif (mol%wfntyp == 3) then
     write (iout,'("wfn type      ","natural orbitals")') 
  endif
  write (iout,'("MOs           ",I6)') mol%nmo
  write (iout,'("primitives    ",I6)') mol%npri
  write (iout,'("charge        ",F6.2)') mol%charge
  write (iout,'("multiplicity  ",I6)') mol%mult

  ! generate the mesh
  mesh = genmesh(mol)
  write (iout,'("mesh size     ",I10)') mesh%n

  ! open scratch files and generate the promolecular density
  call atomin(mol,mesh,qpro)
  write (iout,'("nelec         ",F12.6)') mol%nelec
  write (iout,'("nelec (promol)",F12.6)') qpro

  ! evaluate the density, etc. on the grid
  call evalwfn(mol,mesh)
  write (iout,'("nelec, alpha  ",F12.6)') sum(mesh%w * mesh%rho(:,1))
  write (iout,'("nelec, beta   ",F12.6)') sum(mesh%w * mesh%rho(:,2))
  ntotal = sum(mesh%w * (mesh%rho(:,1)+mesh%rho(:,2)))
  write (iout,'("nelec, total  ",F12.6)') ntotal
  if (abs(ntotal - mol%nelec) > 0.1d0) then
     write (iout,'("WARNING -- inconsistent nelec. I hope you know what you are doing.")')
  endif

  write (iout,'("moments and volumes ")')
  write (iout,'("# i At        <M1^2>             <M2^2>              <M3^2>           Volume              Vfree")')
  do i = 1, mol%n
     write (iout,'(I3,X,A2,X,1p,5(E18.10,X))') i, ptable(mol%z(i)), mol%mm(1,i), &
        mol%mm(2,i), mol%mm(3,i), mol%v(i), frevol(mol%z(i))
  enddo
  write (iout,'("#")')

  write (iout,'("hirshfeld charges ")')
  write (iout,'("# i At        Charge")')
  do i = 1, mol%n
     write (iout,'(I3,X,A2,X,1p,1(E18.10,X))') i, ptable(mol%z(i)), mol%z(i)-mol%q(i)
  enddo
  write (iout,'("#")')

  call edisp(mol,c1br,c2br,egauss,usec9)

  stop 

999 continue
  write (iout,'(/"Usage: postg a1 a2(angstrom) wfnfile [a_hf] [optional]")')
  write (iout,*)
  write (iout,'("  a1,a2 = damping function coefficients. See: http://schooner.chem.dal.ca/wiki/XDM ")')
  write (iout,*)
  write (iout,'("  wfnfile = Gaussian wfn/wfx file ")')
  write (iout,*)
  write (iout,'("  a_hf = functional keyword. One of: ")')
  write (iout,'("         blyp,b3lyp,bhah,bhahlyp,bhandh,bhahlyp,bhandhlyp,camb3lyp,cam-b3lyp, ")')
  write (iout,'("         pbe,pbe0,lcwpbe,lc-wpbe,pw86,pw86pbe,b971,b97-1, or a number between 0.0 ")')
  write (iout,'("         and 1.0 representing the fraction of exact exchange. Default: 0.0.")')
  write (iout,*)
  write (iout,'("  optional = an optional keyword. One of: ")')
  write (iout,'("    c9: calculate the C9 dispersion coefficients (no contribution to the energy).")')
  stop 1

END program postg

