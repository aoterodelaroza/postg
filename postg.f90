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
! Copyright (c) 2013-2025 Alberto Otero de la Roza
! <aoterodelaroza@gmail.com>, Felix Kannemann
! <felix.kannemann@dal.ca>, Erin R. Johnson <erin.johnson@dal.ca>,
! Ross M. Dickson <ross.dickson@dal.ca>, Hartmut Schmider
! <hs7@post.queensu.ca>, Kyle R. Bryenton <kyle.bryenton@gmail.com>,
! and Axel D. Becke <axel.becke@dal.ca>
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
  character*(mline) :: wfnfil, hfword
  logical :: ok, usec9
  integer :: i, narg, narg_i, idx, lp
  real*8 :: qpro, egauss, ntotal
  character*11 :: wfnext

  integer, parameter :: buffer_arr_siz = 10
  character(len=100) :: buffer_arr(buffer_arr_siz)
  real*8 :: xdm_temp, c1br, c2br, z_damp
  integer :: i_code, xdm_damping

  ! init
  chf = 0d0
  call atomic_init()

  ! Initialize Input Buffer
  narg = command_argument_count()
  narg_i = 1
  if (narg < 3) goto 999
  do narg_i = 1, min(narg, 10)
      call getarg(narg_i, buffer_arr(narg_i))
  end do
  narg_i = 1

  ! The first 1-3 lines are BJ/Z Damping Parameters
  ! _Interpretted as_    _As entered_
  ! BJ Damping, a1, a2 = BJ float float
  !                    = float float
  ! Z Damping, z_damp  = Z float
  !                    = float
  read(buffer_arr(narg_i),'(G100.100)', iostat=i_code) xdm_temp
  if (i_code==0) then
      ! Increment and check for a2
      read(buffer_arr(narg_i+1),'(G100.100)', iostat=i_code) xdm_temp
      if (i_code==0) then
          ! Both a1 and a2 were found, must be BJ damping
          xdm_damping = 1
          read(buffer_arr(narg_i),*,iostat=i_code,err=999) c1br
          read(buffer_arr(narg_i+1),*,iostat=i_code,err=999) c2br
          narg_i = narg_i + 2
      else
          ! Only a1 was found, must be Z damping
          xdm_damping = 2
          read(buffer_arr(narg_i),*,iostat=i_code,err=999) z_damp
          narg_i = narg_i + 1
      end if
  else
      ! It's a string, so it should be either "BJ" or "Z"
      select case (trim(lower(buffer_arr(narg_i))))
      case('bj', 'bj_damp') ! Currently the default, so technically redundant
          xdm_damping = 1 ;
          read(buffer_arr(narg_i+1),*,iostat=i_code,err=999) c1br
          read(buffer_arr(narg_i+2),*,iostat=i_code,err=999) c2br
          narg_i = narg_i + 3
      case('z', 'z_damp')
          xdm_damping = 2 ;
          read(buffer_arr(narg_i+1),*,iostat=i_code,err=999) z_damp
          narg_i = narg_i + 2
      case default
          ! Input unrecognized, crash
          write (iout,'("Input Unrecognized: ",A)') buffer_arr(narg_i)
          goto 999
      end select
  end if

  ! Write initial header information
  write(iout,'("* POSTG OUTPUT")')
  if (xdm_damping == 1) then
      write(iout,'("Damping Type  BJ-Damping")')
      write(iout,'("a1        ",F12.6)') c1br
      write(iout,'("a2(ang)   ",F12.6)') c2br
      c2br=c2br/0.52917720859d0
  else if (xdm_damping == 2) then
      write(iout,'("Damping Type  Z-Damping")')
      write(iout,'("z_damp  ",I12)') nint(z_damp)
  end if

  ! The next input is the wavefunction file
  wfnfil = trim(buffer_arr(narg_i))
  inquire(file=wfnfil,exist=ok)
  if (.not.ok) then
      write (iout,'("wfn file not found: ",A)') wfnfil
      stop 1
  end if
  narg_i = narg_i + 1

  ! The next input is the functional (if unknown, read HF mixing fraction from file)
  if (narg_i <= narg) then
      hfword = trim(lower(buffer_arr(narg_i)))  ! KRB 2024-01-28 Made it so it trimmed and lowered only on assignment
      lp = 1
      if (.not.isreal(chf,hfword,lp)) then
          if (hfword == "blyp") then
              chf = chf_blyp
              write(iout,'("a_hf          blyp")')
          elseif (hfword == "b3lyp") then
              chf = chf_b3lyp
              write(iout,'("a_hf          b3lyp")')
          elseif (hfword == "bhandhlyp" .or. hfword == "bhandh"&
              .or. hfword == "bhah" .or. hfword == "bhahlyp") then
              chf = chf_bhahlyp
              write(iout,'("a_hf          bhandhlyp")')
          elseif (hfword == "camb3lyp" .or. hfword == "cam-b3lyp" ) then
              chf = chf_camb3lyp
              write(iout,'("a_hf          camb3lyp")')
          elseif (hfword == "pbe") then
              chf = chf_pbe
              write(iout,'("a_hf          pbe")')
          elseif (hfword == "pbe0") then
              chf = chf_pbe0
              write(iout,'("a_hf          pbe0")')
          elseif (hfword == "lcwpbe" .or. hfword == "lc-wpbe") then
              chf = chf_lcwpbe
              write(iout,'("a_hf          lcwpbe")')
          elseif (hfword == "pw86" .or. hfword == "pw86pbe") then
              chf = chf_pw86
              write(iout,'("a_hf          pw86pbe")')
          elseif (hfword == "b971" .or. hfword == "b97-1") then
              chf = chf_b971
              write(iout,'("a_hf          b971")')
          elseif (hfword == "hf") then
              chf = 1.0d0
              write(iout,'("a_hf          hf")')
          else
              chf = chf_fromfile
              call frevol_read_file(hfword)
          endif
      else
          write(iout,'("a_hf          ",f12.6)') chf
      endif
  else
      write(iout,'("a_hf          ",f12.6)') chf
  end if
  narg_i = narg_i + 1

  ! If any more lines exist, they're optional keywords. 
  ! Currently, only "C9" exists at this point. 
  ! If more are encountered, call an error.
  usec9 = .false.
  do
      if (narg_i > narg) exit
      select case (trim(lower(buffer_arr(narg_i))))
      case('c9')
          usec9 = .true.
          write(iout,'("use_c9        ",L)') usec9
      case default
          call error("postg","unknown keyword " // buffer_arr(narg_i),2)
      end select
      narg_i = narg_i + 1
  end do


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

  write (iout,'("atomic polarizabilities (bohr^3)")')
  write (iout,'("# i At        alpha")')
  do i = 1, mol%n
     write (iout,'(I3,X,A2,X,1p,5(E18.10,X))') i, ptable(mol%z(i)), mol%v(i) * frepol(mol%z(i)) / frevol(mol%z(i))
  enddo
  write (iout,'("#")')

  write (iout,'("hirshfeld charges ")')
  write (iout,'("# i At        Charge")')
  do i = 1, mol%n
     write (iout,'(I3,X,A2,X,1p,1(E18.10,X))') i, ptable(mol%z(i)), mol%z(i)-mol%q(i)
  enddo
  write (iout,'("#")')

  if (xdm_damping == 1) then
      call edisp_bj(mol,c1br,c2br,egauss,usec9)
  else if (xdm_damping == 2) then
      call edisp_z(mol,z_damp,egauss,usec9)
  else
      call error("postg","Unknown xdm_damping value at edisp() call",2)
  end if
  stop 

999 continue
  write (iout,'(/"Usage: postg [damp] wfnfile [a_hf] [optional]")')
  write (iout,*)
  write (iout,'("  damp = Damping function coefficients. One of: ")')
  write (iout,'("         BJ Damping:   BJ a1 a2(ang)")')
  write (iout,'("                       a1 a2(ang)")')
  write (iout,'("         Z Damping:    Z z_damp")')
  write (iout,'("                       z_damp")')
  write (iout,'("         Note: a1, a2(ang), and z_damp are floating point values")')
  write (iout,'("         See: https://erin-r-johnson.github.io/software/ ")')
  write (iout,*)
  write (iout,'("  wfnfile = Gaussian wfn/wfx file ")')
  write (iout,*)
  write (iout,'("  a_hf = functional keyword. One of: ")')
  write (iout,'("         blyp,b3lyp,bhah,bhahlyp,bhandh,bhahlyp,bhandhlyp,camb3lyp,cam-b3lyp, ")')
  write (iout,'("         pbe,pbe0,lcwpbe,lc-wpbe,pw86,pw86pbe,b971,b97-1, or a number between 0.0 ")')
  write (iout,'("         and 1.0 representing the fraction of exact exchange. Default: 0.0.")')
  write (iout,'("         Alternatively, a_hf can also be a text file with rows Z vfree(Z), where Z")')
  write (iout,'("         is the atomic number and vfree(Z) is the free volume for that atom.")')
  write (iout,*)
  write (iout,'("  optional = an optional keyword. One of: ")')
  write (iout,'("    c9: calculate the C9 dispersion coefficients (no contribution to the energy).")')
!  write (iout,'("    xcdm: calculate same- and opposite-spin dynamical correlation contribution")')
!  write (iout,'("          to the exchange-hole dipole moment (contributes to the energy).")')
  stop 1

END program postg

