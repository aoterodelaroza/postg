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
module wfnmod
  use param
  implicit none

  private
  public :: atomin, readwfn, readwfx, readfchk, readtck, readmolden, evalwfn, edisp

  ! double factorials minus one
  integer, parameter :: dfacm1(0:8) = (/1,1,1,2,3,8,15,48,105/) 

  ! double factorials
  integer, parameter :: dfac(0:8) = (/1,1,2,3,8,15,48,105,384/) 

  ! number of angular components for shell type
  !                                    gs fs ds ps  s  p  d   f   g
  integer, parameter :: nshlt(-4:4) = (/9, 7, 5, 3, 1, 3, 6, 10, 15/) 

  ! initial and final types for cartesian shells
  integer, parameter :: jshl0(0:4) = (/1, 2, 5,  11, 21/) ! s, p, d, f, g
  integer, parameter :: jshl1(0:4) = (/1, 4, 10, 20, 35/) ! s, p, d, f, g

  ! named constants for the sph to car matrices below
  real*8, parameter :: s3 = sqrt(3d0)
  real*8, parameter :: s3_4 = sqrt(3d0/4d0)
  real*8, parameter :: s3_8 = sqrt(3d0/8d0)
  real*8, parameter :: s5_4 = sqrt(5d0/4d0)
  real*8, parameter :: s5_8 = sqrt(5d0/8d0)
  real*8, parameter :: s5_16 = sqrt(5d0/16d0)
  real*8, parameter :: s5_28 = sqrt(5d0/28d0)
  real*8, parameter :: s3_40 = sqrt(3d0/40d0)
  real*8, parameter :: s6_5 = sqrt(6d0/5d0)
  real*8, parameter :: s6 = sqrt(6d0)
  real*8, parameter :: s9_7 = sqrt(9d0/7d0)
  real*8, parameter :: s9_8 = sqrt(9d0/8d0)
  real*8, parameter :: s9_20 = sqrt(9d0/20d0)
  real*8, parameter :: s9_56 = sqrt(9d0/56d0)
  real*8, parameter :: s10 = sqrt(10d0)
  real*8, parameter :: s10_7 = sqrt(10d0/7d0)
  real*8, parameter :: s10_8 = sqrt(10d0/8d0)
  real*8, parameter :: s15 = sqrt(15d0)
  real*8, parameter :: s15_4 = sqrt(15d0/4d0)
  real*8, parameter :: s27_16 = sqrt(27d0/16d0)
  real*8, parameter :: s27_28 = sqrt(27d0/28d0)
  real*8, parameter :: s27_35 = sqrt(27d0/35d0)
  real*8, parameter :: s27_560 = sqrt(27d0/560d0)
  real*8, parameter :: s35_4 = sqrt(35d0/4d0)
  real*8, parameter :: s35_8 = sqrt(35d0/8d0)
  real*8, parameter :: s35_64 = sqrt(35d0/64d0)
  real*8, parameter :: s45 = sqrt(45d0)
  real*8, parameter :: s45_4 = sqrt(45d0/4d0)
  real*8, parameter :: s45_8 = sqrt(45d0/8d0)
  real*8, parameter :: s45_56 = sqrt(45d0/56d0)
  real*8, parameter :: s315_8 = sqrt(315d0/8d0)
  real*8, parameter :: s315_16 = sqrt(315d0/16d0)
  real*8, parameter :: d32 = 3d0/2d0
  real*8, parameter :: d34 = 3d0/4d0
  real*8, parameter :: d38 = 3d0/8d0

  ! -- Real solid harmonics r^l * Slm as a function of Cartesian products. -- 

  ! see: Schlegel and Frisch, Int. J. Quantum Chem. 54 (1995) 83,
  !      "Transformation Between Cartesian and Pure Spherical Harmonic Gaussians"
  !      https://doi.org/10.1002/qua.560540202
  ! Use table 1 and note that:
  !   Rlm = sqrt(2) * Im[Y_l^|m|] if m < 0
  !         Y_l^0 if m =0
  !         sqrt(2) * Re[Y_l^|m|] if m > 0
  ! Can also be found in Helgaker, eqs. 9.1.9 to 9.1.12

  ! dsphcar: l = 2 
  !   spherical molden order: m = 0, 1, -1, 2, -2
  !   Cartesian molden order: xx, yy, zz, xy, xz, yz
  real*8 :: dsphcar(5,6) = reshape((/&
     !  0      1     -1      2     -2
     -0.5d0, 0.0d0, 0.0d0,  s3_4, 0.0d0,& ! xx
     -0.5d0, 0.0d0, 0.0d0, -s3_4, 0.0d0,& ! yy
      1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,& ! zz
      0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0,& ! xy
      0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0,& ! xz
      0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0 & ! yz
     /),shape(dsphcar))

  ! fsphcar: l = 3 
  !   spherical molden order: m = 0, 1, -1, 2, -2, 3, -3
  !   Cartesian molden order: xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
  real*8 :: fsphcar(7,10) = reshape((/&
     ! 0        1       -1         2       -2         3       -3 
     0.0d0,   -s3_8,   0.0d0,    0.0d0,   0.0d0,     s5_8,   0.0d0,& ! xxx 
     0.0d0,   0.0d0,   -s3_8,    0.0d0,   0.0d0,    0.0d0,   -s5_8,& ! yyy 
     1.0d0,   0.0d0,   0.0d0,    0.0d0,   0.0d0,    0.0d0,   0.0d0,& ! zzz 
     0.0d0,  -s3_40,   0.0d0,    0.0d0,   0.0d0,    -s9_8,   0.0d0,& ! xyy 
     0.0d0,   0.0d0,  -s3_40,    0.0d0,   0.0d0,    0.0d0,    s9_8,& ! xxy 
    -s9_20,   0.0d0,   0.0d0,     s3_4,   0.0d0,    0.0d0,   0.0d0,& ! xxz 
     0.0d0,    s6_5,   0.0d0,    0.0d0,   0.0d0,    0.0d0,   0.0d0,& ! xzz 
     0.0d0,   0.0d0,    s6_5,    0.0d0,   0.0d0,    0.0d0,   0.0d0,& ! yzz 
    -s9_20,   0.0d0,   0.0d0,    -s3_4,   0.0d0,    0.0d0,   0.0d0,& ! yyz 
     0.0d0,   0.0d0,   0.0d0,    0.0d0,     1d0,    0.0d0,   0.0d0 & ! xyz 
     /),shape(fsphcar))

  ! gsphcar: l = 4
  !   spherical molden order: m = 0, 1, -1, 2, -2, 3, -3, 4, -4
  !   Cartesian molden order: xxxx yyyy zzzz xxxy xxxz xyyy yyyz xzzz yzzz xxyy xxzz yyzz xxyz xyyz xyzz
  real*8 :: gsphcar(9,15) = reshape((/&
  !    0       1      -1       2      -2        3      -3         4      -4
       d38,  0.0d0,  0.0d0, -s5_16,  0.0d0,   0.0d0,  0.0d0,   s35_64,  0.0d0,& ! xxxx
       d38,  0.0d0,  0.0d0,  s5_16,  0.0d0,   0.0d0,  0.0d0,   s35_64,  0.0d0,& ! yyyy
       1d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,   0.0d0,  0.0d0,    0.0d0,  0.0d0,& ! zzzz
     0.0d0,  0.0d0,  0.0d0,  0.0d0, -s5_28,   0.0d0,  0.0d0,    0.0d0,   s5_4,& ! xxxy
     0.0d0,-s45_56,  0.0d0,  0.0d0,  0.0d0,    s5_8,  0.0d0,    0.0d0,  0.0d0,& ! xxxz
     0.0d0,  0.0d0,  0.0d0,  0.0d0, -s5_28,   0.0d0,  0.0d0,    0.0d0,  -s5_4,& ! xyyy
     0.0d0,  0.0d0,-s45_56,  0.0d0,  0.0d0,   0.0d0,  -s5_8,    0.0d0,  0.0d0,& ! yyyz
     0.0d0,  s10_7,  0.0d0,  0.0d0,  0.0d0,   0.0d0,  0.0d0,    0.0d0,  0.0d0,& ! xzzz
     0.0d0,  0.0d0,  s10_7,  0.0d0,  0.0d0,   0.0d0,  0.0d0,    0.0d0,  0.0d0,& ! yzzz
   s27_560,  0.0d0,  0.0d0,  0.0d0,  0.0d0,   0.0d0,  0.0d0,  -s27_16,  0.0d0,& ! xxyy
   -s27_35,  0.0d0,  0.0d0, s27_28,  0.0d0,   0.0d0,  0.0d0,    0.0d0,  0.0d0,& ! xxzz
   -s27_35,  0.0d0,  0.0d0,-s27_28,  0.0d0,   0.0d0,  0.0d0,    0.0d0,  0.0d0,& ! yyzz
     0.0d0,  0.0d0, -s9_56,  0.0d0,  0.0d0,   0.0d0,   s9_8,    0.0d0,  0.0d0,& ! xxyz
     0.0d0, -s9_56,  0.0d0,  0.0d0,  0.0d0,   -s9_8,  0.0d0,    0.0d0,  0.0d0,& ! xyyz
     0.0d0,  0.0d0,  0.0d0,  0.0d0,   s9_7,   0.0d0,  0.0d0,    0.0d0,  0.0d0 & ! xyzz
     /),shape(gsphcar))

  ! gsphcar: l = 4
  !   spherical fchk order: m = 0, 1, -1, 2, -2, 3, -3, 4, -4
  !   Cartesian fchk order:   zzzz yzzz yyzz yyyz yyyy xzzz xyzz xyyz xyyy xxzz xxyz xxyy xxxz xxxy xxxx
  ! This is just a permutation of gsphcar where the rows are arranged in fchk's primitive order
  real*8 :: gsphcar_fchk(9,15) = reshape((/&
  !    0       1      -1       2      -2        3      -3         4      -4
       1d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,   0.0d0,  0.0d0,    0.0d0,  0.0d0,& ! zzzz
     0.0d0,  0.0d0,  s10_7,  0.0d0,  0.0d0,   0.0d0,  0.0d0,    0.0d0,  0.0d0,& ! yzzz
   -s27_35,  0.0d0,  0.0d0,-s27_28,  0.0d0,   0.0d0,  0.0d0,    0.0d0,  0.0d0,& ! yyzz
     0.0d0,  0.0d0,-s45_56,  0.0d0,  0.0d0,   0.0d0,  -s5_8,    0.0d0,  0.0d0,& ! yyyz
       d38,  0.0d0,  0.0d0,  s5_16,  0.0d0,   0.0d0,  0.0d0,   s35_64,  0.0d0,& ! yyyy
     0.0d0,  s10_7,  0.0d0,  0.0d0,  0.0d0,   0.0d0,  0.0d0,    0.0d0,  0.0d0,& ! xzzz
     0.0d0,  0.0d0,  0.0d0,  0.0d0,   s9_7,   0.0d0,  0.0d0,    0.0d0,  0.0d0,& ! xyzz
     0.0d0, -s9_56,  0.0d0,  0.0d0,  0.0d0,   -s9_8,  0.0d0,    0.0d0,  0.0d0,& ! xyyz
     0.0d0,  0.0d0,  0.0d0,  0.0d0, -s5_28,   0.0d0,  0.0d0,    0.0d0,  -s5_4,& ! xyyy
   -s27_35,  0.0d0,  0.0d0, s27_28,  0.0d0,   0.0d0,  0.0d0,    0.0d0,  0.0d0,& ! xxzz
     0.0d0,  0.0d0, -s9_56,  0.0d0,  0.0d0,   0.0d0,   s9_8,    0.0d0,  0.0d0,& ! xxyz
   s27_560,  0.0d0,  0.0d0,  0.0d0,  0.0d0,   0.0d0,  0.0d0,  -s27_16,  0.0d0,& ! xxyy
     0.0d0,-s45_56,  0.0d0,  0.0d0,  0.0d0,    s5_8,  0.0d0,    0.0d0,  0.0d0,& ! xxxz
     0.0d0,  0.0d0,  0.0d0,  0.0d0, -s5_28,   0.0d0,  0.0d0,    0.0d0,   s5_4,& ! xxxy
       d38,  0.0d0,  0.0d0, -s5_16,  0.0d0,   0.0d0,  0.0d0,   s35_64,  0.0d0 & ! xxxx
     /),shape(gsphcar))

contains

  subroutine atomin(m,mesh,qpro)
    use atomicdata
    use tools_math, only: spline
    implicit none

    type(molecule), intent(in) :: m
    type(tmesh), intent(inout) :: mesh
    real*8, intent(out) :: qpro
    
    character*(mline) :: afile
    real*8 :: rmid, rdum(8), h, q, x(3), r, dq, arho1, arho2, arho, rdata
    integer :: intq
    real*8, allocatable, dimension(:) :: a1, b1, c1, f1
    real*8, allocatable, dimension(:) :: promol, hirsh
    integer :: i, j, l, ndata, idum1, idum2, is, istat, kk, nn
    character*2 :: compar
    logical :: ok
    integer :: isenv

    allocate(hirsh(mesh%n),promol(mesh%n),stat=istat)
    if (istat /= 0) call error('atomin','could not allocate hirsh/promol',2)
    hirsh = 0d0
    promol = 0d0

    open(unit=imosa,status='scratch',form='unformatted')
    do i = 1, m%n
       if (m%z(i) < 1) cycle
       rmid = 1d0 / m%z(i)**third

       if (m%z(i) <= 2) then
          ndata = 200
       elseif (m%z(i) <= 10) then
          ndata = 400
       elseif (m%z(i) <= 18) then
          ndata = 600
       elseif (m%z(i) <= 36) then
          ndata = 800
       elseif (m%z(i) <= 54) then
          ndata = 1000
       elseif (m%z(i) <= 86) then
          ndata = 1200
       elseif (m%z(i) <= 94) then
          ndata = 1400
       else
          call error('atomin','atomic number out of range',2)
       endif
       allocate(f1(0:ndata),a1(0:ndata),b1(0:ndata),c1(0:ndata))
       if (istat /= 0) call error('atomin','could not allocate memory for f1, a1, b1, c1',2)
       
       h = 1d0 / (ndata+1)
       f1(0) = 0d0
       do j = 1, ndata
          q = h * j
          rdata = rmid * q / (1.d0-q)
          f1(j) = ftot(j,m%z(i))
       enddo
       call spline(h,f1,a1,b1,c1,ndata,0.d0)

       do kk = 1, mesh%n
          x = mesh%x(:,kk) - m%x(:,i)
          r = sqrt(x(1)**2+x(2)**2+x(3)**2)
          q = r / (r + rmid)
          intq = int((ndata+1) * q)
          dq = q - intq * h
          arho = abs((f1(intq)+dq*(a1(intq)+dq*(b1(intq)+dq*c1(intq)))))/r**2
          hirsh(kk) = arho
          promol(kk) = promol(kk) + arho
       enddo

       write(imosa) (hirsh(j),j=1,mesh%n)
       deallocate(f1,a1,b1,c1)
    enddo
    qpro = sum(mesh%w * promol)

    open(unit=ihrsh,status='scratch',form='unformatted')
    rewind(imosa)
    do i = 1, m%n
       if (m%z(i) < 1) cycle
       read(imosa) (hirsh(j),j=1,mesh%n)
       write(ihrsh) (hirsh(j)/max(promol(j),1d-40),j=1,mesh%n)
    enddo
    deallocate(hirsh,promol)
    close(imosa)

  end subroutine atomin

  !> Read wfn file
  function readwfn(file,egauss) result(m)

    character*(mline), intent(in) :: file
    real*8, intent(out) :: egauss
    type(molecule) :: m

    character*4 :: orbtyp
    character*2 :: dums, elem
    integer :: i, j, istat, imax, icount, ioc, num1, num2
    real*8 :: zreal, ene, ene0
    logical :: isfrac
    integer :: nalpha
    character*8 :: dum1, dum3, dum4
    character*18 :: dum2
    character*1024 :: line

    integer :: idum, idum2

    ! set title
    m%name = file
    m%useecp = .false.

    ! read number of atoms, primitives, orbitals
    open(luwfn,file=file,status='old')
    read (luwfn,*)
    read (luwfn,101) orbtyp, m%nmo, m%npri, m%n

    ! atomic positions and numbers
    allocate(m%x(3,m%n),stat=istat)
    if (istat /= 0) call error('readwfn','could not allocate memory for atomic positions',2)
    allocate(m%z(m%n),stat=istat)
    if (istat /= 0) call error('readwfn','could not allocate memory for atomic numbers',2)
    m%x = 0d0
    m%z = 0
    m%charge = 0d0
    do i = 1, m%n
       read(luwfn,106) elem, m%x(:,i), zreal
       m%charge = m%charge + zreal
       m%z(i) = elem2z(elem)
       if (m%z(i) /= nint(zreal)) m%useecp = .true.
    end do

    ! center assignments, types of primitives
    allocate(m%icenter(m%npri),stat=istat)
    if (istat /= 0) call error('readwfn','could not allocate memory for icenter',2)
    allocate(m%itype(m%npri),stat=istat)
    if (istat /= 0) call error('readwfn','could not allocate memory for itype',2)
    read(luwfn,102) (m%icenter(i),i=1,m%npri)
    read(luwfn,102) (m%itype(i),i=1,m%npri)
    if (any(m%itype(1:m%npri) > 56)) then
       call error("readwfn","primitive type not supported",2)
    endif

    ! primitive exponents
    allocate(m%e(m%npri),stat=istat)
    if (istat /= 0) call error('readwfn','could not allocate memory for exponents',2)
    read(luwfn,103) (m%e(i),i=1,m%npri)
    
    ! deal with ecps
    dums=""
    do while (dums.ne."MO")
       read (luwfn,'(A2)') dums
    enddo
    backspace(luwfn)
    
    ! occupations and orbital coefficients
    allocate(m%occ(m%nmo),stat=istat)
    if (istat /= 0) call error('readwfn','could not allocate memory for occupations',2)
    allocate(m%c(m%nmo,m%npri),stat=istat)
    if (istat /= 0) call error('readwfn','could not allocate memory for orbital coefficients',2)
    m%mult = 1
    isfrac = .false.
    num1 = 0
    num2 = 0
    ene0 = -1d30
    nalpha = -1
    m%nelec = 0d0
    do i = 1, m%nmo
       read(luwfn,104) m%occ(i), ene
       read(luwfn,105) (m%c(i,j),j=1,m%npri)
       m%nelec = m%nelec + m%occ(i)
       ioc = nint(m%occ(i))
       if (abs(ioc-m%occ(i)) > 1d-10) then
          isfrac = .true.
       else if (ioc == 1) then
          num1 = num1 + 1
       else if (ioc == 2) then
          num2 = num2 + 1
       endif
       if (ene < ene0-1d-3) nalpha = i-1
       ene0 = ene
    end do
    read(luwfn,*) dum1
    read(luwfn,'(A80)') line
    line = trim(adjustl(line))
    line = line(index(line,'=')+1:)
    read (line,*) egauss

    ! figure out charge and multiplicity
    ! 0 - closed, 1 - open, 2 - restricted open, 3 - fractional
    m%charge = m%charge - m%nelec
    if (isfrac) then
       m%wfntyp = 3
       call error("readwfn","natural orbital wfn files not supported",2)
    else if (num1 == 0) then
       m%wfntyp = 0
    else if (num2 == 0) then
       m%wfntyp = 1
       if (nalpha == -1) then
          m%mult = nint(m%nelec) + 1
       else
          m%mult = nalpha - (nint(m%nelec) - nalpha) + 1
       endif
       if (m%mult < 0) call error("readwfn","nbeta > nalpha",2)
    else
       m%wfntyp = 2
       m%mult = count(nint(m%occ(1:m%nmo)) == 1) + 1
    endif
    close(luwfn)

101 format (4X,A4,10X,3(I5,15X))
102 format(20X,20I3)
103 format(10X,5E14.7)
104 format(35X,F12.7,15X,F12.6)
105 format(5(E16.8))
106 format(2X,A2,20X,3F12.8,10X,F5.1)

  end function readwfn

  !> Read wfn file
  function readwfx(file,egauss) result(m)

    character*(mline), intent(in) :: file
    real*8, intent(out) :: egauss
    type(molecule) :: m

    integer :: i, j, istat, ncore, kk, lp, idum
    real*8 :: zreal
    character*(mline) :: line, tag
    logical :: keyw(7)

    ! set title
    m%name = file
    m%useecp = .false.

    ! first pass
    open (luwfn,file=file,status='old')
    m%n = 0
    m%nmo = 0
    m%charge = 0
    m%mult = 0
    ncore = 0
    m%npri = 0
    egauss = 0d0
    do while (.true.)
       read(luwfn,'(A)',end=10) line
       line = adjustl(line)
       if (line(1:1) == "<" .and. line(2:2) /= "/") then
          if (trim(line) == "<Number of Nuclei>") then
             read (luwfn,*) m%n
          elseif (trim(line) == "<Number of Occupied Molecular Orbitals>") then
             read (luwfn,*) m%nmo
          elseif (trim(line) == "<Net Charge>") then
             read (luwfn,*) m%charge
          elseif (trim(line) == "<Electronic Spin Multiplicity>") then
             read (luwfn,*) m%mult
          elseif (trim(line) == "<Number of Core Electrons>") then
             read (luwfn,*) ncore
          elseif (trim(line) == "<Number of Primitives>") then
             read(luwfn,*) m%npri
          elseif (trim(line) == "<Energy = T + Vne + Vee + Vnn>") then
             read(luwfn,*) egauss
          endif
       endif
    enddo
10  continue
    
    if (m%n == 0) call error("readwfx","Number of Nuclei tag not found",2)
    if (m%nmo == 0) call error("readwfx","Number of Occupied Molecular Orbitals tag not found",2)
    if (m%mult == 0) call error("readwfx","Electronic Spin Multiplicity tag not found",2)
    if (m%npri == 0) call error("readwfx","Number of Primitives tag not found",2)
    if (ncore > 0) m%useecp = .true.

    ! allocate memory
    allocate(m%x(3,m%n),stat=istat)
    if (istat /= 0) call error('readwfx','could not allocate memory for atomic positions',2)
    allocate(m%z(m%n),stat=istat)
    if (istat /= 0) call error('readwfx','could not allocate memory for atomic numbers',2)
    allocate(m%icenter(m%npri),stat=istat)
    if (istat /= 0) call error('readwfx','could not allocate memory for icenter',2)
    allocate(m%itype(m%npri),stat=istat)
    if (istat /= 0) call error('readwfx','could not allocate memory for itype',2)
    allocate(m%e(m%npri),stat=istat)
    if (istat /= 0) call error('readwfx','could not allocate memory for exponents',2)
    allocate(m%occ(m%nmo),stat=istat)
    if (istat /= 0) call error('readwfx','could not allocate memory for occupations',2)
    allocate(m%c(m%nmo,m%npri),stat=istat)
    if (istat /= 0) call error('readwfx','could not allocate memory for orbital coefficients',2)
    
    ! second pass
    rewind(luwfn)
    keyw = .false.
    do while (.true.)
       read(luwfn,'(A)',end=20) line
       line = adjustl(line)
       if (line(1:1) == "<" .and. line(2:2) /= "/") then
          if (trim(line) == "<Atomic Numbers>") then
             m%z = read_integers(luwfn,m%n)
             keyw(1) = .true.
          elseif (trim(line) == "<Nuclear Cartesian Coordinates>") then
             m%x = reshape(read_reals1(luwfn,3*m%n),shape(m%x))
             keyw(2) = .true.
          elseif (trim(line) == "<Primitive Centers>") then
             m%icenter = read_integers(luwfn,m%npri)
             keyw(3) = .true.
          elseif (trim(line) == "<Primitive Types>") then
             m%itype = read_integers(luwfn,m%npri)
             if (any(m%itype(1:m%npri) > 56)) &
                call error("readwfx","primitive type not supported",2)
             keyw(4) = .true.
          elseif (trim(line) == "<Primitive Exponents>") then
             m%e = read_reals1(luwfn,m%npri)
             keyw(5) = .true.
          elseif (trim(line) == "<Molecular Orbital Occupation Numbers>") then
             m%occ = read_reals1(luwfn,m%nmo)
             m%nelec = sum(m%occ)
             keyw(6) = .true.
          elseif (trim(line) == "<Molecular Orbital Primitive Coefficients>") then
             read(luwfn,*)
             do i = 1, m%nmo
                read(luwfn,*)
                read(luwfn,*)
                m%c(i,:) = read_reals1(luwfn,m%npri)
             enddo
             keyw(7) = .true.
          endif
       endif
    enddo
20  continue
    if (any(.not.keyw)) call error("readwfx","missing array in wfx file",2)

    ! wavefuntion type
    ! 0 - closed, 1 - open, 2 - restricted open, 3 - fractional
    if (m%mult == 1) then
       m%wfntyp = 0
    else
       if (any(m%occ > 1)) then
          m%wfntyp = 2
       else
          m%wfntyp = 1
       endif
    end if

    close(luwfn)

  end function readwfx
  
  !> Read fchk file
  function readfchk(file,egauss) result(m)

    character*(mline), intent(in) :: file
    real*8, intent(out) :: egauss
    type(molecule) :: m

    character*(mline) :: line
    integer :: ityp
    integer :: lp, idum, nalpha, nbeta, nshel, ncshel, lmax, nbassph, nbascar
    integer :: istat, i, j, k, k1, k2, l, ifac
    integer :: acent, nn, nm, nl, nc, ns, ncar, nsph
    logical :: ok, isecp
    logical, allocatable :: icdup(:)
    integer, allocatable :: ishlt(:), ishlpri(:), ishlat(:), itemp(:,:)
    real*8, allocatable :: xat(:), exppri(:), ccontr(:), pccontr(:), motemp(:), mocoef(:,:)
    real*8, allocatable :: cnorm(:), cpri(:), rtemp(:,:)
    real*8 :: acoef, norm

    ! translation between primitive ordering fchk -> postg
    !         1   2 3 4    5  6  7  8  9 10    11  12  13  14  15  16  17  18  19  20
    ! fchk:   s   x y z   xx yy zz xy xz yz   xxx yyy zzz xyy xxy xxz xzz yzz yyz xyz
    ! postg:  s   x y z   xx yy zz xy xz yz   xxx yyy zzz xxy xxz yyz xyy xzz yzz xyz
    !
    !           21   22   23   24   25   26   27   28   29   30   31   32   33   34   35   
    ! fchk:   zzzz yzzz yyzz yyyz yyyy xzzz xyzz xyyz xyyy xxzz xxyz xxyy xxxz xxxy xxxx
    ! postg:  xxxx yyyy zzzz xxxy xxxz xyyy yyyz xzzz yzzz xxyy xxzz yyzz xxyz xyyz xyzz 

    ! reorder the primitives within the same shell
    ! typtrans(fchk) = postg
    integer, parameter :: typtrans(35) = (/&
       !1   2  3  4    5  6  7  8  9  10    11  12  13  14  15  16  17  18  19  20
       1,   2, 3, 4,   5, 6, 7, 8, 9, 10,   11, 12, 13, 17, 14, 15, 18, 19, 16, 20,&
       !21 22  23  24  25  26  27  28  29  30  31  32  33  34  35
       23, 29, 32, 27, 22, 28, 35, 34, 26, 31, 33, 30, 25, 24, 21&
       /)

    ! set title and defaults
    m%name = file
    m%useecp = .false.
    m%wfntyp = 0

    ! first pass: dimensions
    isecp = .false.
    open(luwfn,file=file,status='old')
    do while (.true.)
       read(luwfn,'(A)',end=20) line
       line = adjustl(line)
       lp = 45
       if (line(1:6) == "Charge") then
          ok = isinteger(idum,line,lp)
          m%charge = idum
       elseif (line(1:12) == "Multiplicity") then
          ok = isinteger(m%mult,line,lp)
       elseif (line(1:19) == "Number of electrons") then
          ok = isinteger(idum,line,lp)
          m%nelec = idum
       elseif (line(1:15) == "Number of atoms") then
          ok = isinteger(m%n,line,lp)
       elseif (line(1:25) == "Number of alpha electrons") then
          ok = isinteger(nalpha,line,lp)
       elseif (line(1:25) == "Number of basis functions") then
          ok = isinteger(nbassph,line,lp)
       elseif (line(1:24) == "Number of beta electrons") then
          ok = isinteger(nbeta,line,lp)
       elseif (line(1:27) == "Number of contracted shells") then
          ok = isinteger(ncshel,line,lp)
       elseif (line(1:26) == "Number of primitive shells") then
          ok = isinteger(nshel,line,lp)
       elseif (line(1:24) == "Highest angular momentum") then
          ok = isinteger(lmax,line,lp)
       elseif (line(1:12) == "Total Energy") then
          ok = isreal(egauss,line,lp)
       elseif (line(1:21) == "Beta Orbital Energies") then
          m%wfntyp = 1
       elseif (line(1:8) == "ECP-LMax") then
          isecp = .true.
       endif
    enddo
20  continue
    
    ! ECPs not implemented yet
    if (isecp) call error("readfchk","ECPs not supported.",1)

    ! Count the number of MOs
    if (m%wfntyp == 0) then
       m%nmo = nint(m%nelec) / 2
    else if (m%wfntyp == 1) then
       m%nmo = nint(m%nelec)
    endif

    ! allocate sutff
    allocate(m%occ(m%nmo),stat=istat)
    if (istat /= 0) call error('readfchk','could not allocate memory for occ',2)
    allocate(ishlt(ncshel),ishlpri(ncshel),ishlat(ncshel),stat=istat)
    if (istat /= 0) call error('readfchk','could not allocate memory for shell data',2)
    allocate(exppri(nshel),ccontr(nshel),pccontr(nshel),stat=istat)
    if (istat /= 0) call error('readfchk','could not allocate memory for primitive data',2)

    ! type of wavefunction -> occupations
    if (m%wfntyp == 1) then
       m%occ = 1
    else
       m%occ = 2
    end if

    ! rewind
    rewind(luwfn)

    ! second pass
    allocate(m%x(3,m%n),m%z(m%n),xat(3*m%n),stat=istat)
    if (istat /= 0) call error('readfchk','could not allocate memory for geometry',2)
    allocate(motemp(nbassph*m%nmo),stat=istat)
    if (istat /= 0) call error('readfchk','could not allocate memory for MO coefs',2)
    do while (.true.)
       read(luwfn,'(A)',end=30) line
       line = adjustl(line)
       lp = 45
       if (line(1:11) == "Shell types") then
          do i = 0, (ncshel-1)/6
             read(luwfn,'(6I12)',end=30) (ishlt(6*i+j),j=1,min(6,ncshel-6*i))
          enddo
       elseif (line(1:30) == "Number of primitives per shell") then
          do i = 0, (ncshel-1)/6
             read(luwfn,'(6I12)',end=30) (ishlpri(6*i+j),j=1,min(6,ncshel-6*i))
          enddo
       elseif (line(1:30) == "Shell to atom map") then
          do i = 0, (ncshel-1)/6
             read(luwfn,'(6I12)',end=30) (ishlat(6*i+j),j=1,min(6,ncshel-6*i))
          enddo
       elseif (line(1:29) == "Current cartesian coordinates") then
          do i = 0, (3*m%n-1)/5
             read(luwfn,'(5E16.8)',end=30) (xat(5*i+j),j=1,min(5,3*m%n-5*i))
          enddo
       elseif (line(1:14) == "Atomic numbers") then
          do i = 0, (m%n-1)/6
             read(luwfn,'(6I12)',end=30) (m%z(6*i+j),j=1,min(6,m%n-6*i))
          enddo
       elseif (line(1:19) == "Primitive exponents") then
          do i = 0, (nshel-1)/5
             read(luwfn,'(5E16.8)',end=30) (exppri(5*i+j),j=1,min(5,nshel-5*i))
          enddo
       elseif (line(1:24) == "Contraction coefficients") then
          do i = 0, (nshel-1)/5
             read(luwfn,'(5E16.8)',end=30) (ccontr(5*i+j),j=1,min(5,nshel-5*i))
          enddo
       elseif (line(1:31) == "P(S=P) Contraction coefficients") then
          do i = 0, (nshel-1)/5
             read(luwfn,'(5E16.8)',end=30) (pccontr(5*i+j),j=1,min(5,nshel-5*i))
          enddo
       elseif (line(1:21) == "Alpha MO coefficients") then
          do i = 0, (nalpha*nbassph-1)/5
             read(luwfn,'(5E16.8)',end=30) (motemp(5*i+j),j=1,min(5,nalpha*nbassph-5*i))
          enddo
       elseif (line(1:21) == "Beta MO coefficients") then
          do i = 0, (nbeta*nbassph-1)/5
             read(luwfn,'(5E16.8)',end=30) (motemp(nalpha*nbassph+5*i+j),j=1,min(5,nbeta*nbassph-5*i))
          enddo
       endif
    enddo
30  continue

    ! we are done with the file
    close(luwfn)

    ! unfold sp shells next
    allocate(icdup(ncshel))
    nshel = 0
    ncshel = 0
    do i = 1, size(icdup)
       if (ishlt(i) == -1) then
          icdup(i) = .true.
          ncshel = ncshel + 2
          nshel = nshel + 2 * ishlpri(i)
       else
          icdup(i) = .false.
          ncshel = ncshel + 1
          nshel = nshel + ishlpri(i)
       endif
    end do

    ! transfer the information to the temporary sp-unfolded arrays
    allocate(itemp(ncshel,3),rtemp(nshel,2))
    nn = 0
    nm = 0
    nl = 0
    do i = 1, size(icdup)
       if (icdup(i)) then
          nn = nn + 1
          itemp(nn,1) = 0
          itemp(nn,2) = ishlpri(i)
          itemp(nn,3) = ishlat(i)
          rtemp(nl+1:nl+ishlpri(i),1) = exppri(nm+1:nm+ishlpri(i))
          rtemp(nl+1:nl+ishlpri(i),2) = ccontr(nm+1:nm+ishlpri(i))
          nl = nl + ishlpri(i)

          nn = nn + 1
          itemp(nn,1) = 1
          itemp(nn,2) = ishlpri(i)
          itemp(nn,3) = ishlat(i)
          rtemp(nl+1:nl+ishlpri(i),1) = exppri(nm+1:nm+ishlpri(i))
          rtemp(nl+1:nl+ishlpri(i),2) = pccontr(nm+1:nm+ishlpri(i))
          nl = nl + ishlpri(i)
       else
          nn = nn + 1
          itemp(nn,1) = ishlt(i)
          itemp(nn,2) = ishlpri(i)
          itemp(nn,3) = ishlat(i)
          rtemp(nl+1:nl+ishlpri(i),1) = exppri(nm+1:nm+ishlpri(i))
          rtemp(nl+1:nl+ishlpri(i),2) = ccontr(nm+1:nm+ishlpri(i))
          nl = nl + ishlpri(i)
       endif
       nm = nm + ishlpri(i)
    end do

    ! move the sp-unfolded information back and reallocate
    deallocate(ishlt,ishlpri,ishlat,exppri,ccontr,pccontr,icdup)
    allocate(ishlt(ncshel),ishlpri(ncshel),ishlat(ncshel),exppri(nshel),ccontr(nshel))
    ishlt = itemp(:,1)
    ishlpri = itemp(:,2)
    ishlat = itemp(:,3)
    exppri = rtemp(:,1)
    ccontr = rtemp(:,2)
    deallocate(itemp,rtemp)

    ! atomic coordinates
    m%x = reshape(xat,shape(m%x))
    deallocate(xat)

    ! count the number of primitives and basis functions
    m%npri = 0
    nbascar = 0
    nbassph = 0
    do i = 1, ncshel
       ityp = ishlt(i)
       nbascar = nbascar + nshlt(abs(ityp))
       nbassph = nbassph + nshlt(ityp)
       m%npri = m%npri + nshlt(abs(ishlt(i))) * ishlpri(i)
    enddo

    ! convert spherical basis functions to Cartesian and build the mocoef
    ! deallocate the temporary motemp
    allocate(mocoef(m%nmo,nbascar))
    nc = 0
    ns = 0
    do j = 1, ncshel
       nsph = nshlt(ishlt(j))
       ncar = nshlt(abs(ishlt(j)))
       if (nsph == ncar) then
          do i = 1, m%nmo
             mocoef(i,nc+1:nc+ncar) = motemp((i-1)*nbassph+ns+1:(i-1)*nbassph+ns+nsph)
          end do
       elseif (ishlt(j) == -2) then
          do i = 1, m%nmo
             mocoef(i,nc+1:nc+ncar) = matmul(motemp((i-1)*nbassph+ns+1:(i-1)*nbassph+ns+nsph),dsphcar)
          end do
       elseif (ishlt(j) == -3) then
          do i = 1, m%nmo
             mocoef(i,nc+1:nc+ncar) = matmul(motemp((i-1)*nbassph+ns+1:(i-1)*nbassph+ns+nsph),fsphcar)
          end do
       elseif (ishlt(j) == -4) then
          do i = 1, m%nmo
             mocoef(i,nc+1:nc+ncar) = matmul(motemp((i-1)*nbassph+ns+1:(i-1)*nbassph+ns+nsph),gsphcar_fchk)
          end do
       else
          call error('readfchk','h and higher primitives not supported yet',2)
       endif
       ns = ns + nsph
       nc = nc + ncar
    end do
    deallocate(motemp)

    ! Assign primitive center and type, exponents, etc.
    allocate(m%icenter(m%npri),stat=istat)
    if (istat /= 0) call error('readfchk','could not allocate memory for icenter',2)
    allocate(m%itype(m%npri),stat=istat)
    if (istat /= 0) call error('readfchk','could not allocate memory for itype',2)
    allocate(m%e(m%npri),stat=istat)
    if (istat /= 0) call error('readfchk','could not allocate memory for exponents',2)
    allocate(m%c(m%nmo,m%npri),cpri(m%npri),stat=istat)
    if (istat /= 0) call error('readfchk','could not allocate memory for coeffs',2)

    ! normalize the primitive coefficients without the angular part
    allocate(cnorm(maxval(ishlpri)))
    nm = 0
    nn = 0
    do i = 1, ncshel
       do j = jshl0(abs(ishlt(i))), jshl1(abs(ishlt(i)))
          ityp = typtrans(j)

          ! normalize primitive coefficients
          do k = 1, ishlpri(i)
             cnorm(k) = ccontr(nm+k) * gnorm(ityp,exppri(nm+k)) 
          end do

          ! normalize basis function
          norm = aonorm(ityp,ishlpri(i),exppri(nm+1:nm+ishlpri(i)),cnorm(1:ishlpri(i)))

          ! calculate and assign the normalized primitive coefficients
          do k = 1, ishlpri(i)
             nn = nn + 1
             cpri(nn) = cnorm(k) * norm
          end do
       end do
       nm = nm + ishlpri(i)
    end do
    deallocate(cnorm)

    ! build the wavefunction coefficients for the primitives
    nn = 0
    nm = 0
    nl = 0
    do i = 1, ncshel
       do j = jshl0(abs(ishlt(i))), jshl1(abs(ishlt(i)))
          ityp = typtrans(j)
          nl = nl + 1
          do k = 1, ishlpri(i)
             nn = nn + 1
             m%icenter(nn) = ishlat(i)
             m%itype(nn) = ityp
             m%e(nn) = exppri(nm+k)
             m%c(:,nn) = cpri(nn) * mocoef(:,nl)
          end do
       end do
       nm = nm + ishlpri(i)
    end do

    deallocate(ishlt,ishlpri,ishlat)
    deallocate(exppri,ccontr)
    deallocate(mocoef)


  end function readfchk

  !> Read a molden file. Different programs write molden files in
  !> different ways.  Right now, this routine has only been tested
  !> with psi4 (only versions starting June 2016; previous versions of
  !> psi4 wrote molden files differently) and orca (4.0.2).  Also,
  !> psi4 writes restricted open wavefunctions the same as
  !> unrestricted wavefunctions, so that is not tested either, and
  !> probably won't work.
  function readmolden(file,egauss) result(m)

    character*(mline), intent(in) :: file
    real*8, intent(out) :: egauss
    type(molecule) :: m

    character*(mline) :: line, word, word1, word2, keyword
    integer :: ncshel, nshel, nbascar, nbassph
    integer :: i, j, k, l, idum, idum1, ni, nj, k1, k2
    real*8 :: rdum, norm, cons
    integer, allocatable :: ishlt(:), ishlpri(:), ishlat(:)
    real*8, allocatable :: exppri(:), ccontr(:), motemp(:), cpri(:), mocoef(:,:), cnorm(:)
    integer :: istat, nm, ityp, lnmoa, lnmob, lnmo
    logical :: isang, is5d, is7f, is9g
    integer :: nn, nl, nalpha, nsph, ncar, nc, ns, imoldentype
    logical :: ok, isalpha

    integer, parameter :: mt_psi4 = 1
    integer, parameter :: mt_orca = 2

    ! guide to the variables in this routine:
    !   ncshel = number of contraction shells (# of blocks in [GTO])
    !     -> ishlt(icshel) = angular momentum for this shell (s=0,p=1,d=2,dsph=-2,etc.)
    !     -> ishlpri(icshel) = number of primitives in this shell
    !     -> ishlat(icshel) = atom index where this shell is based
    !
    !   nshel = number of primitives (not counting the angular parts)
    !     -> exppri(ishel) = primitive exponent
    !     -> ccontr(ishel) = contraction coefficient for this primitive
    !
    !   m%npri = number of primitives (including angular parts)
    !
    !   nbassph = number of basis functions (spherical)
    !
    !   nbascar = number of basis functions (Cartesian)
    !
    !   lnmo = number of molecular orbitals, regardless of occupation
    !   m%nmo = number of occupied molecular orbitals
    !     -> m%occ(imo) = occupation of the molecular orbitals
    !
    !   m%nelec = number of electrons
    !   nalpha = number of alpha electrons
    !   m%wfntyp = type of wavefunction (0 = closed, 1 = open)
    !   m%charge = charge
    !   m%mult = multiplicity
    ! 

    ! translation between primitive ordering molden -> postg
    !         1   2 3 4    5  6  7  8  9 10    11  12  13  14  15  16  17  18  19  20
    ! molden: s   x y z   xx yy zz xy xz yz   xxx yyy zzz xyy xxy xxz xzz yzz yyz xyz
    ! postg:  s   x y z   xx yy zz xy xz yz   xxx yyy zzz xxy xxz yyz xyy xzz yzz xyz
    !
    !           21   22   23   24   25   26   27   28   29   30   31   32   33   34   35   
    ! molden: xxxx yyyy zzzz xxxy xxxz xyyy yyyz xzzz yzzz xxyy xxzz yyzz xxyz xyyz xyzz
    ! postg:  xxxx yyyy zzzz xxxy xxxz xyyy yyyz xzzz yzzz xxyy xxzz yyzz xxyz xyyz xyzz 
    !
    ! h primitives and higher not supported in molden format.
    !
    !     typtrans(molden) = postg
    integer, parameter :: typtrans(35) = (/&
       !1   2  3  4    5  6  7  8  9  10    11  12  13  14  15  16  17  18  19  20
       1,   2, 3, 4,   5, 6, 7, 8, 9, 10,   11, 12, 13, 17, 14, 15, 18, 19, 16, 20,&
       !21 22  23  24  25  26  27  28  29  30  31  32  33  34  35
       21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35&
       /)

    line = ""
    m%name = file
    m%useecp = .false.
    is5d = .false.
    is7f = .false.
    is9g = .false.
    imoldentype = mt_psi4

    ! parse the molden file, first pass -> read dimensions prior to allocation
    open(luwfn,file=file,status='old')

    do while(next_keyword())
       if (trim(keyword) == "title") then
          read(luwfn,'(A)',end=40) line
          if (index(line,"created by orca_2mkl") > 0) then
             imoldentype = mt_orca
          end if
       else if (trim(keyword) == "atoms") then
          ! read the number of atoms 
          m%n = 0
          read(luwfn,'(A)',end=40) line
          do while(index(lower(line),"[") == 0 .and. len(trim(line)) > 0)
             m%n = m%n + 1
             read(luwfn,'(A)',end=40) line
          end do
       elseif (trim(keyword) == "gto") then
          ! read the number of shells and primitives
          ncshel = 0
          nshel = 0
          do i = 1, m%n
             read(luwfn,'(A)',end=40) line
             read(luwfn,'(A)',end=40) line
             do while (index(line,".") /= 0)
                read (line,*) word, idum, rdum
                ncshel = ncshel + 1
                nshel = nshel + idum
                do j = 1, idum
                   read(luwfn,'(A)',end=40) line
                end do
                read(luwfn,'(A)',end=40) line
             end do
          end do
       elseif (trim(keyword) == "mo") then
          ! read the number of MOs, number of electrons, and number of alpha electrons
          m%nelec = 0
          nalpha = 0
          m%wfntyp = -1
          do while(.true.)
             read(luwfn,'(A)',end=30) line
             if (index(lower(line),"[") /= 0) exit
             if (index(lower(line),"ene=") /= 0) then
                ! spin
                read(luwfn,'(A)',end=30) line
                read (line,*) word1, word2
                isalpha = (trim(lower(word2)) == "alpha")

                ! occupation
                read(luwfn,'(A)',end=30) line
                read (line,*) word, rdum
                idum = nint(rdum)
                if (abs(rdum-idum) > 1d-6) &
                   call error('readmolden','can not do fractional occupations with molden',2)
                if (idum == 1 .or. idum == 2) then
                   m%nelec = m%nelec + idum
                   if (isalpha) then
                      nalpha = nalpha + idum
                   endif
                   if (m%wfntyp < 0 .and. idum == 2) m%wfntyp = 0
                   if (idum == 1) m%wfntyp = 1
                elseif (idum == 0) then
                   ! do not read MOs with zero occupation
                   continue
                else
                   call error('readmolden','wrong integer occupation',2)
                endif
             end if
          end do
30        continue
       elseif (trim(keyword) == "5d" .or. trim(keyword) == "5d7f") then
          is5d = .true.
          is7f = .true.
          read(luwfn,'(A)',end=40) line
       elseif (trim(keyword) == "5d10f") then
          is5d = .true.
          is7f = .false.
          read(luwfn,'(A)',end=40) line
       elseif (trim(keyword) == "7f") then
          is5d = .false.
          is7f = .true.
          read(luwfn,'(A)',end=40) line
       elseif (trim(keyword) == "9g") then
          is9g = .true.
          read(luwfn,'(A)',end=40) line
       else
          ! didn't catch a keyword, move on
          read(luwfn,'(A)',end=40) line
       end if
    end do

    ! orca always uses spherical
    if (imoldentype == mt_orca) then
       is5d = .true.
       is7f = .true.
       is9g = .true.
    end if

    ! type of wavefunction -> number of MOs
    if (m%wfntyp == 1) then
       m%nmo = nint(m%nelec)
    else
       m%nmo = nint(m%nelec) / 2
       nalpha = nalpha / 2
    end if
       
    ! allocate stuff
    allocate(m%occ(m%nmo),stat=istat)
    if (istat /= 0) call error('readmolden','alloc. memory for occ',2)
    allocate(ishlt(ncshel),ishlpri(ncshel),ishlat(ncshel),stat=istat)
    if (istat /= 0) call error('readmolden','alloc. memory for shell types',2)
    allocate(exppri(nshel),ccontr(nshel),stat=istat)
    if (istat /= 0) call error('readmolden','alloc. memory for prim. shells',2)

    ! type of wavefunction -> occupations
    if (m%wfntyp == 1) then
       m%occ = 1
    else
       m%occ = 2
    end if
       
    ! rewind
    rewind(luwfn)

    ! second pass: populate the arrays
    ! read the geometry
    read(luwfn,'(A)',end=40) line
    do while(.not.trim(lower(line)) == "[molden format]")
       read(luwfn,'(A)',end=40) line
    end do

    ! geometry header -> detect the units for the geometry
    read(luwfn,'(A)',end=40) line
    do while(index(lower(line),"[atoms]") == 0)
       read(luwfn,'(A)',end=40) line
    end do
    read(line,*) word1, word2
    isang = (trim(lower(word2)) == "(angs)".or.trim(lower(word2)) == "(ang)") 

    ! read the atomic numbers and positions
    read(luwfn,'(A)',end=40) line
    allocate(m%x(3,m%n),m%z(m%n),stat=istat)
    if (istat /= 0) call error('readmolden','could not allocate memory for geometry',2)
    do i = 1, m%n
       read(line,*) word1, idum1, m%z(i), m%x(:,i)
       read(luwfn,'(A)',end=40) line
    end do
    if (isang) m%x = m%x / 0.52917720859d0

    ! calculate charge and multiplicity. No energy info in molden files.
    m%charge = 0
    do i = 1, m%n
       m%charge = m%charge + m%z(i)
    end do
    m%charge = m%charge - m%nelec
    m%mult = abs(nalpha - (nint(m%nelec) - nalpha)) + 1
    egauss = 0d0

    ! Basis set specification 
    do while(.not.trim(lower(line)) == "[gto]")
       read(luwfn,'(A)',end=40) line
    end do
    ni = 0
    nj = 0
    do i = 1, m%n
       read(luwfn,'(A)',end=40) line
       read(luwfn,'(A)',end=40) line
       ! atomic blocks are separated by blank lines
       do while (len_trim(line) > 0)
          ni = ni + 1

          ! header -> l, npri, ???
          read(line,*) word, idum 

          ! read exponents and coefficients
          do j = 1, idum
             nj = nj + 1
             read(luwfn,'(A)',end=40) line
             read(line,*) exppri(nj), ccontr(nj)
          end do

          ! assign atomic centers, number of primitives, and shell types.
          ishlat(ni) = i
          ishlpri(ni) = idum
          if (lower(trim(word)) == "s") then
             ishlt(ni) = 0
          else if (lower(trim(word)) == "p") then
             ishlt(ni) = 1
          else if (lower(trim(word)) == "sp") then
             call error("readmolden","can't handle SP in gamess format",2)
          else if (lower(trim(word)) == "d") then
             if (is5d) then
                ishlt(ni) = -2
             else
                ishlt(ni) = 2
             end if
          else if (lower(trim(word)) == "f") then
             if (is7f) then
                ishlt(ni) = -3
             else
                ishlt(ni) = 3
             end if
          else if (lower(trim(word)) == "g") then
             if (is9g) then
                ishlt(ni) = -4
             else
                ishlt(ni) = 4
             end if
          else
             call error("readmolden","basis set type >g not supported in molden files",2)
          endif
          read(luwfn,'(A)',end=40) line
       end do
    end do

    ! Count the number of primitives and basis functions
    m%npri = 0
    nbascar = 0
    nbassph = 0
    do i = 1, ncshel
       ityp = ishlt(i)
       nbascar = nbascar + nshlt(abs(ityp))
       nbassph = nbassph + nshlt(ityp)
       m%npri = m%npri + nshlt(abs(ityp)) * ishlpri(i)
    enddo

    ! Allocate space to read the the molecular orbital information
    allocate(motemp(nbassph*m%nmo),stat=istat)
    if (istat /= 0) call error('readmolden','alloc. memory for MO coefs',2)

    ! advance to the MO coefficients
    do while(index(lower(line),"[mo]") == 0)
       read(luwfn,'(A)',end=40) line
    end do

    ! read the MO coefficients
    lnmoa = 0
    lnmob = nalpha
    lnmo = 0
    do while(.true.)
       read(luwfn,'(A)',end=10) line

       ! is this an alpha electron?
       if (index(lower(line),"spin=") /= 0 .and. m%wfntyp > 0) then
          read (line,*) word1, word2
          isalpha = (trim(lower(word2)) == "alpha") 
       end if

       ! is this orbital occupied?
       ! assumes that the MO coefficiens come right after the "Occup" tag
       if (index(lower(line),"occup=") /= 0) then
          read(line,*) word, rdum
          idum = nint(rdum)
          if (idum > 0) then

             if (m%wfntyp > 0 .and. isalpha) then
                lnmoa = lnmoa + 1
                lnmo = lnmoa
             elseif (m%wfntyp > 0) then
                lnmob = lnmob + 1
                lnmo = lnmob 
             else
                lnmo = lnmo + 1
             endif

             do j = 1, nbassph
                read(luwfn,*,end=40) idum, motemp((lnmo-1)*nbassph+j)
             end do
          end if
       end if
    end do
10  continue

    ! we are done with the file
    close(luwfn)

    ! check the number of MOs is correct
    if (m%wfntyp == 0) then
       if (m%nmo /= lnmo) &
          call error('readmolden','inconsistent number of MOs in the second pass',2)
    else
       if (nalpha /= lnmoa) &
          call error('readmolden','inconsistent number of MOs (nalpha) in the second pass',2)
       if (m%nmo /= lnmob) &
          call error('readmolden','inconsistent number of MOs (total) in the second pass',2)
    endif

    ! If this is an orca molden file, flip the sign of the MO coefficients in these cases:
    ! - abs(m)=3 for the f shells
    ! - abs(m)>=3 for the g shells
    if (imoldentype == mt_orca) then
       do j = 1, ncshel
          nsph = nshlt(ishlt(j))
          if (ishlt(j) == -3) then
             do i = 1, m%nmo
                motemp((i-1)*nbassph+ns+6:(i-1)*nbassph+ns+7) = -motemp((i-1)*nbassph+ns+6:(i-1)*nbassph+ns+7)
             end do
          else if (ishlt(j) == -4) then
             do i = 1, m%nmo
                motemp((i-1)*nbassph+ns+6:(i-1)*nbassph+ns+9) = -motemp((i-1)*nbassph+ns+6:(i-1)*nbassph+ns+9)
             end do
          end if
          ns = ns + nsph
       end do
    end if

    ! convert spherical basis functions to Cartesian and build the mocoef
    ! deallocate the temporary motemp
    allocate(mocoef(m%nmo,nbascar))
    nc = 0
    ns = 0
    do j = 1, ncshel
       nsph = nshlt(ishlt(j))
       ncar = nshlt(abs(ishlt(j)))
       if (nsph == ncar) then
          do i = 1, m%nmo
             mocoef(i,nc+1:nc+ncar) = motemp((i-1)*nbassph+ns+1:(i-1)*nbassph+ns+nsph)
          end do
       elseif (ishlt(j) == -2) then
          do i = 1, m%nmo
             mocoef(i,nc+1:nc+ncar) = matmul(motemp((i-1)*nbassph+ns+1:(i-1)*nbassph+ns+nsph),dsphcar)
          end do
       elseif (ishlt(j) == -3) then
          do i = 1, m%nmo
             mocoef(i,nc+1:nc+ncar) = matmul(motemp((i-1)*nbassph+ns+1:(i-1)*nbassph+ns+nsph),fsphcar)
          end do
       elseif (ishlt(j) == -4) then
          do i = 1, m%nmo
             mocoef(i,nc+1:nc+ncar) = matmul(motemp((i-1)*nbassph+ns+1:(i-1)*nbassph+ns+nsph),gsphcar)
          end do
       else
          call error('readmolden','h and higher primitives not supported in molden format',2)
       endif
       ns = ns + nsph
       nc = nc + ncar
    end do
    deallocate(motemp)

    ! from this point on, only Cartesians
    ishlt = abs(ishlt)

    ! Allocate molecule arrays for primitive center, type, exponents, etc.
    allocate(m%icenter(m%npri),stat=istat)
    if (istat /= 0) call error('readmolden','could not allocate memory for icenter',2)
    allocate(m%itype(m%npri),stat=istat)
    if (istat /= 0) call error('readmolden','could not allocate memory for itype',2)
    allocate(m%e(m%npri),stat=istat)
    if (istat /= 0) call error('readmolden','could not allocate memory for exponents',2)
    allocate(m%c(m%nmo,m%npri),cpri(m%npri),stat=istat)
    if (istat /= 0) call error('readmolden','could not allocate memory for coeffs',2)

    ! normalize the primitive coefficients without the angular part
    allocate(cnorm(maxval(ishlpri)))
    nm = 0
    nn = 0
    do i = 1, ncshel
       do j = jshl0(ishlt(i)), jshl1(ishlt(i))
          ityp = typtrans(j)
          
          ! normalize primitive coefficients
          do k = 1, ishlpri(i)
             cnorm(k) = ccontr(nm+k) * gnorm(ityp,exppri(nm+k))

             ! remove the normalization factor from the orca molden files
             if (imoldentype == mt_orca) then
                cnorm(k) = cnorm(k) / gnorm_orca(ityp,exppri(nm+k))
             end if
          end do

          ! AO normalization factor
          norm = aonorm(ityp,ishlpri(i),exppri(nm+1:nm+ishlpri(i)),cnorm(1:ishlpri(i)))

          ! calculate and assign the normalized primitive coefficients
          do k = 1, ishlpri(i)
             nn = nn + 1
             cpri(nn) = cnorm(k) * norm
          end do
       end do
       nm = nm + ishlpri(i)
    end do
    deallocate(cnorm)

    ! build the wavefunction coefficients for the primitives
    nn = 0
    nm = 0
    nl = 0
    do i = 1, ncshel
       do j = jshl0(ishlt(i)), jshl1(ishlt(i))
          ityp = typtrans(j)
          nl = nl + 1
          do k = 1, ishlpri(i)
             nn = nn + 1
             m%icenter(nn) = ishlat(i)
             m%itype(nn) = ityp
             m%e(nn) = exppri(nm+k)
             m%c(:,nn) = cpri(nn) * mocoef(:,nl)
          end do
       end do
       nm = nm + ishlpri(i)
    end do

    ! clean up and exit
    deallocate(ishlt,ishlpri,ishlat)
    deallocate(exppri,ccontr)
    deallocate(mocoef,cpri)

    return
40  continue
    call error("readmolden","unexpected end of file",2)

  contains

    function next_keyword()

      integer :: istart, iend
      logical :: next_keyword

      keyword = ""
      next_keyword = .false.

      do while(index(lower(line),"[") == 0)
         read(luwfn,'(A)',end=40,err=40) line
      end do
      next_keyword =.true.
      istart = index(lower(line),"[") + 1
      iend = index(lower(line),"]") - 1
      keyword = lower(line(istart:iend))
      return
40    continue
      return

    end function next_keyword

  end function readmolden

  !> Read terachem checkpoint file
  function readtck(file,egauss) result(m)

    character*(mline), intent(in) :: file
    real*8, intent(out) :: egauss
    type(molecule) :: m

    character*(mline) :: line, word, word2
    integer :: lp, idum, nalpha, nbeta, nshel, ncshel, lmax, nbas
    integer :: istat, i, j, k, l, ifac, ii, jj, itypmax
    integer :: acent, nn, nm, nl
    logical :: ok
    integer, allocatable :: ishlt(:), ishlpri(:), ishlat(:), npribas(:), nbaspri(:)
    real*8, allocatable :: xat(:), exppri(:), ccontr(:), mocoef(:)
    real*8 :: aexp, acoef, rdum, norm

    ! set title
    m%name = file
    m%useecp = .false.

    ! parse the tck file, first pass
    open(luwfn,file=file,status='old')
    read(luwfn,*,end=30) m%n
    read(luwfn,*,end=30) idum
    m%nelec = idum
    read(luwfn,*,end=30) idum
    m%charge = idum
    read(luwfn,*,end=30) m%mult
    read(luwfn,*,end=30) nbas
    read(luwfn,*,end=30) m%nmo
    read(luwfn,*,end=30) ncshel
    read(luwfn,*,end=30) nshel
    read(luwfn,'(/)',end=30)

    ! wfntyp -> closed-shell
    m%wfntyp = 0
    if (m%wfntyp == 0) then
       allocate(m%occ(m%nmo),stat=istat)
       if (istat /= 0) call error('readtck','could not allocate memory for occ',2)
       m%occ = 2d0
    endif

    ! allocates
    allocate(m%x(3,m%n),m%z(m%n),xat(3*m%n),stat=istat)
    if (istat /= 0) call error('readtck','could not allocate memory for geometry',2)
    allocate(ishlt(ncshel),ishlpri(ncshel),ishlat(ncshel),stat=istat)
    if (istat /= 0) call error('readtck','could not allocate memory for shell types',2)
    allocate(exppri(nshel),ccontr(nshel),stat=istat)
    if (istat /= 0) call error('readtck','could not allocate memory for prim. shells',2)
    allocate(mocoef(nbas*m%nmo),stat=istat)
    if (istat /= 0) call error('readtck','could not allocate memory for MO coefs',2)

    ! geometry
    do i = 1, m%n
       read(luwfn,*,end=30) word, rdum, m%x(:,i)
       m%z(i) = nint(rdum)
    enddo

    read(luwfn,'(/)',end=30)
    do i = 1,m%n
       read(luwfn,*,end=30)
       nn = 0
       do while(.true.)
          read(luwfn,'(A80)',end=30) word
          lp = 1
          ok = isinteger(idum,word,lp)
          if (ok) then
             nn = nn + 1
             read(word,*,end=30) ii, word2, jj, aexp, acoef
             if (lower(trim(word2)) == "s") then
                ishlt(ii) = 0
             else if (lower(trim(word2)) == "p") then
                ishlt(ii) = 1
             else if (lower(trim(word2)) == "sp") then
                call error("readtck","can't handle SP in gamess format",2)
             else if (lower(trim(word2)) == "d") then
                ishlt(ii) = 2
             else if (lower(trim(word2)) == "f") then
                ishlt(ii) = 3
             else
                call error("readtck","basis set type not supported",2)
             endif
             ishlat(ii) = i
             exppri(jj) = aexp
             ccontr(jj) = acoef
          else
             if (nn > 0) ishlpri(ii) = nn
             nn = 0
             if (trim(word) /= "") then
                exit
             endif
          endif
       enddo
    enddo

    ! energy
    read (luwfn,*) egauss

    ! advance to the MO 
    word = ""
    do while(word /= "Occupied")
       read(luwfn,*,end=30) word
    end do
    do i = 1, nbas
       read(luwfn,*,end=30) (mocoef(j*nbas+i),j=0,m%nmo-1)
    enddo

    ! Count the number of primitives
    m%npri = 0
    do i = 1, ncshel
       if (ishlt(i) == 0) then
          ifac = 1
       else if (ishlt(i) == 1) then
          ifac = 3
       else if (ishlt(i) == -1) then
          ifac = 4
          call error('readtck','ishlt = -1 not supported',2)
       else if (ishlt(i) == 2) then
          ifac = 6
       else if (ishlt(i) == 3) then
          ifac = 10
       endif
       m%npri = m%npri + ifac * ishlpri(i)
    enddo

    ! Normalize the primitive shells
    nm = 0
    do i = 1, ncshel
       l = ishlt(i)
       norm = 0d0
       do j = 1, ishlpri(i)
          do k = 1, ishlpri(i)
             norm = norm + ccontr(nm+j) * ccontr(nm+k) * &
                sqrt((2*exppri(nm+j))**(l+1.5d0)) * sqrt((2*exppri(nm+k))**(l+1.5d0)) / &
                (exppri(nm+j) + exppri(nm+k))**(l+1.5d0)
          end do
       end do
       norm = sqrt(norm)
       do j = 1, ishlpri(i)
          ccontr(nm+j) = ccontr(nm+j) / norm
       end do
       nm = nm + ishlpri(i)
    end do

    ! Assign primitive center and type, exponents, etc.
    allocate(m%icenter(m%npri),stat=istat)
    if (istat /= 0) call error('readtck','could not allocate memory for icenter',2)
    allocate(m%itype(m%npri),stat=istat)
    if (istat /= 0) call error('readtck','could not allocate memory for itype',2)
    allocate(m%e(m%npri),stat=istat)
    if (istat /= 0) call error('readtck','could not allocate memory for exponents',2)
    allocate(m%c(m%nmo,m%npri),stat=istat)
    if (istat /= 0) call error('readtck','could not allocate memory for coeffs',2)
    allocate(npribas(m%npri),nbaspri(nbas),stat=istat)
    if (istat /= 0) call error('readtck','could not allocate memory for pri/bas index',2)
    nn = 0
    nm = 0
    nl = 0
    do i = 1, ncshel
       acent = ishlat(i)
       if (ishlt(i) == 0) then
          nl = nl + 1
          do k = 1, ishlpri(i)
             nn = nn + 1
             npribas(nn) = nl
             nbaspri(nl) = nn
             m%icenter(nn) = acent
             m%itype(nn) = 1
             m%e(nn) = exppri(nm+k)
             m%c(1:m%nmo,nn) = gnorm(m%itype(nn),m%e(nn)) * ccontr(nm+k)
          end do
       else if (ishlt(i) == 1) then
          do j = 2, 4
             nl = nl + 1
             do k = 1, ishlpri(i)
                nn = nn + 1
                npribas(nn) = nl
                nbaspri(nl) = nn
                m%icenter(nn) = acent
                m%itype(nn) = j
                m%e(nn) = exppri(nm+k)
                m%c(1:m%nmo,nn) = gnorm(m%itype(nn),m%e(nn)) * ccontr(nm+k)
             end do
          end do
       else if (ishlt(i) == 2) then
          do j = 5, 10
             nl = nl + 1
             do k = 1, ishlpri(i)
                nn = nn + 1
                npribas(nn) = nl
                nbaspri(nl) = nn
                m%icenter(nn) = acent
                ! terachem sequence is xy,xz,yz,xx,yy,zz
                if (j <= 7) then
                   m%itype(nn) = j+3
                else
                   m%itype(nn) = j-3
                endif
                m%e(nn) = exppri(nm+k)
                m%c(1:m%nmo,nn) = gnorm(m%itype(nn),m%e(nn)) * ccontr(nm+k)
             end do
          end do
       else if (ishlt(i) == 3) then
          do j = 11, 20
             nl = nl + 1
             do k = 1, ishlpri(i)
                nn = nn + 1
                npribas(nn) = nl
                nbaspri(nl) = nn
                m%icenter(nn) = acent
                m%itype(nn) = j
                m%e(nn) = exppri(nm+k)
                m%c(1:m%nmo,nn) = gnorm(m%itype(nn),m%e(nn)) * ccontr(nm+k)
             end do
          end do
       endif
       nm = nm + ishlpri(i)
    end do

    ! re-order MO coefficients in input
    itypmax = maxval(m%itype)
    nl = 0
    do l = 0, 3
       do i = 1, nbas
          ok = (l==0) .and. (m%itype(nbaspri(i))==1)
          ok = ok .or. (l==1) .and. (m%itype(nbaspri(i))>=2.and.m%itype(nbaspri(i))<=4)
          ok = ok .or. (l==2) .and. (m%itype(nbaspri(i))>=5.and.m%itype(nbaspri(i))<=10)
          ok = ok .or. (l==3) .and. (m%itype(nbaspri(i))>=11.and.m%itype(nbaspri(i))<=20)
          if (ok) then
             nl = nl + 1
             do k = 1, m%npri
                if (npribas(k) /= i) cycle
                do j = 1,m%nmo
                   m%c(j,k) = m%c(j,k) * mocoef((j-1)*nbas+nl)
                end do
             end do
          endif
          nn = nn + npribas(i)
       end do
    enddo

    deallocate(ishlt,ishlpri,ishlat)
    deallocate(exppri,ccontr)
    deallocate(mocoef,npribas,nbaspri)
    close(luwfn)

    return
30  continue
    call error("readtck","unexpected end of file",2)

  end function readtck

  !> Calculate the normalization factor of a primitive of a given type
  !> with exponent a.
  !> [ (4a)^(i+j+k) * (2a/pi)^(3/2) / (2i-1)!! / (2j-1)!! / (2k-1)!! ]^(1/2)
  !> (Note: (2l)! = (2l-1)!! * 2^l * l!.)
  function gnorm(type,a) result(N)
    use param, only: pi
    integer, intent(in) :: type
    real*8, intent(in) :: a
    real*8 :: N

    if (type == 1) then
       ! 1
       ! 1
       N = 2**(3d0/4d0) * a**(3d0/4d0) / pi**(3d0/4d0)
    else if (type >= 2 .and. type <= 4) then
       ! 2 3 4
       ! x y z
       N = 2**(7d0/4d0) * a**(5d0/4d0) / pi**(3d0/4d0)
    else if (type >= 5 .and. type <= 7) then
       ! 5  6  7
       ! xx yy zz
       N = 2**(11d0/4d0) * a**(7d0/4d0) / pi**(3d0/4d0) / sqrt(3d0)
    else if (type >= 8 .and. type <= 10) then
       ! 7  8  9
       ! xy xz yz
       N = 2**(11d0/4d0) * a**(7d0/4d0) / pi**(3d0/4d0) 
    else if (type >= 11 .and. type <= 13) then
       ! 11  12  13
       ! xxx yyy zzz
       N = 2**(15d0/4d0) * a**(9d0/4d0) / pi**(3d0/4d0) / sqrt(15d0)
    else if (type >= 14 .and. type <= 19) then
       ! 14  15  16  17  18  19  
       ! xyy xxy xxz xzz yzz yyz 
       N = 2**(15d0/4d0) * a**(9d0/4d0) / pi**(3d0/4d0) / sqrt(3d0)
    else if (type == 20) then
       ! 20
       ! xyz
       N = 2**(15d0/4d0) * a**(9d0/4d0) / pi**(3d0/4d0)
    else if (type >= 21 .and. type <= 23) then
       ! 21   22   23   
       ! xxxx yyyy zzzz 
       N = 2**(19d0/4d0) * a**(11d0/4d0) / pi**(3d0/4d0) / sqrt(105d0)
    else if (type >= 24 .and. type <= 29) then
       ! 24   25   26   27   28   29  
       ! xxxy xxxz xyyy yyyz xzzz yzzz
       N = 2**(19d0/4d0) * a**(11d0/4d0) / pi**(3d0/4d0) / sqrt(15d0)
    else if (type >= 30 .and. type <= 32) then
       ! 30   31   32  
       ! xxyy xxzz yyzz
       N = 2**(19d0/4d0) * a**(11d0/4d0) / pi**(3d0/4d0) / sqrt(9d0)
    else if (type >= 33 .and. type <= 35) then
       ! 33   34   35   
       ! xxyz xyyz xyzz 
       N = 2**(19d0/4d0) * a**(11d0/4d0) / pi**(3d0/4d0) / sqrt(3d0)
    else
       call error("gnorm","fixme: primitive type not supported",2)
    endif

  endfunction gnorm

  !> Calculate the normalization factor of a primitive of a given
  !> angular momentum l with exponent a, orca convention. This is used
  !> to "de-normalize" the orca primitive coefficients.
  function gnorm_orca(type,a) result(N)
    use param, only: pi
    integer, intent(in) :: type
    real*8, intent(in) :: a
    real*8 :: N

    if (type == 1) then
       N = 2**(3d0/4d0) * a**(3d0/4d0) / pi**(3d0/4d0)
    else if (type >= 2 .and. type <= 4) then
       N = 2**(7d0/4d0) * a**(5d0/4d0) / pi**(3d0/4d0)
    else if (type >= 5 .and. type <= 10) then
       N = 2**(11d0/4d0) * a**(7d0/4d0) / pi**(3d0/4d0) / sqrt(3d0)
    else if (type >= 11 .and. type <= 20) then
       N = 2**(15d0/4d0) * a**(9d0/4d0) / pi**(3d0/4d0) / sqrt(15d0)
    else if (type >= 21 .and. type <= 35) then
       N = 2**(19d0/4d0) * a**(11d0/4d0) / pi**(3d0/4d0) / sqrt(105d0)
    else
       call error("gnorm","fixme: primitive type not supported",2)
    endif

  endfunction gnorm_orca

  !> Calculate the normalization factor of an atomic orbital composed
  !> of n primitives of type type. The exponents are a and the
  !> coefficients are c.
  function aonorm(type,n,a,c)
    use param, only: pi
    integer, intent(in) :: type, n
    real*8, intent(in) :: a(n), c(n)
    real*8 :: aonorm

    real*8 :: dfac
    integer :: l
    integer :: i, j

    if (type == 1) then
       ! 1
       ! 1
       dfac = 1d0
       l = 0
    else if (type >= 2 .and. type <= 4) then
       ! 2 3 4
       ! x y z
       dfac = 1d0
       l = 1
    else if (type >= 5 .and. type <= 7) then
       ! 5  6  7
       ! xx yy zz
       dfac = 3d0
       l = 2
    else if (type >= 8 .and. type <= 10) then
       ! 7  8  9
       ! xy xz yz
       dfac = 1d0
       l = 2
    else if (type >= 11 .and. type <= 13) then
       ! 11  12  13
       ! xxx yyy zzz
       dfac = 15d0
       l = 3
    else if (type >= 14 .and. type <= 19) then
       ! 14  15  16  17  18  19  
       ! xyy xxy xxz xzz yzz yyz 
       dfac = 3d0
       l = 3
    else if (type == 20) then
       ! 20
       ! xyz
       dfac = 1d0
       l = 3
    else if (type >= 21 .and. type <= 23) then
       ! 21   22   23   
       ! xxxx yyyy zzzz 
       dfac = 105d0
       l = 4
    else if (type >= 24 .and. type <= 29) then
       ! 24   25   26   27   28   29  
       ! xxxy xxxz xyyy yyyz xzzz yzzz
       dfac = 15d0
       l = 4
    else if (type >= 30 .and. type <= 32) then
       ! 30   31   32  
       ! xxyy xxzz yyzz
       dfac = 9d0
       l = 4
    else if (type >= 33 .and. type <= 35) then
       ! 33   34   35   
       ! xxyz xyyz xyzz 
       dfac = 3d0
       l = 4
    else
       call error("aonorm","fixme: primitive type not supported",2)
    endif

    aonorm = 0d0
    do i = 1, n
       do j = 1, n
          aonorm = aonorm + c(i) * c(j) / (a(i)+a(j))**(l+3d0/2d0)
       end do
    end do
    aonorm = aonorm * dfac * pi**(3d0/2d0) / 2**l
    aonorm = 1d0 / sqrt(aonorm)

  endfunction aonorm

  subroutine evalwfn(m,mesh)
    use param
    !     
    !     HARTMUT SCHMIDER March 2005
    !     Produces various properties on a Gaussian type cube grid
    !
    !     Adapted by Erin R. Johnson and Axel D. Becke, March 2005
    !     Uses grid generated by numol and returns
    !     properties and wavefunctions at all grid points
    !
    type(molecule), intent(inout) :: m
    type(tmesh), intent(inout) :: mesh

    integer :: i, istat, nmo1, nmo2
    integer :: inuc
    real*8 :: quads, dsigs, hirsh(mesh%n), r, r1, r2

    istat = 0
    if (.not.allocated(mesh%rho)) allocate(mesh%rho(mesh%n,2),stat=istat)
    if (istat /= 0) call error('evalwfn','could not allocate memory for rho',2)
    if (.not.allocated(mesh%b)) allocate(mesh%b(mesh%n,2),stat=istat)
    if (istat /= 0) call error('evalwfn','could not allocate memory for tau',2)
    
    ! 0 - closed, 1 - open, 2 - restricted open, 3 - fractional
    call propts(m,mesh%n,mesh%x,mesh%rho,mesh%b)

    if (.not.allocated(m%mm)) allocate(m%mm(3,m%n),stat=istat)
    if (istat /= 0) call error('evalwfn','could not allocate memory for moments',2)
    if (.not.allocated(m%v)) allocate(m%v(m%n),stat=istat)
    if (istat /= 0) call error('evalwfn','could not allocate memory for volumes',2)
    if (.not.allocated(m%q)) allocate(m%q(m%n),stat=istat)
    if (istat /= 0) call error('evalwfn','could not allocate memory for charges',2)
    m%mm = 0d0
    m%v = 0d0
    m%q = 0d0

    rewind(ihrsh)
    do inuc = 1, m%n
       if (m%z(inuc) < 1) cycle
       read(ihrsh)(hirsh(i),i=1,mesh%n)
       m%mm(:,inuc) = 0d0
       m%v(inuc) = 0d0
       m%q(inuc) = 0d0

       ! calculate hole dipole and moments
       do i = 1, mesh%n
          r = sqrt(dot_product(mesh%x(:,i)-m%x(:,inuc),mesh%x(:,i)-m%x(:,inuc)))
          r1 = max(0.d0,r-mesh%b(i,1))
          r2 = max(0.d0,r-mesh%b(i,2))

          m%mm(1,inuc)=m%mm(1,inuc) + mesh%w(i) * hirsh(i) * &
             (mesh%rho(i,1) * (r-r1)**2 + mesh%rho(i,2) * (r-r2)**2)
          m%mm(2,inuc)=m%mm(2,inuc) + mesh%w(i) * hirsh(i) * &
             (mesh%rho(i,1) * (r**2-r1**2)**2 + mesh%rho(i,2) * (r**2-r2**2)**2)
          m%mm(3,inuc)=m%mm(3,inuc) + mesh%w(i) * hirsh(i) * &
             (mesh%rho(i,1) * (r**3-r1**3)**2 + mesh%rho(i,2) * (r**3-r2**3)**2)
          m%v(inuc) = m%v(inuc) + mesh%w(i) * hirsh(i) * &
             (mesh%rho(i,1)+mesh%rho(i,2)) * r**3
          m%q(inuc) = m%q(inuc) + mesh%w(i) * hirsh(i) * (mesh%rho(i,1)+mesh%rho(i,2))
       enddo
    enddo

  end subroutine evalwfn

 subroutine propts(m,nr,r,rho,b)
   
   type(molecule), intent(in) :: m
   integer, intent(in) :: nr
   real*8, intent(in) :: r(3,nr)
   real*8, intent(out) :: rho(nr,2), b(nr,2)

   integer :: i, j, nn, ityp, ipri, iat, ipria, ix, l(3)
   integer :: imo, nspin, n0(2), n1(2), nmo1
   real*8 :: al, x0(3), ex, xl(3,0:2), xl2
   real*8 :: chi(m%npri,10), maxc(m%npri), dd(3,m%n), d2(m%n)
   real*8 :: phi(m%nmo,10), gg(3), hh(3), quads, drho2, d2rho, taup
   real*8 :: dsigs, aocc, prho(2), pb(2)
   logical :: ldopri(m%npri,10)
   
   real*8, parameter :: cutoff_pri = 1d-15
   real*8, parameter :: small = 1d-10

   integer, parameter :: li(3,56) = reshape((/&
      0,0,0, & ! s
      1,0,0, 0,1,0, 0,0,1, & ! p
      2,0,0, 0,2,0, 0,0,2, 1,1,0, 1,0,1, 0,1,1, & !d
      3,0,0, 0,3,0, 0,0,3, 2,1,0, 2,0,1, 0,2,1, &
      1,2,0, 1,0,2, 0,1,2, 1,1,1,& ! f
      4,0,0, 0,4,0, 0,0,4, 3,1,0, 3,0,1, 1,3,0, 0,3,1, 1,0,3,&
      0,1,3, 2,2,0, 2,0,2, 0,2,2, 2,1,1, 1,2,1, 1,1,2,& ! g
      0,0,5, 0,1,4, 0,2,3, 0,3,2, 0,4,1, 0,5,0, 1,0,4, 1,1,3,&
      1,2,2, 1,3,1, 1,4,0, 2,0,3, 2,1,2, 2,2,1, 2,3,0, 3,0,2,&
      3,1,1, 3,2,0, 4,0,1, 4,1,0, 5,0,0/),shape(li)) ! h

   ! identify the max coefficients
   maxc = 0d0
   do imo = 1, m%nmo
      do ipri = 1, m%npri
         maxc(ipri) = max(maxc(ipri),abs(m%c(imo,ipri)))
      enddo
   enddo

   !$omp parallel do private(dd,d2,ityp,iat,al,ex,l,xl,chi,ldopri,phi,prho,pb,&
   !$omp taup,gg,hh,aocc,drho2,d2rho,dsigs,quads,nmo1) schedule(dynamic)
   do i = 1, nr
      ! calculate distances
      do iat = 1, m%n
         dd(:,iat) = r(:,i) - m%x(:,iat)
         d2(iat) = dd(1,iat)*dd(1,iat)+dd(2,iat)*dd(2,iat)+dd(3,iat)*dd(3,iat)
      enddo

      do ipri = 1, m%npri
         ityp = m%itype(ipri)
         iat = m%icenter(ipri)
         al = m%e(ipri)
         ex = exp(-al * d2(iat))

         l = li(1:3,ityp)
         do ix = 1, 3
            if (l(ix) == 0) then
               xl(ix,0) = 1d0
               xl(ix,1) = 0d0
               xl(ix,2) = 0d0
            else if (l(ix) == 1) then
               xl(ix,0) = dd(ix,iat)
               xl(ix,1) = 1d0
               xl(ix,2) = 0d0
            else if (l(ix) == 2) then
               xl(ix,0) = dd(ix,iat) * dd(ix,iat)
               xl(ix,1) = 2d0 * dd(ix,iat)
               xl(ix,2) = 2d0
            else if (l(ix) == 3) then
               xl(ix,0) = dd(ix,iat) * dd(ix,iat) * dd(ix,iat)
               xl(ix,1) = 3d0 * dd(ix,iat) * dd(ix,iat)
               xl(ix,2) = 6d0 * dd(ix,iat)
            else if (l(ix) == 4) then
               xl2 = dd(ix,iat) * dd(ix,iat)
               xl(ix,0) = xl2 * xl2
               xl(ix,1) = 4d0 * xl2 * dd(ix,iat)
               xl(ix,2) = 12d0 * xl2
            else if (l(ix) == 5) then
               xl2 = dd(ix,iat) * dd(ix,iat)
               xl(ix,0) = xl2 * xl2 * dd(ix,iat)
               xl(ix,1) = 5d0 * xl2 * xl2
               xl(ix,2) = 20d0 * xl2 * dd(ix,iat)
            else
               call error('pri012','power of L not supported',2)
            end if
         end do

         chi(ipri,1) = xl(1,0)*xl(2,0)*xl(3,0)*ex
         chi(ipri,2) = (xl(1,1)-2*al*dd(1,iat)**(l(1)+1))*xl(2,0)*xl(3,0)*ex
         chi(ipri,3) = (xl(2,1)-2*al*dd(2,iat)**(l(2)+1))*xl(1,0)*xl(3,0)*ex
         chi(ipri,4) = (xl(3,1)-2*al*dd(3,iat)**(l(3)+1))*xl(1,0)*xl(2,0)*ex
         chi(ipri,5) = (xl(1,2)-2*al*(2*l(1)+1)*xl(1,0)&
            +4*al*al*dd(1,iat)**(l(1)+2))*xl(2,0)*xl(3,0)*ex
         chi(ipri,6) = (xl(2,2)-2*al*(2*l(2)+1)*xl(2,0)&
            +4*al*al*dd(2,iat)**(l(2)+2))*xl(3,0)*xl(1,0)*ex
         chi(ipri,7) = (xl(3,2)-2*al*(2*l(3)+1)*xl(3,0)&
            +4*al*al*dd(3,iat)**(l(3)+2))*xl(1,0)*xl(2,0)*ex
         chi(ipri,8) = (xl(1,1)-2*al*dd(1,iat)**(l(1)+1))*&
            (xl(2,1)-2*al*dd(2,iat)**(l(2)+1))*xl(3,0)*ex
         chi(ipri,9) = (xl(1,1)-2*al*dd(1,iat)**(l(1)+1))*&
            (xl(3,1)-2*al*dd(3,iat)**(l(3)+1))*xl(2,0)*ex
         chi(ipri,10)= (xl(3,1)-2*al*dd(3,iat)**(l(3)+1))*&
            (xl(2,1)-2*al*dd(2,iat)**(l(2)+1))*xl(1,0)*ex

         do ix = 1, 10
            ldopri(ipri,ix) = (abs(chi(ipri,ix))*maxc(ipri) > cutoff_pri)
         enddo
      enddo ! ipri = 1, npri

      ! build the MO avlues at the point
      phi = 0d0
      do ix = 1, 10
         do ipri = 1, m%npri
            if (.not.ldopri(ipri,ix)) cycle
            do imo = 1, m%nmo
               phi(imo,ix) = phi(imo,ix) + m%c(imo,ipri)*chi(ipri,ix)
            enddo
         enddo
      enddo

      ! contribution to the density, etc.
      prho = 0d0
      pb = 0d0
      taup = 0d0
      gg = 0d0
      hh = 0d0
      if (m%wfntyp == 0) then
         do imo = 1, m%nmo
            aocc = m%occ(imo) * 0.5d0
            prho(1) = prho(1) + aocc * phi(imo,1) * phi(imo,1)
            gg = gg + 2 * aocc * phi(imo,1) * phi(imo,2:4) 
            hh = hh + 2 * aocc * (phi(imo,1)*phi(imo,5:7)+phi(imo,2:4)**2)
            taup = taup + aocc * (phi(imo,2)*phi(imo,2)+phi(imo,3)*phi(imo,3)+phi(imo,4)*phi(imo,4))
         enddo
         prho(2) = prho(1)
         if (prho(1) > small) then
            drho2 = gg(1)*gg(1)+gg(2)*gg(2)+gg(3)*gg(3)
            d2rho = hh(1)+hh(2)+hh(3)
            dsigs = taup - 0.25d0 * drho2 / max(prho(1),1d-30)
            quads = (d2rho - 2d0 * dsigs) / 6d0
            call bhole(prho(1),quads,1d0,pb(1))
            pb(2) = pb(1)
         endif
      else if (m%wfntyp == 1) then
         nmo1 = (m%nmo + m%mult - 1)/2
         do imo = 1, nmo1
            aocc = m%occ(imo)
            prho(1) = prho(1) + aocc * phi(imo,1) * phi(imo,1)
            gg = gg + 2 * aocc * phi(imo,1) * phi(imo,2:4) 
            hh = hh + 2 * aocc * (phi(imo,1)*phi(imo,5:7)+phi(imo,2:4)**2)
            taup = taup + aocc * (phi(imo,2)*phi(imo,2)+phi(imo,3)*phi(imo,3)+phi(imo,4)*phi(imo,4))
         enddo
         if (prho(1) > small) then
            drho2 = gg(1)*gg(1)+gg(2)*gg(2)+gg(3)*gg(3)
            d2rho = hh(1)+hh(2)+hh(3)
            dsigs = taup - 0.25d0 * drho2 / max(prho(1),1d-30)
            quads = (d2rho - 2d0 * dsigs) / 6d0
            call bhole(prho(1),quads,1d0,pb(1))
         endif
         taup = 0d0
         gg = 0d0
         hh = 0d0
         do imo = nmo1+1, m%nmo
            aocc = m%occ(imo)
            prho(2) = prho(2) + aocc * phi(imo,1) * phi(imo,1)
            gg = gg + 2 * aocc * phi(imo,1) * phi(imo,2:4) 
            hh = hh + 2 * aocc * (phi(imo,1)*phi(imo,5:7)+phi(imo,2:4)**2)
            taup = taup + aocc * (phi(imo,2)*phi(imo,2)+phi(imo,3)*phi(imo,3)+phi(imo,4)*phi(imo,4))
         enddo
         if (prho(2) > small) then
            drho2 = gg(1)*gg(1)+gg(2)*gg(2)+gg(3)*gg(3)
            d2rho = hh(1)+hh(2)+hh(3)
            dsigs = taup - 0.25d0 * drho2 / max(prho(2),1d-30)
            quads = (d2rho - 2d0 * dsigs) / 6d0
            call bhole(prho(2),quads,1d0,pb(2))
         endif
      else if (m%wfntyp == 2) then
         nmo1 = m%nmo - m%mult + 1
         do imo = 1, nmo1
            aocc = m%occ(imo) * 0.5d0
            prho(2) = prho(2) + aocc * phi(imo,1) * phi(imo,1)
            gg = gg + 2 * aocc * phi(imo,1) * phi(imo,2:4) 
            hh = hh + 2 * aocc * (phi(imo,1)*phi(imo,5:7)+phi(imo,2:4)**2)
            taup = taup + aocc * (phi(imo,2)*phi(imo,2)+phi(imo,3)*phi(imo,3)+phi(imo,4)*phi(imo,4))
         enddo
         if (prho(2) > small) then
            drho2 = gg(1)*gg(1)+gg(2)*gg(2)+gg(3)*gg(3)
            d2rho = hh(1)+hh(2)+hh(3)
            dsigs = taup - 0.25d0 * drho2 / max(prho(2),1d-30)
            quads = (d2rho - 2d0 * dsigs) / 6d0
            call bhole(prho(2),quads,1d0,pb(2))
         endif
         prho(1) = prho(2)
         do imo = nmo1+1, m%nmo
            aocc = m%occ(imo)
            prho(1) = prho(1) + aocc * phi(imo,1) * phi(imo,1)
            gg = gg + 2 * aocc * phi(imo,1) * phi(imo,2:4) 
            hh = hh + 2 * aocc * (phi(imo,1)*phi(imo,5:7)+phi(imo,2:4)**2)
            taup = taup + aocc * (phi(imo,2)*phi(imo,2)+phi(imo,3)*phi(imo,3)+phi(imo,4)*phi(imo,4))
         enddo
         if (prho(1) > small) then
            drho2 = gg(1)*gg(1)+gg(2)*gg(2)+gg(3)*gg(3)
            d2rho = hh(1)+hh(2)+hh(3)
            dsigs = taup - 0.25d0 * drho2 / max(prho(1),1d-30)
            quads = (d2rho - 2d0 * dsigs) / 6d0
            call bhole(prho(1),quads,1d0,pb(1))
         endif
      else
         call error("evalwfn","wfn type not implemented",2)
      endif
      !$omp critical (write)
      rho(i,:) = prho(:)
      b(i,:) = pb(:)
      !$omp end critical (write)
   enddo ! i = 1, nr
   !$omp end parallel do

 end subroutine propts

 subroutine bhole(rho,quad,hnorm,b)
   use param

   real*8, intent(in) :: rho, quad, hnorm
   real*8, intent(out) :: b 

   real*8 :: rhs, x0, shift, x, x1, expo, prefac, alf, f, df
   integer :: i

   rhs=third2*(pi*rho/hnorm)**third2*rho/quad
   x0=2.d0
   shift=1.d0
   if(rhs.lt.0.d0)go to 10
   if(rhs.gt.0.d0)go to 20
10 do i=1,16
      x=x0-shift
      call xfuncs(x,rhs,f,df)
      if(f.lt.0.d0)go to 88
      shift=0.1d0*shift
   enddo
   write(iout,1002)
   stop
20 do i=1,16
      x=x0+shift
      call xfuncs(x,rhs,f,df)
      if(f.gt.0.d0)go to 88
      shift=0.1d0*shift
   enddo
   write(iout,1002)
   stop
88 continue
   do i=1,100
      call xfuncs(x,rhs,f,df)
      x1=x-f/df
      if(dabs(x1-x).lt.1.d-10)go to 111
      x=x1
   enddo
   write(iout,1001)
   stop
111 x=x1
   expo=dexp(-x)
   prefac=rho/expo
   alf=(8.d0*pi*prefac/hnorm)**third
   b=x/alf
   return
1001 format(' ','bhole: newton algorithm fails to converge!')
1002 format(' ','bhole: newton algorithm fails to initialize!')
 end subroutine bhole


 subroutine xfuncs(x,rhs,f,df)
   real*8, intent(in) :: x, rhs
   real*8, intent(out) :: f, df

   real*8 :: expo23

   expo23=dexp(-2.d0/3.d0*x)
   f = x*expo23/(x-2.d0) - rhs
   df=2.d0/3.d0*(2.d0*x-x**2-3.d0)/(x-2.d0)**2*expo23
 end subroutine xfuncs

 subroutine edisp(m,a1,a2,egauss,usec9)
   type(molecule), intent(in) :: m
   real*8, intent(in) :: a1, a2, egauss
   logical, intent(in) :: usec9

   integer :: i, j, k, k1, k2
   real*8 :: d, atpol(m%n), fac, rvdw, c6, c8, c10, rc
   real*8 :: c6com, c8com, c10com, xij(3), ifac
   real*8 :: e, f(3,m%n), q(3,m%n,3,m%n), qfac
   real*8 :: c9, qi, qj, qk

   do i = 1, m%n
      if (m%z(i) < 1) cycle
      atpol(i) = m%v(i) * frepol(m%z(i)) / frevol(m%z(i))
      ! write (iout,'("MOMENT",X,1p,2(E18.10,X))') m%mm(1,i), atpol(i)
   enddo

   write (iout,'("coefficients and distances (a.u.)")')
   write (iout,'("# i  j       dij            C6               C8               C10              Rc           Rvdw")') 
   e = 0d0
   f = 0d0
   q = 0d0
   do i = 1, m%n
      if (m%z(i) < 1) cycle
      do j = i, m%n
         if (m%z(j) < 1) cycle
         xij = m%x(:,j)-m%x(:,i)
         d = sqrt(dot_product(xij,xij))
         fac = atpol(i)*atpol(j)/(m%mm(1,i)*atpol(j)+m%mm(1,j)*atpol(i))
         c6 = fac*m%mm(1,i)*m%mm(1,j)
         c8 = 1.5d0*fac*(m%mm(1,i)*m%mm(2,j)+m%mm(2,i)*m%mm(1,j))
         c10 = 2.d0*fac*(m%mm(1,i)*m%mm(3,j)+m%mm(3,i)*m%mm(1,j))&
            +4.2d0*fac*m%mm(2,i)*m%mm(2,j)
         rc = (sqrt(c8/c6) + sqrt(sqrt(c10/c6)) +&
            sqrt(c10/c8)) / 3.D0
         rvdw = a1 * rc + a2
         if (d > 1d-5) then
            e = e - c6 / (rvdw**6 + d**6) - c8 / (rvdw**8+d**8) - &
               c10 / (rvdw**10 + d**10)
            c6com = 6.d0*c6*d**4/(rvdw**6+d**6)**2
            c8com = 8.d0*c8*d**6/(rvdw**8+d**8)**2
            c10com = 10.d0*c10*d**8/(rvdw**10+d**10)**2
            f(:,i) = f(:,i) + (c6com+c8com+c10com) * xij
            f(:,j) = f(:,j) - (c6com+c8com+c10com) * xij
            do k1 = 1, 3
               do k2 = 1, 3
                  if (k1 == k2) then
                     ifac = 1d0
                  else
                     ifac = 0d0
                  endif
                  qfac = &
                     c6com  * (-ifac - 4*xij(k1)*xij(k2)/d**2 + 12*xij(k1)*xij(k2)*d**4/(rvdw**6+d**6)) + &
                     c8com  * (-ifac - 6*xij(k1)*xij(k2)/d**2 + 16*xij(k1)*xij(k2)*d**6/(rvdw**8+d**8)) + &
                     c10com * (-ifac - 8*xij(k1)*xij(k2)/d**2 + 20*xij(k1)*xij(k2)*d**8/(rvdw**10+d**10)) 
                  q(k1,i,k2,j) = qfac
                  q(k2,j,k1,i) = qfac
               enddo
            enddo
         endif
         write (iout,'(I3,X,I3,1p,E14.6,X,3(E16.9,X),2(E13.6,X))') &
            i, j, d, c6, c8, c10, rc, rvdw
      end do
   end do
   write (iout,'("#")')

   if (usec9) then
      write (iout,'("three-body dispersion coefficients (a.u.)")')
      write (iout,'("# i  j  k            C9")')
      do i = 1, m%n
         if (m%z(i) < 1) cycle
         do j = i, m%n
            if (m%z(j) < 1) cycle
            do k = j, m%n
               if (m%z(k) < 1) cycle

               qi = m%mm(1,i) / atpol(i)
               qj = m%mm(1,j) / atpol(j)
               qk = m%mm(1,k) / atpol(k)

               c9 = m%mm(1,i)*m%mm(1,j)*m%mm(1,k)* (qi+qj+qk) / ((qi+qj)*(qi+qk)*(qj+qk))
               write (iout,'(3(I3,X)1p,E16.9)') i, j, k, c9

            end do
         end do
      end do
      write (iout,'("#")')
   end if

   ! sum rules for the second derivatives
   do i = 1, m%n
      do k1 = 1, 3
         do k2 = 1, 3
            q(k1,i,k2,i) = 0d0
            do j = 1, m%n
               if (j == i) cycle
               q(k1,i,k2,i) = q(k1,i,k2,i) - q(k1,i,k2,j)
            enddo
         enddo
      enddo
   enddo

   write (iout,'("dispersion energy ",1p,E20.12)') e
   write (iout,'("scf energy ",1p,E20.12)') egauss
   write (iout,'("total energy (SCF+XDM) ",1p,E20.12)') egauss+e
   write (iout,'("dispersion forces ")')
   write (iout,'("# i          Fx                   Fy                   Fz")')
   do i = 1, m%n
      write (iout,'(I3,X,1p,3(E20.12,X))') i, f(:,i)
   enddo
   write (iout,'("#")')
   write (iout,'("dispersion force constant matrix ")')
   write (iout,'("# i  xyz   j   xyz    Exixj ")')
   do i = 1, m%n
      do k1 = 1, 3
         do j = 1, i-1
            do k2 = 1, 3
               write (iout,'(4(I3,X),1p,E20.12)') i, k1, j, k2, q(k1,i,k2,j)
            enddo
         enddo
         do k2 = 1, k1
            write (iout,'(4(I3,X),1p,E20.12)') i, k1, i, k2, q(k1,i,k2,i)
         enddo
      enddo
   enddo
   write (iout,'("#"/)')
   
 end subroutine edisp

 function read_integers(lu,n) result(x)
   integer, intent(in) :: lu, n
   integer :: x(n)

   integer :: kk, lp, idum
   character*(mline) :: line

   kk = 0
   lp = 1
   read(lu,'(A)',end=999) line
   do while(.true.)
      if (.not.isinteger(idum,line,lp)) then
         line = adjustl(line(lp:))
         if (line(1:1) == "<") exit
         lp = 1
         read(lu,'(A)',end=999) line
         line = adjustl(line)
         if (line(1:1) == "<") exit
      else
         kk = kk + 1
         if (kk > n) call error("read_integers","exceeded size of the array",2)
         x(kk) = idum
      endif
   enddo
   
   return
999 call error("read_integers","unexpected end of file",2)

 endfunction read_integers

 function read_reals1(lu,n) result(x)
   integer, intent(in) :: lu, n
   real*8 :: x(n)

   integer :: kk, lp
   real*8 :: rdum
   character*(mline) :: line

   kk = 0
   lp = 1
   read(lu,'(A)',end=999) line
   do while(.true.)
      if (.not.isreal(rdum,line,lp)) then
         line = adjustl(line(lp:))
         if (line(1:1) == "<") exit
         lp = 1
         read(lu,'(A)',end=999) line
         line = adjustl(line)
         if (line(1:1) == "<") exit
      else
         kk = kk + 1
         if (kk > n) call error("read_reals1","exceeded size of the array",2)
         x(kk) = rdum
      endif
   enddo
   
   return
999 call error("read_reals1","unexpected end of file",2)

 endfunction read_reals1

end module wfnmod
