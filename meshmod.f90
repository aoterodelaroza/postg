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
module meshmod
  use param

  private
  public :: genmesh

contains

  function genmesh(m) result(mesh)
    use tools_math
    implicit none
    ! MESH POINTS, NUCLEAR WEIGHTS, AND INTEGRATION WEIGHTS.

    type(molecule), intent(in) :: m
    type(tmesh) :: mesh

    real*8 :: rr(m%n,m%n), rmid, r, r1, r2, hypr, vp0, vpsum, vpi
    integer :: i, j, k, kk
    real*8, allocatable :: rads(:), wrads(:), xang(:), yang(:), zang(:), wang(:)
    integer :: nr, nang, ir, il, istat, mang, mr
    real*8 :: cutoff(m%n,m%n), x(3)
    real*8, allocatable :: meshrl(:,:,:), meshx(:,:,:,:)

    ! interatomic distances
    rr = 0d0
    do i = 1, m%n
       do j = i+1, m%n
          rr(i,j) = sqrt((m%x(1,i)-m%x(1,j))**2+(m%x(2,i)-m%x(2,j))**2+(m%x(3,i)-m%x(3,j))**2)
          rr(j,i) = rr(i,j)
       enddo
    enddo

    ! allocate space for the mesh
    mesh%n = 0
    do i = 1, m%n
       if (m%z(i) < 1) cycle
       mesh%n = mesh%n + z2nr(m%z(i)) * z2nang(m%z(i))
    enddo
    allocate(mesh%w(mesh%n),mesh%x(3,mesh%n),stat=istat)

    ! allocate work arrays
    mr = -1
    mang = -1
    do i = 1, m%n
       if (m%z(i) < 1) cycle
       mang = max(mang,z2nang(m%z(i)))
       mr = max(mr,z2nr(m%z(i)))
    end do
    allocate(meshrl(mang,mr,m%n),meshx(3,mang,mr,m%n))
    allocate(rads(mr),wrads(mr),stat=istat)
    if (istat /= 0) call error('genmesh','could not allocate memory for radial meshes',2)
    allocate(xang(mang),yang(mang),zang(mang),wang(mang),stat=istat)
    if (istat /= 0) call error('readwfn','could not allocate memory for angular meshes',2)
    
    ! Precompute the mesh weights with multiple threads. The job has to be
    ! split in two because the nodes have to be positioned in the array in 
    ! the correct order 
    !$omp parallel do private(nr,nang,rmid,ir,r,il,x,j,k,r1,r2,hypr,&
    !$omp cutoff,vp0,vpsum,vpi) firstprivate(rads,wrads,xang,yang,zang,wang)
    do i = 1, m%n
       if (m%z(i) < 1) cycle

       ! radial mesh
       nr = z2nr(m%z(i))
       nang = z2nang(m%z(i))
       rmid = 1d0/real(m%z(i),8)**third
       call rmesh(nr,rmid,rads,wrads)

       ! angular mesh

       if (nang == 6) then
          call ld0006(xang,yang,zang,wang,nang)
       else if (nang == 14) then
          call ld0014(xang,yang,zang,wang,nang)
       else if (nang == 26) then
          call ld0026(xang,yang,zang,wang,nang)
       else if (nang == 38) then
          call ld0038(xang,yang,zang,wang,nang)
       else if (nang == 50) then
          call ld0050(xang,yang,zang,wang,nang)
       else if (nang == 74) then
          call ld0074(xang,yang,zang,wang,nang)
       else if (nang == 86) then
          call ld0086(xang,yang,zang,wang,nang)
       else if (nang == 110) then
          call ld0110(xang,yang,zang,wang,nang)
       else if (nang == 146) then
          call ld0146(xang,yang,zang,wang,nang)
       else if (nang == 170) then
          call ld0170(xang,yang,zang,wang,nang)
       else if (nang == 194) then
          call ld0194(xang,yang,zang,wang,nang)
       else if (nang == 230) then
          call ld0230(xang,yang,zang,wang,nang)
       else if (nang == 266) then
          call ld0266(xang,yang,zang,wang,nang)
       else if (nang == 302) then
          call ld0302(xang,yang,zang,wang,nang)
       else if (nang == 350) then
          call ld0350(xang,yang,zang,wang,nang)
       else if (nang == 434) then
          call ld0434(xang,yang,zang,wang,nang)
       else if (nang == 590) then
          call ld0590(xang,yang,zang,wang,nang)
       else if (nang == 770) then
          call ld0770(xang,yang,zang,wang,nang)
       else if (nang == 974) then
          call ld0974(xang,yang,zang,wang,nang)
       else if (nang == 1202) then
          call ld1202(xang,yang,zang,wang,nang)
       else if (nang == 1454) then
          call ld1454(xang,yang,zang,wang,nang)
       else if (nang == 1730) then
          call ld1730(xang,yang,zang,wang,nang)
       else if (nang == 2030) then
          call ld2030(xang,yang,zang,wang,nang)
       else if (nang == 2354) then
          call ld2354(xang,yang,zang,wang,nang)
       else if (nang == 2702) then
          call ld2702(xang,yang,zang,wang,nang)
       else if (nang == 3074) then
          call ld3074(xang,yang,zang,wang,nang)
       else if (nang == 3470) then
          call ld3470(xang,yang,zang,wang,nang)
       else if (nang == 3890) then
          call ld3890(xang,yang,zang,wang,nang)
       else if (nang == 4334) then
          call ld4334(xang,yang,zang,wang,nang)
       else if (nang == 4802) then
          call ld4802(xang,yang,zang,wang,nang)
       else if (nang == 5294) then
          call ld5294(xang,yang,zang,wang,nang)
       else if (nang == 5810) then
          call ld5810(xang,yang,zang,wang,nang)
       else
          call error("genmesh","unknown value of nang",2)
       end if

       ! 3d mesh, do not parallelize to get the nodes in order
       do ir = 1, nr
          r = rads(ir)
          do il = 1, nang
             x = m%x(:,i) + r * (/xang(il),yang(il),zang(il)/)
             do j = 2, m%n
                if (m%z(j) < 1) cycle
                do k = 1, j-1
                   if (m%z(k) < 1) cycle
                   r1 = sqrt((x(1)-m%x(1,j))**2+(x(2)-m%x(2,j))**2+(x(3)-m%x(3,j))**2)
                   r2 = sqrt((x(1)-m%x(1,k))**2+(x(2)-m%x(2,k))**2+(x(3)-m%x(3,k))**2)
                   hypr = (r1-r2) / rr(j,k)
                   hypr = 1.5d0*hypr-0.5d0*hypr**3
                   hypr = 1.5d0*hypr-0.5d0*hypr**3
                   hypr = 1.5d0*hypr-0.5d0*hypr**3
                   hypr = 1.5d0*hypr-0.5d0*hypr**3
                   cutoff(j,k) = (1d0-hypr) / 2d0
                   cutoff(k,j) = (1d0+hypr) / 2d0
                enddo
                cutoff(j,j) = 1d0
             enddo
             cutoff(1,1) = 1d0
             vp0 = 1d0
             vpsum = 0d0
             do j = 1, m%n
                if (m%z(j) < 1) cycle
                vp0=vp0*cutoff(i,j)
                vpi=1d0
                do k = 1, m%n
                   if (m%z(k) < 1) cycle
                   vpi = vpi * cutoff(j,k)
                enddo
                vpsum = vpsum + vpi
             enddo
             !$omp critical (mmesh)
             meshrl(il,ir,i) = vp0/vpsum * wrads(ir) * wang(il)
             meshx(:,il,ir,i) = x
             !$omp end critical (mmesh)
          enddo
       enddo
    end do
    !$omp end parallel do

    ! clean up
    if (allocated(rads)) deallocate(rads)
    if (allocated(wrads)) deallocate(wrads)
    if (allocated(xang)) deallocate(xang)
    if (allocated(yang)) deallocate(yang)
    if (allocated(zang)) deallocate(zang)
    if (allocated(wang)) deallocate(wang)

    ! fill the 3d mesh
    kk = 0
    do i = 1, m%n
       if (m%z(i) < 1) cycle
       nr = z2nr(m%z(i))
       nang = z2nang(m%z(i))
       do ir = 1, nr
          do il = 1, nang
             kk = kk + 1
             mesh%w(kk) = meshrl(il,ir,i)
             mesh%x(:,kk) = meshx(:,il,ir,i)
          enddo
       enddo
    enddo

  end function genmesh

  subroutine rmesh(n,rmid,r,wintr)
    implicit none
    !
    !     RADIAL MESH AND INTEGRATION WEIGHTS,
    !     DERIVATIVES OF VARIABLE R WITH RESPECT TO VARIABLE Q.
    !
    !     THE Q-MESH IS UNIFORM ON THE INTERVAL (0,+1). TRANSFORMATION IS
    !
    !                        Q
    !         R  =  RMID ---------
    !                    ( 1 - Q )
    !
    !     ALSO,
    !     3 AND 7-POINT FINITE DIFFERENCE MATRICES FOR D2(WRT)R ON Q-MESH.
    !

    integer, intent(in) :: n
    real*8, intent(in) :: rmid
    real*8, intent(out) :: r(n), wintr(n)

    real*8 :: h, q
    integer :: i

    h = 1d0/real(n+1,8)
    do i=1,n
       q=h*i
       r(i) = rmid * q / (1.d0-q)
       wintr(i) = fourpi * h * r(i)**2 * rmid / (1.d0-q)**2
    enddo

  end subroutine rmesh

end module meshmod
