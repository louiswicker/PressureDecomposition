! =======================================================================
! 
! 

  subroutine pdcomp2024(ni,nj,nk,wbc,ebc,sbc,nbc,dx,dy,atri,ctri,btri,rhs3d,sol3d)

  USE SINGLETON

  implicit none

  integer, intent(in)                       :: ni,nj,nk
  integer, intent(in)                       :: wbc,ebc,sbc,nbc
  real,    intent(in)                       :: dx,dy
  real,    intent(in),  dimension(nk)       :: atri,ctri, btri
  real,    intent(in),  dimension(ni,nj,nk) :: rhs3d
  real,    intent(out), dimension(ni,nj,nk) :: sol3d

  real*8 :: pavgin

!-----------------------------------------------------------------------
!
!  pdcomp - a fortran90 subroutine to retrieve pressure perturbations
!           from cloud-scale numerical model output.  Three terms 
!           are retrieved:
!              pb  = buoyancy pressure
!              pdn = nonlinear dynamic pressure
!              pdl = linear dynamic pressure
!              ptdn = total dynamic pressure  (JT)
!
!  Version 1.03                           Last modified:  20 February 2013
!
!  Author:  George H. Bryan
!           Mesoscale and Microscale Meteorology Division
!           National Center for Atmospheric Research
!           Boulder, Colorado, USA
!           gbryan@ucar.edu
!
!  Disclaimer:  This code is made available WITHOUT WARRANTY.
!
!  References:  Rotunno and Klemp, 1982, MWR, p. 150
!               Weisman and Rotunno, 2000, JAS, p. 1461
!
!-----------------------------------------------------------------------
!
! Input:
!   Integer variables:
!     ni  = number of grid points in x
!     nj  = number of grid points in y
!     nk  = number of grid points in z
!
!     wbc = west boundary condition (see below)
!     ebc = east boundary condition (see below)
!     sbc = south boundary condition (see below)
!     nbc = north boundary condition (see below)
!
!  for boundary conditions:   1 = periodic ;   2 = open ;   3 = rigid wall
!       
!
!   Real variables:
!     dx  = grid spacing in x (m)
!     dy  = grid spacing in y (m)  (must be same as dx, for now!)
!
!   Real one-dimensional arrays:
!     zh  (nk)      = height of model's half levels (scalar points) (m)
!     rho0(nk)      = base-state density (kg/m^3)
!     th0 (nk)      = base-state potential temperature (K)
!     qv0 (nk)      = base-state water vapor mixing ratio (kg/kg)
!     pi0 (nk)      = base-state nondimensional pressure (dimensionless)
!      u0 (nk)      = base-state wind profile in x-direction (m/s)
!      v0 (nk)      = base-state wind profile in y-direction (m/s)
!
! Output:
!
!   Real three-dimensional arrays:
!     pb (ni,nj,nk) = buoyancy pressure
!     pdn(ni,nj,nk) = nonlinear dynamic pressure
!     pdl(ni,nj,nk) = linear dynamic pressure
!
!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!            No need to modify anything below here:
!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! JT notes:
! 1.  uses 6th order approximations to 1st derivative
!-----------------------------------------------------------------------

  integer :: i,j,k,nloop,ipb,ipe,jpb,jpe,kpb,kpe,imirror,jmirror
  real :: rdx,rdy
  real*8 :: dpi,pavg,frac
  real, dimension(0:nk+1) :: thr0
  real, dimension(:,:,:), allocatable :: dum1,dum2,dum3,divx,uten,vten,wten,buoy

  real, dimension(0:nk+1) :: r1,rf0,rr0,mh,mf,zf

  complex, dimension(:,:), allocatable :: rhs,trans
  complex, dimension(:,:,:), allocatable :: deft
  complex, dimension(0:nk+1) :: diag,lgbth,lgbph
  real,    dimension(0:nk+1) :: cfa,cfc

  real, dimension(ni,nj,nk) :: bb 

!-----------------------------------------------------------------------

  real, parameter :: g     = 9.81
  real, parameter :: rd    = 287.04
  real, parameter :: rv    = 461.5
  real, parameter :: cp    = 1005.7
  real, parameter :: pi    = 3.14159265358979323
  real, parameter :: p00   = 100000.0

  real, parameter :: eps   = rd/rv
  real, parameter :: reps  = 1.0/eps
  real, parameter :: rp00  = 1.0/p00
  real, parameter :: rddcp = rd/cp

!-----------------------------------------------------------------------

  print *,'PDCOMP2:  ni, nj, nz', ni, nj, nk

  if(wbc.eq.1.and.ebc.ne.1)then
    print *
    print *,'  Inconsistent boundary conditions'
    print *
    stop 1111
  endif
  if(ebc.eq.1.and.wbc.ne.1)then
    print *
    print *,'  Inconsistent boundary conditions'
    print *
    stop 1112
  endif
  if(sbc.eq.1.and.nbc.ne.1)then
    print *
    print *,'  Inconsistent boundary conditions'
    print *
    stop 1113
  endif
  if(nbc.eq.1.and.sbc.ne.1)then
    print *
    print *,'  Inconsistent boundary conditions'
    print *
    stop 1114
  endif
  if( (ebc.lt.1.or.ebc.gt.3) .or.   &
      (wbc.lt.1.or.wbc.gt.3) .or.   &
      (sbc.lt.1.or.sbc.gt.3) .or.   &
      (nbc.lt.1.or.nbc.gt.3) )then
    print *
    print *,'  Invalid setting for boundary conditions'
    print *
    stop 1115
  endif

!-----------------------------------------------------------------------

  ipb=1
  ipe=ni

  jpb=1
  jpe=nj

  imirror = 0
  jmirror = 0

  if( (wbc.eq.2.or.wbc.eq.3).or.(ebc.eq.2.or.ebc.eq.3) )then

    imirror = 1
    ipe = ni*2

  endif

  if( (sbc.eq.2.or.sbc.eq.3).or.(nbc.eq.2.or.nbc.eq.3) )then

    jmirror = 1
    jpe = nj*2

  endif

  kpb=0
  kpe=nk+1

  print *,'  ipb,ipe         = ',ipb,ipe
  print *,'  jpb,jpe         = ',jpb,jpe
  print *,'  kpb,kpe         = ',kpb,kpe
  print *,'  imirror,jmirror = ',imirror,jmirror
  allocate(  deft(ipb:ipe,jpb:jpe,kpb:kpe) )

!----- constants -----

  rdx = 1.0/dx
  rdy = 1.0/dy

  dpi = 4.0d0*datan(1.0d0)

!----- retrieve pressure -----

    print *,'PDCOMP2  '

    deft(:,:,:) = 0.0

!-----------------------------------------------------------------------
!  forcing for buoyancy pressure

    do k=1,nk
    do j=1,nj
    do i=1,ni
      deft(i,j,k) = rhs3d(i,j,k)  
    enddo
    enddo
    enddo

    DO j=jpb,jpe
    DO i=ipb,ipe
      deft(i,j, 1) = rhs3d(i,j, 1) + atri( 1)*rhs3d(i,j, 1)  ! reflective bc 
      deft(i,j,nk) = rhs3d(i,j,nk) + ctri(nk)*rhs3d(i,j,nk)
    ENDDO
    ENDDO

!-----------------------------------------------------------------------
!  p solver

  write(*,*) '  ipb,ipe,jpb,jpe,kpb,kpe = ',ipb,ipe,jpb,jpe,kpb,kpe
  write(*,*) '   alloc 1 '
  allocate(   rhs(ipb:ipe,jpb:jpe) )
  write(*,*) '   alloc 2 '
  allocate( trans(ipb:ipe,jpb:jpe) )

  DO k=1,nk

    do j=1,nj
    do i=1,ni
      rhs(i,j)=deft(i,j,k)
    enddo
    enddo

    if(imirror.eq.1)then

      do j=1,nj
      do i=1,ni
        rhs(ipe+1-i,j)=rhs(i,j)
      enddo
      enddo

    endif

    if(jmirror.eq.1)then

      do j=1,nj
      do i=1,ni
        rhs(i,jpe+1-j)=rhs(i,j)
      enddo
      enddo

    endif

    if(imirror.eq.1.and.jmirror.eq.1)then

      do j=1,nj
      do i=1,ni
        rhs(ipe+1-i,jpe+1-j)=rhs(i,j)
      enddo
      enddo

    endif

    trans=fft(rhs)

    DO j=jpb,jpe
    DO i=ipb,ipe
      deft(i,j,k)=trans(i,j)
    ENDDO
    ENDDO

  ENDDO

  DO j=jpb,jpe
  DO i=ipb,ipe

    do k=1,nk
      diag(k)=2.0d0*( dcos(2.0d0*dpi*dble(i-1)/dble(ipe))          &
                    +dcos(2.0d0*dpi*dble(j-1)/dble(jpe))          &
                    -2.0d0)/(dx*dx) - atri(k) - ctri(k)
    enddo

    if(i.eq.1.and.j.eq.1)then   ! know the solution at the top corner p_nh = 0.0

      r1(nk+1) = 0.0
      r1(nk)   = 0.0

      do k=nk,2,-1
        r1(k-1)=(deft(i,j,k)-ctri(k)*r1(k+1)-diag(k)*r1(k))/atri(k)
      enddo
      do k=1,nk
        deft(i,j,k)=cmplx( r1(k) , 0.0 )
      enddo

    else    ! Standard Thomas algorithm

      lgbth(1)=-ctri(1)/diag(1)
      lgbph(1)= deft(i,j,1)/diag(1)

      do k=2,nk
        lgbth(k)=-ctri(k)/(atri(k)*lgbth(k-1)+diag(k))
        lgbph(k)=(deft(i,j,k)-atri(k)*lgbph(k-1))/(atri(k)*lgbth(k-1)+diag(k))
      enddo
      deft(i,j,nk)=lgbph(nk)
      do k=nk-1,1,-1
        deft(i,j,k)=lgbth(k)*deft(i,j,k+1)+lgbph(k)
      enddo
    endif

  ENDDO
  ENDDO

  DO k=1,nk

    do j=jpb,jpe
    do i=ipb,ipe
      rhs(i,j)=deft(i,j,k)
    enddo
    enddo

    trans=fft(rhs,inv=.true.)

    do j=1,nj
    do i=1,ni
      deft(i,j,k)=real(trans(i,j))
    enddo
    enddo

  ENDDO

  deallocate(   rhs )
  deallocate( trans )

!----- adjust mean pressure --------------------------------------------
!  for pb, adjust mean pi along upper boundary to match total
!  for pd, adjust mean pi along upper boundary to be zero

! pavg = 0.0d0
! frac = 0.0d0
!
! do j=1,nj
! do i=1,ni
!   frac = frac + real(deft(i,j,nk))
! enddo
! enddo
!
! frac = frac / (ni*nj)
!
! pavg = pavgin
!
! print *
! print *,'  frac,pavg = ',frac,pavg
!
! frac = pavg - frac
! pavg = 0.0d0
!
! offset solution by mean pressure
!
! do k=1,nk
! do j=1,nj
! do i=1,ni
!   deft(i,j,k) = deft(i,j,k) + frac
!   if(k.eq.nk) pavg = pavg + deft(i,j,k)
! enddo
! enddo
! enddo
!
! pavg = pavg / (ni*nj)
!
! print *,'  pavg      = ',pavg
! print *

!---------------------------------------------------

  do k=1,nk
  do j=1,nj
  do i=1,ni
      sol3d(i,j,k)=real(deft(i,j,k))
  enddo
  enddo
  enddo

  deallocate(deft)

  RETURN
  END SUBROUTINE PDCOMP2024

