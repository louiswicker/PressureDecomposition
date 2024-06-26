! =======================================================================
! 
! 

  SUBROUTINE SOLVE_ELLIP(nx,ny,nz,wbc,ebc,sbc,nbc,dx,dy,atri,ctri,btri,rhs3d,sol3d)

  USE SINGLETON

  implicit none

  integer, intent(in)                       :: nx, ny, nz
  integer, intent(in)                       :: wbc, ebc, sbc, nbc
  real,    intent(in)                       :: dx,dy
  real,    intent(in),  dimension(nx,ny,nz) :: atri, ctri, btri
  real,    intent(in),  dimension(nx,ny,nz) :: rhs3d
  real,    intent(out), dimension(nx,ny,nz) :: sol3d

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
!     nx  = number of grid points in x
!     ny  = number of grid points in y
!     nz  = number of grid points in z
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
! Output:
!
!   Real three-dimensional arrays:
!     soln (nx,ny,nz) 
!
!-----------------------------------------------------------------------

  integer :: i,j,k,nloop,ipb,ipe,jpb,jpe,kpb,kpe,imirror,jmirror
  real    :: rdx,rdy
  real*8 :: dpi,pavg,frac
  real, dimension(0:nz+1) :: thr0
  real, dimension(:,:,:), allocatable :: dum1,dum2,dum3,divx,uten,vten,wten,buoy

  real, dimension(0:nz+1) :: r1,rf0,rr0,mh,mf,zf

  complex, dimension(:,:),   allocatable :: rhs,trans
  complex, dimension(:,:,:), allocatable :: deft
  complex, dimension(0:nz+1)             :: diag,lgbth,lgbph
  real,    dimension(0:nz+1)             :: cfa,cfc

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

  write(6,*) 
  write(6,*) ' ---> SOLVE_ELLIPTICAL:  nx, ny, nz', nx, ny, nz

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
  ipe=nx

  jpb=1
  jpe=ny

  imirror = 0
  jmirror = 0

  if( (wbc.eq.2.or.wbc.eq.3).or.(ebc.eq.2.or.ebc.eq.3) )then

    imirror = 1
    ipe = nx*2

  endif

  if( (sbc.eq.2.or.sbc.eq.3).or.(nbc.eq.2.or.nbc.eq.3) )then

    jmirror = 1
    jpe = ny*2

  endif

  kpb=0
  kpe=nz+1

! print *,'  ipb,ipe         = ',ipb,ipe
! print *,'  jpb,jpe         = ',jpb,jpe
! print *,'  kpb,kpe         = ',kpb,kpe
! print *,'  imirror,jmirror = ',imirror,jmirror

  allocate(  deft(ipb:ipe,jpb:jpe,kpb:kpe) )

!----- constants -----

  rdx = 1.0/dx
  rdy = 1.0/dy

  dpi = 4.0d0*datan(1.0d0)

!----- retrieve pressure -----

    deft(:,:,:) = 0.0

!-----------------------------------------------------------------------
!  forcing for buoyancy pressure

    do k=1,nz
    do j=1,ny
    do i=1,nx
      deft(i,j,k) = rhs3d(i,j,k)  
    enddo
    enddo
    enddo

!-----------------------------------------------------------------------
!  p solver


! write(6,*) '  ipb,ipe,jpb,jpe,kpb,kpe = ',ipb,ipe,jpb,jpe,kpb,kpe
  write(6,*) 

  allocate(   rhs(ipb:ipe,jpb:jpe) )
  allocate( trans(ipb:ipe,jpb:jpe) )

  DO k=1,nz

    do j=1,ny
    do i=1,nx
      rhs(i,j)=deft(i,j,k)
    enddo
    enddo

    if(imirror.eq.1)then

      do j=1,ny
      do i=1,nx
        rhs(ipe+1-i,j)=rhs(i,j)
      enddo
      enddo

    endif

    if(jmirror.eq.1)then

      do j=1,ny
      do i=1,nx
        rhs(i,jpe+1-j)=rhs(i,j)
      enddo
      enddo

    endif

    if(imirror.eq.1.and.jmirror.eq.1)then

      do j=1,ny
      do i=1,nx
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

    DO k = 1,nz

      diag(k) = 2.0d0*( dcos(2.0d0*dpi*dble(i-1)/dble(ipe))                     &
                      + dcos(2.0d0*dpi*dble(j-1)/dble(jpe)) - 2.0d0 ) / (dx*dx) &
                      + btri(i,j,k)

    ENDDO

    IF( i .eq. 1 .and. j .eq. 1 ) THEN  

      r1(nz+1) = 0.0
      r1(nz)   = 0.0

      DO k = nz,2,-1

        r1(k-1) = (deft(i,j,k) - ctri(i,j,k)*r1(k+1) - diag(k)*r1(k)) / atri(i,j,k)

      ENDDO

      DO k = 1,nz

        deft(i,j,k) = cmplx( r1(k), 0.0 )

      ENDDO

    ELSE    ! Standard Thomas algorithm

      lgbth(1) = -ctri(i,j,1) / diag(1)
      lgbph(1) =  deft(i,j,1) / diag(1)

      DO k = 2,nz

        lgbth(k) = -ctri(i,j,k) / (atri(i,j,k)*lgbth(k-1) + diag(k))
        lgbph(k) = (deft(i,j,k) - atri(i,j,k)*lgbph(k-1)) / (atri(i,j,k)*lgbth(k-1) + diag(k))

      ENDDO

      deft(i,j,nz) = lgbph(nz)

      DO k = nz-1,1,-1

        deft(i,j,k) = lgbth(k)*deft(i,j,k+1) + lgbph(k)

      ENDDO

    ENDIF

  ENDDO
  ENDDO

  DO k=1,nz

    do j=jpb,jpe
    do i=ipb,ipe
      rhs(i,j)=deft(i,j,k)
    enddo
    enddo

    trans=fft(rhs,inv=.true.)

    do j=1,ny
    do i=1,nx
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
! do j=1,ny
! do i=1,nx
!   frac = frac + real(deft(i,j,nz))
! enddo
! enddo
!
! frac = frac / (nx*ny)
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
! do k=1,nz
! do j=1,ny
! do i=1,nx
!   deft(i,j,k) = deft(i,j,k) + frac
!   if(k.eq.nz) pavg = pavg + deft(i,j,k)
! enddo
! enddo
! enddo
!
! pavg = pavg / (nx*ny)
!
! print *,'  pavg      = ',pavg
! print *

!---------------------------------------------------

  do k=1,nz
  do j=1,ny
  do i=1,nx
      sol3d(i,j,k)=real(deft(i,j,k))
  enddo
  enddo
  enddo

  deallocate(deft)

  RETURN
  END SUBROUTINE SOLVE_ELLIP

