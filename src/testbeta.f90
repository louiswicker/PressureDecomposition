
    PROGRAM TEST_BETA

    use netcdf

    implicit none

!-----------------------------------------------------------------------
!
!  080919:  Program to read CM1 GrADS-format output and calculate
!           pressure perturbations from the subroutine "pdcomp"
!
!  Last modified:  20 February 2013
!
!-----------------------------------------------------------------------

  ! NOTE:  this code (as well as pdcomp) assumes NO TERRAIN
  ! (i.e., perfectly flat lower boundary must be used)

  ! (NOTE:  can't use this code with horizontal grid stretching)
  ! (dx and dy must be constant)

  !  for boundary conditions:   1 = periodic ;   2 = open ;   3 = rigid wall
  !       wbc = west boundary condition
  !       ebc = east boundary condition
  !       sbc = south boundary condition
  !       nbc = north boundary condition

  ! BTYPE:  0 - conventional buoyancy
  !         1 - Davies-Jones (2003) formulation (Horiz Laplacian of density)

!-----------------------------------------------------------------------
! Jeevanjee and Romps perturbations

    real, parameter :: radh = 1000.
    real, parameter :: radz = 500.
    real, parameter :: drho = 300.
    real, parameter :: dthea= 3.0 / 300.
    real, parameter :: xc   = 0.
    real, parameter :: yc   = 0.
    real, parameter :: zc   = 0.

! Grid
    
    integer, parameter :: nx = 21
    integer, parameter :: ny = 21
    integer, parameter :: nz = 40
    integer, parameter :: nv = 2

    integer, parameter :: pow = 4

    real,    parameter :: dx = 250.
    real,    parameter :: dy = 250.
    real,    parameter :: dz = 30.

    real, dimension(:,:,:,:),allocatable :: rhs

    integer  :: i,j,k,n

    real, dimension(:),allocatable :: xh, yh, zf, zh, mfc, mfe, atri, btri, ctri, tmpz
    real, dimension(:),allocatable :: prs0, u0, v0, rho0, pi0, th0, qv0, thv0, rhoE

    real, dimension(:,:,:),allocatable :: u, v, w, tmp

    real, dimension(:,:,:),allocatable :: qvpert,thpert,thv,prspert

    real, dimension(:,:,:),allocatable :: den, rho, th, qv, prs

    real zfac, hradius, zradius, rho_avg
    real tv0, tv1, pavgin
    integer, parameter :: wbc = 1
    integer, parameter :: ebc = 1
    integer, parameter :: sbc = 1
    integer, parameter :: nbc = 1

    character*12, parameter :: outfile = "test_out.nc"
    character*10, parameter :: testfile = "beta_in.nc   "

    character*10, dimension(nv) :: var_names

    real :: pii 

! Declarations for creating local base state
!-----------------------------------------------------------------------

! WK thermodynamic sounding parameters

    real, parameter :: ztr  = 12000. 
    real, parameter :: thtr = 343.   
    real, parameter :: ttr  = 213.   
    real, parameter :: tsfc = 300.   
    real, parameter :: psfc = 100000.


    real, parameter :: p00    = 100000.0
    real, parameter :: rp00   = 1.0/p00
    real, parameter :: rd     = 287.04
    real, parameter :: cp     = 1005.7
    real, parameter :: rv     = 461.5
    real, parameter :: cv     = cp-rd
    real, parameter :: g      = 9.81
    real, parameter :: eps    = rd/rv
    real, parameter :: reps   = 1.0 / eps
    real, parameter :: repsm1 = eps - 1.0
    real, parameter :: rovcp  = rd/cp
    real, parameter :: rdocp  = rd/cp
    real, parameter :: rdocv  = rd/cv
    real, parameter :: rdorv  = rd/rv
    real, parameter :: cvdrd  = cv/rd
    real, parameter :: cpdrd  = cp/rd

    var_names(1) = "IC_Beta   "
    var_names(2) = "Soln_Beta "

    pii = 4.0*atan(1.0)

!-----------------Read netCDF--------------------

    write(*,*) ' ---> TestPP'

    allocate( xh(nx) )
    allocate( yh(ny) )
    allocate( zh(nz) )
    allocate( prs0(nz) )
    allocate( u0(nz) )
    allocate( v0(nz) )
    allocate( rho0(nz) )
    allocate( pi0(nz) )
    allocate( th0(nz) )
    allocate( thv0(nz) )
    allocate( qv0(nz) )
    allocate( zf(nz+1) )
    allocate( rhoE(nz+1) )
    allocate( mfc(nz) )
    allocate( mfe(nz+1) )

! Tridiagonal coefficients

    allocate( atri(nz+1) )
    allocate( btri(nz) )
    allocate( ctri(nz+1) )
    allocate( tmpz(nz+1) )

    allocate( u(nx,ny,nz) )
    allocate( v(nx,ny,nz) )
    allocate( w(nx,ny,nz) )
    allocate( tmp(nx,ny,nz) )
    allocate( thpert(nx,ny,nz) )
    allocate( thv(nx,ny,nz) )
    allocate( th(nx,ny,nz) )
    allocate( rho(nx,ny,nz) )
    allocate( den(nx,ny,nz) )
    allocate( prs(nx,ny,nz) )

! alloc space for IC and solution

    allocate( rhs(nx,ny,nz,nv) )

    rhs(:,:,:,:) = 0.0

! Simple WK sounding

    DO k = 1,nz

      zh(k) = (float(k) - 0.5) * dz
      zf(k) = (float(k) - 1.0) * dz

      IF( zh(k) .le. ztr ) THEN
          zfac   = ( zh(k) / ztr ) ** 1.25
          th0(k) = tsfc + ( thtr - tsfc ) * zfac
          qv0(k) = 1. - 0.75 * zfac
      ELSE
          th0(k) = thtr * exp (g * (zh(k) - ztr) / (cp * ttr) )
          qv0(k) = 0.25
      ENDIF

    ENDDO
    zf(nz+1) = (float(k)) * dz

    mfc(:) = 1.0   ! these are the mapping factors which are 1.0 for constant grid
    mfe(:) = 1.0

! PI0 via vertical integration of hydrostatic equation

    tv0    = th0(1)

    pi0(1) = (psfc/p00)**rdocp-0.5*g/(dz*tv0*cp)

    DO k = 2,nz

       tv1    = th0(k)
       pi0(k) = pi0(k-1) - 2.0*g*dz / ((tv0+tv1)*cp)
       tv0    = tv1

    ENDDO

! Get base state rho

    DO k = 1,nz

      prs0(k) = p00 * pi0(k)**(cp/rd)

      rho0(k) = prs0(k) / (rd * th0(k) * pi0(k))

      write(*,*) k, zh(k), pi0(k), prs0(k), th0(k), rho0(k)

    ENDDO

! Horizontal grid

    DO i = 1,nx
     xh(i) = -(nx/2)*dx + float(i-1) * dx
    ENDDO

    DO j = 1,ny
     yh(j) = -(ny/2)*dy + float(j-1) * dy
    ENDDO

! Compute tridiagonal coefficients for pressure/density formulation

    rhoE(:) = 0.0

    DO k = 2,nz
       rhoE(k) = 0.5*(rho0(k)+rho0(k-1))
    ENDDO

    rhoE(1)    = rho0(1)  ! These depend on boundary condition
    rhoE(nz+1) = rho0(nz)

    DO k = 1,nz

!      atri(k) = mfc(k)*mfe(k)*rhoE(k) / (dz*dz*rho0(k))
!      ctri(k) = mfc(k)*mfe(k+1)*rhoE(k+1) / (dz*dz*rho0(k))
!      btri(k) = - atri(k) - ctri(k)

       atri(k) = mfc(k)*mfe(k) / (dz*dz)

       ctri(k) = mfc(k)*mfe(k+1) / (dz*dz)

       btri(k) = - atri(k) - ctri(k)

    ENDDO
    
!------- Compute 3D state using J & R 2015

    write(*,*) ' ---> TEST_BETA: computing full state density'

! Calculate buoyancy

    DO k=1,nz

    DO j=1,ny
    DO i=1,nx

        hradius = sqrt((xh(i) - xc)**2 + (yh(j) - yc)**2 )
        zradius = sqrt((zh(k) - zc)**2)

        rhs(i,j,k,1) = rho0(k) - rho0(1)/drho * exp( -(hradius/radh)**pow - (zradius/radz)**pow)  ! J&R pp. 3202

    ENDDO
    ENDDO

    ENDDO
    
! Write out field so we can test the retrieve program...

    call writenc2(testfile, nx, ny, nz, 1, xh, yh, zh, rhs, 'den')

! Call Horizontal laplacian operator

    call DELSQH(rhs(1,1,1,1), tmp, dx, dy, nx, ny, nz, 'DENSITY')
    rhs(:,:,:,2) = -g * tmp(:,:,:) 

! Set boundary conditions for beta=0 at ground, this is a reflective bc

    btri( 1) = btri( 1) - atri( 1) 
    btri(nz) = btri(nz) - ctri(nz) 

! Solve elliptic system for Beta

    call writemxmn(rhs(1,1,1,1), nx, ny, nz, var_names(1))
    call pdcomp2024(nx, ny, nz, wbc, ebc, sbc, nbc, dx, dy, atri, ctri, btri, rhs(1,1,1,2), tmp)
    rhs(:,:,:,2) = tmp(:,:,:)

    write(*,*) ' ---> TEST_BETA:  Finished computing BETAs'

    call writemxmn(rhs(1,1,1,2), nx, ny, nz, var_names(2))

    write(6,FMT='(" ------------------------------------------------------------")')
    write(*,*)

! Write data to netCDF-format file:

    write(*,*) ' ---> TEST_BETA:  Writing netCDF4 file'

    call writenc2(outfile, nx, ny, nz, 2, xh, yh, zh, rhs, var_names)

    write(*,*) ' ---> TEST_BETA:  Wrote netCDF4 file'

    END PROGRAM TEST_BETA
