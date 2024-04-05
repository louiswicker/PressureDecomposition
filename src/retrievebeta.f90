
    PROGRAM RETRIEVE_BETA

    use netcdf

    implicit none

!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------

  ! This code currently only works with no terrain
  ! (i.e., perfectly flat lower boundary must be used)

  ! (NOTE:  can't use this code with horizontal grid stretching)
  ! (dx and dy must be constant)

  !  for boundary conditions:   1 = periodic ;   2 = open ;   3 = rigid wall
  !       wbc = west boundary condition
  !       ebc = east boundary condition
  !       sbc = south boundary condition
  !       nbc = north boundary condition

  ! The code computes effective buoyancy using the Davies-Jones (2003)
  ! Jeevanjee & Romps (2015) formulation

!-----------------------------------------------------------------------
    integer, parameter :: nv = 4

    real, dimension(:,:,:,:),allocatable :: rhs

    integer  :: i,j,k,n

    real, dimension(:),allocatable :: xh, yh, zf, zh, mfc, mfe

    real, dimension(:,:,:),allocatable :: tmp 
    real, dimension(:),    allocatable :: atri, btri, ctri, tmpz

    integer :: wbc, ebc, sbc, nbc 
    integer :: nx, ny, nz

    character(len=100) :: infile, outfile

    real :: dx, dy

    real, parameter :: g = 9.81

    character*10, dimension(nv) :: var_names

    namelist /inputs/ infile,outfile,nx,ny,nz,wbc,ebc,sbc,nbc

!----------------- 3D field names for output netCDF4 file --------------------!

    var_names(1) = "IC_Beta   "
    var_names(2) = "Soln_Beta "

!----------------- Read namelist in --------------------!

    write(*,*) ' ---> Retrieve_Beta: Reading in namelist'

    open(8,file='input.nml')
    read(8,nml=inputs)

    write(*,*)
    write(*,*) ' ---> Retrieve_Beta:  Namelist variables:'
    write(*,*) '   infile      = ',infile
    write(*,*) '   outfile     = ',outfile
    write(*,*) '   nx          = ',nx
    write(*,*) '   ny          = ',ny
    write(*,*) '   nz          = ',nz
    write(*,*) '   wbc         = ',wbc
    write(*,*) '   ebc         = ',ebc
    write(*,*) '   sbc         = ',sbc
    write(*,*) '   nbc         = ',nbc

!----------------- Allocate Storage --------------------!

    allocate( xh(nx) )
    allocate( yh(ny) )
    allocate( zh(nz) )
    allocate( mfc(nz) )
    allocate( zf(nz+1) )
    allocate( mfe(nz+1) )

! Tridiagonal coefficients

    allocate( atri(nz+1) )
    allocate( btri(nz) )
    allocate( ctri(nz+1) )
    allocate( tmpz(nz+1) )

! 3D arrays for solution

    allocate( tmp(nx,ny,nz) )

    allocate( rhs(nx,ny,nz,nv) )

    tmp(:,:,:)   = 0.0
    rhs(:,:,:,:) = 0.0

!----------------- Read 3D density from input netCDF4 file --------------------!

    write(*,*) ' ---> Retrieve_Beta: Reading in 3D density'

    CALL READNC2( infile, nx, ny, nz, xh, yh, zh, tmp ) 

    rhs(:,:,:,1) = tmp

    write(6,*) xh

    write(*,*) ' ---> Retrieve_Beta: Read in 3D density'

    zf(1) = 0.0

    DO k=2,nz+1
      zf(k) = 2.0*zh(k-1) - zf(k-1)
    ENDDO

    DO k=1,nz
      mfc(k) = 1.0 / (zf(k+1) - zf(k))
    ENDDO

    DO k=2,nz
      mfe(k) = 1.0 / (zh(k) - zh(k-1))
    ENDDO

    mfe(1)    = mfe(2)
    mfe(nz+1) = mfe(nz)

    dx = xh(2) - xh(1)   ! assume constant grid in horizontal
    dy = yh(2) - yh(1)

! Compute tridiagonal coefficients for pressure/density formulation

    DO k = 1,nz

       atri(k) = mfc(k)*mfe(k) 

       ctri(k) = mfc(k)*mfe(k+1) 

       btri(k) = - atri(k) - ctri(k)

    ENDDO
    
! Call Horizontal laplacian operator

    write(*,*) ' ---> Retrieve_Beta: computing laplacian of density'

    call DELSQH(rhs(1,1,1,1), tmp, dx, dy, nx, ny, nz)

    call writemxmn(tmp, nx, ny, nz, 'DEL^2 DENSITY')

    rhs(:,:,:,2) = -g * tmp(:,:,:)

! Set boundary conditions for beta=0 at ground, this is a reflective bc

    btri( 1) = btri( 1) - atri( 1) 
    btri(nz) = btri(nz) - ctri(nz) 

! Solve elliptic system for Beta

    write(*,*) ' ---> Retrieve_Beta: Ready to solve elliptic system'

    call pdcomp2024(nx, ny, nz, wbc, ebc, sbc, nbc, dx, dy, atri, ctri, btri, rhs(1,1,1,2), tmp)

    rhs(:,:,:,2) = tmp(:,:,:)

    write(*,*) ' ---> Retrieve_Beta:  Finished computing BETA'

    call writemxmn(rhs(1,1,1,2), nx, ny, nz, var_names(2))

    write(6,FMT='(" ------------------------------------------------------------")')
    write(*,*)

! Write data to netCDF-format file:

    write(*,*) ' ---> Retrieve_Beta::  Writing netCDF4 to outfile'

    call writenc2(outfile, nx, ny, nz, 2, xh, yh, zh, rhs, var_names)

    write(*,*) ' ---> Retrieve_Beta::  Wrote netCDF4 to outfile'

    END PROGRAM RETRIEVE_BETA
