
    PROGRAM getpp

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
    
    integer, parameter :: nx = 101
    integer, parameter :: ny = 101
    integer, parameter :: nz = 60
    integer, parameter :: nv = 2

    integer, parameter :: pow = 2

    real,    parameter :: dx = 50.
    real,    parameter :: dy = 50.
    real,    parameter :: dz = 25.

    real, dimension(:,:,:,:),allocatable :: rhs

    integer  :: i,j,k,n

    real, dimension(:),allocatable :: zf,zh,prs0,u0,v0,rho0,pi0,th0,qv0, xh, yh
    real, dimension(:),allocatable :: thv0

    real, dimension(:,:,:),allocatable :: u
    real, dimension(:,:,:),allocatable :: v
    real, dimension(:,:,:),allocatable :: w, tmp

    real, dimension(:,:,:),allocatable :: qvpert,thpert,thv,prspert
    real, dimension(:,:,:),allocatable   :: unc
    real, dimension(:,:,:),allocatable   :: vnc
    real, dimension(:,:,:),allocatable   :: wnc

    real, dimension(:,:,:),allocatable :: den, rho, th,qv,prs

    real zfac, hradius, zradius
    real tv0, tv1, pavgin
    integer, parameter :: wbc = 1
    integer, parameter :: ebc = 1
    integer, parameter :: sbc = 1
    integer, parameter :: nbc = 1
    character*7, parameter :: outfile = "test.nc"

    character*10, dimension(4) :: var_names

    

! Declarations for creating local base state
!-----------------------------------------------------------------------

! WK thermodynamic sounding parameters

    real, parameter :: ztr  = 12000. 
    real, parameter :: thtr = 343.   
    real, parameter :: ttr  = 213.   
    real, parameter :: tsfc = 300.   
    real, parameter :: psfc = 100000.

! Jeevanjee and Romps perturbations

    real, parameter :: radh = 1000.
    real, parameter :: radz = 500.
    real, parameter :: drho = 300.
    real, parameter :: xc   = 0.
    real, parameter :: yc   = 0.
    real, parameter :: zc   = 0.

! Other

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
    var_names(2) = "IC_Dens   "
    var_names(3) = "Soln_Beta "
    var_names(4) = "Soln_Dens "

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

    allocate( rhs(nx,ny,nz,2*nv) )

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

!------- Compute 3D state using J & R 2015

    write(*,*) ' ---> GetPP: computing full state..'

! Calculate buoyancy

    DO k=1,nz

    DO j=1,ny
    DO i=1,nx

        hradius = sqrt((xh(i) - xc)**2 + (yh(j) - yc)**2 )
        zradius = sqrt((zh(k) - zc)**2)

        rhs(i,j,k,1) = rho0(k) + rho0(1)/drho * exp( -(hradius/radh)**pow - (zradius/radz)**pow)
        rhs(i,j,k,2) = rhs(i,j,k,1) - rho0(k)

    ENDDO
    ENDDO
    ENDDO

! Call Horizontal laplacian operator

    call DELSQH(rhs(1,1,1,1), tmp, dx, dy, nx, ny, nz, 'DENSITY')
    rhs(:,:,:,3) = -g*tmp(:,:,:)

    write(*,*) ' ---> GetPP:  Finished computing BETAs'

! Compute del(rho0*B) / del_Z

    DO j=1,ny
    DO i=1,nx

      tmp(i,j,1) = 0.5*(rho0(1)*rhs(i,j,1,2) + rho0(1)*rhs(i,j,1,2)) 
      DO k=2,nz
        tmp(i,j,k) = 0.5*(rho0(k-1)*rhs(i,j,k-1,2) + rho0(k)*rhs(i,j,k,2)) 
      ENDDO

      DO k=1,nz-1
        rhs(i,j,k,4) = g*(tmp(i,j,k+1) - tmp(i,j,k)) / dz
      ENDDO
      rhs(i,j,nz,4) = g*rhs(i,j,nz-1,4)

    ENDDO
    ENDDO

    write(6,FMT='(" ------------------------------------------------------------")')

    DO n = 1,nv
      call writemxmn(rhs(1,1,1,n), nx, ny, nz, var_names(n))
    ENDDO

    write(6,FMT='(" ------------------------------------------------------------")')
    write(*,*)

    write(*,*) ' ---> GetPP --> computed buoyancy term'

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Get new pressure perturbation: 

    write(*,*) ' ---> GetPP:  Calling pdcomp'
    write(*,*) ' ---> GetPP:  dx,dy = ',dx,dy
    write(*,*) ' ---> GetPP:  nx,ny,nz = ',nx,ny,nz
    write(*,*)

    DO n = 1,nv 
      call writemxmn(rhs(1,1,1,n), nx, ny, nz, var_names(n))
      call pdcomp2(nx, ny, nz, wbc, ebc, sbc, nbc, dx, dy, zh, rhs(1,1,1,n), rhs(1,1,1,2+n))
    ENDDO

    write(*,*) ' ---> GetPP:  Finished  pdcomp'

! Write data to netCDF-format file:

    write(*,*) 'GETPP:  Before writenc2'

    call writenc2(outfile, nx, ny, nz, 2*nv, rhs, var_names)

    write(*,*) 'GETPP:  After writenc2'

    end program getpp

!=========================================================
!
!
! Del^2 - horiz
!
!
!=========================================================
    SUBROUTINE DELSQH(input, output, dx, dy, nx, ny, nz, label)

    implicit none

    integer, intent(in) :: nx, ny, nz

    real, intent(in) :: dx, dy 

    real, dimension(nx,ny,nz), intent(in)  :: input

    real, dimension(nx,ny,nz), intent(out) :: output

    character(len=*), intent(in) :: label

    integer :: i,j,k

    output(:,:,:) = 0.0

    DO k=1,nz   ! outer loop

      DO j=2,ny-1
      DO i=2,nx-1

        output(i,j,k) = (input(i-1,j,k) - 2.0*input(i,j,k) + input(i+1,j,k)) / (dx**2) &
                      + (input(i,j-1,k) - 2.0*input(i,j,k) + input(i,j+1,k)) / (dy**2)

      ENDDO
      ENDDO

      DO j = 2,ny-1
        output(1,j,k)  = (input(nx,  j,  k) - 2.0*input(1, j,k) + input(2, j,  k)) / (dx**2) &
                       + (input(1,   j-1,k) - 2.0*input(1, j,k) + input(1, j+1,k)) / (dy**2)
        output(nx,j,k) = (input(nx-1,j,  k) - 2.0*input(nx,j,k) + input(1, j,  k)) / (dx**2) &
                       + (input(nx,  j-1,k) - 2.0*input(nx,j,k) + input(nx,j+1,k)) / (dy**2)
      ENDDO

      DO i = 2,nx-1
        output(i,1,k)  = (input(i-1,   1,k) - 2.0*input(i, 1,k) + input(i+1, 1,k)) / (dx**2) &
                       + (input(i,  ny-1,k) - 2.0*input(i, 1,k) + input(i,   2,k)) / (dy**2)
        output(i,ny,k) = (input(i-1,ny,  k) - 2.0*input(i,ny,k) + input(i+1,ny,k)) / (dx**2) &
                       + (input(i,  ny-1,k) - 2.0*input(i,ny,k) + input(i,   1,k)) / (dy**2)
      ENDDO

      print *, label, ': ', k, maxval(output(:,:,k)), minval(output(:,:,k))

    ENDDO

    RETURN

    END SUBROUTINE DELSQH


!=================== Write netCDF ==========================

    subroutine writenc2(filename, nx, ny, nz, nv, vars, labels)

    USE netcdf

    implicit none
 
    character(len=*),          intent(in) :: filename
    integer,                   intent(in) :: nx, ny, nz, nv

    real, dimension(nx,ny,nz,nv), intent(in) :: vars

    character(len=10), dimension(nv), intent(in) :: labels

! Local declarations

    integer :: n

    integer :: ncid,status,niDimID,njDimID,nkDimID,timeDimID

    integer, dimension(nv) :: VarID

    real, dimension(nx,ny,nz) :: tmp

!----------------- Create and open netCDF -----------------

    status = nf90_create(trim(filename),NF90_64BIT_OFFSET,ncid)

    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

!------------ Define dimensions and variables -------------

    status = nf90_def_dim(ncid,"nx",nx,niDimID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)
    
    status = nf90_def_dim(ncid,"ny",ny,njDimID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_dim(ncid,"nz",nz,nkDimID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_dim(ncid,"time",1,timeDimID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    DO n = 1,nv

      status = nf90_def_var(ncid,labels(n),nf90_float, &
                           (/niDimID,njDimID,nkDimID/),VarID(n))
                       !   (/niDimID,njDimID,nkDimID,timeDimID/),VarID(n))
      if(status /= nf90_NoErr) write(*,*) labels(n), nf90_strerror(status)

    ENDDO

    status = nf90_enddef(ncid)

    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

!------------ Write dimensions and variables -------------
    DO n = 1,nv

      tmp(:,:,:) = vars(:,:,:,n)

      status = nf90_put_var(ncid, VarID(n), tmp)
      if(status /= nf90_NoErr) write(*,*) labels(n), nf90_strerror(status)

    ENDDO

    status = nf90_close(ncid)

    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    return
    end
!=========================================================

subroutine writemxmn(array, nx, ny, nz, label)

  implicit none
  integer, intent(in)                   :: nx,ny,nz
  real, dimension(nx,ny,nz), intent(in) :: array
  character(len=*), intent(in)          :: label

  real thesum, theavg, stddev

  thesum = SUM(array)

  theavg = thesum / float(nx*ny*nz)

  stddev = SUM( (array - theavg)**2 )

  stddev=sqrt(stddev/float(nx*ny*nz))

  write(6,*)

  write(6,FMT='(" ---> VAR: ",a,2x,"MAX: ",g10.2,2x,"MIN: ",g10.2,2x,"STDDEV: " g10.2)') label, &
                maxval(array), minval(array), stddev

  write(6,*)

  return
  end
