
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
  ! Jeevanjee & Romps (2015) formulation from the 3D/4D density field.


  ! This code has been modified to solve for CM1 and FV3-SOLO.  This
  ! requires (due to vertical grid in FV3) a 4D height array. This
  ! complicates a few things.  However, I have choosen to keep the
  ! rest of the vertical coordinates 3D to keep it readable as the
  ! coefficients are only needed as temporarys within the time loop.

!-----------------------------------------------------------------------

    real, dimension(:,:,:,:),allocatable :: rhs, soln

    integer  :: i,j,k,n

    real, dimension(:),       allocatable :: xh, yh, time
    real, dimension(:,:,:,:), allocatable :: zh

    real, dimension(:,:,:),   allocatable :: zf, mfc, mfe, tmp 
    real, dimension(:,:,:),   allocatable :: atri, btri, ctri

    integer :: wbc, ebc, sbc, nbc 
    integer :: nt, nx, ny, nz

! Command line args

    character(len=256) :: arg
    integer            :: narg

    character(len=256) :: infile, outfile, model

    real               :: dx, dy

    character*10       :: var_names

    real, parameter :: grav = 9.806

!----------------- 3D field names for output netCDF4 file --------------------!

    var_names = "Soln_Beta "

!----------------- Read namelist in --------------------!

    write(*,*) ' ---> Retrieve_Beta: Parsing Command Line arguments'

    narg = 0

    do while ( narg  <  command_argument_count() )

        narg = narg + 1

        call get_command_argument(narg, arg)
        print '(2a, /)', arg

        select case (arg)
            case ('-i', '--input')
                call get_command_argument(narg+1, arg)
                infile = arg
                print '(2a, /)', infile

            case ('-o', '--output')
                call get_command_argument(narg+1, arg)
                outfile = arg
                print '(2a, /)', outfile
        end select
    end do

    write(*,*) ' ---> Retrieve_Beta: Reading dims and attrs'

    CALL READ_NC4_DIMS( infile, nt, nx, ny, nz )

    CALL READ_NC4_ATT( infile, wbc, ebc, sbc, nbc, model )

    write(*,*)
    write(*,*) ' ---> Retrieve_Beta:  Input parameters\n'
    write(*,*)
    write(*,*) '  source model = ',model
    write(*,*)
    write(*,*) '  infile       = ',infile
    write(*,*) '  outfile      = ',outfile
    write(*,*) '  ntimes       = ',nt
    write(*,*) '  nx           = ',nx
    write(*,*) '  ny           = ',ny
    write(*,*) '  nz           = ',nz
    write(*,*) '  wbc          = ',wbc
    write(*,*) '  ebc          = ',ebc
    write(*,*) '  sbc          = ',sbc
    write(*,*) '  nbc          = ',nbc

!----------------- Allocate Storage --------------------!

    allocate( time(nt) )
    allocate( xh(nx) )
    allocate( yh(ny) )
    allocate( zh(nx,ny,nz,nt) )
    allocate( mfc(nx,ny,nz) )
    allocate( zf(nx,ny,nz+1) )
    allocate( mfe(nx,ny,nz+1) )

! Tridiagonal coefficients

    allocate( atri(nx,ny,nz+1) )
    allocate( btri(nx,ny,nz  ) )
    allocate( ctri(nx,ny,nz+1) )

! 4D arrays for forcing, temporary, and solution

    allocate( tmp (nx,ny,nz) )
    allocate( rhs (nx,ny,nz,nt) )
    allocate( soln(nx,ny,nz,nt) )

    tmp (:,:,:)   = 0.0
    rhs (:,:,:,:) = 0.0
    soln(:,:,:,:) = 0.0

!----------------- Read 3D density from input netCDF4 file --------------------!

    write(*,*) ' ---> Retrieve_Beta: Reading in 3D density'

    CALL READ_NC4_FILE( infile, nt, nx, ny, nz, xh, yh, zh, time, rhs ) 

    write(*,*) ' ---> Retrieve_Beta: Read in 3D density'

! assume constant grid in horizontal

    dx = (xh(2) - xh(1))
    dy = (yh(2) - yh(1))

! DO a time loop....

    DO n = 1,nt

      call writemxmn(rhs(1,1,1,n), nx, ny, nz, 'DENSITY')

      zf(:,:,:) = 0.0

      DO j = 1,ny
      DO i = 1,nx

        DO k=2,nz+1
          zf(i,j,k) = 2.0*zh(i,j,k-1,n) - zf(i,j,k-1)
        ENDDO

        DO k=1,nz
          mfc(i,j,k) = 1.0 / (zf(i,j,k+1) - zf(i,j,k))
        ENDDO

        DO k=2,nz
          mfe(i,j,k) = 1.0 / (zh(i,j,k,n) - zh(i,j,k-1,n))
        ENDDO

        mfe(i,j,1)    = mfe(i,j,2)
        mfe(i,j,nz+1) = mfe(i,j,nz)

! Compute tridiagonal coefficients for pressure/density formulation

        DO k = 1,nz

           atri(i,j,k) = mfc(i,j,k)*mfe(i,j,k) 
    
           ctri(i,j,k) = mfc(i,j,k)*mfe(i,j,k+1) 

           btri(i,j,k) = - atri(i,j,k) - ctri(i,j,k)

        ENDDO

! Set boundary conditions for beta=0 at ground, this is a reflective bc

        btri(i,j, 1) = btri(i,j, 1) - atri(i,j, 1) 
        btri(i,j,nz) = btri(i,j,nz) - ctri(i,j,nz) 

      ENDDO
      ENDDO

! Call Horizontal laplacian operator

      write(*,*) ' ---> Retrieve_Beta: Computing laplacian of density'

      call DELSQH(rhs(1,1,1,n), tmp, dx, dy, nx, ny, nz)

      call writemxmn(tmp, nx, ny, nz, 'LAPLACIAN of DENSITY')

      tmp(:,:,:) = -grav * tmp(:,:,:)

! Solve elliptic system for Beta

      write(*,*) ' ---> Retrieve_Beta: Ready to solve elliptic system for Time: ',time(n)

      call SOLVE_ELLIP(nx, ny, nz, wbc, ebc, sbc, nbc, dx, dy, atri, ctri, btri, tmp, soln(1,1,1,n))

      write(*,*) ' ---> Retrieve_Beta:  Finished computing BETA for Time: ',time(n)

      call writemxmn(soln(1,1,1,n), nx, ny, nz, var_names)

    ENDDO  ! time loop

    write(6,FMT='(" ------------------------------------------------------------")')
    write(*,*)

! Write data to netCDF-format file:

    write(*,*) ' ---> Retrieve_Beta::  Writing netCDF4 to outfile'

    call WRITE_NC4_FILE(outfile, nt, nx, ny, nz, xh, yh, zh, time, soln, var_names)

    write(*,*) ' ---> Retrieve_Beta::  Wrote netCDF4 to outfile'

    END PROGRAM RETRIEVE_BETA
