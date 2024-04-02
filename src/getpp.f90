
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
    
    integer      :: nx,ny,nz,btype,ptype,wbc,ebc,sbc,nbc
    real         :: dx,dy
    integer      :: nk
    integer      :: i,j,k,imax,jmax,ntim
    real(kind=8) :: pavg, pavgin

    real, dimension(:),allocatable :: zf,zh,prs0,u0,v0,rho0,pi0,th0,qv0
    real, dimension(:),allocatable :: thv0

    real, dimension(:,:,:),allocatable :: u
    real, dimension(:,:,:),allocatable :: v
    real, dimension(:,:,:),allocatable :: w
    real, dimension(:,:,:),allocatable :: beta, betaT, betaP, pb,pdn,pdl,ppi, tmp

    real, dimension(:,:,:),allocatable :: ptdn,fpb,fptdn,fpdn
    real, dimension(:,:,:),allocatable :: fpdlx,fpdly,fpdlz
    real, dimension(:,:,:),allocatable :: fpdnx,fpdny,fpdnz
    real, dimension(:,:,:),allocatable :: fpbx,fpby,fpbz
    real, dimension(:,:,:),allocatable :: fptdnx,fptdny,fptdnz
    real, dimension(:,:,:),allocatable :: fpdl
    real, dimension(:,:,:),allocatable :: pex,psh,fpex,fpsh
    real, dimension(:,:,:),allocatable :: psum,fpsum
    real, dimension(:,:,:),allocatable :: qvpert,thpert,thv,prspert
    real, dimension(:,:,:),allocatable   :: unc
    real, dimension(:,:,:),allocatable   :: vnc
    real, dimension(:,:,:),allocatable   :: wnc

    real, dimension(:,:,:),allocatable :: den, rho, th,qv,prs
    real, dimension(:,:,:),allocatable :: qc,qr,qi,qs,qg,qhl

    character(len=100)  :: infile,outfile

    real, dimension(:,:),allocatable    :: qtot

! Declarations for creating local base state

    real :: sumng,sumprs0,sumth0,sumqv0,sumu0,sumv0

    integer, parameter :: avgstart = 1250
    integer, parameter :: avglen = 200

!-----------------------------------------------------------------------

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

    namelist /inputs/ infile,outfile,nx,ny,nz,dx,dy,btype,ptype,wbc,ebc,sbc,nbc

!-----------------Read netCDF--------------------

    write(*,*) ' ---> GetPP: Reading in namelist'

    open(8,file='pdcomp.input')
    read(8,nml=inputs)

    write(*,*)
    write(*,*) ' ---> GetPP:  Namelist variables:'
    write(*,*) '   infile      = ',infile
    write(*,*) '   outfile     = ',outfile
    write(*,*) '   nx          = ',nx
    write(*,*) '   ny          = ',ny
    write(*,*) '   nz          = ',nz
    write(*,*) '   dx          = ',dx
    write(*,*) '   dy          = ',dy
    write(*,*) '   ptype       = ',ptype
    write(*,*) '   btype       = ',btype
    write(*,*) '   wbc         = ',wbc
    write(*,*) '   ebc         = ',ebc
    write(*,*) '   sbc         = ',sbc
    write(*,*) '   nbc         = ',nbc

    IF ( ptype .ne. 27 .and. ptype .ne. 5 .and. ptype .ne. 1 .and. ptype .ne. 0 .and. ptype .ne. -1 ) THEN

      write(*,*) ' ---> GetPP:  Invalid microphysics (ptype): Can only use ptype = (0, 1, 5, or 27)'
      IF ( ptype .eq.  0) write(*,*) ' ---> GetPP:  PTYPE == 0, will only read in QV field '
      IF ( ptype .eq. -1) write(*,*) ' ---> GetPP:  PTYPE == -1, will only read in THETA field'
      STOP

    ENDIF

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
    allocate( u(-2:nx+4,-2:ny+3,nz) )
    allocate( v(-2:nx+3,-2:ny+4,nz) )
    allocate( w(-2:nx+3,-2:ny+3,nz+1) )
    allocate( tmp(nx,ny,nz) )
    allocate( beta(nx,ny,nz) )
    allocate( betaT(nx,ny,nz) )
    allocate( betaP(nx,ny,nz) )
    allocate( pb(nx,ny,nz) )
    allocate( pdn(nx,ny,nz) )
    allocate( pdl(nx,ny,nz) )
    allocate( ptdn(nx,ny,nz) )
    allocate( ppi(nx,ny,nz) )
    allocate( fpdlx(nx,ny,nz) )
    allocate( fpdly(nx,ny,nz) )
    allocate( fpdlz(nx,ny,nz) )
    allocate( fpdnx(nx,ny,nz) )
    allocate( fpdny(nx,ny,nz) )
    allocate( fpdnz(nx,ny,nz) )
    allocate( fpbx(nx,ny,nz) )
    allocate( fpby(nx,ny,nz) )
    allocate( fpbz(nx,ny,nz) )
    allocate( fptdnx(nx,ny,nz) )
    allocate( fptdny(nx,ny,nz) )
    allocate( fptdnz(nx,ny,nz) )
    allocate( pex(nx,ny,nz) )
    allocate( psh(nx,ny,nz) )
    allocate( fpex(nx,ny,nz) )
    allocate( fpsh(nx,ny,nz) )
    allocate( psum(nx,ny,nz) )
    allocate( fpsum(nx,ny,nz) )
    allocate( qvpert(nx,ny,nz) )
    allocate( thpert(nx,ny,nz) )
    allocate( thv(nx,ny,nz) )
    allocate( prspert(nx,ny,nz) )
    allocate( th(nx,ny,nz) )
    allocate( rho(nx,ny,nz) )
    allocate( den(nx,ny,nz) )
    allocate( qv(nx,ny,nz) )
    allocate( prs(nx,ny,nz) )
    allocate( qc(nx,ny,nz) )
    allocate( qr(nx,ny,nz) )
    allocate( qi(nx,ny,nz) )
    allocate( qs(nx,ny,nz) )
    allocate( qg(nx,ny,nz) )
    allocate( qhl(nx,ny,nz) )
    allocate( unc(nx+1,ny,nz) )
    allocate( vnc(nx,ny+1,nz) )
    allocate( wnc(nx,ny,nz+1) )
    allocate( qtot(nx,ny) )

    allocate( fpb(nx,ny,nz) )
    allocate( fpdl(nx,ny,nz) )
    allocate( fpdn(nx,ny,nz) )
    allocate( fptdn(nx,ny,nz) )

!-------Calculate values needed by pdcomp--------

    call readnc(infile, nx, ny, nz,                &
                zh,zf,th0,qv0,pi0,u0,v0,prs0,rho0, &
                rho, th, qv, qc, qr, qi, qs, qg, qhl,   &
                prs, unc, vnc, wnc, thpert, qvpert, prspert, &
                ptype )

    write(*,*) ' ---> GetPP: Read in full state'

! Calculate buoyancy

    DO k=1,nz

      qtot(:,:)     = 0.0
      den(:,:,k)    = 0.0
      beta(:,:,k)   = 0.0
      betaT(:,:,k)  = 0.0
      betaP(:,:,k)  = 0.0
      pb(:,:,k)     = 0.0

    DO j=1,ny
    DO i=1,nx

      IF ( ptype .eq. 27 ) THEN

        qtot(i,j) = qtot(i,j)+qc(i,j,k)+qr(i,j,k)+qi(i,j,k)+qs(i,j,k)+qg(i,j,k)+qhl(i,j,k)

      ELSEIF ( ptype .eq. 5 ) THEN

        qtot(i,j) = qtot(i,j)+qc(i,j,k)+qr(i,j,k)+qi(i,j,k) +qs(i,j,k)+qhl(i,j,k)

      ELSEIF ( ptype .eq. 1 ) THEN

        qtot(i,j) = qtot(i,j)+qc(i,j,k)+qr(i,j,k)

      ENDIF

      ppi(i,j,k) = (prs(i,j,k)*rp00)**rovcp - pi0(k)  ! compute pi-prime

      IF ( ptype < 0 ) THEN

        pb(i,j,k)  = g*(th(i,j,k)-th0(k))/th0(k)

      ELSEIF ( ptype == 0 ) THEN

        pb(i,j,k)  = g*((th(i,j,k)-th0(k))/th0(k) + repsm1*(qv(i,j,k)-qv0(k)))

!       pb(i,j,k)  = g*((thv(i,j,k)-thv0(k))/thv0(k))

      ELSEIF ( ptype > 0 ) THEN

! full theta_rho buoyancy 

        thv0(k)      = th0(k) * (1.0 + reps*qv0(k)) / (1.0 + qv(i,j,k))
        thv(i,j,k)   = th(i,j,k)*(1+reps*qv(i,j,k)) / (1.0 + qv(i,j,k) + qtot(i,j))
        pb(i,j,k)    = g*(thv(i,j,k)/thv0(k) - 1.0)

! use the Dawson approximation for splitting beta...

        thv0(k)      = th0(k)   *(1.0 + eps*qv0(k)) 
        thv(i,j,k)   = th(i,j,k)*(1.0 + eps*qv(i,j,k)) - qtot(i,j)
        betaT(i,j,k) = -(g * rho0(k) / thv0(k)) * ( thv(i,j,k) - thv0(k) )
        betaP(i,j,k) =  (g * (1.0-rdocp) * rho0(k) / prs0(k)) * (prs(i,j,k) - prs0(k))

! full beta from DJ 2003

        beta(i,j,k)  = g*rho(i,j,k) * (1.0 + qv(i,j,k) + qtot(i,j))

      ENDIF

    ENDDO
    ENDDO
    ENDDO

! Call orizontal laplacian operator

    call DELSQH(betaT, tmp, dx, dy, nx, ny, nz, 'BETA_T')
    betaT(:,:,:) = tmp(:,:,:)

    call DELSQH(betaP, tmp, dx, dy, nx, ny, nz, 'BETA_P')
    betaP(:,:,:) = tmp(:,:,:)

    call DELSQH(beta, tmp, dx, dy, nx, ny, nz, 'BETA_DJ')
    beta(:,:,:) = tmp(:,:,:)

    write(*,*) ' ---> GetPP:  Finished computing BETAs'

    write(6,FMT='(" ------------------------------------------------------------")')
    call writemxmn(beta, nx, ny, nz, 'BETA DJ')
    call writemxmn(betaP, nx, ny, nz, 'BETA P')
    call writemxmn(betaT, nx, ny, nz, 'BETA T')
    call writemxmn(pb, nx, ny, nz, 'BUOY')
    write(6,FMT='(" ------------------------------------------------------------")')
    write(*,*)


    write(*,*) ' ---> GetPP --> computed buoyancy term'

    do k=1,nz
    do j=1,ny+1
    do i=1,nx
      v(i,j,k) = vnc(i,j,k)
    enddo
    enddo
    do j=1,ny
    do i=1,nx+1
      u(i,j,k) = unc(i,j,k)
    enddo
    enddo
    enddo

    do k=1,nz+1
    do j=1,ny
    do i=1,nx
      w(i,j,k) = wnc(i,j,k)
    enddo
    enddo
    enddo

    write(*,*) ' ---> GetPP --> finished copying velocity field'
    write(*,*) ''

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    pavg = 0.0d0

    do j=1,ny
    do i=1,nx
      pavg = pavg + ppi(i,j,nz)
    enddo
    enddo

    pavg = pavg / (nx*ny)
    pavgin = pavg
    call writemxmn(ppi, nx, ny, nz, 'PII-PERT')
    write(*,*) ' ---> PII-AVG: ',pavg

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Get new pressure perturbation: 

    write(*,*) ' ---> GetPP:  Calling pdcomp'
    write(*,*) ' ---> GetPP:  dx,dy = ',dx,dy
    write(*,*) ' ---> GetPP:  nx,ny,nz = ',nx,ny,nz
    write(*,*)

!   call pdcomp(nx,ny,nz,wbc,ebc,sbc,nbc,dx,dy,          &
!               zh,rho0,th0,qv0,pi0,u0,v0,u,v,w,pavgin,  &
!               betaT, betaP, beta, pdn, pdl, ptdn, fpb, fptdn,fpdn,          &
!               pex,psh,psum,fpex,fpsh,fpsum,fpdl,fpdlx, &
!               fpdly,fpdlz,fpdnx,fpdny,fpdnz,fpbx,fpby, &
!               fpbz,fptdnx,fptdny,fptdnz)

    call pdcomp2(nx,ny,nz,wbc,ebc,sbc,nbc,dx,dy,          &
                 zh,rho0,th0,qv0,pi0,u0,v0,u,v,w,pavgin,  &
                 beta, betaP, betaT, pb, fpb)

    tmp = betaP+betaT

    write(*,*) ' ---> GetPP:  Finished  pdcomp'

! Write data to netCDF-format file:

    write(*,*) 'GETPP:  Before writenc'

!   call writenc(outfile,nx,ny,nz, betaT, betaP, pb, fpb, pdl, pdn, fptdnx,  &
!                fptdny,fptdnz,fpbx,fpby,fpbz,fpdnx,  &
!                fpdny,fpdnz,fpdlx,fpdly,fpdlz)

    call writenc2(outfile,nx,ny,nz, betaP, betaT, beta, tmp, pb, fpb)

    write(*,*) 'GETPP:  After writenc'

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
!=========================================================
!
! Things to add:  pass ptype in to read correct vars
!
!     call readnc(infile, nx, ny, nz,                &
!               zh,zf,th0,qv0,pi0,u0,v0,prs0,rho0, &
!               rho, th, qv, qc, qr, qi, qs, qg, qhl,   &
!               prs, unc, vnc, wnc, thpert, qvpert, prspert, &
!               ptype, btype )
!
!
!=================== Read NetCDF ==========================
    SUBROUTINE readnc(filename, nx, ny, nz,                        &
                      zh,zf,th0,qv0,pi0,u0,v0,prs0,rho0,           &
                      rho, th, qv, qc, qr, qi, qs, qg, qhl,        &
                      prs, unc, vnc, wnc, thpert, qvpert, prspert, &
                      ptype ) 
    use netcdf

    implicit none

    integer, intent(in) :: nx, ny, nz, ptype
    character(len=100), intent(in) :: filename

    real, dimension(nz), intent(out) :: th0,qv0,zh,zf,pi0,u0,v0,prs0
    real, dimension(nz), intent(out) :: rho0

    real, dimension(nx,ny,nz),   intent(out) :: rho,wnc,th,qv,prs,qc,qr,qi,qs,qg,qhl
    real, dimension(nx,ny,nz),   intent(out) :: thpert,qvpert,prspert
    real, dimension(nx+1,ny,nz), intent(out) :: unc
    real, dimension(nx,ny+1,nz), intent(out) :: vnc

    real    :: t0
    integer :: k
    integer :: varid, ncid, status

    real, parameter :: p00   = 100000.0
    real, parameter :: rp00  = 1.0/p00
    real, parameter :: rd    = 287.04
    real, parameter :: cp    = 1005.7
    real, parameter :: rv     = 461.5
    real, parameter :: cv     = cp-rd
    real, parameter :: g     = 9.81
    real, parameter :: reps   = rv/rd
    real, parameter :: repsm1 = rv/rd-1.0
    real, parameter :: rovcp  = rd/cp
    real, parameter :: rddcp  = rd/cp
    real, parameter :: rddcv  = rd/cv
    real, parameter :: rddrv  = rd/rv
    real, parameter :: cvdrd  = cv/rd
    real, parameter :: cpdrd  = cp/rd

    logical rho_not_here

!------------------Open netCDF-----------------------------

    status = nf90_open(trim(filename),nf90_nowrite,ncid)

    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

!----------Get 1D variables needed from netcdf----------------

    status = nf90_inq_varid(ncid,"zh",varid)
    status = nf90_get_var(ncid,varid,zh,start=(/1/),count=(/nz/))

    status = nf90_inq_varid(ncid,"zf",varid)
    status = nf90_get_var(ncid,varid,zf,start=(/1/),count=(/nz+1/))

    status = nf90_inq_varid(ncid,"th0",varid)
    status = nf90_get_var(ncid,varid,th0,start=(/1,1,1,1/),count=(/1,1,nz,1/))

    status = nf90_inq_varid(ncid,"qv0",varid)
    status = nf90_get_var(ncid,varid,qv0,start=(/1,1,1,1/),count=(/1,1,nz,1/))

    status = nf90_inq_varid(ncid,"prs0",varid)
    status = nf90_get_var(ncid,varid,prs0,start=(/1,1,1,1/),count=(/1,1,nz,1/))

    status = nf90_inq_varid(ncid,"u0",varid)
    status = nf90_get_var(ncid,varid,u0,start=(/1,1,1,1/),count=(/1,1,nz,1/))

    status = nf90_inq_varid(ncid,"v0",varid)
    status = nf90_get_var(ncid,varid,v0,start=(/1,1,1,1/),count=(/1,1,nz,1/))

!----------Get variables needed from netcdf----------------

    status = nf90_inq_varid(ncid,"rho",varid)

    IF( status /= nf90_NoErr) THEN
        write(*,*) 'No 3D density in file'
        rho_not_here = .true.
    ELSE
      status = nf90_get_var(ncid,varid,rho,start=(/1,1,1,1/),count=(/nx,ny,nz,1/))
        rho_not_here = .false.
    ENDIF

    status = nf90_inq_varid(ncid,"th",varid)
    status = nf90_get_var(ncid,varid,th,start=(/1,1,1,1/),count=(/nx,ny,nz,1/))

    status = nf90_inq_varid(ncid,"qv",varid)
    status = nf90_get_var(ncid,varid,qv,start=(/1,1,1,1/),count=(/nx,ny,nz,1/))

    status = nf90_inq_varid(ncid,"qc",varid)
    status = nf90_get_var(ncid,varid,qc,start=(/1,1,1,1/),count=(/nx,ny,nz,1/))

    status = nf90_inq_varid(ncid,"qr",varid)
    status = nf90_get_var(ncid,varid,qr,start=(/1,1,1,1/),count=(/nx,ny,nz,1/))

    IF( ptype .eq. 5 .or. ptype .eq. 27 ) THEN

        status = nf90_inq_varid(ncid,"qi",varid)
        status = nf90_get_var(ncid,varid,qi,start=(/1,1,1,1/),count=(/nx,ny,nz,1/))

        status = nf90_inq_varid(ncid,"qs",varid)
        status = nf90_get_var(ncid,varid,qs,start=(/1,1,1,1/),count=(/nx,ny,nz,1/))

        status = nf90_inq_varid(ncid,"qg",varid)
        status = nf90_get_var(ncid,varid,qg,start=(/1,1,1,1/),count=(/nx,ny,nz,1/))

        status = nf90_inq_varid(ncid,"qhl",varid)
        status = nf90_get_var(ncid,varid,qhl,start=(/1,1,1,1/),count=(/nx,ny,nz,1/))

    ENDIF

    status = nf90_inq_varid(ncid,"prs",varid)
    status = nf90_get_var(ncid,varid,prs,start=(/1,1,1,1/),count=(/nx,ny,nz,1/))

    status = nf90_inq_varid(ncid,"u",varid)
    status = nf90_get_var(ncid,varid,unc,start=(/1,1,1,1/),count=(/nx+1,ny,nz,1/))

    status = nf90_inq_varid(ncid,"v",varid)
    status = nf90_get_var(ncid,varid,vnc,start=(/1,1,1,1/),count=(/nx,ny+1,nz,1/))

    status = nf90_inq_varid(ncid,"w",varid)
    status = nf90_get_var(ncid,varid,wnc,start=(/1,1,1,1/),count=(/nx,ny,nz+1,1/))

! Compute rho0
    DO k=1,nz
      pi0(k)  = (prs0(k)/p00)**(rd/cp)
      t0      = th0(k)*(prs0(k)/p00)**(0.2854)
      rho0(k) = prs0(k)/(rd*t0)  ! dry base state density
      zh(k)   = zh(k)*1000.0
      zf(k)   = zf(k)*1000.0
!     thv0(k) = th0(k)*(1.0+repsm1*qv0(k))
    enddo

!------------------Close netCDF----------------------------

    status = nf90_close(ncid)

    write(*,*) ' ---> GetPP:  Finished reading netCDF'

    write(6,FMT='(" ------------------------------------------------------------")') 
    call writemxmn(prs, nx, ny, nz, 'PRES')
    call writemxmn(th, nx, ny, nz, 'TH')
    call writemxmn(qv, nx, ny, nz, 'QV')
    call writemxmn(qc, nx, ny, nz, 'QC')
    call writemxmn(qr, nx, ny, nz, 'QR')
    call writemxmn(unc, nx+1, ny, nz, 'U')
    call writemxmn(vnc, nx, ny+1, nz, 'V')
    call writemxmn(wnc, nx, ny, nz+1, 'W')
    write(6,FMT='(" ------------------------------------------------------------")') 
    write(*,*)


    RETURN
    END SUBROUTINE
!=========================================================

subroutine writemxmn(array, nx, ny, nz, label)

  implicit none
  integer, intent(in)                   :: nx,ny,nz
  real, dimension(nx,ny,nz), intent(in) :: array
  character(len=*), intent(in)          :: label

  write(6,*)

  write(6,FMT='(" ---> VAR: ", a, 2x, "MAX: ", g10.2, 2x, "MIN: ", g10.2)') label, &
                maxval(array), minval(array)

  write(6,*)

  return
  end

!=================== Write netCDF ==========================
!   call writenc(outfile,nx,ny,nz, bb, pb, beta, fpb, pdl, pdn, fptdnx,  &
!                fptdny,fptdnz,fpbx,fpby,fpbz,fpdnx,  &
!                fpdny,fpdnz,fpdlx,fpdly,fpdlz)

    subroutine writenc(filename, nx, ny, nz,                       &
               betaT, betaP, pb, fpb,                              & 
               pdl, pdn, fptdnx, fptdny, fptdnz, fpbx, fpby, fpbz, &
               fpdnx,fpdny,fpdnz,fpdlx,fpdly, fpdlz)

    USE netcdf

    implicit none
 
    character(len=100),        intent(inout) :: filename

    integer,                   intent(in) :: nx,ny,nz

    real, dimension(nx,ny,nz), intent(in) :: fpdlx,fpdly,fpdlz
    real, dimension(nx,ny,nz), intent(in) :: fpdnx,fpdny,fpdnz
    real, dimension(nx,ny,nz), intent(in) :: fpbx,fpby,fpbz
    real, dimension(nx,ny,nz), intent(in) :: fptdnx,fptdny,fptdnz
    real, dimension(nx,ny,nz), intent(in) :: fpb, betaT, betaP, pb, pdl, pdn

    integer :: ncid,status,niDimID,njDimID,nkDimID,timeDimID
    integer :: fptdnVarID,fpbVarID,fpdnVarID,fpdlVarID
    integer :: fptdnxVarID,fptdnyVarID,fptdnzVarID
    integer :: fpbxVarID,fpbyVarID,fpbzVarID
    integer :: fpdlxVarID,fpdlyVarID,fpdlzVarID
    integer :: fpdnxVarID,fpdnyVarID,fpdnzVarID
    integer :: bdpdzVarID, betaTVarID, pbVarID, pdlVarID, pdnVarID, betaPVarID

!----------------- Create and open netCDF -----------------

    status = nf90_create(trim(filename),NF90_64BIT_OFFSET,ncid)

    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

!------------ Define dimensions and variables -------------

    status = nf90_def_dim(ncid,"ni",nx,niDimID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)
    
    status = nf90_def_dim(ncid,"nj",ny,njDimID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_dim(ncid,"nk",nz,nkDimID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_dim(ncid,"time",1,timeDimID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"fptdnx",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),fptdnxVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"fptdny",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),fptdnyVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"fptdnz",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),fptdnzVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"fpb",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),fpbVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"fpbx",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),fpbxVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"fpby",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),fpbyVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"fpbz",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),fpbzVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"fpdnx",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),fpdnxVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"fpdny",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),fpdnyVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

     status = nf90_def_var(ncid,"fpdnz",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),fpdnzVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"fpdlx",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),fpdlxVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"fpdly",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),fpdlyVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"fpdlz",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),fpdlzVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"pb",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),pbVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"betaT",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),betaTVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"betaP",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),betaPVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"pdl",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),pdlVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"pdn",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),pdnVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

!   status = nf90_def_var(ncid,"bdpdz",nf90_float, &
!                         (/niDimID,njDimID,nkDimID,timeDimID/),bdpdzVarID)
!   if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

!    status = nf90_def_var(ncid,"fpsh",nf90_float, &
!                         (/niDimID,njDimID,nkDimID,timeDimID/),fpshrVarID)
!    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_enddef(ncid)

    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

!------------ Write dimensions and variables -------------

    status = nf90_put_var(ncid,fptdnxVarID,fptdnx)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_put_var(ncid,fptdnyVarID,fptdny)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_put_var(ncid,fptdnzVarID,fptdnz)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_put_var(ncid,fpbVarID,fpb)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_put_var(ncid,fpbxVarID,fpbx)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_put_var(ncid,fpbyVarID,fpby)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_put_var(ncid,fpbzVarID,fpbz)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_put_var(ncid,fpdnxVarID,fpdnx)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_put_var(ncid,fpdnyVarID,fpdny)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_put_var(ncid,fpdnzVarID,fpdnz)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_put_var(ncid,fpdlxVarID,fpdlx)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_put_var(ncid,fpdlyVarID,fpdly)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_put_var(ncid,fpdlzVarID,fpdlz)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_put_var(ncid,pbVarID,pb)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_put_var(ncid,betaTVarID,betaT)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_put_var(ncid,betaPVarID,betaP)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_put_var(ncid,pdlVarID,pdl)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_put_var(ncid,pdnVarID,pdn)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

!   status = nf90_put_var(ncid,bdpdzVarID,bdpdz)
!   if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_close(ncid)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    return
    end

!=================== Write netCDF ==========================
!   call writenc(outfile,nx,ny,nz, beta, betaP, betaT, pb, 
!                dpbdz)

    subroutine writenc2(filename, nx, ny, nz,                       &
                       beta, betaT, betaP, betasum, pb, dpbdz)

    USE netcdf

    implicit none
 
    character(len=100),        intent(inout) :: filename

    integer,                   intent(in) :: nx,ny,nz

    real, dimension(nx,ny,nz), intent(in) :: beta, betaT, betaP, pb, dpbdz, betasum

    integer :: ncid,status,niDimID,njDimID,nkDimID,timeDimID

    integer :: fptdnVarID,fpbVarID,fpdnVarID,fpdlVarID
    integer :: fptdnxVarID,fptdnyVarID,fptdnzVarID
    integer :: fpbxVarID,fpbyVarID,fpbzVarID
    integer :: fpdlxVarID,fpdlyVarID,fpdlzVarID
    integer :: fpdnxVarID,fpdnyVarID,fpdnzVarID
    integer :: dpbdzVarID, betaVarID, betaTVarID, pbVarID, betaPVarID, bsumVarID

!----------------- Create and open netCDF -----------------

    status = nf90_create(trim(filename),NF90_64BIT_OFFSET,ncid)

    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

!------------ Define dimensions and variables -------------

    status = nf90_def_dim(ncid,"ni",nx,niDimID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)
    
    status = nf90_def_dim(ncid,"nj",ny,njDimID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_dim(ncid,"nk",nz,nkDimID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_dim(ncid,"time",1,timeDimID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"dpbdz",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),DpbDzVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"pb",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),pbVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"beta",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),betaVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"betaT",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),betaTVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"betaP",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),betaPVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"betaSum",nf90_float, &
                         (/niDimID,njDimID,nkDimID,timeDimID/),bsumVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_enddef(ncid)

    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

!------------ Write dimensions and variables -------------

    status = nf90_put_var(ncid,DpbDzVarID,dpbdz)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_put_var(ncid,pbVarID,pb)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_put_var(ncid,betaTVarID,betaT)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_put_var(ncid,betaPVarID,betaP)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_put_var(ncid,betaVarID,beta)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_put_var(ncid,bsumVarID,betasum)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_close(ncid)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    return
    end

! Move to a subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Write data to GrADS-format file:
!    datfile = trim(fl) // trim(datext)
!    open(unit=45,file=datfile,status='unknown',   &
!         form='unformatted',access='direct',recl=4*nx*ny)
!    orec = 1

!    do n=1,5
!    do n=1,9
!    do n=1,4
!    do k=1,nz
!      if(n.eq.1) write(45,rec=orec) (( bb(i,j,k),i=1,nx),j=1,ny)
!      if(n.eq.2) write(45,rec=orec) ((ppi(i,j,k),i=1,nx),j=1,ny)
!      if(n.eq.3) write(45,rec=orec) (( pb(i,j,k),i=1,nx),j=1,ny)
!      if(n.eq.4) write(45,rec=orec) ((pdn(i,j,k),i=1,nx),j=1,ny)
!      if(n.eq.5) write(45,rec=orec) ((pdl(i,j,k),i=1,nx),j=1,ny)
!      if(n.eq.6) write(45,rec=orec) ((ptdn(i,j,k),i=1,nx),j=1,ny)
!      if(n.eq.7) write(45,rec=orec) ((pex(i,j,k),i=1,nx),j=1,ny)
!      if(n.eq.8) write(45,rec=orec) ((psh(i,j,k),i=1,nx),j=1,ny)
!      if(n.eq.9) write(45,rec=orec) ((psum(i,j,k),i=1,nx),j=1,ny)
!      if(n.eq.1) write(45,rec=orec) ((fptdn(i,j,k),i=1,nx),j=1,ny)
!      if(n.eq.2) write(45,rec=orec) ((fpb(i,j,k),i=1,nx),j=1,ny)
!      if(n.eq.3) write(45,rec=orec) ((fpdn(i,j,k),i=1,nx),j=1,ny)
!      if(n.eq.4) write(45,rec=orec) ((fpex(i,j,k),i=1,nx),j=1,ny)
!      if(n.eq.5) write(45,rec=orec) ((fpsh(i,j,k),i=1,nx),j=1,ny)
!      if(n.eq.4) write(45,rec=orec) ((fpdl(i,j,k),i=1,nx),j=1,ny)

!      orec=orec+1
!    enddo
!    enddo
!    close(unit=45)

!--------------------------------------------
!  grads descriptor file for diagnostic output:
!    ctlfile = trim(fl) // trim(ctlext)
!    open(unit=30,file=ctlfile,status='unknown')
!    write(30,201)
!!!    write(30,212)
!    write(30,202)
!    write(30,203)
!    write(30,204) nx,0.001*0.5*dx,0.001*dx
!    write(30,205) ny,0.001*0.5*dy,0.001*dy
!    write(30,206) nz
!    do k=1,nz
!      write(30,211) 0.001*zh(k)
!211   format(2x,f12.6)
!    enddo
!    write(30,207)
!    write(30,208) 5
!    write(30,208) 9
!    write(30,208) 4
!    write(30,209) 'b       ',nz,'buoyancy from model                               '
!    write(30,209) 'ppi     ',nz,'actual pi^prime from model                        '
!    write(30,209) 'pb      ',nz,'diagnosed pi^prime:  buoyant component            '
!    write(30,209) 'pdn     ',nz,'diagnosed pi^prime:  nonlinear dynamic component  '
!    write(30,209) 'pdl     ',nz,'diagnosed pi^prime:  linear dynamic component     '
!    write(30,209) 'ptdn    ',nz,'diagnosed pi^prime:  total dynamic component      '
!    write(30,209) 'pex     ',nz,'diagnosed pi^prime:  exten dynamic component      '
!    write(30,209) 'psh     ',nz,'diagnosed pi^prime:  shear dynamic component      '
!    write(30,209) 'psum    ',nz,'diagnosed pi^prime:  sum of shr + ext compnt      '
!    write(30,209) 'fptdn   ',nz,'total dynamic pressure force                      '
!    write(30,209) 'fpb     ',nz,'buoyancy pressure force                           '
!    write(30,209) 'fpdn    ',nz,'nonlinear dynamic pressure force                  '
!    write(30,209) 'fpex    ',nz,'nonlinear dynamic pressure ext                    '
!    write(30,209) 'fpsh    ',nz,'nonlinear dynamic pressure shr                    '
!    write(30,209) 'fpsum   ',nz,'nonlinear dynamic pressure sum                    '
!    write(30,209) 'fpdl    ',nz,'linear dynamic pressure                           '
!    write(30,210)
!    close(unit=30)

!201 format('dset ^pdcomp.dat')
!212 format('options template')
!202 format('title CM1 output')
!203 format('undef -99999999.')
!204 format('xdef ',i6,' linear ',f12.6,1x,f12.6)
!205 format('ydef ',i6,' linear ',f12.6,1x,f12.6)
!206 format('zdef ',i6,' levels')
!207 format('tdef 1000 linear 00Z03JUL0001 1YR')
!208 format('vars ',i3)
!209 format(a8,1x,i6,' 99 ',a50)
!210 format('endvars')
