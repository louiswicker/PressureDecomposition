
!===========================================================
!
!   WRITENC2 (simple netCDF4 write utility)
!
!===========================================================


    SUBROUTINE WRITENC2(filename, nx, ny, nz, nv, x, y, z, vars, labels)

    USE netcdf

    implicit none
 
    character(len=*),          intent(in) :: filename
    integer,                   intent(in) :: nx, ny, nz, nv

    real, dimension(nx,ny,nz,nv), intent(in) :: vars
    real, dimension(nx), intent(in)          :: x
    real, dimension(ny), intent(in)          :: y
    real, dimension(nz), intent(in)          :: z

    character(len=10), dimension(nv), intent(in) :: labels

! Local declarations

    integer :: n

    integer :: ncid, status
    integer :: nxDimID, nyDimID, nzDimID, ntDimID
    integer ::  xVarID,  yVarID,  zVarID,  tVarID

    integer, dimension(nv) :: VarID

    real, dimension(nx,ny,nz) :: tmp

!----------------- Create and open netCDF -----------------

    status = nf90_create(trim(filename),NF90_64BIT_OFFSET,ncid)

    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

!------------ Define dimensions and variables -------------

    status = nf90_def_dim(ncid,"nx",nx,nxDimID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)
    
    status = nf90_def_dim(ncid,"ny",ny,nyDimID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_dim(ncid,"nz",nz,nzDimID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_dim(ncid,"time",1,ntDimID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"xh",nf90_float,(/nxDimID/), xVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"yh",nf90_float,(/nyDimID/), yVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    status = nf90_def_var(ncid,"zh",nf90_float,(/nzDimID/), zVarID)
    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    DO n = 1,nv

      status = nf90_def_var(ncid,labels(n),nf90_float, &
                           (/nxDimID,nyDimID,nzDimID/),VarID(n))
      if(status /= nf90_NoErr) write(*,*) labels(n), nf90_strerror(status)

    ENDDO

    status = nf90_enddef(ncid)

    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

!------------ Write coordinates and variables -------------

    status = nf90_put_var(ncid, xVarID, x)
    if(status /= nf90_NoErr) write(*,*) 'X-COORD: ', nf90_strerror(status)

    status = nf90_put_var(ncid, yVarID, y)
    if(status /= nf90_NoErr) write(*,*) 'Y-COORD: ', nf90_strerror(status)

    status = nf90_put_var(ncid, zVarID, z)
    if(status /= nf90_NoErr) write(*,*) 'Z-COORD: ', nf90_strerror(status)

    DO n = 1,nv

      tmp(:,:,:) = vars(:,:,:,n)

      status = nf90_put_var(ncid, VarID(n), tmp)
      if(status /= nf90_NoErr) write(*,*) labels(n), nf90_strerror(status)

    ENDDO

    status = nf90_close(ncid)

    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

    RETURN
    END SUBROUTINE WRITENC2

!===========================================================
!
!   READNC2 (simple netCDF4 read utility)
!
!===========================================================

    SUBROUTINE READNC2( filename, nt, nx, ny, nz, xc, yc, zc, den ) 

    use netcdf

    implicit none

    integer, intent(in) :: nt, nx, ny, nz
    character(len=100), intent(in) :: filename

    real, dimension(nx), intent(out) :: xc
    real, dimension(ny), intent(out) :: yc

    real, dimension(nt,nx,ny,nz), intent(out) :: den, zc

    integer :: k
    integer :: varid, ncid, status

!------------------Open netCDF-----------------------------

    status = nf90_open(trim(filename),nf90_nowrite,ncid)

    if(status /= nf90_NoErr) write(*,*) nf90_strerror(status)

!----------Get 1D variables needed from netcdf----------------

    status = nf90_inq_varid(ncid,"time",varid)
    status = nf90_get_var(ncid,varid,time,start=(/1/),count=(/nt/))

    status = nf90_inq_varid(ncid,"zh",varid)
    status = nf90_get_var(ncid,varid,zc,start=(/1/),count=(/nx,ny,nz,nt/))

    status = nf90_inq_varid(ncid,"yh",varid)
    status = nf90_get_var(ncid,varid,yc,start=(/1/),count=(/ny/))

    status = nf90_inq_varid(ncid,"xh",varid)
    status = nf90_get_var(ncid,varid,xc,start=(/1/),count=(/nx/))

!----------Get 3D variables needed from netcdf----------------

    status = nf90_inq_varid(ncid,"den",varid)

    IF( status /= nf90_NoErr) THEN
        write(*,*) ' ----> Retrieve_Beta/READNC2: No 3D density in file, stopping'
        stop 999
    ELSE
      status = nf90_get_var(ncid,varid,den,start=(/1,1,1,1/),count=(/nx,ny,nz,nt/))
    ENDIF

!------------------Close netCDF----------------------------

    status = nf90_close(ncid)

    RETURN
    END SUBROUTINE READNC2
