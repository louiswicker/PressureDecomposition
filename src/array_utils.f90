!=========================================================
!
!   WRITEMXMN:  simple check for array max/min/sttdev
!
!=========================================================
  SUBROUTINE WRITEMXMN(array, nx, ny, nz, label)

    implicit none
    integer, intent(in)                   :: nx,ny,nz
    real, dimension(nx,ny,nz), intent(in) :: array
    character(len=*), intent(in)          :: label

    real thesum, theavg, stddev

    thesum = SUM(array)

    theavg = thesum / float(nx*ny*nz)

    stddev = SUM( (array - theavg)**2 )

    stddev=sqrt(stddev/float(nx*ny*nz))

    write(6,FMT='(" ------------------------------------------------------------")')

    write(6,*)

    write(6,FMT='(" ---> VAR: ",a,2x,"MAX: ",g10.2,2x,"MIN: ",g10.2,2x,"STDDEV: " g10.2)') label, &
                  maxval(array), minval(array), stddev

    write(6,*)

    write(6,FMT='(" ------------------------------------------------------------")')
    write(6,*)

  RETURN
  END SUBROUTINE WRITEMXMN

!=========================================================
!
! DELSQH:  A horizontal laplacian operator 
!
!=========================================================
  SUBROUTINE DELSQH(input, output, dx, dy, nx, ny, nz)

    implicit none

    integer, intent(in) :: nx, ny, nz

    real, intent(in) :: dx, dy

    real, dimension(nx,ny,nz), intent(in)  :: input

    real, dimension(nx,ny,nz), intent(out) :: output

    integer :: i,j,k

    output(:,:,:) = 0.0

    call writemxmn(input, nx, ny, nz, 'DEL^2')
    write(6,*) dx, dy

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

    ENDDO

  RETURN

  END SUBROUTINE DELSQH
