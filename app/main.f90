module uniplot
    implicit none
    private
    public :: fig_t
    type :: fig_t
        private
        integer :: nx, ny
        integer, allocatable :: array(:, :)
  contains
        procedure :: init     
        procedure :: point 
        procedure :: show  
        procedure :: line0
        procedure :: line
    end type fig_t
contains
    subroutine init(fig, nx, ny)
        class(fig_t), intent(out) :: fig 
        integer, intent(in) :: nx, ny
        fig%nx = nx
        fig%ny = ny 
        allocate(fig%array(0:(nx+1)/2, 0:(ny+3)/4) )
    end subroutine init  

  subroutine point(fig, ix, iy)
      class(fig_t), intent(in out) :: fig 
      integer, intent(in) :: ix, iy
      integer :: iax, iay
      iax = ix / 2
      iay = iy / 4
      ! clipping
      if (0<=ix .and. ix<fig%nx .and. 0<=iy .and. iy<fig%ny) then
          fig%array(iax, iay) = ior(fig%array(iax, iay), icode(mod(ix, 2), mod(iy, 4)))
      end if     
  end subroutine point

  pure elemental integer function icode(kx, ky)
      integer, intent(in) :: kx, ky
      if (ky == 3) then
          icode = 64 + 64 * kx
      else ! 0, 1, 2
          icode = 2**(ky + 3*kx)  
      end if
  end function icode

  subroutine line0(fig, ix0, iy0, ix1, iy1)
      class(fig_t), intent(in out) :: fig 
      integer, intent(in) :: ix0, iy0, ix1, iy1
      integer :: i, ix, iy, nx, ny
      real :: d
      nx = ix1 - ix0
      ny = iy1 - iy0
      if (nx == 0 .and. ny ==0) then
          call fig%point(ix, iy) 
      else if (abs(nx) < abs(ny)) then
          d = nx / real(ny)
          do i = 0, abs(ny)
              ix = nint(ix0 + d * sign(i, ny))
              iy = iy0 + sign(i, ny)
              call fig%point(ix, iy)
            end do  
      else
          d = ny / real(nx)
          do i = 0, abs(nx)
              iy = nint(iy0 + d * sign(i, nx))
              ix = ix0 + sign(i, nx)
              call fig%point(ix, iy)
          end do  
      end if  
  end subroutine line0

  subroutine show(fig)
      class(fig_t), intent(in) :: fig 
      integer :: iy
      do iy = 0, ubound(fig%array, 2)
          print '(*(a))', reverse_endian(shift_code(fig%array(:, iy)))
      end do
  end subroutine show    

  pure elemental integer function shift_code(k)
      integer, intent(in) :: k
      integer, parameter :: n0 = 14852224 !Z'E2A080' !14852224
      shift_code = n0 + 256 * (k /64) + mod(k, 64)  !E2A180, E2A280, E2A380      
  end function shift_code    
  
  pure elemental character(len = 4) function reverse_endian(i)
      integer, intent(in) :: i
      character:: tmp(4)
      tmp = transfer(i, ' ', size = 4)
      reverse_endian = transfer(tmp(4:1:-1), '    ')  !array 4 to len 4
  end function reverse_endian
  
  
  subroutine line(fig, x, y, ipen)
      class(fig_t), intent(in out) :: fig
      real, intent(in) :: x, y
      integer, intent(in) :: ipen
      integer, save :: ix0 = 0, iy0 = 0 
      integer :: ix, iy
      real, parameter :: xn = 80.0, yn = 100.0, fx = 1.0, fy = 0.85
      ix = nint( fx * x + xn)      
      iy = nint(-fy * y + yn)
      if (ipen == 1) call fig%line0(ix0, iy0, ix, iy)
      ix0 = ix
      iy0 = iy
  end subroutine line

end module uniplot
program uniplot_main
    use :: uniplot
    implicit none
    real, parameter :: pi = 4 * atan(1.0)
    real, allocatable :: x(:), h(:), theta(:)
    integer :: n
    n = 100
    allocate(x(n), h(n), theta(n))
    call random_seed()
    call random_number(x)
    call random_number(h)
    call random_number(theta)

plot: block 
        type(fig_t) :: fig1
        integer :: i, ix0, iy0, ix1, iy1, k
        k = 100
         
        print *
        print *, 'Buffon''s needle: estimated pi =', real(2 * n) / count(abs(h + sin(2 * pi * theta) - 0.5) > 0.5) 
        call fig1%init(3 * k, k)
        ! draw lines
        call fig1%line0(0, k/4, k + k-1, k/4)
        call fig1%line0(0, 2*k/4, k + k-1, 2*k/4)
        
        do i = 1, n 
            ix0 = k/4 + int(1.5 * k * x(i))
            iy0 = int(k/4 * h(i)) + k/4
            ix1 = ix0 + k/4 * cos(2 * pi * theta(i))
            iy1 = iy0 + k/4 * sin(2 * pi * theta(i))
            call fig1%line0(ix0, iy0, ix1, iy1)
        end do 
        call fig1%show()
    end block plot
  end program uniplot_main