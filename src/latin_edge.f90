subroutine latin_edge ( dim_num, point_num, seed, x )

!*****************************************************************************80
!
!! LATIN_EDGE returns edge points in a Latin square.
!
!  Discussion:
!
!    In each spatial dimension, there will be exactly one
!    point with the coordinate value
!
!      ( 0, 1, 2, ..., point_num-1 ) / ( point_num - 1 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points, which should
!    be at least 2!
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator,
!    needed if the portable UNIFORM routine is being used.
!
!    Output, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) :: base = 1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) perm(point_num)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(dim_num,point_num)
!f2py intent(in) dim_num
!f2py intent(in) point_num
!f2py intent(out), depend(dim_num, point_num), dimension(dim_num, point_num) :: x

  if ( point_num == 1 ) then

    x(1:dim_num,1) = 0.5D+00

  else

    do i = 1, dim_num

      call perm_uniform ( point_num, base, seed, perm )

      do j = 1, point_num
        x(i,j) = real (   perm(j) - 1, kind = 8 ) &
               / real ( point_num - 1, kind = 8 )
      end do

    end do

  end if

  return
end