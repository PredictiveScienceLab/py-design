subroutine lambert1 ( n, eta )

!*****************************************************************************80
!
!! LAMBERT1 computes the Lambert sequence in 1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    J P Lambert,
!    Quasi-Random Sequences for Optimization and Numerical Integration,
!    in Numerical Integration,
!    edited by P Keast and G Fairweather,
!    D Reidel, 1987, pages 193-203.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of the sequence
!    to compute.
!
!    Output, real ( kind = 8 ) ETA(1,N), the elements of the sequence.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) eta(1,n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) t
  real ( kind = 8 ) x
!f2py integer intent(in) :: n
!f2py real*8 intent(out),depend(n),dimension(1,n) :: eta

  eta(1,1) = 0.0D+00

  do j = 2, n

    t = 1.0D+00

    do

      t = t / 2.0D+00

      if ( x + t < 1.0D+00 ) then
        exit
      end if

    end do

    x = x + t - 1.0D+00

    if ( x < 0.0D+00 ) then
      x = x + 2.0D+00 * t
    end if

    eta(1,j) = x

  end do

  return
end
subroutine lambert2 ( n, eta )

!*****************************************************************************80
!
!! LAMBERT2 computes the Lambert sequence in 2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    J P Lambert,
!    Quasi-Random Sequences for Optimization and Numerical Integration,
!    in Numerical Integration,
!    edited by P Keast and G Fairweather,
!    D Reidel, 1987, pages 193-203.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of the sequence
!    to compute.
!
!    Output, real ( kind = 8 ) ETA(2,N), the elements of the sequence.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) eta(2,n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!f2py integer intent(in) :: n
!f2py real*8 intent(out),depend(n),dimension(2,n) :: eta

  eta(1:2,1) = 0.0D+00

  do j = 2, n

    t = 1.0D+00

    do

      t = t / 2.0D+00
      u = 1.0D+00 - t

      if ( x < u .or. t <= y ) then
        exit
      end if

    end do

    x = x - u

    if ( x < 0.0D+00 ) then
      x = x + 2.0D+00 * t
      if ( y < t ) then
        y = y + t
      else
        y = y - t
      end if
    end if

    eta(1:2,j) = (/ x, y /)

  end do

  return
end
subroutine lambert3 ( n, eta )

!*****************************************************************************80
!
!! LAMBERT3 computes the Lambert sequence in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    J P Lambert,
!    Quasi-Random Sequences for Optimization and Numerical Integration,
!    in Numerical Integration,
!    edited by P Keast and G Fairweather,
!    D Reidel, 1987, pages 193-203.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of the sequence
!    to compute.
!
!    Output, real ( kind = 8 ) ETA(3,N), the elements of the sequence.
!
  implicit none

  integer n

  real ( kind = 8 ) eta(3,n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z
!f2py integer intent(in) :: n
!f2py real*8 intent(out),depend(n),dimension(3,n) :: eta

  eta(1:3,1) = 0.0D+00

  do j = 2, n

    t = 1.0D+00

    do

      t = t / 2.0D+00
      u = 1.0D+00 - t

      if ( x < u .or. y < u .or. t <= z ) then
        exit
      end if

    end do

    x = x - u
    y = y - u
    z = z - t

    if ( x < 0.0D+00 ) then

      x = x + 2.0D+00 * t
      if ( y < 0.0D+00 ) then
        y = y + 2.0D+00 * t
      end if
      if ( z < 0.0D+00 ) then
        z = z + 2.0D+00 * t
      end if

    else

      if ( y < 0.0D+00 ) then
        if ( z < 0.0D+00 ) then
          y = y + t
          z = z + 2.0D+00 * t
        else if ( y < 0.0D+00 ) then
          y = y + 2.0D+00 * t
          z = z + t
        end if
      else if ( 0.0D+00 <= y ) then
        y = y + t
      end if

    end if

    eta(1:3,j) = (/ x, y, z /)

  end do

  return
end
subroutine lambert4 ( n, eta )

!*****************************************************************************80
!
!! LAMBERT4 computes the Lambert sequence in 4D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    J P Lambert,
!    Quasi-Random Sequences for Optimization and Numerical Integration,
!    in Numerical Integration,
!    edited by P Keast and G Fairweather,
!    D Reidel, 1987, pages 193-203.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of the sequence
!    to compute.
!
!    Output, real ( kind = 8 ) ETA(4,N), the elements of the sequence.
!
  implicit none

  integer n

  real ( kind = 8 ) eta(4,n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z
!f2py integer intent(in) :: n
!f2py real*8 intent(out),depend(n),dimension(4,n) :: eta

  eta(1:4,1) = 0.0D+00

  do j = 2, n

    t = 1.0D+00

    do

      t = t / 2.0D+00
      u = 1.0D+00 - t

      if ( x < u .or. y < u .or. z < u .or. t <= w ) then
        exit
      end if

    end do

    x = x - u
    y = y - u
    z = z - u
    w = w - t

    if ( x < 0.0D+00 ) then

      x = x + 2.0D+00 * t

      if ( y < 0.0D+00 ) then
        y = y + 2.0D+00 * t
      end if
      if ( z < 0.0D+00 ) then
        z = z + 2.0D+00 * t
      end if
      if ( w < 0.0D+00 ) then
        w = w + 2.0D+00 * t
      end if

    else if ( y < 0.0D+00 ) then

      if ( z < 0.0D+00 ) then

        if ( w < 0.0D+00 ) then
          y = y + 2.0D+00 * t
          z = z + 2.0D+00 * t
          w = w + t
        else if ( 0.0D+00 <= w ) then
          y = y + t
          z = z + 2.0D+00 * t
        end if

      else if ( 0.0D+00 <= z ) then

        if ( w < 0.0D+00 ) then
          y = y + 2.0D+00 * t
          z = z + t
          w = w + 2.0D+00 * t
        else if ( 0.0D+00 <= w ) then
          y = y + 2.0D+00 * t
          w = w + t
        end if

      end if

    else if ( z < 0.0D+00 ) then

      if ( w < 0.0D+00 ) then
        z = z + t
        w = w + 2.0D+00 * t
      else if ( 0.0D+00 <= w ) then
        z = z + 2.0D+00 * t
        w = w + t
      end if

    else

      y = y + t

    end if

    eta(1:4,j) = (/ x, y, z, w /)

  end do

  return
end
