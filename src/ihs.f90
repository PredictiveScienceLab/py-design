subroutine covariance ( dim_num, n, x, average, std, covc )

!*****************************************************************************80
!
!! COVARIANCE does a covariance calculation for IHS solutions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 June 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points to be generated.
!
!    Input, integer ( kind = 4 ) X(DIM_NUM,N), the points.
!
!    Output, real ( kind = 8 ) AVERAGE, the average minimum distance.
!
!    Output, real ( kind = 8 ) STD, the standard deviation of the
!    minimum distances.
!
!    Output, real ( kind = 8 ) COVC, the covariance of the minimum distances.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  real ( kind = 8 ) average
  real ( kind = 8 ) covc
  real ( kind = 8 ) dist
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) mindist(n)
  real ( kind = 8 ), parameter :: r8_huge = 1.0D+30
  real ( kind = 8 ) std
  real ( kind = 8 ) vec(dim_num)
  integer ( kind = 4 ) x(dim_num,n)
!
!  Set up the distance matrix.
!
  do i = 1, n
    mindist(i) = r8_huge
    do j = 1, n
      if ( i /= j ) then
        vec(1:dim_num) = real ( x(1:dim_num,i) - x(1:dim_num,j), kind = 8 )
        dist = sqrt ( dot_product ( vec(1:dim_num), vec(1:dim_num) ) )
        mindist(i) = min ( mindist(i), dist )
      end if
    end do
  end do
!
!  Find the average minimum distance.
!
  average = sum ( mindist(1:n) ) / real ( n, kind = 8 )
!
!  Compute the standard deviation of the distances.
!
  call r8vec_std ( n, mindist, std )
!
!  Compute the covariance.
!
  covc = std / average

  return
end
subroutine ihs ( dim_num, n, duplication, seed, x )

!*****************************************************************************80
!
!! IHS implements the improved distributed hypercube sampling algorithm.
!
!  Discussion:
!
!    N Points in an DIM_NUM dimensional Latin hypercube are to be selected.
!    Each of the coordinate dimensions is discretized to the values
!    1 through N.  The points are to be chosen in such a way that
!    no two points have any coordinate value in common.  This is
!    a standard Latin hypercube requirement, and there are many
!    solutions.
!
!    This algorithm differs in that it tries to pick a solution
!    which has the property that the points are "spread out"
!    as evenly as possible.  It does this by determining an optimal
!    even spacing, and using the DUPLICATION factor to allow it
!    to choose the best of the various options available to it.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Brian Beachkofski, Ramana Grandhi,
!    Improved Distributed Hypercube Sampling,
!    American Institute of Aeronautics and Astronautics Paper 2002-1274.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points to be generated.
!
!    Input, integer ( kind = 4 ) DUPLICATION, the duplication factor.  This must
!    be at least 1.  A value of 5 is reasonable.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) X(DIM_NUM,N), the points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) duplication
  integer ( kind = 4 ) n

  integer ( kind = 4 ) avail(dim_num,n)
  integer ( kind = 4 ) best
  integer ( kind = 4 ) count
  real ( kind = 8 ) dist
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) list(duplication*n)
  real ( kind = 8 ) min_all
  real ( kind = 8 ) min_can
  real ( kind = 8 ) opt
  integer ( kind = 4 ) point(dim_num,duplication*n)
  integer ( kind = 4 ) point_index
  real ( kind = 8 ), parameter :: r8_huge = 1.0D+30
  integer ( kind = 4 ) seed
  real ( kind = 8 ) vec(dim_num)
  integer ( kind = 4 ) x(dim_num,n)
!f2py integer intent(in) :: dim_num
!f2py integer intent(in) :: n
!f2py integer optional,intent(in) :: duplication=5
!f2py integer intent(in) :: seed
!f2py integer intent(out),depend(dim_num,n),dimension(dim_num,n) :: x

  opt = real ( n, kind = 8 ) / &
    ( real ( n, kind = 8 ) )**( 1.0D+00 / real ( dim_num, kind = 8 ) )
!
!  Pick the first point.
!
  call i4vec_uniform ( dim_num, 1, n, seed, x(1:dim_num,n) )
!
!  Initialize AVAIL,
!  and set an entry in a random row of each column of AVAIL to N.
!
  do j = 1, n
    avail(1:dim_num,j) = j
  end do

  do i = 1, dim_num
    avail(i,x(i,n)) = n
  end do
!
!  Main loop:
!  Assign a value to X(1:DIM_NUM,COUNT) for COUNT = N-1 down to 2.
!
  do count = n-1, 2, -1
!
!  Generate valid points.
!
    do i = 1, dim_num

      do k = 1, duplication
        list(count*(k-1)+1:k*count) = avail(i,1:count)
      end do

      do k = count*duplication, 1, -1
        point_index = i4_uniform ( 1, k, seed )
        point(i,k) = list(point_index)
        list(point_index) = list(k)
      end do

    end do
!
!  For each candidate, determine the distance to all the
!  points that have already been selected, and save the minimum value.
!
    min_all = r8_huge
    best = 0

    do k = 1, duplication*count

      min_can = r8_huge

      do j = count+1, n
        vec(1:dim_num) = real ( point(1:dim_num,k) - x(1:dim_num,j), kind = 8 )
        dist = sqrt ( dot_product ( vec(1:dim_num), vec(1:dim_num) ) )
        min_can = min ( min_can, dist )
      end do

      if ( abs ( min_can - opt ) < min_all ) then
        min_all = abs ( min_can - opt )
        best = k
      end if

    end do

    x(1:dim_num,count) = point(1:dim_num,best)
!
!  Having chosen X(*,COUNT), update AVAIL.
!
    do i = 1, dim_num

      do j = 1, n
        if ( avail(i,j) == x(i,count) ) then
          avail(i,j) = avail(i,count)
        end if
      end do

    end do

  end do
!
!  For the last point, there's only one choice.
!
  x(1:dim_num,1) = avail(1:dim_num,1)

  return
end
subroutine r8vec_std ( n, a, std )

!*****************************************************************************80
!
!! R8VEC_STD returns the standard deviation of a real vector.
!
!  Discussion:
!
!    The standard deviation of a vector X of length N is defined as
!
!      mean ( X(1:n) ) = sum ( X(1:n) ) / n
!
!      std ( X(1:n) ) = sqrt ( sum ( ( X(1:n) - mean )^2 ) / ( n - 1 ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!    N should be at least 2.
!
!    Input, real ( kind = 8 ) A(N), the vector.
!
!    Output, real ( kind = 8 ) STD, the standard deviation of the vector.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) mean
  real ( kind = 8 ) std

  if ( n < 2 ) then

    std = 0.0D+00

  else

    mean = sum ( a(1:n) ) / real ( n, kind = 8 )

    std = sum ( ( a(1:n) - mean )**2 )

    std = sqrt ( std / real ( n - 1, kind = 8 ) )

  end if

  return
end
subroutine i4vec_uniform ( n, a, b, seed, x )

!*****************************************************************************80
!
!! I4VEC_UNIFORM returns a scaled pseudorandom I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer ( kind = 4 ) values.
!
!    The pseudorandom numbers should be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) X(N), a vector of numbers between A and B.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value
  integer ( kind = 4 ) x(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
    r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) &
      +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
    value = nint ( r, kind = 4 )

    value = max ( value, min ( a, b ) )
    value = min ( value, max ( a, b ) )

    x(i) = value

  end do

  return
end
