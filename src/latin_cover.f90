function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of I4 division.
!
!  Discussion:
!
!    If
!      NREM = I4_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!        I     J     MOD I4_MODP    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number to be divided.
!
!    Input, integer ( kind = 4 ) J, the number that divides I.
!
!    Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder when I is
!    divided by J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) j
  integer ( kind = 4 ) value
!f2py intent(in) i
!f2py intent(in) j
!f2py intent(out) i4_modp

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
    stop
  end if

  value = mod ( i, j )

  if ( value < 0 ) then
    value = value + abs ( j )
  end if

  i4_modp = value

  return
end
function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! I4_WRAP forces an I4 to lie between given limits by wrapping.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    There appears to be a bug in the GFORTRAN compiler which can lead to
!    erroneous results when the first argument of I4_WRAP is an expression.
!    In particular:
!
!    do i = 1, 3
!      if ( test ) then
!        i4 = i4_wrap ( i + 1, 1, 3 )
!      end if
!    end do
!
!    was, when I = 3, returning I4 = 3.  So I had to replace this with
!
!    do i = 1, 3
!      if ( test ) then
!        i4 = i + 1
!        i4 = i4_wrap ( i4, 1, 3 )
!      end if
!    end do
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  Value
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, a value.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds.
!
!    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of the value.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) value
  integer ( kind = 4 ) wide
!f2py intent(in) ival
!f2py intent(in) ilo
!f2py intent(in) ihi
!f2py intent(out) i4_wrap

  jlo = min ( ilo, ihi )
  jhi = max ( ilo, ihi )

  wide = jhi - jlo + 1

  if ( wide == 1 ) then
    value = jlo
  else
    value = jlo + i4_modp ( ival - jlo, wide )
  end if

  i4_wrap = value

  return
end
subroutine latin_cover ( n, p, a )

!*****************************************************************************80
!
!! LATIN_COVER returns a 2D Latin Square Covering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, integer ( kind = 4 ) P(N), a permutation which describes the
!    first Latin square.
!
!    Output, integer ( kind = 4 ) A(N,N), the Latin cover.  A(I,J) = K
!    means that (I,J) is one element of the K-th Latin square.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) k
  integer ( kind = 4 ) p(n)
!f2py intent(in) p
!f2py optional, intent(in), check(n<=len(p)) :: n = len(p)
!f2py intent(out), depend(n), dimensions(n,n) :: a
  call perm_check ( n, p )

  do i = 1, n
    do k = 1, n
      ik = i4_wrap ( i + k - 1, 1, n )
      a(i,p(ik)) = k
    end do
  end do

  return
end
subroutine latin_cover_2d ( n, p, a )

!*****************************************************************************80
!
!! LATIN_COVER_2D returns a 2D Latin Square Covering.
!
!  Discussion:
!
!    This procedure has a chance of being extended to M dimensions.
!
!    A basic solution is computed, and the user is permitted to permute
!    both the I and J coordinates.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, integer ( kind = 4 ) P(2,N), permutations to be applied
!    to the spatial dimensions.
!
!    Output, integer ( kind = 4 ) A(N,N), the Latin cover.  A(I,J) = K
!    means that (I,J) is one element of the K-th Latin square.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n,n)
  integer ( kind = 4 ) b(n,n)
  integer ( kind = 4 ) :: base = 1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) p(2,n)

  call perm_check ( n, p(1,1:n) )
  call perm_check ( n, p(2,1:n) )
!
!  Set up the basic solution.
!
  do i = 1, n
    do j = 1, n
      a(i,j) = i4_wrap ( i - j + base, 0 + base, n - 1 + base )
    end do
  end do
!
!  Apply permutation to dimension I.
!
  do i = 1, n
    b(p(1,i),1:n) = a(i,1:n)
  end do
!
!  Apply permutation to dimension J.
!
  do j = 1, n
    a(1:n,p(2,j)) = b(1:n,j)
  end do

  return
end
subroutine latin_cover_3d ( n, p, a )

!*****************************************************************************80
!
!! LATIN_COVER_3D returns a 3D Latin Square Covering.
!
!  Discussion:
!
!    A basic solution is computed, and the user is permitted to permute
!    I, J and K coordinates.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, integer ( kind = 4 ) P(3,N), permutations to be applied
!    to the spatial dimensions.
!
!    Output, integer ( kind = 4 ) A(N,N,N), the Latin cover.  A(I,J,K) = L
!    means that (I,J,K) is one element of the L-th Latin square.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n,n,n)
  integer ( kind = 4 ) b(n,n,n)
  integer ( kind = 4 ) :: base = 1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) k
  integer ( kind = 4 ) p(3,n)

  call perm_check ( n, p(1,1:n) )
  call perm_check ( n, p(2,1:n) )
  call perm_check ( n, p(3,1:n) )
!
!  Set up the basic solution.
!
  do i = 1, n
    do j = 1, n
      do k = 1, n
        ik = i4_wrap ( i + 1 - k, 1, n )
        jk = i4_wrap ( j + 1 - k, 1, n )
        b(i,j,k) = ik + ( jk - 1 ) * n
      end do
    end do
  end do
!
!  Apply permutation to dimension I.
!
  do i = 1, n
    a(p(1,i),1:n,1:n) = b(i,1:n,1:n)
  end do
!
!  Apply permutation to dimension J.
!
  do j = 1, n
    b(1:n,p(2,j),1:n) = a(1:n,j,1:n)
  end do
!
!  Apply permutation to dimension K.
!
  do k = 1, n
    a(1:n,1:n,p(3,k)) = b(1:n,1:n,k)
  end do

  return
end
subroutine perm_check ( n, p )

!*****************************************************************************80
!
!! PERM_CHECK checks that a vector represents a permutation.
!
!  Discussion:
!
!    The routine verifies that each of the integers from 1
!    to N occurs among the N entries of the permutation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries.
!
!    Input, integer ( kind = 4 ) P(N), the permutation, in standard index form.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ifind
  integer ( kind = 4 ) iseek
  integer ( kind = 4 ) p(n)

  ierror = 0

  do iseek = 1, n

    ierror = iseek

    do ifind = 1, n
      if ( p(ifind) == iseek ) then
        ierror = 0
        exit
      end if
    end do

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PERM_CHECK - Fatal error!'
      write ( *, '(a)' ) '  The input array does not represent'
      write ( *, '(a)' ) '  a proper permutation.  In particular, the'
      write ( *, '(a,i8)' ) '  array is missing the value ', ierror
      stop
    end if

  end do

  return
end