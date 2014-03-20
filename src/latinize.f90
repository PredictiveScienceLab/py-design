subroutine latinize ( m, n, table )

!*****************************************************************************80
!
!! R8MAT_LATINIZE "Latinizes" an R8MAT.
!
!  Discussion:
!
!    On output, each row of the table will have the properties that:
!    1) the minimum and maximum row values are the same as on input;
!    2) the row contains N evenly spaced values between the
!       minimum and maximum;
!    3) in each row, the elements retain their ordering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of columns.
!
!    Input/output, real ( kind = 8 ) TABLE(M,N).  On input, the dataset to
!    be "Latinized".  On output, the Latinized dataset.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) table(m,n)
!f2py intent(in,out), dimension(m, n) :: table
!f2py intent(in), depend(table), check(m <= shape(table, 0)) :: m = shape(table, 0)
!f2py intent(in), depend(table), check(n <= shape(table, 1)) :: n = shape(table, 1)
  real ( kind = 8 ) v_max
  real ( kind = 8 ) v_min

  if ( n <= 2 ) then
    return
  end if

  do i = 1, m

    v_min = minval ( table(i,1:n) )
    v_max = maxval ( table(i,1:n) )

    call r8vec_sort_heap_index_a ( n, table(i,1:n), indx )

    do j = 1, n
      table(i,indx(j)) =  ( real ( n - j,     kind = 8 ) * v_min   &
                          + real (     j - 1, kind = 8 ) * v_max ) &
                          / real ( n     - 1, kind = 8 )
    end do

  end do

  return
end
subroutine r8vec_sort_heap_index_a ( n, a, indx )

!*****************************************************************************80
!
!! R8VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of a real vector.
!
!  Discussion:
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(INDX(I)), I = 1 to N is sorted,
!
!    or explicitly, by the call
!
!      call R8VEC_PERMUTE ( N, A, INDX )
!
!    after which A(I), I = 1 to N is sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), an array to be index-sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), the sort index.  The
!    I-th element of the sorted array is A(INDX(I)).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) aval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l

  if ( n < 1 ) then
    return
  end if

  if ( n == 1 ) then
    indx(1) = 1
    return
  end if

  do i = 1, n
    indx(i) = i
  end do

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      aval = a(indxt)

    else

      indxt = indx(ir)
      aval = a(indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if ( a(indx(j)) < a(indx(j+1)) ) then
          j = j + 1
        end if
      end if

      if ( aval < a(indx(j)) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end