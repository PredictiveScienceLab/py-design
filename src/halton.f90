function arc_cosine ( c )

!*****************************************************************************80
!
!! ARC_COSINE computes the arc cosine function, with argument truncation.
!
!  Discussion:
!
!    If you call your system ACOS routine with an input argument that is
!    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
!    This routine truncates arguments outside the range.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) C, the argument.
!
!    Output, real ( kind = 8 ) ARC_COSINE, an angle whose cosine is C.
!
  implicit none

  real ( kind = 8 ) arc_cosine
  real ( kind = 8 ) c
  real ( kind = 8 ) c2

  c2 = c
  c2 = max ( c2, real ( -1.0D+00, kind = 8 ) )
  c2 = min ( c2, real (  1.0D+00, kind = 8 ) )

  arc_cosine = acos ( c2 )

  return
end
function atan4 ( y, x )

!*****************************************************************************80
!
!! ATAN4 computes the inverse tangent of the ratio Y / X.
!
!  Discussion:
!
!    ATAN4 returns an angle whose tangent is ( Y / X ), a job which
!    the built in functions ATAN and ATAN2 already do.
!
!    However:
!
!    * ATAN4 always returns a positive angle, between 0 and 2 PI,
!      while ATAN and ATAN2 return angles in the interval [-PI/2,+PI/2]
!      and [-PI,+PI] respectively;
!
!    * ATAN4 accounts for the signs of X and Y, (as does ATAN2).  The ATAN
!     function by contrast always returns an angle in the first or fourth
!     quadrants.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) Y, X, two quantities which represent the
!    tangent of an angle.  If Y is not zero, then the tangent is (Y/X).
!
!    Output, real ( kind = 8 ) ATAN4, an angle between 0 and 2 * PI,
!    whose tangent is (Y/X), and which lies in the appropriate quadrant
!    so that the signs of its cosine and sine match those of X and Y.
!
  implicit none

  real ( kind = 8 ) abs_x
  real ( kind = 8 ) abs_y
  real ( kind = 8 ) atan4
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta
  real ( kind = 8 ) theta_0
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  Special cases:
!
  if ( x == 0.0D+00 ) then

    if ( 0.0D+00 < y ) then
      theta = pi / real ( 2.0D+00, kind = 8 )
    else if ( y < 0.0D+00 ) then
      theta = real ( 3.0D+00, kind = 8 ) * pi / real ( 2.0D+00, kind = 8 )
    else if ( y == 0.0D+00 ) then
      theta = 0.0D+00
    end if

  else if ( y == 0.0D+00 ) then

    if ( 0.0D+00 < x ) then
      theta = 0.0D+00
    else if ( x < 0.0D+00 ) then
      theta = pi
    end if
!
!  We assume that ATAN2 is correct when both arguments are positive.
!
  else

    abs_y = abs ( y )
    abs_x = abs ( x )

    theta_0 = atan2 ( abs_y, abs_x )

    if ( 0.0D+00 < x .and. 0.0D+00 < y ) then
      theta = theta_0
    else if ( x < 0.0D+00 .and. 0.0D+00 < y ) then
      theta = pi - theta_0
    else if ( x < 0.0D+00 .and. y < 0.0D+00 ) then
      theta = pi + theta_0
    else if ( 0.0D+00 < x .and. y < 0.0D+00 ) then
      theta = real ( 2.0D+00, kind = 8 ) * pi - theta_0
    end if

  end if

  atan4 = theta

  return
end
function halham_leap_check ( dim_num, leap )

!*****************************************************************************80
!
!! HALHAM_LEAP_CHECK checks LEAP for a Halton or Hammersley sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEAP(DIM_NUM), the leap vector.
!
!    Output, logical, HALHAM_LEAP_CHECK, is true if LEAP is legal.
!
  implicit none

  integer ( kind = 4 ) dim_num

  logical halham_leap_check
  integer ( kind = 4 ) leap(dim_num)

  if ( any ( leap(1:dim_num) < 1 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HALHAM_LEAP_CHECK - Fatal error!'
    write ( *, '(a)' ) '  Some entry of LEAP < 1!'
    write ( *, '(a)' ) ' '
    call i4vec_transpose_print ( dim_num, leap, 'LEAP:  ' )
    halham_leap_check = .false.
  else
    halham_leap_check = .true.
  end if

  return
end
function halham_n_check ( n )

!*****************************************************************************80
!
!! HALHAM_N_CHECK checks N for a Halton or Hammersley sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, logical HALHAM_N_CHECK, is true if N is legal.
!
  implicit none

  logical halham_n_check
  integer ( kind = 4 ) n

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HALHAM_N_CHECK - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    write ( *, '(a,i12)' ) '  N = ', n
    halham_n_check = .false.
  else
    halham_n_check = .true.
  end if

  return
end
function halham_dim_num_check ( dim_num )

!*****************************************************************************80
!
!! HALHAM_DIM_NUM_CHECK checks DIM_NUM for a Halton or Hammersley sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, logical HALHAM_DIM_NUM_CHECK, is true if DIM_NUM is legal.
!
  implicit none

  integer ( kind = 4 ) dim_num
  logical halham_dim_num_check

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HALHAM_DIM_NUM_CHECK - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
    halham_dim_num_check = .false.
  else
    halham_dim_num_check = .true.
  end if

  return
end
function halham_seed_check ( dim_num, seed )

!*****************************************************************************80
!
!! HALHAM_SEED_CHECK checks SEED for a Halton or Hammersley sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SEED(DIM_NUM), the seed vector.
!
!    Output, logical, HALHAM_SEED_CHECK, is true if SEED is legal.
!
  implicit none

  integer ( kind = 4 ) dim_num

  logical halham_seed_check
  integer ( kind = 4 ) seed(dim_num)

  if ( any ( seed(1:dim_num) < 0 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HALHAM_SEED_CHECK - Fatal error!'
    write ( *, '(a)' ) '  Some entry of SEED < 0!'
    write ( *, '(a)' ) ' '
    call i4vec_transpose_print ( dim_num, seed, 'SEED:  ' )
    halham_seed_check = .false.
  else
    halham_seed_check = .true.
  end if

  return
end
function halham_step_check ( step )

!*****************************************************************************80
!
!! HALHAM_STEP_CHECK checks STEP for a Halton or Hammersley sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) STEP, the index of the subsequence element.
!
!    Output, logical HALHAM_STEP_CHECK, is true if STEP is legal.
!
  implicit none

  logical halham_step_check
  integer ( kind = 4 ) step

  if ( step < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HALHAM_STEP_CHECK - Fatal error!'
    write ( *, '(a)' ) '  STEP < 0.'
    write ( *, '(a,i12)' ) '  STEP = ', step
    halham_step_check = .false.
  else
    halham_step_check = .true.
  end if

  return
end
subroutine halham_write ( dim_num, n, step, seed, leap, base, r, file_out_name )

!*****************************************************************************80
!
!! HALHAM_WRITE writes a Halton or Hammersley subsequence to a file.
!
!  Discussion:
!
!    The initial lines of the file are comments, which begin with a
!    '#' character.
!
!    Thereafter, each line of the file contains the DIM_NUM-dimensional
!    components of the next entry of the dataset.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 August 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of elements in the subsequence.
!
!    Input, integer ( kind = 4 ) STEP, the index of the subsequence element.
!    0 <= STEP is required.
!
!    Input, integer ( kind = 4 ) SEED(DIM_NUM), the sequence index for STEP = 0.
!
!    Input, integer ( kind = 4 ) LEAP(DIM_NUM), the successive jumps in
!    the sequence.
!
!    Input, integer ( kind = 4 ) BASE(DIM_NUM), the bases.
!
!    Input, real ( kind = 8 ) R(DIM_NUM,N), the points.
!
!    Input, character ( len = * ) FILE_OUT_NAME, the output file name.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  integer ( kind = 4 ) base(dim_num)
  logical, parameter ::  comment = .false.
  character ( len = * ) file_out_name
  integer ( kind = 4 ) file_out_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  integer ( kind = 4 ) leap(dim_num)
  integer ( kind = 4 ) mhi
  integer ( kind = 4 ) mlo
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) step
  character ( len = 40 ) string

  call get_unit ( file_out_unit )

  open ( unit = file_out_unit, file = file_out_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HALHAM_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file:'
    write ( *, '(a)' ) '  "' // trim ( file_out_name ) // '".'
    stop
  end if

  if ( comment ) then

    write ( file_out_unit, '(a)'       ) '#  ' // trim ( file_out_name )
    write ( file_out_unit, '(a)'       ) '#  created by HALHAM_WRITE.F90'
    write ( file_out_unit, '(a)'       ) '#'
    write ( file_out_unit, '(a)'       ) '#'
    write ( file_out_unit, '(a,i12)'    ) '#  DIM_NUM = ', dim_num
    write ( file_out_unit, '(a,i12)'    ) '#  N =    ', n
    write ( file_out_unit, '(a,i12)'    ) '#  STEP = ', step

    do mlo = 1, dim_num, 5
      mhi = min ( mlo + 5 - 1, dim_num )
      if ( mlo == 1 ) then
        write ( file_out_unit, '(a,5i12)' ) '#  SEED = ', seed(mlo:mhi)
      else
        write ( file_out_unit, '(a,5i12)' ) '#         ', seed(mlo:mhi)
      end if
    end do
    do mlo = 1, dim_num, 5
      mhi = min ( mlo + 5 - 1, dim_num )
      if ( mlo == 1 ) then
        write ( file_out_unit, '(a,5i12)' ) '#  LEAP = ', leap(mlo:mhi)
      else
        write ( file_out_unit, '(a,5i12)' ) '#         ', leap(mlo:mhi)
      end if
    end do
    do mlo = 1, dim_num, 5
      mhi = min ( mlo + 5 - 1, dim_num )
      if ( mlo == 1 ) then
        write ( file_out_unit, '(a,5i12)' ) '#  BASE = ', base(mlo:mhi)
      else
        write ( file_out_unit, '(a,5i12)' ) '#         ', base(mlo:mhi)
      end if
    end do
    write ( file_out_unit, '(a,g14.6)' ) '#  EPSILON (unit roundoff ) = ', &
      epsilon ( r(1,1) )
    write ( file_out_unit, '(a)'      ) '#'

  end if

  write ( string, '(a,i3,a)' ) '(', dim_num, '(2x,f10.6))'

  do j = 1, n
    write ( file_out_unit, string ) r(1:dim_num,j)
  end do

  close ( unit = file_out_unit )

  return
end
subroutine halton ( dim_num, r )

!*****************************************************************************80
!
!! HALTON computes the next element in a leaped Halton subsequence.
!
!  Discussion:
!
!    The DIM_NUM-dimensional Halton sequence is really DIM_NUM separate
!    sequences, each generated by a particular base.
!
!    This routine selects elements of a "leaped" subsequence of the
!    Halton sequence.  The subsequence elements are indexed by a
!    quantity called STEP, which starts at 0.  The STEP-th subsequence
!    element is simply element
!
!      SEED(1:DIM_NUM) + STEP * LEAP(1:DIM_NUM)
!
!    of the original Halton sequence.
!
!
!    This routine "hides" a number of input arguments.  To specify these
!    arguments explicitly, use I4_TO_HALTON instead.
!
!    All the arguments have default values.  However, if you want to
!    examine or change them, you may call the appropriate routine first.
!
!    * DIM_NUM, the spatial dimension,
!      Default: DIM_NUM = 1;
!      Required: 1 <= DIM_NUM is required.
!
!    * STEP, the subsequence index.
!      Default: STEP = 0.
!      Required: 0 <= STEP.
!
!    * SEED(1:DIM_NUM), the Halton sequence element corresponding to STEP = 0.
!      Default SEED = (0, 0, ... 0).
!      Required: 0 <= SEED(1:DIM_NUM).
!
!    * LEAP(1:DIM_NUM), the succesive jumps in the Halton sequence.
!      Default: LEAP = (1, 1, ..., 1).
!      Required: 1 <= LEAP(1:DIM_NUM).
!
!    * BASE(1:DIM_NUM), the Halton bases.
!      Default: BASE = (2, 3, 5, 7, 11, ... ).
!      Required: 1 < BASE(1:DIM_NUM).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John Halton,
!    On the efficiency of certain quasi-random sequences of points
!    in evaluating multi-dimensional integrals,
!    Numerische Mathematik,
!    Volume 2, 1960, pages 84-90.
!
!    John Halton and G B Smith,
!    Algorithm 247: Radical-Inverse Quasi-Random Point Sequence,
!    Communications of the ACM,
!    Volume 7, 1964, pages 701-702.
!
!    Ladislav Kocis and William Whiten,
!    Computational Investigations of Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 2, 1997, pages 266-294.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) R(DIM_NUM), the next element of the
!    leaped Halton subsequence.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) base(dim_num)
  integer ( kind = 4 ) leap(dim_num)
  real    ( kind = 8 ) r(dim_num)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) step
  integer ( kind = 4 ) value(1)

  value(1) = dim_num
  call halton_memory ( 'SET', 'DIM_NUM', 1, value )
  call halton_memory ( 'GET', 'STEP', 1, value )
  step = value(1)
  call halton_memory ( 'GET', 'SEED', dim_num, seed )
  call halton_memory ( 'GET', 'LEAP', dim_num, leap )
  call halton_memory ( 'GET', 'BASE', dim_num, base )

  call i4_to_halton ( dim_num, step, seed, leap, base, r )

  value(1) = 1
  call halton_memory ( 'INC', 'STEP', 1, value )

  return
end
subroutine halton_base_get ( base )

!*****************************************************************************80
!
!! HALTON_BASE_GET gets the base vector for a leaped Halton subsequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) BASE(DIM_NUM), the Halton bases.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) base(*)
  integer ( kind = 4 ) value(1)

  call halton_memory ( 'GET', 'DIM_NUM', 1, value )
  dim_num = value(1)

  call halton_memory ( 'GET', 'BASE', dim_num, base )

  return
end
function halton_base_check ( dim_num, base )

!*****************************************************************************80
!
!! HALTON_BASE_CHECK checks BASE for a Halton sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) BASE(DIM_NUM), the bases.
!
!    Output, logical, HALTON_BASE_CHECK, is true if BASE is legal.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) base(dim_num)
  logical halton_base_check

  if ( any ( base(1:dim_num) <= 1 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HALTON_BASE_CHECK - Fatal error!'
    write ( *, '(a)' ) '  Some entry of BASE is <= 1!'
    write ( *, '(a)' ) ' '
    call i4vec_transpose_print ( dim_num, base, 'BASE:  ' )
    halton_base_check = .false.
  else
    halton_base_check = .true.
  end if

  return
end
subroutine halton_base_set ( base )

!*****************************************************************************80
!
!! HALTON_BASE_SET sets the base vector for a leaped Halton subsequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) BASE(DIM_NUM), the Halton bases.
!
  implicit none

  integer ( kind = 4 ) base(*)
  logical halton_base_check
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) value(1)

  call halton_memory ( 'GET', 'DIM_NUM', 1, value )
  dim_num = value(1)

  if ( .not. halton_base_check ( dim_num, base ) ) then
    stop
  end if

  call halton_memory ( 'SET', 'BASE', dim_num, base )

  return
end
subroutine halton_leap_get ( leap )

!*****************************************************************************80
!
!! HALTON_LEAP_GET gets the leap vector for a leaped Halton subsequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) LEAP(DIM_NUM), the successive jumps in
!    the Halton sequence.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) leap(*)
  integer ( kind = 4 ) value(1)

  call halton_memory ( 'GET', 'DIM_NUM', 1, value )
  dim_num = value(1)

  call halton_memory ( 'GET', 'LEAP', dim_num, leap )

  return
end
subroutine halton_leap_set ( leap )

!*****************************************************************************80
!
!! HALTON_LEAP_SET sets the leap vector for a leaped Halton subsequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LEAP(DIM_NUM), the successive jumps in
!    the Halton sequence.
!
  implicit none

  integer ( kind = 4 ) dim_num
  logical halham_leap_check
  integer ( kind = 4 ) leap(*)
  integer ( kind = 4 ) value(1)

  call halton_memory ( 'GET', 'DIM_NUM', 1, value )
  dim_num = value(1)

  if ( .not. halham_leap_check ( dim_num, leap ) ) then
    stop
  end if

  call halton_memory ( 'SET', 'LEAP', dim_num, leap )

  return
end
subroutine halton_memory ( action, name, dim_num, value )

!*****************************************************************************80
!
!! HALTON_MEMORY holds data associated with a leaped Halton subsequence.
!
!  Discussion:
!
!    If you're going to define a new problem, it's important that
!    you set the value of DIM_NUM before setting the values of BASE,
!    LEAP or SEED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, the desired action.
!    'GET' means get the value of a particular quantity.
!    'SET' means set the value of a particular quantity.
!    'INC' means increment the value of a particular quantity.
!          (Only SEED and STEP can be incremented.)
!
!    Input, character ( len = * ) NAME, the name of the quantity.
!    'BASE' means the Halton base vector.
!    'LEAP' means the Halton leap vector.
!    'DIM_NUM' means the spatial dimension.
!    'SEED' means the Halton seed vector.
!    'STEP' means the Halton step.
!
!    Input/output, integer ( kind = 4 ) DIM_NUM, the dimension of the quantity.
!    If ACTION is 'SET' and NAME is 'BASE', then DIM_NUM is input, and
!    is the number of entries in VALUE to be put into BASE.
!
!    Input/output, integer ( kind = 4 ) VALUE(DIM_NUM), contains a value.
!    If ACTION is 'SET', then on input, VALUE contains values to be assigned
!    to the internal variable.
!    If ACTION is 'GET', then on output, VALUE contains the values of
!    the specified internal variable.
!    If ACTION is 'INC', then on input, VALUE contains the increment to
!    be added to the specified internal variable.
!
  implicit none

  character ( len = * ) action
  integer ( kind = 4 ), allocatable, save, dimension ( : ) :: base
  logical, save :: first_call = .true.
  integer ( kind = 4 ) i
  integer ( kind = 4 ), allocatable, save, dimension ( : ) :: leap
  character ( len = * ) name
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ), save :: dim_num_save = 0
  integer ( kind = 4 ) prime
  integer ( kind = 4 ), allocatable, save, dimension ( : ) :: seed
  integer ( kind = 4 ), save :: step = 0
  integer ( kind = 4 ) value(*)

  if ( first_call ) then
    dim_num_save = 1
    allocate ( base(dim_num_save) )
    allocate ( leap(dim_num_save) )
    allocate ( seed(dim_num_save) )
    base(1) = 2
    leap(1) = 1
    seed(1) = 0
    step = 0
    first_call = .false.
  end if
!
!  If this is a SET DIM_NUM call, and the input value of DIM_NUM
!  differs from the internal value, discard all old information.
!
  if ( action(1:1) == 'S' .or. action(1:1) == 's') then
    if ( name == 'DIM_NUM' .or. name == 'dim_num' ) then
      if ( dim_num_save /= value(1) ) then
        deallocate ( base )
        deallocate ( leap )
        deallocate ( seed )
        dim_num_save = value(1)
        allocate ( base(dim_num_save) )
        allocate ( leap(dim_num_save) )
        allocate ( seed(dim_num_save) )
        do i = 1, dim_num_save
          base(i) = prime ( i )
        end do
        leap(1:dim_num_save) = 1
        seed(1:dim_num_save) = 0
      end if
    end if
  end if
!
!  Set
!
  if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == 'BASE' .or. name == 'base' ) then

      if ( dim_num_save /= dim_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HALTON_MEMORY - Fatal error!'
        write ( *, '(a)' ) '  Internal and input values of DIM_NUM disagree'
        write ( *, '(a)' ) '  while setting BASE.'
        stop
      end if

      base(1:dim_num) = value(1:dim_num)

    else if ( name == 'LEAP' .or. name == 'leap' ) then

      if ( dim_num_save /= dim_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HALTON_MEMORY - Fatal error!'
        write ( *, '(a)' ) '  Internal and input values of DIM_NUM disagree'
        write ( *, '(a)' ) '  while setting LEAP.'
        stop
      end if

      leap(1:dim_num) = value(1:dim_num)

    else if ( name == 'DIM_NUM' .or. name == 'dim_num' ) then

      dim_num_save = value(1)

    else if ( name == 'SEED' .or. name == 'seed' ) then

      if ( dim_num_save /= dim_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HALTON_MEMORY - Fatal error!'
        write ( *, '(a)' ) '  Internal and input values of DIM_NUM disagree'
        write ( *, '(a)' ) '  while setting SEED.'
        stop
      end if

      seed(1:dim_num) = value(1:dim_num)

    else if ( name == 'STEP' .or. name == 'step' ) then

      if ( value(1) < 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HALTON_MEMORY - Fatal error!'
        write ( *, '(a)' ) '  Input value of STEP < 0.'
        stop
      end if

      step = value(1)

    end if
!
!  Get
!
  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == 'BASE' .or. name == 'base' ) then

      value(1:dim_num_save) = base(1:dim_num_save)

    else if ( name == 'LEAP' .or. name == 'leap' ) then

      value(1:dim_num_save) = leap(1:dim_num_save)

    else if ( name == 'DIM_NUM' .or. name == 'dim_num' ) then

      value(1) = dim_num_save

    else if ( name == 'SEED' .or. name == 'seed' ) then

      value(1:dim_num_save) = seed(1:dim_num_save)

    else if ( name == 'STEP' .or. name == 'step' ) then

      value(1) = step

    end if
!
!  Increment
!
  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == 'SEED' .or. name == 'seed' ) then
      if ( dim_num == 1 ) then
        seed(1:dim_num_save) = seed(1:dim_num_save) + value(1)
      else
        seed(1:dim_num_save) = seed(1:dim_num_save) + value(1:dim_num_save)
      end if
    else if ( name == 'STEP' .or. name == 'step' ) then
      step = step + value(1)
    end if

  end if

  return
end
subroutine halton_dim_num_get ( dim_num )

!*****************************************************************************80
!
!! HALTON_DIM_NUM_GET: spatial dimension for a leaped Halton subsequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) value(1)

  call halton_memory ( 'GET', 'DIM_NUM', 1, value )
  dim_num = value(1)

  return
end
subroutine halton_dim_num_set ( dim_num )

!*****************************************************************************80
!
!! HALTON_DIM_NUM_SET sets the spatial dimension for leaped Halton subsequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!    1 <= DIM_NUM is required.
!
  implicit none

  integer ( kind = 4 ) dim_num
  logical halham_dim_num_check
  integer ( kind = 4 ) value(1)

  if ( .not. halham_dim_num_check ( dim_num ) ) then
    stop
  end if

  value(1) = dim_num
  call halton_memory ( 'SET', 'DIM_NUM', 1, value )

  return
end
subroutine halton_seed_get ( seed )

!*****************************************************************************80
!
!! HALTON_SEED_GET gets the seed vector for a leaped Halton subsequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) SEED(DIM_NUM), the Halton sequence index
!    corresponding to STEP = 0.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) seed(*)
  integer ( kind = 4 ) value(1)

  call halton_memory ( 'GET', 'DIM_NUM', 1, value )
  dim_num = value(1)

  call halton_memory ( 'GET', 'SEED', dim_num, seed )

  return
end
subroutine halton_seed_set ( seed )

!*****************************************************************************80
!
!! HALTON_SEED_SET sets the seed vector for a leaped Halton subsequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEED(DIM_NUM), the Halton sequence index
!    corresponding to STEP = 0.
!
  implicit none

  logical halham_seed_check
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) seed(*)
  integer ( kind = 4 ) value(1)

  call halton_memory ( 'GET', 'DIM_NUM', 1, value )
  dim_num = value(1)

  if ( .not. halham_seed_check ( dim_num, seed ) ) then
    stop
  end if

  call halton_memory ( 'SET', 'SEED', dim_num, seed )

  return
end
subroutine halton_sequence ( dim_num, n, r )

!*****************************************************************************80
!
!! HALTON_SEQUENCE computes N elements of a leaped Halton subsequence.
!
!  Discussion:
!
!    The DIM_NUM-dimensional Halton sequence is really DIM_NUM separate
!    sequences, each generated by a particular base.
!
!    This routine selects elements of a "leaped" subsequence of the
!    Halton sequence.  The subsequence elements are indexed by a
!    quantity called STEP, which starts at 0.  The STEP-th subsequence
!    element is simply element
!
!      SEED(1:DIM_NUM) + STEP * LEAP(1:DIM_NUM)
!
!    of the original Halton sequence.
!
!
!    This routine "hides" a number of input arguments.  To specify these
!    arguments explicitly, use I4_TO_HALTON_SEQUENCE instead.
!
!    All the arguments have default values.  However, if you want to
!    examine or change them, you may call the appropriate routine first.
!
!    The arguments that the user may set include:
!
!    * DIM_NUM, the spatial dimension,
!      Default: DIM_NUM = 1;
!      Required: 1 <= DIM_NUM is required.
!
!    * STEP, the subsequence index.
!      Default: STEP = 0.
!      Required: 0 <= STEP.
!
!    * SEED(1:DIM_NUM), the Halton sequence element corresponding to STEP = 0.
!      Default SEED = (0, 0, ... 0).
!      Required: 0 <= SEED(1:DIM_NUM).
!
!    * LEAP(1:DIM_NUM), the succesive jumps in the Halton sequence.
!      Default: LEAP = (1, 1, ..., 1).
!      Required: 1 <= LEAP(1:DIM_NUM).
!
!    * BASE(1:DIM_NUM), the Halton bases.
!      Default: BASE = (2, 3, 5, 7, 11, ... ).
!      Required: 1 < BASE(1:DIM_NUM).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John Halton,
!    On the efficiency of certain quasi-random sequences of points
!    in evaluating multi-dimensional integrals,
!    Numerische Mathematik,
!    Volume 2, 1960, pages 84-90.
!
!    John Halton and G B Smith,
!    Algorithm 247: Radical-Inverse Quasi-Random Point Sequence,
!    Communications of the ACM,
!    Volume 7, 1964, pages 701-702.
!
!    Ladislav Kocis and William Whiten,
!    Computational Investigations of Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 2, 1997, pages 266-294.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of elements desired.
!
!    Output, real ( kind = 8 ) R(DIM_NUM,N), the next N elements of the
!    leaped Halton subsequence.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  integer ( kind = 4 ) base(dim_num)
  integer ( kind = 4 ) leap(dim_num)
  real    ( kind = 8 ) r(dim_num,n)
!f2py integer intent(in) :: dim_num
!f2py integer intent(in) :: n
!f2py real*8 intent(out),depend(dim_num,n),dimension(dim_num,n) :: r
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) step
  integer ( kind = 4 ) value(1)

  value(1) = dim_num
  call halton_memory ( 'SET', 'DIM_NUM', 1, value )
  call halton_memory ( 'GET', 'STEP', 1, value )
  step = value(1)
  call halton_memory ( 'GET', 'SEED', dim_num, seed )
  call halton_memory ( 'GET', 'LEAP', dim_num, leap )
  call halton_memory ( 'GET', 'BASE', dim_num, base )

  call i4_to_halton_sequence ( dim_num, n, step, seed, leap, base, r )

  value(1) = n
  call halton_memory ( 'INC', 'STEP', 1, value )

  return
end
subroutine halton_step_get ( step )

!*****************************************************************************80
!
!! HALTON_STEP_GET gets the "step" for a leaped Halton subsequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) STEP, the index of the subsequence element.
!
  implicit none

  integer ( kind = 4 ) step
  integer ( kind = 4 ) value(1)

  call halton_memory ( 'GET', 'STEP', 1, value )
  step = value(1)

  return
end
subroutine halton_step_set ( step )

!*****************************************************************************80
!
!! HALTON_STEP_SET sets the "step" for a leaped Halton subsequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) STEP, the index of the subsequence element.
!    0 <= STEP is required.
!
  implicit none

  logical halham_step_check
  integer ( kind = 4 ) step
  integer ( kind = 4 ) value(1)

  if ( .not. halham_step_check ( step ) ) then
    stop
  end if

  value(1) = step
  call halton_memory ( 'SET', 'STEP', 1, value )

  return
end
subroutine i4_to_halton ( dim_num, step, seed, leap, base, r )

!*****************************************************************************80
!
!! I4_TO_HALTON computes one element of a leaped Halton subsequence.
!
!  Discussion:
!
!    The DIM_NUM-dimensional Halton sequence is really DIM_NUM separate
!    sequences, each generated by a particular base.
!
!    This routine selects elements of a "leaped" subsequence of the
!    Halton sequence.  The subsequence elements are indexed by a
!    quantity called STEP, which starts at 0.  The STEP-th subsequence
!    element is simply element
!
!      SEED(1:DIM_NUM) + STEP * LEAP(1:DIM_NUM)
!
!    of the original Halton sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John Halton,
!    On the efficiency of certain quasi-random sequences of points
!    in evaluating multi-dimensional integrals,
!    Numerische Mathematik,
!    Volume 2, 1960, pages 84-90.
!
!    John Halton and G B Smith,
!    Algorithm 247: Radical-Inverse Quasi-Random Point Sequence,
!    Communications of the ACM,
!    Volume 7, 1964, pages 701-702.
!
!    Ladislav Kocis and William Whiten,
!    Computational Investigations of Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 2, 1997, pages 266-294.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!    1 <= DIM_NUM is required.
!
!    Input, integer ( kind = 4 ) STEP, the index of the subsequence element.
!    0 <= STEP is required.
!
!    Input, integer ( kind = 4 ) SEED(DIM_NUM), the Halton sequence index
!    corresponding to STEP = 0.
!    0 <= SEED(1:DIM_NUM) is required.
!
!    Input, integer ( kind = 4 ) LEAP(DIM_NUM), the successive jumps in
!    the Halton sequence.
!    1 <= LEAP(1:DIM_NUM) is required.
!
!    Input, integer ( kind = 4 ) BASE(DIM_NUM), the Halton bases.
!    1 < BASE(1:DIM_NUM) is required.
!
!    Output, real ( kind = 8 ) R(DIM_NUM), the STEP-th element of the leaped
!    Halton subsequence.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) base(dim_num)
  real    ( kind = 8 ) base_inv
  integer ( kind = 4 ) digit
  logical halham_leap_check
  logical halham_dim_num_check
  logical halham_seed_check
  logical halham_step_check
  logical halton_base_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) leap(dim_num)
  real    ( kind = 8 ) r(dim_num)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) seed2
  integer ( kind = 4 ) step
!
!  Check the input.
!
  if ( .not. halham_dim_num_check ( dim_num ) ) then
    stop
  end if

  if ( .not. halham_step_check ( step ) ) then
    stop
  end if

  if ( .not. halham_seed_check ( dim_num, seed ) ) then
    stop
  end if

  if ( .not. halham_leap_check ( dim_num, leap ) ) then
    stop
  end if

  if ( .not. halton_base_check ( dim_num, base ) ) then
    stop
  end if
!
!  Calculate the data.
!
  do i = 1, dim_num

    seed2 = seed(i) + step * leap(i)

    r(i) = 0.0D+00

    base_inv = real ( 1.0D+00, kind = 8 ) / real ( base(i), kind = 8 )

    do while ( seed2 /= 0 )
      digit = mod ( seed2, base(i) )
      r(i) = r(i) + real ( digit, kind = 8 ) * base_inv
      base_inv = base_inv / real ( base(i), kind = 8 )
      seed2 = seed2 / base(i)
    end do

  end do

  return
end
subroutine i4_to_halton_sequence ( dim_num, n, step, seed, leap, base, r )

!*****************************************************************************80
!
!! I4_TO_HALTON_SEQUENCE computes N elements of a leaped Halton subsequence.
!
!  Discussion:
!
!    The DIM_NUM-dimensional Halton sequence is really DIM_NUM separate
!    sequences, each generated by a particular base.
!
!    This routine selects elements of a "leaped" subsequence of the
!    Halton sequence.  The subsequence elements are indexed by a
!    quantity called STEP, which starts at 0.  The STEP-th subsequence
!    element is simply element
!
!      SEED(1:DIM_NUM) + STEP * LEAP(1:DIM_NUM)
!
!    of the original Halton sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John Halton,
!    On the efficiency of certain quasi-random sequences of points
!    in evaluating multi-dimensional integrals,
!    Numerische Mathematik,
!    Volume 2, 1960, pages 84-90.
!
!    John Halton and G B Smith,
!    Algorithm 247: Radical-Inverse Quasi-Random Point Sequence,
!    Communications of the ACM,
!    Volume 7, 1964, pages 701-702.
!
!    Ladislav Kocis and William Whiten,
!    Computational Investigations of Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 2, 1997, pages 266-294.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!    1 <= DIM_NUM is required.
!
!    Input, integer ( kind = 4 ) N, the number of elements of the sequence.
!
!    Input, integer ( kind = 4 ) STEP, the index of the subsequence element.
!    0 <= STEP is required.
!
!    Input, integer ( kind = 4 ) SEED(DIM_NUM), the Halton sequence index
!    corresponding to STEP = 0.
!
!    Input, integer ( kind = 4 ) LEAP(DIM_NUM), the succesive jumps in the
!    Halton sequence.
!
!    Input, integer ( kind = 4 ) BASE(DIM_NUM), the Halton bases.
!
!    Output, real ( kind = 8 ) R(DIM_NUM,N), the next N elements of the
!    leaped Halton subsequence, beginning with element STEP.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  integer ( kind = 4 ) base(dim_num)
  real    ( kind = 8 ) base_inv
  integer ( kind = 4 ) digit(n)
  logical halham_leap_check
  logical halham_n_check
  logical halham_dim_num_check
  logical halham_seed_check
  logical halham_step_check
  logical halton_base_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) leap(dim_num)
  real    ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) seed2(n)
  integer ( kind = 4 ) step
!
!  Check the input.
!
  if ( .not. halham_dim_num_check ( dim_num ) ) then
    stop
  end if

  if ( .not. halham_n_check ( n ) ) then
    stop
  end if

  if ( .not. halham_step_check ( step ) ) then
    stop
  end if

  if ( .not. halham_seed_check ( dim_num, seed ) ) then
    stop
  end if

  if ( .not. halham_leap_check ( dim_num, leap ) ) then
    stop
  end if

  if ( .not. halton_base_check ( dim_num, base ) ) then
    stop
  end if
!
!  Calculate the data.
!
  r(1:dim_num,1:n) = 0.0D+00

  do i = 1, dim_num

    do j = 1, n
      seed2(j) = seed(i) + ( step + j - 1 ) * leap(i)
    end do

    base_inv = real ( 1.0D+00, kind = 8 ) / real ( base(i), kind = 8 )

    do while ( any ( seed2(1:n) /= 0 ) )
      digit(1:n) = mod ( seed2(1:n), base(i) )
      r(i,1:n) = r(i,1:n) + real ( digit(1:n), kind = 8 ) * base_inv
      base_inv = base_inv / real ( base(i), kind = 8 )
      seed2(1:n) = seed2(1:n) / base(i)
    end do

  end do

  return
end
subroutine i4vec_transpose_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_TRANSPOSE_PRINT prints an I4VEC "transposed".
!
!  Example:
!
!    A = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 /)
!    TITLE = 'My vector:  '
!
!    My vector:      1    2    3    4    5
!                    6    7    8    9   10
!                   11
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  character ( len = 11 ) string
  character ( len = * ) title
  integer ( kind = 4 ) title_len

  if ( 0 < len ( title ) ) then

    title_len = len ( title )

    write ( string, '(a,i3,a)' ) '(', title_len, 'x,5i12)'

    do ilo = 1, n, 5
      ihi = min ( ilo + 5 - 1, n )
      if ( ilo == 1 ) then
        write ( *, '(a, 5i12)' ) title, a(ilo:ihi)
      else
        write ( *, string      )        a(ilo:ihi)
      end if
    end do

  else

    do ilo = 1, n, 5
      ihi = min ( ilo + 5 - 1, n )
      write ( *, '(5i12)' ) a(ilo:ihi)
    end do

  end if

  return
end
subroutine u1_to_sphere_unit_2d ( u, x )

!*****************************************************************************80
!
!! U1_TO_SPHERE_UNIT_2D maps a point in the unit interval to the unit circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) U, a point in the unit interval.
!
!    Output, real ( kind = 8 ) X(2), the corresponding point on the circle.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) u
  real ( kind = 8 ) x(2)

  angle = real ( 2.0D+00, kind = 8 ) * pi * u

  x(1) = cos ( angle )
  x(2) = sin ( angle )

  return
end
subroutine u2_to_ball_unit_2d ( u, x )

!*****************************************************************************80
!
!! U2_TO_BALL_UNIT_2D maps points from the unit box to the unit ball in 2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) U(2), a point in the unit square.
!
!    Output, real ( kind = 8 ) X(2), the corresponding point in the
!    unit ball.
!
  implicit none

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) theta
  real ( kind = 8 ) u(2)
  real ( kind = 8 ) x(2)

  r = sqrt ( u(1) )
  theta = real ( 2.0D+00, kind = 8 ) * pi * u(2)

  x(1) = r * cos ( theta )
  x(2) = r * sin ( theta )

  return
end
subroutine u2_to_sphere_unit_3d ( u, x )

!*****************************************************************************80
!
!! U2_TO_SPHERE_UNIT_3D maps a point in the unit box onto the unit sphere in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) U(2), the point in the unit box.
!
!    Output, real ( kind = 8 ) X(3), the corresponding point on the unit sphere.
!
  implicit none

  real ( kind = 8 ) arc_cosine
  real ( kind = 8 ) phi
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta
  real ( kind = 8 ) u(2)
  real ( kind = 8 ) vdot
  real ( kind = 8 ) x(3)
!
!  Pick a uniformly random VDOT, which must be between -1 and 1.
!  This represents the dot product of the random vector with the Z unit vector.
!
!  Note: this works because the surface area of the sphere between
!  Z and Z + dZ is independent of Z.  So choosing Z uniformly chooses
!  a patch of area uniformly.
!
  vdot = real ( 2.0D+00, kind = 8 ) * u(1) - real ( 1.0D+00, kind = 8 )

  phi = arc_cosine ( vdot )
!
!  Pick a uniformly random rotation between 0 and 2 Pi around the
!  axis of the Z vector.
!
  theta = real ( 2.0D+00, kind = 8 ) * pi * u(2)

  x(1) = cos ( theta ) * sin ( phi )
  x(2) = sin ( theta ) * sin ( phi )
  x(3) = cos ( phi )

  return
end
subroutine u3_to_ball_unit_3d ( u, x )

!*****************************************************************************80
!
!! U3_TO_BALL_UNIT_3D maps points from the unit box to the unit ball in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) U(3), a point in the unit box in 3D.
!
!    Output, real ( kind = 8 ) X(3), the corresponding point in the
!    unit ball in 3D.
!
  implicit none

  real ( kind = 8 ) arc_cosine
  real ( kind = 8 ) phi
  real ( kind = 8 ) r
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta
  real ( kind = 8 ) u(3)
  real ( kind = 8 ) vdot
  real ( kind = 8 ) x(3)
!
!  Pick a uniformly random VDOT, which must be between -1 and 1.
!  This represents the dot product of the random vector with the Z unit vector.
!
!  Note: this works because the surface area of the sphere between
!  Z and Z + dZ is independent of Z.  So choosing Z uniformly chooses
!  a patch of area uniformly.
!
  vdot = real ( 2.0D+00, kind = 8 ) * u(1) - real ( 1.0D+00, kind = 8 )

  phi = arc_cosine ( vdot )
!
!  Pick a uniformly random rotation between 0 and 2 Pi around the
!  axis of the Z vector.
!
  theta = real ( 2.0D+00, kind = 8 ) * pi * u(2)
!
!  Pick a random radius R.
!
  r = u(3)**( real ( 1.0D+00, kind = 8 ) / real ( 3.0D+00, kind = 8 ) )

  x(1) = r * cos ( theta ) * sin ( phi )
  x(2) = r * sin ( theta ) * sin ( phi )
  x(3) = r * cos ( phi )

  return
end

subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
