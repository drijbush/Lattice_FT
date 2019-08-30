Module distributed_matrix_module

  !! This module implements the operations on both real and complex distributed matrices
  !  Once compilers become more mature submodule would be a good way to keep this one under control
  !  once for real, one for complex

  Use numbers_module       , Only : wp
  Use matrix_mapping_module, Only : matrix_mapping 
  
  Implicit None

  Integer, Parameter :: distributed_matrix_INVALID = -1
  Integer, Parameter :: distributed_matrix_NOT_ME  = -2

  Type, Abstract, Public :: distributed_matrix
     !! An abstract for the base type. This deal with dimensions, mappings and transposes
     !! but contains no data
     Type( matrix_mapping )              , Private :: matrix_map             !! The mapping of the matrix onto the processes
     Integer, Dimension( : ), Allocatable, Private :: global_to_local_rows   !! Map the global row index to a local one
     Integer, Dimension( : ), Allocatable, Private :: global_to_local_cols   !! Map the global column index to a local one
     Integer, Dimension( : ), Allocatable, Private :: local_to_global_rows   !! Map the local row index to a global one
     Integer, Dimension( : ), Allocatable, Private :: local_to_global_cols   !! Map the local column index to a global one
     Logical                             , Private :: daggered = .False.     !! If true use the matrix in daggered form
   Contains
     ! Public methods that are NOT overridden
     Procedure, Public :: get_maps               => matrix_get_maps            !! Get all the mapping arrays
     Procedure, Public :: global_to_local        => matrix_global_to_local     !! Get an array for mapping global indices to local  ones
     Procedure, Public :: local_to_global        => matrix_local_to_global     !! Get an array for mapping local  indices to global ones
     Procedure, Public :: size                   => matrix_size                !! Get the dimensions of the matrix
     Procedure, Public :: get_comm               => matrix_communicator        !! Get the communicator containing the processes holding the matrix
     Generic  , Public :: Operator( * )          => multiply                   !! Multiply two matrices together
     Generic  , Public :: Operator( * )          => rscal_multiply             !! Pre -scale by a real scalar
     Generic  , Public :: Operator( * )          => multiply_rscal             !! Post-scale by a real scalar
     Generic  , Public :: Operator( + )          => plus                       !! Unary plus operation
     Generic  , Public :: Operator( + )          => add                        !! Add two matrices together
     Generic  , Public :: Operator( + )          => add_diagonal               !! Add a general matrix to a diagonal matrix
     Generic  , Public :: Operator( + )          => diagonal_add               !! Add a general matrix to a diagonal matrix
     Generic  , Public :: Operator( - )          => minus                      !! Unary minus operation
     Generic  , Public :: Operator( - )          => subtract                   !! Subtract two matrices 
     Generic  , Public :: Operator( - )          => subtract_diagonal          !! Subtract a diagonal matrix from a generla matrix
     Generic  , Public :: Operator( - )          => diagonal_subtract          !! Subtract a general matrix from a diagonal matrix
     Generic  , Public :: Operator( .Dagger.   ) => matrix_dagger              !! Apply the dagger operator to the matrix
     Generic  , Public :: Operator( .Choleski. ) => choleski                   !! choleski decompose a matrix
     Generic  , Public :: Operator( .Trinv.    ) => tr_inv                     !! invert a traingular matrix
     Generic  , Public :: set_by_global          => set_global_real, set_global_complex !! Set a matrix using global indexing
     Generic  , Public :: get_by_global          => get_global_real, get_global_complex !! Get from a matrix using global indexing
     ! Public methods that are overridden
     Procedure( create     ), Deferred, Public :: create                     !! Create storage for the data of the matrix 
     Procedure( local_size ), Deferred, Public :: local_size                 !! Get the dimensions of the local part of the matrix
     Procedure( remap_op   ), Deferred, Public :: remap                      !! Remap the data to another distribution
     Procedure( diag_op    ), Deferred, Public :: diag                       !! Diagonalise the (assumed Hermitian) matrix
     Procedure( extract_op ), Deferred, Public :: extract                    !! Extract a patch of one matrix into another matrix
     ! Private implementations
     Procedure,                                             Private :: matrix_dagger           !! Apply the dagger operator to the matrix
     Procedure(  set_global_real    ), Deferred,            Private :: set_global_real         !! Set values with a real    array using global indexing
     Procedure(  set_global_complex ), Deferred,            Private :: set_global_complex      !! Set values with a complex array using global indexing
     Procedure(  get_global_real    ), Deferred,            Private :: get_global_real         !! Get values from a real    array using global indexing
     Procedure(  get_global_complex ), Deferred,            Private :: get_global_complex      !! Get values from a complex array using global indexing
     Procedure(           binary_op ), Deferred,            Private :: multiply
     Procedure(      real_binary_op ), Deferred, Pass( B ), Private :: real_multiply
     Procedure(   complex_binary_op ), Deferred, Pass( B ), Private :: complex_multiply
     Procedure(        pre_rscal_op ), Deferred, Pass( A ), Private :: rscal_multiply
     Procedure(       post_rscal_op ), Deferred,            Private :: multiply_rscal
     Procedure(            unary_op ), Deferred,            Private :: plus
     Procedure(           binary_op ), Deferred,            Private :: add
     Procedure(    post_diagonal_op ), Deferred,            Private :: add_diagonal
     Procedure(     pre_diagonal_op ), Deferred, Pass( A ), Private :: diagonal_add
     Procedure(      real_binary_op ), Deferred, Pass( B ), Private :: real_add
     Procedure(   complex_binary_op ), Deferred, Pass( B ), Private :: complex_add
     Procedure(            unary_op ), Deferred,            Private :: minus
     Procedure(           binary_op ), Deferred,            Private :: subtract
     Procedure(      real_binary_op ), Deferred, Pass( B ), Private :: real_subtract
     Procedure(   complex_binary_op ), Deferred, Pass( B ), Private :: complex_subtract
     Procedure(    post_diagonal_op ), Deferred,            Private :: subtract_diagonal
     Procedure(     pre_diagonal_op ), Deferred, Pass( A ), Private :: diagonal_subtract
     Procedure(        real_diag_op ), Deferred, Pass( Q ), Private :: real_diag
     Procedure(     complex_diag_op ), Deferred, Pass( Q ), Private :: complex_diag
     Procedure(            unary_op ), Deferred,            Private :: choleski
     Procedure(            unary_op ), Deferred,            Private :: tr_inv
     Procedure(       real_remap_op ), Deferred, Pass( B ), Private :: real_remap
     Procedure(    complex_remap_op ), Deferred, Pass( B ), Private :: complex_remap
  End type distributed_matrix

  Type, Extends( distributed_matrix ), Public :: real_distributed_matrix
     !! An instance of a distributed matrix that holds real data
     Real( wp ), Dimension( :, : ), Allocatable, Private    :: data          
   Contains
     ! Public methods
     Procedure, Public :: create        => matrix_create_real              !! Create storage for the data of the matrix 
     Procedure, Public :: local_size    => matrix_local_size_real          !! Get the dimensions of the local part of the matrix
     Procedure, Public :: remap         => real_remap                      !! Remap the data to another distribution
     Procedure, Public :: diag          => real_diag                       !! Diagonalise the (assumed symmetric) matrix
     Procedure, Public :: extract       => extract_real
     ! Private implementations
     Procedure,            Private :: set_global_real    => real_matrix_set_global_real
     Procedure,            Private :: set_global_complex => real_matrix_set_global_complex
     Procedure,            Private :: get_global_real    => real_matrix_get_global_real
     Procedure,            Private :: get_global_complex => real_matrix_get_global_complex
     Procedure,            Private :: multiply           => real_multiply
     Procedure, Pass( B ), Private :: real_multiply      => real_multiply_real
     Procedure, Pass( B ), Private :: complex_multiply   => complex_multiply_real
     Procedure, Pass( A ), Private :: rscal_multiply     => rscal_multiply_real
     Procedure,            Private :: multiply_rscal     => real_multiply_rscal
     Procedure,            Private :: plus               => plus_real
     Procedure,            Private :: add                => real_add
     Procedure, Pass( B ), Private :: real_add           => real_add_real
     Procedure, Pass( B ), Private :: complex_add        => complex_add_real
     Procedure,            Private :: add_diagonal       => real_add_diagonal
     Procedure, Pass( A ), Private :: diagonal_add       => diagonal_add_real
     Procedure,            Private :: minus              => minus_real
     Procedure,            Private :: subtract           => real_subtract
     Procedure, Pass( B ), Private :: real_subtract      => real_subtract_real
     Procedure, Pass( B ), Private :: complex_subtract   => complex_subtract_real
     Procedure,            Private :: subtract_diagonal  => real_subtract_diagonal
     Procedure, Pass( A ), Private :: diagonal_subtract  => diagonal_subtract_real
     Procedure, Pass( Q ), Private :: real_diag          => real_diag_real
     Procedure, Pass( Q ), Private :: complex_diag       => complex_diag_real
     Procedure           , Private :: choleski           => choleski_real
     Procedure           , Private :: tr_inv             => tr_inv_real
     Procedure, Pass( B ), Private :: real_remap         => real_remap_real
     Procedure, Pass( B ), Private :: complex_remap      => complex_remap_real
  End type real_distributed_matrix

  Type, Extends( distributed_matrix ), Public :: complex_distributed_matrix
     !! An instance of a distributed matrix that holds complex data
     Complex( wp ), Dimension( :, : ), Allocatable, Private :: data          
   Contains
     ! Public methods
     Procedure, Public :: create     => matrix_create_complex                 !! Create storage for the data of the matrix 
     Procedure, Public :: local_size => matrix_local_size_complex             !! Get the dimensions of the local part of the matrix
     Procedure, Public :: remap      => complex_remap                         !! Remap the data to another distribution
     Procedure, Public :: diag       => complex_diag                          !! Diagonalise the (assumed Hermitian) matrix
     Procedure, Public :: extract    => extract_complex
     ! Private implementations
     Procedure,            Private :: set_global_real    => complex_matrix_set_global_real
     Procedure,            Private :: set_global_complex => complex_matrix_set_global_complex
     Procedure,            Private :: get_global_real    => complex_matrix_get_global_real
     Procedure,            Private :: get_global_complex => complex_matrix_get_global_complex
     Procedure,            Private :: multiply           => complex_multiply
     Procedure, Pass( B ), Private :: real_multiply      => real_multiply_complex
     Procedure, Pass( B ), Private :: complex_multiply   => complex_multiply_complex
     Procedure, Pass( A ), Private :: rscal_multiply     => rscal_multiply_complex
     Procedure,            Private :: multiply_rscal     => complex_multiply_rscal
     Procedure,            Private :: plus               => plus_complex
     Procedure,            Private :: add                => complex_add
     Procedure, Pass( B ), Private :: real_add           => real_add_complex
     Procedure, Pass( B ), Private :: complex_add        => complex_add_complex
     Procedure,            Private :: add_diagonal       => complex_add_diagonal
     Procedure, Pass( A ), Private :: diagonal_add       => diagonal_add_complex
     Procedure,            Private :: minus              => minus_complex
     Procedure,            Private :: subtract           => complex_subtract
     Procedure, Pass( B ), Private :: real_subtract      => real_subtract_complex
     Procedure, Pass( B ), Private :: complex_subtract   => complex_subtract_complex
     Procedure,            Private :: subtract_diagonal  => complex_subtract_diagonal
     Procedure, Pass( A ), Private :: diagonal_subtract  => diagonal_subtract_complex
     Procedure, Pass( Q ), Private :: real_diag          => real_diag_complex
     Procedure, Pass( Q ), Private :: complex_diag       => complex_diag_complex
     Procedure           , Private :: choleski           => choleski_complex
     Procedure           , Private :: tr_inv             => tr_inv_complex
     Procedure, Pass( B ), Private :: real_remap         => real_remap_complex
     Procedure, Pass( B ), Private :: complex_remap      => complex_remap_complex
  End type complex_distributed_matrix

  Public :: distributed_matrix_init
  Public :: distributed_matrix_comm_to_base
  Public :: distributed_matrix_finalise
  Public :: distributed_matrix_set_default_blocking
  
  Private

  Integer, Parameter, Private :: diag_work_size_fiddle_factor = 4 ! From experience Scalapack sometimes returns too small a work size
  
  Integer, Parameter, Private :: default_block_fac = 96
  Integer,            Private :: block_fac = default_block_fac

  Abstract Interface

     Subroutine create( A, m, n, source_matrix )
       !! Create storage for the data of the matrix
       Import :: distributed_matrix
       Implicit None
       Class( distributed_matrix ), Intent(   Out ) :: A
       Integer                    , Intent( In    ) :: m
       Integer                    , Intent( In    ) :: n
       Class( distributed_matrix ), Intent( In    ) :: source_matrix
     End Subroutine create

     Function local_size( A, dim ) Result( n )
       !! Get the dimensions of the local part of the matrix
       Import :: distributed_matrix
       Implicit None
       Integer                                   :: n
       Class( distributed_matrix ), Intent( In ) :: A
       Integer                    , Intent( In ) :: dim
     End Function local_size

     Subroutine set_global_real( A, m, n, p, q, data )
       !! Set values with a real    array using global indexing
       Import :: wp
       Import :: distributed_matrix
       Implicit None
       Class( distributed_matrix )    , Intent( InOut ) :: A
       Integer                        , Intent( In    ) :: m
       Integer                        , Intent( In    ) :: n
       Integer                        , Intent( In    ) :: p
       Integer                        , Intent( In    ) :: q
       Real( wp ), Dimension( m:, p: ), Intent( In    ) :: data
     End Subroutine set_global_real
     Subroutine set_global_complex( A, m, n, p, q, data )
       !! Set values with a complex array using global indexing
       Import :: wp
       Import :: distributed_matrix
       Implicit None
       Class( distributed_matrix )       , Intent( InOut ) :: A
       Integer                           , Intent( In    ) :: m
       Integer                           , Intent( In    ) :: n
       Integer                           , Intent( In    ) :: p
       Integer                           , Intent( In    ) :: q
       Complex( wp ), Dimension( m:, p: ), Intent( In    ) :: data
     End Subroutine set_global_complex

     Subroutine get_global_real( A, m, n, p, q, data )
       !! Get values from a real array using global indexing
       Import :: wp
       Import :: distributed_matrix
       Implicit None
       Class( distributed_matrix )    , Intent( In    ) :: A
       Integer                        , Intent( In    ) :: m
       Integer                        , Intent( In    ) :: n
       Integer                        , Intent( In    ) :: p
       Integer                        , Intent( In    ) :: q
       Real( wp ), Dimension( m:, p: ), Intent(   Out ) :: data
     End Subroutine get_global_real
     Subroutine get_global_complex( A, m, n, p, q, data )
       !! Get values from a complex array using global indexing
       Import :: wp
       Import :: distributed_matrix
       Implicit None
       Class( distributed_matrix )       , Intent( In    ) :: A
       Integer                           , Intent( In    ) :: m
       Integer                           , Intent( In    ) :: n
       Integer                           , Intent( In    ) :: p
       Integer                           , Intent( In    ) :: q
       Complex( wp ), Dimension( m:, p: ), Intent(   Out ) :: data
     End Subroutine get_global_complex

     Function pre_rscal_op( s, A ) Result( B )
       !! Pre-scale a matrix with a scalar
       Import :: wp
       Import :: distributed_matrix
       Implicit None
       Class( distributed_matrix ), Allocatable  :: B
       Real( wp )                 , Intent( In ) :: s
       Class( distributed_matrix ), Intent( In ) :: A
     End Function pre_rscal_op
     Function post_rscal_op( A, s ) Result( B )
       !! Post-scale a matrix with a scalar
       Import :: wp
       Import :: distributed_matrix
       Implicit None
       Class( distributed_matrix ), Allocatable  :: B
       Class( distributed_matrix ), Intent( In ) :: A
       Real( wp )                 , Intent( In ) :: s
     End Function post_rscal_op
     
     Function pre_diagonal_op( d, A ) Result( B )
       !! Pre-apply a digonal matrix to a matrix
       Import :: wp
       Import :: distributed_matrix
       Implicit None
       Class( distributed_matrix ), Allocatable  :: B
       Real( wp ), Dimension( : ) , Intent( In ) :: d
       Class( distributed_matrix ), Intent( In ) :: A
     End Function pre_diagonal_op
     Function post_diagonal_op( A, d ) Result( B )
       !! Post-apply a digonal matrix to a matrix
       Import :: wp
       Import :: distributed_matrix
       Implicit None
       Class( distributed_matrix ), Allocatable  :: B
       Class( distributed_matrix ), Intent( In ) :: A
       Real( wp ), Dimension( : ) , Intent( In ) :: d
     End Function post_diagonal_op

     Function unary_op( A ) Result( C )
       !! A unary operation on a base class object
       Import :: distributed_matrix
       Implicit None
       Class( distributed_matrix ), Allocatable  :: C
       Class( distributed_matrix ), Intent( In ) :: A
     End Function unary_op

     Function binary_op( A, B ) Result( C )
       !! A binary operation between two base class objects
       Import :: distributed_matrix
       Implicit None
       Class( distributed_matrix ), Allocatable  :: C
       Class( distributed_matrix ), Intent( In ) :: A
       Class( distributed_matrix ), Intent( In ) :: B
     End Function binary_op
     Function real_binary_op( A, B ) Result( C )
       !! A binary operation between a real matrix and the base class
       Import ::      distributed_matrix
       Import :: real_distributed_matrix
       Implicit None
       Class(      distributed_matrix ), Allocatable  :: C
       Class( real_distributed_matrix ), Intent( In ) :: A
       Class(      distributed_matrix ), Intent( In ) :: B
     End Function real_binary_op
     Function complex_binary_op( A, B ) Result( C )
       !! A binary operation between a complex matrix and the base class
       Import ::         distributed_matrix
       Import :: complex_distributed_matrix
       Implicit None
       Class(         distributed_matrix ), Allocatable  :: C
       Class( complex_distributed_matrix ), Intent( In ) :: A
       Class(         distributed_matrix ), Intent( In ) :: B
     End Function complex_binary_op

     Subroutine diag_op( A, Q, E )
       !! Diagonalising a base class matrix with vectors returned as base class
       Import :: wp
       Import :: distributed_matrix
       Implicit None
       Class( distributed_matrix ),             Intent( In    ) :: A
       Class( distributed_matrix ),             Intent(   Out ) :: Q
       Real( wp ), Dimension( : ), Allocatable, Intent(   Out ) :: E
     End Subroutine diag_op
     Subroutine real_diag_op( A, Q, E )
       !! Diagonalising a base class matrix with vectors returned as base class
       Import :: wp
       Import :: distributed_matrix
       Import :: real_distributed_matrix
       Implicit None
       Class( real_distributed_matrix ),        Intent( In    ) :: A
       Class(      distributed_matrix ),        Intent(   Out ) :: Q
       Real( wp ), Dimension( : ), Allocatable, Intent(   Out ) :: E
     End Subroutine real_diag_op
     Subroutine complex_diag_op( A, Q, E )
       !! Diagonalising a base class matrix with vectors returned as base class
       Import :: wp
       Import :: distributed_matrix
       Import :: complex_distributed_matrix
       Implicit None
       Class( complex_distributed_matrix ),     Intent( In    ) :: A
       Class(         distributed_matrix ),     Intent(   Out ) :: Q
       Real( wp ), Dimension( : ), Allocatable, Intent(   Out ) :: E
     End Subroutine complex_diag_op

     Subroutine remap_op( A, is_A_dummy, parent_comm, B, is_B_dummy ) 
       !! A remap operation between two base class objects
       Import :: distributed_matrix
       Implicit None
       Class( distributed_matrix ), Intent( In    ) :: A
       Logical                    , Intent( In    ) :: is_A_dummy
       Integer                    , Intent( In    ) :: parent_comm
       Class( distributed_matrix ), Intent( InOut ) :: B
       Logical                    , Intent( In    ) :: is_B_dummy
     End Subroutine remap_op
     Subroutine real_remap_op( A, is_A_dummy, parent_comm, B, is_B_dummy ) 
       !! A remap operation between the base class and a real matrix
       Import ::      distributed_matrix
       Import :: real_distributed_matrix
       Implicit None
       Class( real_distributed_matrix ), Intent( In    ) :: A
       Logical                         , Intent( In    ) :: is_A_dummy
       Integer                         , Intent( In    ) :: parent_comm
       Class(      distributed_matrix ), Intent( InOut ) :: B
       Logical                         , Intent( In    ) :: is_B_dummy
     End Subroutine real_remap_op
     Subroutine complex_remap_op( A, is_A_dummy, parent_comm, B, is_B_dummy ) 
       !! A remap operation between the base class and a real matrix
       Import ::         distributed_matrix
       Import :: complex_distributed_matrix
       Implicit None
       Class( complex_distributed_matrix ), Intent( In    ) :: A
       Logical                         , Intent( In    ) :: is_A_dummy
       Integer                         , Intent( In    ) :: parent_comm
       Class(      distributed_matrix ), Intent( InOut ) :: B
       Logical                         , Intent( In    ) :: is_B_dummy
     End Subroutine complex_remap_op

     Function extract_op( A, m, n, p, q ) Result( B )
       !! Extract a pacth from one matrix to form another matrix
       Import :: wp
       Import :: distributed_matrix
       Implicit None
       Class( distributed_matrix )    , Allocatable  :: B
       Class( distributed_matrix )    , Intent( In ) :: A
       Integer                        , Intent( In ) :: m
       Integer                        , Intent( In ) :: n
       Integer                        , Intent( In ) :: p
       Integer                        , Intent( In ) :: q
     End Function extract_op

  End Interface
  
Contains

  !####################################################################
  ! Non-type bound procedures

  Subroutine distributed_matrix_init

    Use matrix_mapping_module, Only : matrix_mapping_init

    !! Initialise the matrix system

    Call matrix_mapping_init
    
  End Subroutine distributed_matrix_init

  Subroutine distributed_matrix_comm_to_base( comm, base_matrix )

    !! Converts an MPI communicator into the data structures
    !! required to describe a matrix mapped onto it
    
    Use matrix_mapping_module, Only : matrix_mapping_comm_to_base

    Integer                      , Intent( In    ) :: comm
    Class  ( distributed_matrix ), Intent(   Out ) :: base_matrix 

    Type( matrix_mapping ) :: base_matrix_mapping
    
    Call matrix_mapping_comm_to_base( comm, base_matrix_mapping )

    base_matrix%matrix_map = base_matrix_mapping
    base_matrix%global_to_local_rows = [ distributed_matrix_INVALID ]
    base_matrix%global_to_local_cols = [ distributed_matrix_INVALID ]
    base_matrix%local_to_global_rows = [ distributed_matrix_INVALID ]
    base_matrix%local_to_global_cols = [ distributed_matrix_INVALID ]
    
  End Subroutine distributed_matrix_comm_to_base

  Subroutine distributed_matrix_finalise

    !! Finalise the matrix system

    Use matrix_mapping_module, Only : matrix_mapping_finalise

    ! Reset the blocking factor to the default
    Call distributed_matrix_set_default_blocking( default_block_fac )

    Call matrix_mapping_finalise
    
  End Subroutine distributed_matrix_finalise

  Subroutine distributed_matrix_set_default_blocking( bfac )

    !! Set the default blocking factor
    
    Integer, Intent( In ) :: bfac

    block_fac = bfac
    
  End Subroutine distributed_matrix_set_default_blocking

  !###################################################################################
  ! Methods implemented on the base type and inherited by the extended types

  Subroutine matrix_get_maps( A, gl_rows, gl_cols, lg_rows, lg_cols )

    !! For the given matrix get all the mapping arrays

    Class( distributed_matrix )         , Intent( In    ) :: A
    Integer, Dimension( : ), Allocatable, Intent(   Out ) :: gl_rows
    Integer, Dimension( : ), Allocatable, Intent(   Out ) :: gl_cols
    Integer, Dimension( : ), Allocatable, Intent(   Out ) :: lg_rows
    Integer, Dimension( : ), Allocatable, Intent(   Out ) :: lg_cols

    ! Note using allocate on set
    gl_rows = A%global_to_local_rows
    gl_cols = A%global_to_local_cols
    lg_rows = A%local_to_global_rows
    lg_cols = A%local_to_global_cols
    
  End Subroutine matrix_get_maps

  Function matrix_global_to_local( A, what ) Result( gl_indexing )

    !! For the given matrix get an array for mapping global indices to local  ones
    !! Get row mapping arrays if WHAT='C' or 'c', row ones if it is equal to 'R' or 'r'

    Integer, Dimension( : ), Allocatable :: gl_indexing
    
    Class( distributed_matrix ), Intent( In ) :: A
    Character( Len = * )       , Intent( In ) :: what

    Select Case( what )
    Case Default
       Stop "Illegal WHAT in global_to_local"
    Case( 'R', 'r' )
       gl_indexing = A%global_to_local_rows
    Case( 'C', 'c' )
       gl_indexing = A%global_to_local_cols
    End Select

  End Function matrix_global_to_local
  
  Function matrix_local_to_global( A, what ) Result( lg_indexing )

    !! For the given matrix get an array for mapping local indices to global  ones
    !! Get row mapping arrays if WHAT='C' or 'c', row ones if it is equal to 'R' or 'r'

    Integer, Dimension( : ), Allocatable :: lg_indexing
    
    Class( distributed_matrix ), Intent( In ) :: A
    Character( Len = * )       , Intent( In ) :: what

    Select Case( what )
    Case Default
       Stop "Illegal WHAT in local_to_global"
    Case( 'R', 'r' )
       lg_indexing = A%local_to_global_rows
    Case( 'C', 'c' )
       lg_indexing = A%local_to_global_cols
    End Select

  End Function matrix_local_to_global

  Function matrix_size( A, dim ) Result( n )

    !! Get the dimensions of the matrix

    Integer :: n

    Class( distributed_matrix ), Intent( In ) :: A
    Integer                    , Intent( In ) :: dim

    If( dim <= 2 ) Then

       If( .Not. A%daggered ) Then
          Select Case( dim )
          Case Default
             Stop "Ilegal dimension in matrix_size"
          Case( 1 )
             Call A%matrix_map%get_data( m = n )
          Case( 2 )
             Call A%matrix_map%get_data( n = n )
          End Select
       Else
          Select Case( dim )
          Case Default
             Stop "Ilegal dimension in matrix_size"
          Case( 1 )
             Call A%matrix_map%get_data( n = n )
          Case( 2 )
             Call A%matrix_map%get_data( m = n )
          End Select
       End If
       
    Else

       Stop "Ilegal dimension in matrix_size"
       
    End If

  End Function matrix_size

  Function matrix_communicator( A ) Result( c )

    !! Get the communicator containing the processes holding the matrix

    Integer :: c

    Class( distributed_matrix ), Intent( In ) :: A

    c = A%matrix_map%get_comm()
    
  End Function matrix_communicator

  Pure Function matrix_dagger( A ) Result( tm )

    !! Apply the dagger operator to the matrix
    !  Note this does NOT do any communication. Instead we simply set the flag
    !  saying that until further notice the matrix should be used in tranposed form
    
    Class( distributed_matrix ), Allocatable :: tm

    Class( distributed_matrix ), Intent( In ) :: A

    Allocate( tm, Source = A )
    tm%daggered = .Not. tm%daggered
    
  End Function matrix_dagger

  !##########################################################################################
  ! Over-ridding routines

  Subroutine matrix_create_real( A, m, n, source_matrix )

    !! Create the data for a real MxN matrix. The distribution is the same as the provided source matrix.

    Use, Intrinsic :: ieee_arithmetic, Only : ieee_value, ieee_support_nan, ieee_signaling_nan, &
         ieee_support_halting, ieee_get_halting_mode, ieee_set_halting_mode,                    &
         ieee_get_flag, ieee_set_flag, ieee_invalid

    Use Scalapack_interfaces , Only : numroc

    Class( real_distributed_matrix ), Intent(   Out ) :: A
    Integer                         , Intent( In    ) :: m
    Integer                         , Intent( In    ) :: n
    Class(      distributed_matrix ), Intent( In    ) :: source_matrix

    Integer :: nprow, myprow, mb, lda
    Integer :: npcol, mypcol, nb, sda
    Integer :: ctxt

    Logical :: is_invalid, is_halting

    ! NEED TO FIX IF N, M SMALLER THAN BLOCKING FAC
    mb = block_fac
    nb = block_fac
    mb = Min( mb, nb )
    nb = mb

    Call source_matrix%matrix_map%get_data( nprow = nprow, myprow = myprow )
    lda = numroc( m, mb, myprow, 0, nprow )

    Call source_matrix%matrix_map%get_data( npcol = npcol, mypcol = mypcol )
    sda = numroc( n, nb, mypcol, 0, npcol )

    Call source_matrix%matrix_map%get_data( ctxt = ctxt )

    Call A%matrix_map%set( source_matrix%matrix_map%proc_mapping, ctxt, m, n, mb, nb, 0, 0, lda )

    Call set_local_to_global( A%local_to_global_rows, m, mb, myprow, nprow, lda )
    Call set_local_to_global( A%local_to_global_cols, n, nb, mypcol, npcol, sda )

    Call set_global_to_local( A%global_to_local_rows, m, mb, myprow, nprow )
    Call set_global_to_local( A%global_to_local_cols, n, nb, mypcol, npcol )

    Allocate( A%data( 1:lda, 1:sda  ) )
    ! Try to initialise with signalling NANs
    If( ieee_support_nan( A%data ) .And. ieee_support_halting( ieee_invalid ) ) Then
       ! This processor supports ieee maths and control of halting - Use this as carefully as possible to initialise
       ! matrix type objects to a value that can help detect their use when unitilised - namely a ignalling NaN
       ! First get the current halting mode for ieee invalid 
       Call ieee_get_halting_mode( ieee_invalid, is_halting )
       ! Now deliberately turn halting off
       Call ieee_set_halting_mode( ieee_invalid, .False. )
       ! Get the current value of the invalid flag to avoid spurious signalling caused by the below
       Call ieee_get_flag( ieee_invalid, is_invalid )
       A%data = ieee_value( A%data, ieee_signaling_nan )
       ! Reset the invalid flag to what it was before we deliberatley used a signalling NaN to avoid missing ones later
       Call ieee_set_flag( ieee_invalid, is_invalid )
       ! And reset the halting mode
       Call ieee_set_halting_mode( ieee_invalid, is_halting )
    Else
       ! The processor doesn't support ieee or doesn't support halting control under ieee
       ! Simply initialise to a very big value
       A%data = Huge( A%data )
    End If

    ! Indicate that A is not tranposed
    A%daggered = .False.

  End Subroutine matrix_create_real

  Subroutine matrix_create_complex( A, m, n, source_matrix )

    !! Create the data for a complex MxN matrix. The distribution is the same as the provided source matrix.

    Use, Intrinsic :: ieee_arithmetic, Only : ieee_value, ieee_support_nan, ieee_signaling_nan, &
         ieee_support_halting, ieee_get_halting_mode, ieee_set_halting_mode,                    &
         ieee_get_flag, ieee_set_flag, ieee_invalid

    Use Scalapack_interfaces , Only : numroc

    Class( complex_distributed_matrix ), Intent(   Out ) :: A
    Integer                            , Intent( In    ) :: m
    Integer                            , Intent( In    ) :: n
    Class(         distributed_matrix ), Intent( In    ) :: source_matrix

    Integer :: nprow, myprow, mb, lda
    Integer :: npcol, mypcol, nb, sda
    Integer :: ctxt

    Logical :: is_invalid, is_halting

    ! Need to fix if n, m smaller than blocking fac
    mb = block_fac
    nb = block_fac
    mb = Min( mb, nb )
    nb = mb

    Call source_matrix%matrix_map%get_data( nprow = nprow, myprow = myprow )
    lda = numroc( m, mb, myprow, 0, nprow )

    Call source_matrix%matrix_map%get_data( npcol = npcol, mypcol = mypcol )
    sda = numroc( n, nb, mypcol, 0, npcol )

    Call source_matrix%matrix_map%get_data( ctxt = ctxt )

    Call A%matrix_map%set( source_matrix%matrix_map%proc_mapping, ctxt, m, n, mb, nb, 0, 0, lda )

    Call set_local_to_global( A%local_to_global_rows, m, mb, myprow, nprow, lda )
    Call set_local_to_global( A%local_to_global_cols, n, nb, mypcol, npcol, sda )

    Call set_global_to_local( A%global_to_local_rows, m, mb, myprow, nprow )
    Call set_global_to_local( A%global_to_local_cols, n, nb, mypcol, npcol )

    Allocate( A%data( 1:lda, 1:sda  ) )
    ! Try to initialise with signalling NANs
    If( ieee_support_nan( Real( A%data ) ) .And. ieee_support_halting( ieee_invalid ) ) Then
       ! This processor supports ieee maths and control of halting - Use this as carefully as possible to initialise
       ! matrix type objects to a value that can help detect their use when unitilised - namely a ignalling NaN
       ! First get the current halting mode for ieee invalid 
       Call ieee_get_halting_mode( ieee_invalid, is_halting )
       ! Now deliberately turn halting off
       Call ieee_set_halting_mode( ieee_invalid, .False. )
       ! Get the current value of the invalid flag to avoid spurious signalling caused by the below
       Call ieee_get_flag( ieee_invalid, is_invalid )
       A%data = Cmplx( ieee_value( 0.0_wp, ieee_signaling_nan ), &
                       ieee_value( 0.0_wp, ieee_signaling_nan ), wp )
       ! Reset the invalid flag to what it was before we deliberatley used a signalling NaN to avoid missing ones later
       Call ieee_set_flag( ieee_invalid, is_invalid )
       ! And reset the halting mode
       Call ieee_set_halting_mode( ieee_invalid, is_halting )
    Else
       ! The processor doesn't support ieee or doesn't support halting control under ieee
       ! Simply initialise to a very big value
       A%data = Huge( Real( A%data, wp ) )
    End If

    ! Indicate that A is not tranposed
    A%daggered = .False.
        
  End Subroutine matrix_create_complex
  
  Function matrix_local_size_real( A, dim ) Result( n )

    !! Get the dimensions of the local part of the matrix

    Integer :: n

    Class( real_distributed_matrix ), Intent( In ) :: A
    Integer                    , Intent( In ) :: dim

    If( dim <= 2 ) Then

       n = Size( A%data, Dim = dim )

    Else

       Stop "Ilegal dimension in matrix_local_size"

    End If
       
  End Function matrix_local_size_real

  Function matrix_local_size_complex( A, dim ) Result( n )

    !! Get the dimensions of the local part of the matrix

    Integer :: n

    Class( complex_distributed_matrix ), Intent( In ) :: A
    Integer                    , Intent( In ) :: dim

    If( dim <= 2 ) Then

       n = Size( A%data, Dim = dim )

    Else

       Stop "Ilegal dimension in matrix_local_size"

    End If
       
  End Function matrix_local_size_complex

  Subroutine real_matrix_set_global_real( A, m, n, p, q, data )

    !! Sets the data ( m:n, p:q ) in the global matrix

    Class( real_distributed_matrix ), Intent( InOut ) :: A
    Integer                         , Intent( In    ) :: m
    Integer                         , Intent( In    ) :: n
    Integer                         , Intent( In    ) :: p
    Integer                         , Intent( In    ) :: q
    Real( wp ), Dimension( m:, p: ) , Intent( In    ) :: data

    Integer :: i_glob, j_glob
    Integer :: i_loc , j_loc

    ! Could be optimised by introducing ranges for the mapping arryas

    If( .Not. A%daggered ) Then
       Do j_glob = p, q
          j_loc = A%global_to_local_cols( j_glob )
          If( j_loc == distributed_matrix_NOT_ME ) Cycle
          Do i_glob = m, n
             i_loc = A%global_to_local_rows( i_glob )
             If( i_loc == distributed_matrix_NOT_ME ) Cycle
             A%data( i_loc, j_loc ) = data( i_glob, j_glob )
          End Do
       End Do
    Else
       Do j_glob = p, q
          j_loc = A%global_to_local_rows( j_glob )
          If( j_loc == distributed_matrix_NOT_ME ) Cycle
          Do i_glob = m, n
             i_loc = A%global_to_local_cols( i_glob )
             If( i_loc == distributed_matrix_NOT_ME ) Cycle
             A%data( j_loc, i_loc ) = data( i_glob, j_glob )
          End Do
       End Do
    End If
       
  End Subroutine real_matrix_set_global_real

  Subroutine real_matrix_set_global_complex( A, m, n, p, q, data )

    !! Sets the data ( m:n, p:q ) in the global matrix

    Class( real_distributed_matrix )  , Intent( InOut ) :: A
    Integer                           , Intent( In    ) :: m
    Integer                           , Intent( In    ) :: n
    Integer                           , Intent( In    ) :: p
    Integer                           , Intent( In    ) :: q
    Complex( wp ), Dimension( m:, p: ), Intent( In    ) :: data

    Stop "Trying to set real matrix with complex data in real_matrix_set_global_complex"
    ! Shut the compiler up about unused vars
    Write( *, * ) A%data, m, n, p, q, data

  End Subroutine real_matrix_set_global_complex

  Subroutine complex_matrix_set_global_real( A, m, n, p, q, data )

    !! Sets the data ( m:n, p:q ) in the global matrix

    Class( complex_distributed_matrix ), Intent( InOut ) :: A
    Integer                            , Intent( In    ) :: m
    Integer                            , Intent( In    ) :: n
    Integer                            , Intent( In    ) :: p
    Integer                            , Intent( In    ) :: q
    Real( wp ), Dimension( m:, p: )    , Intent( In    ) :: data

    Stop "Trying to set complex matrix with real data in complex_matrix_set_global_real"
    ! Shut the compiler up about unused vars
    Write( *, * ) A%data, m, n, p, q, data

  End Subroutine complex_matrix_set_global_real

  Subroutine complex_matrix_set_global_complex( A, m, n, p, q, data )

    !! Sets the data ( m:n, p:q ) in the global matrix

    Class( complex_distributed_matrix ), Intent( InOut ) :: A
    Integer                            , Intent( In    ) :: m
    Integer                            , Intent( In    ) :: n
    Integer                            , Intent( In    ) :: p
    Integer                            , Intent( In    ) :: q
    Complex( wp ), Dimension( m:, p: ) , Intent( In    ) :: data

    Integer :: i_glob, j_glob
    Integer :: i_loc , j_loc

    ! Could be optimised by introducing ranges for the mapping arrays
    
    If( .Not. A%daggered ) Then
       Do j_glob = p, q
          j_loc = A%global_to_local_cols( j_glob )
          If( j_loc == distributed_matrix_NOT_ME ) Cycle
          Do i_glob = m, n
             i_loc = A%global_to_local_rows( i_glob )
             If( i_loc == distributed_matrix_NOT_ME ) Cycle
             A%data( i_loc, j_loc ) = data( i_glob, j_glob )
          End Do
       End Do
    Else
       Do j_glob = p, q
          j_loc = A%global_to_local_rows( j_glob )
          If( j_loc == distributed_matrix_NOT_ME ) Cycle
          Do i_glob = m, n
             i_loc = A%global_to_local_cols( i_glob )
             If( i_loc == distributed_matrix_NOT_ME ) Cycle
             A%data( j_loc, i_loc ) = Conjg( data( i_glob, j_glob ) )
          End Do
       End Do
    End If
       
  End Subroutine complex_matrix_set_global_complex

  Subroutine real_matrix_get_global_real( A, m, n, p, q, data )

    Use mpi, Only : MPI_Type_create_f90_real, MPI_Sizeof, MPI_Type_match_size, MPI_Allreduce, MPI_In_place, MPI_Sum, &
         MPI_Double_precision, MPI_Typeclass_real
    
    !! Gets the data ( m:n, p:q ) in the global A

    Class( real_distributed_matrix ), Intent( In    ) :: A
    Integer                         , Intent( In    ) :: m
    Integer                         , Intent( In    ) :: n
    Integer                         , Intent( In    ) :: p
    Integer                         , Intent( In    ) :: q
    Real( wp ), Dimension( m:, p: ) , Intent(   Out ) :: data
    
    Real( wp ) :: rdum

    Integer :: i_glob, j_glob
    Integer :: i_loc , j_loc
    Integer :: handle
    Integer :: rsize, error
    
    data = 0.0_wp
    If( .Not. A%daggered ) Then
       Do j_glob = p, q
          j_loc = A%global_to_local_cols( j_glob )
          If( j_loc == distributed_matrix_NOT_ME ) Cycle
          Do i_glob = m, n
             i_loc = A%global_to_local_rows( i_glob )
             If( i_loc == distributed_matrix_NOT_ME ) Cycle
             data( i_glob, j_glob ) = A%data( i_loc, j_loc )
          End Do
       End Do
    Else
       Do j_glob = p, q
          j_loc = A%global_to_local_rows( j_glob )
          If( j_loc == distributed_matrix_NOT_ME ) Cycle
          Do i_glob = m, n
             i_loc = A%global_to_local_cols( i_glob )
             If( i_loc == distributed_matrix_NOT_ME ) Cycle
             data( i_glob, j_glob ) = A%data( j_loc, i_loc )
          End Do
       End Do
    End If
    ! Generate a portable MPI data type handle from the variable to be communicated
    Call MPI_Type_create_f90_real( Precision( data ), Range( data ), handle, error )
    ! Replicate the data
!!!!HACK TO WORK AROUND BUG IN MVAPICH2
!!!!    Call MPI_Allreduce( MPI_IN_PLACE, data, Size( data ), handle, MPI_SUM, A%matrix_map%get_comm(), error )
    Call MPI_Sizeof( rdum, rsize, error )
    ! Note MPI_Type_match_size does NOT create a new handle, it returns the value of an exisiting one. Hence no need to free
    Call MPI_Type_match_size( MPI_TYPECLASS_REAL, rsize, handle, error )
    Call MPI_Allreduce( MPI_IN_PLACE, data, Size( data ), handle, MPI_SUM, A%matrix_map%get_comm(), error )
       
  End Subroutine real_matrix_get_global_real

  Subroutine real_matrix_get_global_complex( A, m, n, p, q, data )

    !! Gets the data ( m:n, p:q ) in the global A

    Class( real_distributed_matrix )   , Intent( In    ) :: A
    Integer                            , Intent( In    ) :: m
    Integer                            , Intent( In    ) :: n
    Integer                            , Intent( In    ) :: p
    Integer                            , Intent( In    ) :: q
    Complex( wp ), Dimension( m:, p: ) , Intent(   Out ) :: data

    Stop "Trying to get from real matrix into complex data in real_matrix_get_global_complex"
    ! Shut the compiler up about unused vars
    data = 0.0_wp
    Write( *, * ) A%data, m, n, p, q

  End Subroutine real_matrix_get_global_complex

  Subroutine complex_matrix_get_global_real( A, m, n, p, q, data )
    
    !! Gets the data ( m:n, p:q ) in the global A

    Class( complex_distributed_matrix ), Intent( In    ) :: A
    Integer                            , Intent( In    ) :: m
    Integer                            , Intent( In    ) :: n
    Integer                            , Intent( In    ) :: p
    Integer                            , Intent( In    ) :: q
    Real( wp ), Dimension( m:, p: )    , Intent(   Out ) :: data

    Stop "Trying to get from complex matrix into real data in complex_matrix_get_global_real"
    ! Shut the compiler up about unused vars
    data = 0.0_wp
    Write( *, * ) A%data, m, n, p, q

  End Subroutine complex_matrix_get_global_real

  Subroutine complex_matrix_get_global_complex( A, m, n, p, q, data )

    Use mpi, Only : MPI_Type_create_f90_complex, MPI_Sizeof, MPI_Type_match_size, MPI_Allreduce, MPI_In_place, MPI_Sum, &
         MPI_Double_complex, MPI_Typeclass_complex

    !! Gets the data ( m:n, p:q ) in the global matrix

    Class( complex_distributed_matrix ), Intent( In    ) :: A
    Integer                            , Intent( In    ) :: m
    Integer                            , Intent( In    ) :: n
    Integer                            , Intent( In    ) :: p
    Integer                            , Intent( In    ) :: q
    Complex( wp ), Dimension( m:, p: ) , Intent(   Out ) :: data

    Complex( wp ) :: cdum

    Integer :: i_glob, j_glob
    Integer :: i_loc , j_loc
    Integer :: csize, handle
    Integer :: error
    
    data = 0.0_wp
    If( .Not. A%daggered ) Then
       Do j_glob = p, q
          j_loc = A%global_to_local_cols( j_glob )
          If( j_loc == distributed_matrix_NOT_ME ) Cycle
          Do i_glob = m, n
             i_loc = A%global_to_local_rows( i_glob )
             If( i_loc == distributed_matrix_NOT_ME ) Cycle
             data( i_glob, j_glob ) = A%data( i_loc, j_loc )
          End Do
       End Do
    Else
       Do j_glob = p, q
          j_loc = A%global_to_local_rows( j_glob )
          If( j_loc == distributed_matrix_NOT_ME ) Cycle
          Do i_glob = m, n
             i_loc = A%global_to_local_cols( i_glob )
             If( i_loc == distributed_matrix_NOT_ME ) Cycle
             data( i_glob, j_glob ) = Conjg( A%data( j_loc, i_loc ) )
          End Do
       End Do
    End If
    ! Generate a portable MPI data type handle from the variable to be communicated
    Call MPI_Type_create_f90_complex( Precision( data ), Range( data ), handle, error )
    ! Replicate the data
!!!!HACK TO WORK AROUND BUG IN MVAPICH2
!!!!    Call MPI_Allreduce( MPI_IN_PLACE, data, Size( data ), handle, MPI_SUM, A%matrix_map%get_comm(), error )
    Call MPI_sizeof( cdum, csize, error )
    ! Note MPI_Type_match_size does NOT create a new handle, it returns the value of an exisiting one. Hence no need to free
    Call MPI_type_match_size( MPI_Typeclass_complex, csize, handle, error )
    Call MPI_Allreduce( MPI_In_place, data, Size( data ), handle, MPI_Sum, A%matrix_map%get_comm(), error )
       
  End Subroutine complex_matrix_get_global_complex

  ! Multiplication routines

  Function real_multiply( A, B ) Result( C )

    !! Multiply a real matrix by something yet to be determined
  
    Class(      distributed_matrix ), Allocatable :: C

    Class( real_distributed_matrix ), Intent( In ) :: A
    Class(      distributed_matrix ), Intent( In ) :: B

    ! We know what A is. Now use dispatch on B to pick out what it is
    C = B%real_multiply( A )
    
  End Function real_multiply

  Function complex_multiply( A, B ) Result( C )
  
    !! Multiply a complex matrix by something yet to be determined

    Class(         distributed_matrix ), Allocatable :: C

    Class( complex_distributed_matrix ), Intent( In ) :: A
    Class(         distributed_matrix ), Intent( In ) :: B

    ! We know what A is. Now use dispatch on B to pick out what it is
    C = B%complex_multiply( A )
    
  End Function complex_multiply

  Function real_multiply_real( A, B ) Result( C )

    !! Multiply a real matrix by a real matrix
    
    Use Scalapack_interfaces , Only : pdgemm

    Class(      distributed_matrix ), Allocatable :: C

    Class( real_distributed_matrix ), Intent( In ) :: A
    Class( real_distributed_matrix ), Intent( In ) :: B

    Type( real_distributed_matrix ) :: T

    Integer :: ma, na
    Integer :: mb, nb
    Integer :: m, n, k

    Character :: t1, t2

    ! Work out the dimensions of the result
    t1 = Merge( 'T', 'N', A%daggered )
    t2 = Merge( 'T', 'N', B%daggered )
    
    Call A%matrix_map%get_data( m = ma, n = na )
    Call B%matrix_map%get_data( m = mb, n = nb )
    
    If( t1 == 'N' .And. t2 == 'N' ) Then
       m = ma
       n = nb
       k = na
    Else If( t1 == 'T' .And. t2 == 'N' ) Then
       m = na
       n = nb
       k = ma
    Else If( t1 == 'N' .And. t2 == 'T' ) Then
       m = ma
       n = mb
       k = na
    Else If( t1 == 'T' .And. t2 == 'T' ) Then
       m = na
       n = mb
       k = ma
    Else
       Stop 'How did we get here in real_multiply_real???'
    End If

    ! Create a result matrix of the right type
    Call T%create( m, n, A )

    ! Do the multiplication
    Call pdgemm( t1, t2, m, n, k, 1.0_wp, A%data, 1, 1, A%matrix_map%get_descriptor(), &
                                          B%data, 1, 1, B%matrix_map%get_descriptor(), &
                                  0.0_wp, T%data, 1, 1, T%matrix_map%get_descriptor() )

    ! And store the result
    C = T
    
  End Function real_multiply_real
     
  Function real_multiply_complex( A, B ) Result( C )

    !! Multiply a real matrix by complex matrix - ILLEGAL

    Class(         distributed_matrix ), Allocatable :: C

    Class(    real_distributed_matrix ), Intent( In ) :: A
    Class( complex_distributed_matrix ), Intent( In ) :: B

    Stop "Illegal combination of arguments in real_multiply_complex"
    Deallocate( C )
    Write( *, * ) A%data, B%data

  End Function real_multiply_complex

  Function complex_multiply_real( A, B ) Result( C )

    !! Multiply a complex matrix by real matrix - ILLEGAL

    Class(         distributed_matrix ), Allocatable :: C

    Class( complex_distributed_matrix ), Intent( In ) :: A
    Class(    real_distributed_matrix ), Intent( In ) :: B

    Stop "Illegal combination of arguments in complex_multiply_real"
    Deallocate( C )
    Write( *, * ) A%data, B%data
    
  End Function complex_multiply_real

  Function complex_multiply_complex( A, B ) Result( C )

    !! Multiply a complex matrix by complex matrix

    Use Scalapack_interfaces , Only : pzgemm

    Class(         distributed_matrix ), Allocatable :: C

    Class( complex_distributed_matrix ), Intent( In ) :: A
    Class( complex_distributed_matrix ), Intent( In ) :: B

    Type( complex_distributed_matrix ) :: T
    
    Integer :: ma, na
    Integer :: mb, nb
    Integer :: m, n, k

    Character :: t1, t2

    ! Work out the shape of the result
    t1 = Merge( 'C', 'N', A%daggered )
    t2 = Merge( 'C', 'N', B%daggered )
       
    Call A%matrix_map%get_data( m = ma, n = na )
    Call B%matrix_map%get_data( m = mb, n = nb )

    If( t1 == 'N' .And. t2 == 'N' ) Then
       m = ma
       n = nb
       k = na
    Else If( t1 == 'C' .And. t2 == 'N' ) Then
       m = na
       n = nb
       k = ma
    Else If( t1 == 'N' .And. t2 == 'C' ) Then
       m = ma
       n = mb
       k = na
    Else If( t1 == 'C' .And. t2 == 'C' ) Then
       m = na
       n = mb
       k = ma
    Else
       Stop 'How did we get here in matrix_multiply_complex???'
    End If
    
    Call T%create( m, n, A )

    Call pzgemm( t1, t2, m, n, k, ( 1.0_wp, 0.0_wp ), A%data, 1, 1, A%matrix_map%get_descriptor(), &
                                                      B%data, 1, 1, B%matrix_map%get_descriptor(), &
                                  ( 0.0_wp, 0.0_wp ), T%data, 1, 1, T%matrix_map%get_descriptor() )

    C = T
    
  End Function complex_multiply_complex

  Function real_multiply_rscal( A, s ) Result( B )

    !! Post-scale a real matrix

    Class(      distributed_matrix ), Allocatable :: B

    Class( real_distributed_matrix ), Intent( In ) :: A
    Real( wp )                      , Intent( In ) :: s

    Type( real_distributed_matrix ) :: T

    T = A
    T%data = s * T%data
    B = T
    
  End Function real_multiply_rscal

  Function rscal_multiply_real( s, A ) Result( B )

    !! Pre-scale a real matrix

    Class(      distributed_matrix ), Allocatable :: B

    Real( wp )                      , Intent( In ) :: s
    Class( real_distributed_matrix ), Intent( In ) :: A

    Type( real_distributed_matrix ) :: T

    T = A
    T%data = s * A%data
    B = T
    
  End Function rscal_multiply_real

  Function  complex_multiply_rscal( A, s ) Result( B )

    !! Post-scale a complex matrix

    Class(         distributed_matrix ), Allocatable :: B

    Class( complex_distributed_matrix ), Intent( In ) :: A
    Real( wp )                         , Intent( In ) :: s

    Type( complex_distributed_matrix ) :: T

    T = A
    T%data = s * A%data
    B = T
    
  End Function complex_multiply_rscal

  Function  rscal_multiply_complex( s, A ) Result( B )

    !! Pre-scale a complex matrix

    Class(         distributed_matrix ), Allocatable :: B

    Real( wp )                         , Intent( In ) :: s
    Class( complex_distributed_matrix ), Intent( In ) :: A

    Type( complex_distributed_matrix ) :: T

    T = A
    T%data = s * A%data
    B = T
    
  End Function rscal_multiply_complex

  ! Addition routines

  Function real_add( A, B ) Result( C )

    !! Add a real matrix to something yet to be determined
  
    Class(      distributed_matrix ), Allocatable :: C

    Class( real_distributed_matrix ), Intent( In ) :: A
    Class(      distributed_matrix ), Intent( In ) :: B

    ! We know what A is. Now use dispatch on B to pick out what it is
    C = B%real_add( A )
    
  End Function real_add

  Function complex_add( A, B ) Result( C )
  
    !! Add a complex matrix to something yet to be determined

    Class(         distributed_matrix ), Allocatable :: C

    Class( complex_distributed_matrix ), Intent( In ) :: A
    Class(         distributed_matrix ), Intent( In ) :: B

    ! We know what A is. Now use dispatch on B to pick out what it is
    C = B%complex_add( A )
    
  End Function complex_add

  Function real_add_real( A, B ) Result( C )

    !! Adds two real matrices A and B, returning the result in C
    ! Slightly complicated by the tranpose options on A and B, i.e. A and B may be flagged that they have
    ! previously been tranposed, and so we have to add the correct forms together. Do this in a way
    ! to avoid communication whereever possible, but as it gets a bit fiddly we will do
    ! it very slowly and explicitly
    ! The actual addition will be by pdgeadd( op, ... A, ... C ) which does C -> C + op( A ), so by careful
    ! use of the transposes we can keep things sane
    
    Use Scalapack_interfaces, Only : pdgeadd

    Class(      distributed_matrix ), Allocatable :: C

    Class( real_distributed_matrix ), Intent( In ) :: A
    Class( real_distributed_matrix ), Intent( In ) :: B

    Type( real_distributed_matrix ) :: T

    Character :: tA, TB

    Integer :: mA, nA
    Integer :: mB, nB
    Integer :: mT, nT
    
    ! Work out the tranposes
    tA = Merge( 'T', 'N', A%daggered )
    tB = Merge( 'T', 'N', B%daggered )

    ! Get the shapes of the input matrices
    Call A%matrix_map%get_data( m = mA, n = nA )
    Call B%matrix_map%get_data( m = mB, n = nB )
    
    ! Consider each in turn
    ! We will generate the result in T which is a real_distributed_matrix
    ! and then once we have it return the result via C

    If     ( tA == 'N' .And. tB == 'N' ) Then
       ! Neither matrix to be used in tranposed form
       ! Perform A -> A + B
       ! First check the matrices are compatible
       If( mA /= mB .Or. nA /= nB ) Then
          Stop "Incompatible sizes in real_add_real" 
       End If
       ! Shape of the result
       mT = mA
       nT = nA
       ! Now generate the result in T
       T = A
       Call pdgeadd( 'N', mT, nT, 1.0_wp, B%data, 1, 1, B%matrix_map%get_descriptor(), &
                                  1.0_wp, T%data, 1, 1, T%matrix_map%get_descriptor() )
       
    Else If( tA == 'N' .And. tB == 'T' ) Then
       ! A not transposed, B to be used in transposed form
       ! Perform A -> A + Tranpose( B )
       ! First check the matrices are compatible
       If( mA /= nB .Or. nA /= mB ) Then
          Stop "Incompatible sizes in real_add_real" 
       End If
       ! Shape of the result
       mT = mA
       nT = nA
       ! Now generate the result in T
       T = A
       Call pdgeadd( 'T', mT, nT, 1.0_wp, B%data, 1, 1, B%matrix_map%get_descriptor(), &
                                  1.0_wp, T%data, 1, 1, T%matrix_map%get_descriptor() )

    Else If( tA == 'T' .And. tB == 'N' ) Then
       ! A tranposed, B not transposed
       ! Perform B -> B + Tranpose( A )
       ! First check the matrices are compatible
       If( mA /= nB .Or. nA /= mB ) Then
          Stop "Incompatible sizes in real_add_real" 
       End If
       ! Shape of the result
       mT = mB
       nT = nB
       ! Now generate the result in T
       T = B
       Call pdgeadd( 'T', mT, nT, 1.0_wp, A%data, 1, 1, A%matrix_map%get_descriptor(), &
                                  1.0_wp, T%data, 1, 1, T%matrix_map%get_descriptor() )

    Else If( tA == 'T' .And. tB == 'T' ) Then
       ! Both matrices in transposed form
       ! THIS IS THE TRICKSY ONE
       ! Perform A -> A + B AND THEN indicate the matrix is returned in tranposed form to avoid comms
       ! First check the matrices are compatible
       If( mA /= mB .Or. nA /= nB ) Then
          Stop "Incompatible sizes in real_add_real" 
       End If
       ! Shape of the result
       mT = mA
       nT = nA
       ! Now generate the result in T
       T = A
       Call pdgeadd( 'N', mT, nT, 1.0_wp, B%data, 1, 1, B%matrix_map%get_descriptor(), &
                                  1.0_wp, T%data, 1, 1, T%matrix_map%get_descriptor() )
       ! Now "tranpose" the result
       T%daggered = .True.

    Else
       Stop "How did we get here in real_add_real?"
    End If

    ! Return the result
    C = T

  End Function real_add_real

  Function real_add_complex( A, B ) Result( C )

    !! Add a real matrix to a complex matrix - ILLEGAL

    Class(         distributed_matrix ), Allocatable :: C

    Class(    real_distributed_matrix ), Intent( In ) :: A
    Class( complex_distributed_matrix ), Intent( In ) :: B

    Stop "Illegal combination of arguments in real_add_complex"
    Deallocate( C )
    Write( *, * ) A%data, B%data

  End Function real_add_complex

  Function complex_add_real( A, B ) Result( C )

    !! Add a complex matrix to a real matrix - ILLEGAL

    Class(         distributed_matrix ), Allocatable :: C

    Class( complex_distributed_matrix ), Intent( In ) :: A
    Class(    real_distributed_matrix ), Intent( In ) :: B

    Stop "Illegal combination of arguments in complex_add_real"
    Deallocate( C )
    Write( *, * ) A%data, B%data
    
  End Function complex_add_real
  
  Function complex_add_complex( A, B ) Result( C )

    !! Adds two real matrices A and B, returning the result in C
    ! Slightly complicated by the tranpose options on A and B, i.e. A and B may be flagged that they have
    ! previously been tranposed, and so we have to add the correct forms together. Do this in a way
    ! to avoid communication whereever possible, but as it gets a bit fiddly we will do
    ! it very slowly and explicitly
    ! The actual addition will be by pdgeadd( op, ... A, ... C ) which does C -> C + op( A ), so by careful
    ! use of the transposes we can keep things sane
    
    Use Scalapack_interfaces, Only : pzgeadd

    Class(         distributed_matrix ), Allocatable :: C

    Class( complex_distributed_matrix ), Intent( In ) :: A
    Class( complex_distributed_matrix ), Intent( In ) :: B

    Type( complex_distributed_matrix ) :: T

    Character :: tA, TB

    Integer :: mA, nA
    Integer :: mB, nB
    Integer :: mT, nT
    
    ! Work out the tranposes
    tA = Merge( 'C', 'N', A%daggered )
    tB = Merge( 'C', 'N', B%daggered )

    ! Get the shapes of the input matrices
    Call A%matrix_map%get_data( m = mA, n = nA )
    Call B%matrix_map%get_data( m = mB, n = nB )
    
    ! Consider each in turn
    ! We will generate the result in T which is a real_distributed_matrix
    ! and then once we have it return the result via C

    If     ( tA == 'N' .And. tB == 'N' ) Then
       ! Neither matrix to be used in tranposed form
       ! Perform A -> A + B
       ! First check the matrices are compatible
       If( mA /= mB .Or. nA /= nB ) Then
          Stop "Incompatible sizes in real_add_real" 
       End If
       ! Shape of the result
       mT = mA
       nT = nA
       ! Now generate the result in T
       T = A
       Call pzgeadd( 'N', mT, nT, ( 1.0_wp, 0.0_wp ), B%data, 1, 1, B%matrix_map%get_descriptor(), &
                                  ( 1.0_wp, 0.0_wp ), T%data, 1, 1, T%matrix_map%get_descriptor() )
       
    Else If( tA == 'N' .And. tB == 'C' ) Then
       ! A not transposed, B to be used in transposed form
       ! Perform A -> A + Tranpose( B )
       ! First check the matrices are compatible
       If( mA /= nB .Or. nA /= mB ) Then
          Stop "Incompatible sizes in real_add_real" 
       End If
       ! Shape of the result
       mT = mA
       nT = nA
       ! Now generate the result in T
       T = A
       Call pzgeadd( 'C', mT, nT, ( 1.0_wp, 0.0_wp ), B%data, 1, 1, B%matrix_map%get_descriptor(), &
                                  ( 1.0_wp, 0.0_wp ), T%data, 1, 1, T%matrix_map%get_descriptor() )

    Else If( tA == 'C' .And. tB == 'N' ) Then
       ! A tranposed, B not transposed
       ! Perform B -> B + Tranpose( A )
       ! First check the matrices are compatible
       If( mA /= nB .Or. nA /= mB ) Then
          Stop "Incompatible sizes in real_add_real" 
       End If
       ! Shape of the result
       mT = mB
       nT = nB
       ! Now generate the result in T
       T = B
       Call pzgeadd( 'C', mT, nT, ( 1.0_wp, 0.0_wp ), A%data, 1, 1, A%matrix_map%get_descriptor(), &
                                  ( 1.0_wp, 0.0_wp ), T%data, 1, 1, T%matrix_map%get_descriptor() )

    Else If( tA == 'C' .And. tB == 'C' ) Then
       ! Both matrices in transposed form
       ! THIS IS THE TRICKSY ONE
       ! Perform A -> A + B AND THEN indicate the matrix is returned in tranposed form to avoid comms
       ! First check the matrices are compatible
       If( mA /= mB .Or. nA /= nB ) Then
          Stop "Incompatible sizes in real_add_real" 
       End If
       ! Shape of the result
       mT = mA
       nT = nA
       ! Now generate the result in T
       T = A
       Call pzgeadd( 'N', mT, nT, ( 1.0_wp, 0.0_wp ), B%data, 1, 1, B%matrix_map%get_descriptor(), &
                                  ( 1.0_wp, 0.0_wp ), T%data, 1, 1, T%matrix_map%get_descriptor() )
       ! Now "tranpose" the result
       T%daggered = .True.

    Else
       Stop "How did we get here in real_add_real?"
    End If

    ! Return the result
    C = T

  End Function complex_add_complex

  Function real_add_diagonal( A, d ) Result( B )

    !! Add a real matrix to a diagonal matrix

    Class(      distributed_matrix ), Allocatable :: B

    Class( real_distributed_matrix ),                 Intent( In ) :: A
    Real( wp )                      , Dimension( : ), Intent( In ) :: d

    Type( real_distributed_matrix ) :: T
    
    Integer :: m, n
    Integer :: i_glob
    Integer :: i_loc, j_loc
    
    Call A%matrix_map%get_data( m = m, n = n )

    If( m == n .And. Size( d ) == n ) Then

       T = A
       Do i_glob = 1, n
          i_loc = A%global_to_local_rows( i_glob )
          j_loc = A%global_to_local_cols( i_glob )
          If(  i_loc /= distributed_matrix_NOT_ME .And. &
               j_loc /= distributed_matrix_NOT_ME ) Then
             T%data( i_loc, j_loc ) = A%data( i_loc, j_loc ) + d( i_glob )
          End If
       End Do

    Else

       Stop "Inconsistent matrix dimensions in real_add_diagonal "

    End If

    B = T
    
  End Function real_add_diagonal

  Function complex_add_diagonal( A, d ) Result( B )

    !! Add a complex matrix to a diagonal matrix

    Class(         distributed_matrix ), Allocatable :: B

    Class( complex_distributed_matrix ),                 Intent( In ) :: A
    Real( wp )                         , Dimension( : ), Intent( In ) :: d

    Type( complex_distributed_matrix ) :: T

    Integer :: m, n
    Integer :: i_glob
    Integer :: i_loc, j_loc
    
    Call A%matrix_map%get_data( m = m, n = n )

    If( m == n .And. Size( d ) == n ) Then

       T = A
       Do i_glob = 1, n
          i_loc = A%global_to_local_rows( i_glob )
          j_loc = A%global_to_local_cols( i_glob )
          If(  i_loc /= distributed_matrix_NOT_ME .And. &
               j_loc /= distributed_matrix_NOT_ME ) Then
             T%data( i_loc, j_loc ) = A%data( i_loc, j_loc ) + d( i_glob )
          End If
       End Do

    Else

       Stop "Inconsistent matrix dimensions in complex_add_diagonal "

    End If

    B = T
    
  End Function complex_add_diagonal

  Function diagonal_add_real( d, A ) Result( B )

    !! Add a diagonal matrix to a real matrix 

    Class(      distributed_matrix ), Allocatable :: B

    Real( wp )                      , Dimension( : ), Intent( In ) :: d
    Class( real_distributed_matrix ),                 Intent( In ) :: A

    B = A + d
    
  End Function diagonal_add_real

  Function diagonal_add_complex( d, A ) Result( B )

    !! Add a diagonal matrix to a complex matrix

    Class(         distributed_matrix ), Allocatable :: B

    Real( wp )                         , Dimension( : ), Intent( In ) :: d
    Class( complex_distributed_matrix ),                 Intent( In ) :: A

    B = A + d

  End Function diagonal_add_complex

  ! Subtraction Routines
  
  Function real_subtract( A, B ) Result( C )

    !! Subtract from a real matrix something yet to be determined
  
    Class(      distributed_matrix ), Allocatable :: C

    Class( real_distributed_matrix ), Intent( In ) :: A
    Class(      distributed_matrix ), Intent( In ) :: B

    ! We know what A is. Now use dispatch on B to pick out what it is
    C = B%real_subtract( A )
    
  End Function real_subtract

  Function complex_subtract( A, B ) Result( C )
  
    !! Subtract from a complex matrix something yet to be determined

    Class(         distributed_matrix ), Allocatable :: C

    Class( complex_distributed_matrix ), Intent( In ) :: A
    Class(         distributed_matrix ), Intent( In ) :: B

    ! We know what A is. Now use dispatch on B to pick out what it is
    C = B%complex_subtract( A )
    
  End Function complex_subtract

  Function real_subtract_real( A, B ) Result( C )

    !! Forms A - B, returning the result in C
    ! Slightly complicated by the tranpose options on A and B, i.e. A and B may be flagged that they have
    ! previously been tranposed, and so we have to subtract the correct forms together. Do this in a way
    ! to avoid communication whereever possible, but as it gets a bit fiddly we will do
    ! it very slowly and explicitly
    ! The actual subtraction will be by pdgeadd( op, ... A, ... C ) which does C -> C + op( A ), so by careful
    ! use of the transposes we can keep things sane
    
    Use Scalapack_interfaces, Only : pdgeadd

    Class(      distributed_matrix ), Allocatable :: C

    Class( real_distributed_matrix ), Intent( In ) :: A
    Class( real_distributed_matrix ), Intent( In ) :: B

    Type( real_distributed_matrix ) :: T

    Character :: tA, TB

    Integer :: mA, nA
    Integer :: mB, nB
    Integer :: mT, nT
    
    ! Work out the tranposes
    tA = Merge( 'T', 'N', A%daggered )
    tB = Merge( 'T', 'N', B%daggered )

    ! Get the shapes of the input matrices
    Call A%matrix_map%get_data( m = mA, n = nA )
    Call B%matrix_map%get_data( m = mB, n = nB )
    
    ! Consider each in turn
    ! We will generate the result in T which is a real_distributed_matrix
    ! and then once we have it return the result via C

    If     ( tA == 'N' .And. tB == 'N' ) Then
       ! Neither matrix to be used in tranposed form
       ! Perform A -> -B + A
       ! First check the matrices are compatible
       If( mA /= mB .Or. nA /= nB ) Then
          Stop "Incompatible sizes in real_subtract_real" 
       End If
       ! Shape of the result
       mT = mA
       nT = nA
       ! Now generate the result in T
       T = A
       Call pdgeadd( 'N', mT, nT, -1.0_wp, B%data, 1, 1, B%matrix_map%get_descriptor(), &
                                   1.0_wp, T%data, 1, 1, T%matrix_map%get_descriptor() )
       
    Else If( tA == 'N' .And. tB == 'T' ) Then
       ! A not transposed, B to be used in transposed form
       ! Perform A -> -Tranpose( B ) + A
       ! First check the matrices are compatible
       If( mA /= nB .Or. nA /= mB ) Then
          Stop "Incompatible sizes in real_subtract_real" 
       End If
       ! Shape of the result
       mT = mA
       nT = nA
       ! Now generate the result in T
       T = A
       Call pdgeadd( 'T', mT, nT, -1.0_wp, B%data, 1, 1, B%matrix_map%get_descriptor(), &
                                   1.0_wp, T%data, 1, 1, T%matrix_map%get_descriptor() )

    Else If( tA == 'T' .And. tB == 'N' ) Then
       ! A tranposed, B not transposed
       ! Perform B -> Tranpose( A ) - B
       ! First check the matrices are compatible
       If( mA /= nB .Or. nA /= mB ) Then
          Stop "Incompatible sizes in real_subtract_real" 
       End If
       ! Shape of the result
       mT = mB
       nT = nB
       ! Now generate the result in T
       T = B
       Call pdgeadd( 'T', mT, nT,  1.0_wp, A%data, 1, 1, A%matrix_map%get_descriptor(), &
                                  -1.0_wp, T%data, 1, 1, T%matrix_map%get_descriptor() )

    Else If( tA == 'T' .And. tB == 'T' ) Then
       ! Both matrices in transposed form
       ! THIS IS THE TRICKSY ONE
       ! Perform A -> -B + A AND THEN indicate the matrix is returned in tranposed form to avoid comms
       ! First check the matrices are compatible
       If( mA /= mB .Or. nA /= nB ) Then
          Stop "Incompatible sizes in real_subtract_real" 
       End If
       ! Shape of the result
       mT = mA
       nT = nA
       ! Now generate the result in T
       T = A
       Call pdgeadd( 'N', mT, nT, -1.0_wp, B%data, 1, 1, B%matrix_map%get_descriptor(), &
                                   1.0_wp, T%data, 1, 1, T%matrix_map%get_descriptor() )
       ! Now "tranpose" the result
       T%daggered = .True.

    Else
       Stop "How did we get here in real_subtract_real?"
    End If

    ! Return the result
    C = T

  End Function real_subtract_real

  Function real_subtract_complex( A, B ) Result( C )

    !! Subtract from a real matrix a complex matrix - ILLEGAL

    Class(         distributed_matrix ), Allocatable :: C

    Class(    real_distributed_matrix ), Intent( In ) :: A
    Class( complex_distributed_matrix ), Intent( In ) :: B

    Stop "Illegal combination of arguments in real_subtract_complex"
    ! Silly lines to stop compiler warning about unused vars
    Deallocate( C )
    Write( *, * ) A%data, B%data

  End Function real_subtract_complex

  Function complex_subtract_real( A, B ) Result( C )

    !! Subtract from a complex matrix a real matrix - ILLEGAL

    Class(         distributed_matrix ), Allocatable :: C

    Class( complex_distributed_matrix ), Intent( In ) :: A
    Class(    real_distributed_matrix ), Intent( In ) :: B

    Stop "Illegal combination of arguments in complex_subtract_real"
    ! Silly lines to stop compiler warning about unused vars
    Deallocate( C )
    Write( *, * ) A%data, B%data
    
  End Function complex_subtract_real
  
  Function complex_subtract_complex( A, B ) Result( C )

    !! Froms A - B, returning the result in C
    ! Slightly complicated by the tranpose options on A and B, i.e. A and B may be flagged that they have
    ! previously been tranposed, and so we have to subtract the correct forms together. Do this in a way
    ! to avoid communication whereever possible, but as it gets a bit fiddly we will do
    ! it very slowly and explicitly
    ! The actual subtraction will be by pzgeadd( op, ... A, ... C ) which does C -> C + op( A ), so by careful
    ! use of the transposes we can keep things sane
    
    Use Scalapack_interfaces, Only : pzgeadd

    Class(         distributed_matrix ), Allocatable :: C

    Class( complex_distributed_matrix ), Intent( In ) :: A
    Class( complex_distributed_matrix ), Intent( In ) :: B

    Type( complex_distributed_matrix ) :: T

    Character :: tA, TB

    Integer :: mA, nA
    Integer :: mB, nB
    Integer :: mT, nT
    
    ! Work out the tranposes
    tA = Merge( 'C', 'N', A%daggered )
    tB = Merge( 'C', 'N', B%daggered )

    ! Get the shapes of the input matrices
    Call A%matrix_map%get_data( m = mA, n = nA )
    Call B%matrix_map%get_data( m = mB, n = nB )
    
    ! Consider each in turn
    ! We will generate the result in T which is a real_distributed_matrix
    ! and then once we have it return the result via C

    If     ( tA == 'N' .And. tB == 'N' ) Then
       ! Neither matrix to be used in tranposed form
       ! Perform A -> -B + A
       ! First check the matrices are compatible
       If( mA /= mB .Or. nA /= nB ) Then
          Stop "Incompatible sizes in real_subtract_real" 
       End If
       ! Shape of the result
       mT = mA
       nT = nA
       ! Now generate the result in T
       T = A
       Call pzgeadd( 'N', mT, nT, - ( 1.0_wp, 0.0_wp ), B%data, 1, 1, B%matrix_map%get_descriptor(), &
                                    ( 1.0_wp, 0.0_wp ), T%data, 1, 1, T%matrix_map%get_descriptor() )
       
    Else If( tA == 'N' .And. tB == 'C' ) Then
       ! A not transposed, B to be used in transposed form
       ! Perform A -> -Tranpose( B ) + B
       ! First check the matrices are compatible
       If( mA /= nB .Or. nA /= mB ) Then
          Stop "Incompatible sizes in real_subtract_real" 
       End If
       ! Shape of the result
       mT = mA
       nT = nA
       ! Now generate the result in T
       T = A
       Call pzgeadd( 'C', mT, nT, - ( 1.0_wp, 0.0_wp ), B%data, 1, 1, B%matrix_map%get_descriptor(), &
                                    ( 1.0_wp, 0.0_wp ), T%data, 1, 1, T%matrix_map%get_descriptor() )

    Else If( tA == 'C' .And. tB == 'N' ) Then
       ! A tranposed, B not transposed
       ! Perform B -> Tranpose( A ) - B
       ! First check the matrices are compatible
       If( mA /= nB .Or. nA /= mB ) Then
          Stop "Incompatible sizes in real_subtract_real" 
       End If
       ! Shape of the result
       mT = mB
       nT = nB
       ! Now generate the result in T
       T = B
       Call pzgeadd( 'C', mT, nT,   ( 1.0_wp, 0.0_wp ), A%data, 1, 1, A%matrix_map%get_descriptor(), &
                                  - ( 1.0_wp, 0.0_wp ), T%data, 1, 1, T%matrix_map%get_descriptor() )

    Else If( tA == 'C' .And. tB == 'C' ) Then
       ! Both matrices in transposed form
       ! THIS IS THE TRICKSY ONE
       ! Perform A -> -B + A AND THEN indicate the matrix is returned in tranposed form to avoid comms
       ! First check the matrices are compatible
       If( mA /= mB .Or. nA /= nB ) Then
          Stop "Incompatible sizes in real_subtract_real" 
       End If
       ! Shape of the result
       mT = mA
       nT = nA
       ! Now generate the result in T
       T = A
       Call pzgeadd( 'N', mT, nT, - ( 1.0_wp, 0.0_wp ), B%data, 1, 1, B%matrix_map%get_descriptor(), &
                                    ( 1.0_wp, 0.0_wp ), T%data, 1, 1, T%matrix_map%get_descriptor() )
       ! Now "tranpose" the result
       T%daggered = .True.

    Else
       Stop "How did we get here in real_subtract_real?"
    End If

    ! Return the result
    C = T

  End Function complex_subtract_complex

  Function real_subtract_diagonal( A, d ) Result( B )

    !! Subtract a real digaonal matrix from a real one

    Class(      distributed_matrix ), Allocatable :: B

    Class( real_distributed_matrix ),                 Intent( In ) :: A
    Real( wp )                      , Dimension( : ), Intent( In ) :: d

    Type( real_distributed_matrix ) :: T
    
    Integer :: m, n
    Integer :: i_glob
    Integer :: i_loc, j_loc
    
    Call A%matrix_map%get_data( m = m, n = n )

    If( m == n .And. Size( d ) == n ) Then

       T = A
       Do i_glob = 1, n
          i_loc = A%global_to_local_rows( i_glob )
          j_loc = A%global_to_local_cols( i_glob )
          If(  i_loc /= distributed_matrix_NOT_ME .And. &
               j_loc /= distributed_matrix_NOT_ME ) Then
             T%data( i_loc, j_loc ) = A%data( i_loc, j_loc ) - d( i_glob )
          End If
       End Do

    Else

       Stop "Inconsistent matrix dimensions in real_subtract_diagonal "

    End If

    B = T
    
  End Function real_subtract_diagonal

  Function complex_subtract_diagonal( A, d ) Result( B )

    !! Subtract a real digaonal matrix from a complex one

    Class(         distributed_matrix ), Allocatable :: B

    Class( complex_distributed_matrix ),                 Intent( In ) :: A
    Real( wp )                         , Dimension( : ), Intent( In ) :: d

    Type( complex_distributed_matrix ) :: T

    Integer :: m, n
    Integer :: i_glob
    Integer :: i_loc, j_loc
    
    Call A%matrix_map%get_data( m = m, n = n )

    If( m == n .And. Size( d ) == n ) Then

       T = A
       Do i_glob = 1, n
          i_loc = A%global_to_local_rows( i_glob )
          j_loc = A%global_to_local_cols( i_glob )
          If(  i_loc /= distributed_matrix_NOT_ME .And. &
               j_loc /= distributed_matrix_NOT_ME ) Then
             T%data( i_loc, j_loc ) = A%data( i_loc, j_loc ) - d( i_glob )
          End If
       End Do

    Else

       Stop "Inconsistent matrix dimensions in complex_subtract_diagonal "

    End If

    B = T
    
  End Function complex_subtract_diagonal

  Function diagonal_subtract_real( d, A ) Result( B )

    !! Subtract a real matrix from a digaonal one

    Class(      distributed_matrix ), Allocatable :: B

    Real( wp )                      , Dimension( : ), Intent( In ) :: d
    Class( real_distributed_matrix ),                 Intent( In ) :: A

    Type( real_distributed_matrix ) :: T
    
    Integer :: m, n
    Integer :: i_glob
    Integer :: i_loc, j_loc
    
    Call A%matrix_map%get_data( m = m, n = n )

    If( m == n .And. Size( d ) == n ) Then

       T = A
       T%data = -T%data
       Do i_glob = 1, n
          i_loc = A%global_to_local_rows( i_glob )
          j_loc = A%global_to_local_cols( i_glob )
          If(  i_loc /= distributed_matrix_NOT_ME .And. &
               j_loc /= distributed_matrix_NOT_ME ) Then
             T%data( i_loc, j_loc ) = - A%data( i_loc, j_loc ) + d( i_glob )
          End If
       End Do

    Else

       Stop "Inconsistent matrix dimensions in diagonal_subtract_real"

    End If

    B = T
    
  End Function diagonal_subtract_real

  Function diagonal_subtract_complex( d, A ) Result( B )

    !! Subtract a complex matrix from a digaonal one

    Class(         distributed_matrix ), Allocatable :: B

    Real( wp )                         , Dimension( : ), Intent( In ) :: d
    Class( complex_distributed_matrix ),                 Intent( In ) :: A

    Type( complex_distributed_matrix ) :: T

    Integer :: m, n
    Integer :: i_glob
    Integer :: i_loc, j_loc
    
    Call A%matrix_map%get_data( m = m, n = n )

    If( m == n .And. Size( d ) == n ) Then

       T = A
       T%data = -T%data
       Do i_glob = 1, n
          i_loc = A%global_to_local_rows( i_glob )
          j_loc = A%global_to_local_cols( i_glob )
          If(  i_loc /= distributed_matrix_NOT_ME .And. &
               j_loc /= distributed_matrix_NOT_ME ) Then
             T%data( i_loc, j_loc ) = - A%data( i_loc, j_loc ) + d( i_glob )
          End If
       End Do

    Else

       Stop "Inconsistent matrix dimensions in diagonal_subtract_complex"

    End If

    B = T
    
  End Function diagonal_subtract_complex

  ! Diagonalisation routines

  Subroutine real_diag( A, Q, E )

    !! Diagonalise a real matrix
    
    Implicit None

    Class( real_distributed_matrix ),              Intent( In    ) :: A
    Class(      distributed_matrix ),              Intent(   Out ) :: Q
    Real( wp ), Dimension( : )      , Allocatable, Intent(   Out ) :: E

    Call Q%real_diag( A, E )

  End Subroutine real_diag

  Subroutine complex_diag( A, Q, E )

    !! Diagonalise a complex matrix
    
    Implicit None

    Class( complex_distributed_matrix ),              Intent( In    ) :: A
    Class(         distributed_matrix ),              Intent(   Out ) :: Q
    Real( wp ), Dimension( : )         , Allocatable, Intent(   Out ) :: E

    Call Q%complex_diag( A, E )

  End Subroutine complex_diag

  Subroutine real_diag_real( A, Q, E )

    !! Diagonalise a real symmetric matrix

    Use Scalapack_interfaces, Only : pdsyevd

    Implicit None

    Class( real_distributed_matrix ),              Intent( In    ) :: A
    Class( real_distributed_matrix ),              Intent(   Out ) :: Q
    Real( wp ), Dimension( : )      , Allocatable, Intent(   Out ) :: E

    Real( wp ), Dimension( :, : ), Allocatable :: tmp_a

    Real( wp ), Dimension( : ), Allocatable :: work

    Integer, Dimension( : ), Allocatable :: iwork
    
    Integer :: nwork
    Integer :: npcol
    Integer :: m
    Integer :: info

    ! Give Q the same mapping as A
    Call A%matrix_map%get_data( m = m, npcol = npcol )
    Call Q%create( m, m, A )

    Allocate( E( 1:m ) )
    
    ! The diag overwrites the matrix. Horrible so use a temporary
    tmp_A = A%data
    
    ! Workspace size enquiry
    Allocate( work( 1:1 ), iwork( 1:1 ) )
    Call pdsyevd( 'V', 'U', m, tmp_A, 1, 1, A%matrix_map%get_descriptor(), E, Q%data, 1, 1, Q%matrix_map%get_descriptor(), &
         work, -1, iwork, 0, info )
    nwork = Nint( work( 1 ) )
    nwork = nwork * diag_work_size_fiddle_factor ! From experience ...
    Deallocate( work, iwork )
    Allocate(  work( 1:nwork ) )
    ! Scalapack recipe is behind the strange numbers
    Allocate( iwork( 1:7 * m + 8 * npcol + 2 ) )
    ! Do the diag
    Call pdsyevd( 'V', 'U', m, tmp_A, 1, 1, A%matrix_map%get_descriptor(), E, Q%data, 1, 1, Q%matrix_map%get_descriptor(), &
         work, Size( work ), iwork, Size( iwork ), info )

    If( info /= 0 ) Then
       Deallocate( E )
    End If

  End Subroutine real_diag_real

  Subroutine real_diag_complex( A, Q, E )

    !! Diagonalise a real matrix returning complex vectors - ILLEGAL

    Implicit None

    Class(    real_distributed_matrix ),              Intent( In    ) :: A
    Class( complex_distributed_matrix ),              Intent(   Out ) :: Q
    Real( wp ), Dimension( : )         , Allocatable, Intent(   Out ) :: E

    Stop "Illegal combination of arguments in complex_multiply_real"
    ! Silly lines to stop compiler warning about unused vars
    Deallocate( E )
    Write( *, * ) A%data, Q%data

  End Subroutine real_diag_complex

  Subroutine complex_diag_real( A, Q, E )

    !! Diagonalise a complex matrix returning real vectors - ILLEGAL

    Implicit None

    Class( complex_distributed_matrix ),              Intent( In    ) :: A
    Class(    real_distributed_matrix ),              Intent(   Out ) :: Q
    Real( wp ), Dimension( : )         , Allocatable, Intent(   Out ) :: E

    Stop "Illegal combination of arguments in complex_multiply_real"
    ! Silly lines to stop compiler warning about unused vars
    Deallocate( E )
    Write( *, * ) A%data, Q%data

  End Subroutine complex_diag_real

  Subroutine complex_diag_complex( A, Q, E )

    !! Diagonalise a complex Hermitian matrix

    Use Scalapack_interfaces, Only : pzheevd

    Implicit None

    Class( complex_distributed_matrix ),              Intent( In    ) :: A
    Class( complex_distributed_matrix ),              Intent(   Out ) :: Q
    Real( wp ), Dimension( : )         , Allocatable, Intent(   Out ) :: E

    Complex( wp ), Dimension( :, : ), Allocatable :: tmp_a

    Complex( wp ), Dimension( : ), Allocatable :: cwork

    Real( wp ), Dimension( : ), Allocatable :: rwork

    Integer, Dimension( : ), Allocatable :: iwork
    
    Integer :: ncwork, nrwork
    Integer :: npcol
    Integer :: m
    Integer :: info

    ! Give Q the same mapping as A
    Call A%matrix_map%get_data( m = m, npcol = npcol )
    Call Q%create( m, m, A )
    
    Allocate( E( 1:m ) )

    ! The diag overwrites the matrix. Horrible so use a temporary
    tmp_A = A%data
       
    ! Workspace size enquiry
    Allocate( cwork( 1:1 ), rwork( 1:1 ), iwork( 1:1 ) )
    Call pzheevd( 'V', 'U', m, tmp_A, 1, 1, A%matrix_map%get_descriptor(), E, Q%data, 1, 1, Q%matrix_map%get_descriptor(), &
         cwork, -1, rwork, -1, iwork, 0, info )
    ncwork = Nint( Real( cwork( 1 ), wp ) )
    ncwork = ncwork * diag_work_size_fiddle_factor ! From experience ...
    nrwork = Nint( rwork( 1 ) )
    nrwork = nrwork * diag_work_size_fiddle_factor ! From experience ...
    Deallocate( cwork, rwork, iwork )
    Allocate( cwork( 1:ncwork ) )
    Allocate( rwork( 1:nrwork ) )
    ! Scalapack recipe is behind the strange numbers
    Allocate( iwork( 1:7 * m + 8 * npcol + 2 ) )
    ! Do the diag
    Call pzheevd( 'V', 'U', m, tmp_A, 1, 1, A%matrix_map%get_descriptor(), E, Q%data, 1, 1, Q%matrix_map%get_descriptor(), &
            cwork, Size( cwork ), rwork, Size( rwork ), iwork, Size( iwork ), info )

    If( info /= 0 ) Then
       Deallocate( E )
    End If

  End Subroutine complex_diag_complex

  ! Choleski Routines
  
  Function choleski_real( A ) Result( C )

    !! Choleski decompose into lower triangular factors a real symmetric positive definite matrix
    !! Return value is deallocated on error

    Use Scalapack_interfaces, Only : pdpotrf
    
    Class(      distributed_matrix ), Allocatable :: C

    Class( real_distributed_matrix ), Intent( In ) :: A

    Type( real_distributed_matrix ) :: T
    
    Integer :: m
    Integer :: i_glob, j_glob
    Integer :: i, j
    Integer :: error

    T = A
    
    ! Zero Upper half of T
    Do j = 1, Size( T%data, Dim = 2 )
       j_glob = T%local_to_global_cols( j )
       Do i = 1, Size( T%data, Dim = 1 )
          i_glob = T%local_to_global_rows( i )
          If( j_glob > i_glob ) Then
             T%data( i, j ) = 0.0_wp
          End If
       End Do
    End Do

    Call T%matrix_map%get_data( m = m )
    Call pdpotrf( 'L', m, T%data, 1, 1, T%matrix_map%get_descriptor(), error )
    If( error == 0 ) Then
       C = T
    End If
    
  End Function choleski_real

  Function choleski_complex( A ) Result( C )

    !! Choleski decompose into lower triangular factors a Hermitian positive definite matrix
    !! Return value is deallocated on error

    Use Scalapack_interfaces, Only : pzpotrf

    Class(         distributed_matrix ), Allocatable :: C

    Class( complex_distributed_matrix ), Intent( In ) :: A

    Type( complex_distributed_matrix ) :: T

    Integer :: m
    Integer :: i_glob, j_glob
    Integer :: i, j
    Integer :: error

    T = A
    
    ! Zero Upper half of T
    Do j = 1, Size( T%data, Dim = 2 )
       j_glob = T%local_to_global_cols( j )
       Do i = 1, Size( T%data, Dim = 1 )
          i_glob = T%local_to_global_rows( i )
          If( j_glob > i_glob ) Then
             T%data( i, j ) = 0.0_wp
          End If
       End Do
    End Do

    Call T%matrix_map%get_data( m = m )
    Call pzpotrf( 'L', m, T%data, 1, 1, T%matrix_map%get_descriptor(), error )
    If( error == 0 ) Then
       C = T
    End If
    
  End Function choleski_complex

  ! Triangular invert Routines
  
  Function tr_inv_real( A ) Result( C )

    !! Invert a lower triangular real matrix
    !! Return value is deallocated on error

    Use Scalapack_interfaces, Only : pdtrtri
    
    Class(      distributed_matrix ), Allocatable :: C

    Class( real_distributed_matrix ), Intent( In ) :: A

    Type( real_distributed_matrix ) :: T
    
    Integer :: m
    Integer :: error

    T = A
    
    Call T%matrix_map%get_data( m = m )
    Call pdtrtri( 'L', 'N', m, T%data, 1, 1, T%matrix_map%get_descriptor(), error )
    If( error == 0 ) Then
       C = T
    End If
    
  End Function tr_inv_real
  
  Function tr_inv_complex( A ) Result( C )

    !! Invert a lower triangular complex matrix
    !! Return value is deallocated on error

    Use Scalapack_interfaces, Only : pztrtri

    Class(         distributed_matrix ), Allocatable :: C

    Class( complex_distributed_matrix ), Intent( In ) :: A

    Type( complex_distributed_matrix ) :: T

    Integer :: m
    Integer :: error

    T = A
    
    Call T%matrix_map%get_data( m = m )
    Call pztrtri( 'L', 'N', m, T%data, 1, 1, T%matrix_map%get_descriptor(), error )
    If( error == 0 ) Then
       C = T
    End If
    
  End Function tr_inv_complex

  ! Unary plus/minus routines
  
  Function plus_real( A ) Result( C )

    !! Unary minus operation
    
    Class(      distributed_matrix ), Allocatable :: C

    Class( real_distributed_matrix ), Intent( In ) :: A

    Type( real_distributed_matrix ) :: T
    
    T = A
    T%data = + T%data
    C = T
    
  End Function plus_real
  
  Function plus_complex( A ) Result( C )

    !! Unary minus operation
    
    Class(         distributed_matrix ), Allocatable :: C

    Class( complex_distributed_matrix ), Intent( In ) :: A

    Type( complex_distributed_matrix ) :: T
    
    T = A
    T%data = + T%data
    C = T
    
  End Function plus_complex
  
  Function minus_real( A ) Result( C )

    !! Unary minus operation
    
    Class(      distributed_matrix ), Allocatable :: C

    Class( real_distributed_matrix ), Intent( In ) :: A

    Type( real_distributed_matrix ) :: T
    
    T = A
    T%data = - T%data
    C = T
    
  End Function minus_real
  
  Function minus_complex( A ) Result( C )

    !! Unary minus operation
    
    Class(         distributed_matrix ), Allocatable :: C

    Class( complex_distributed_matrix ), Intent( In ) :: A

    Type( complex_distributed_matrix ) :: T
    
    T = A
    T%data = - T%data
    C = T
    
  End Function minus_complex

  ! Matrix extract routines
  
  Function extract_real( A, m, n, p, q ) Result( B )

    !! Extract a pach of one matrix into another matrix

    Use scalapack_interfaces, Only : pdgemr2d, pdgeadd

    Class( distributed_matrix ), Allocatable :: B
    
    Class( real_distributed_matrix ), Intent( In    ) :: A
    Integer                         , Intent( In    ) :: m
    Integer                         , Intent( In    ) :: n
    Integer                         , Intent( In    ) :: p 
    Integer                         , Intent( In    ) :: q

    Type( real_distributed_matrix ) :: T1, T2

    Integer :: mb
    Integer :: nb
    Integer :: a_ctxt

    mb = n - m + 1
    nb = q - p + 1

    Call T1%create( mb, nb, A )
    Call A%matrix_map%get_data( ctxt = a_ctxt )

    If( .Not. A%daggered ) Then

       ! Non-tranposed matrix - just extract the patch
       Call pdgemr2d( mb, nb,  A%data, m, p,  A%matrix_map%get_descriptor(), &
                              T1%data, 1, 1, T1%matrix_map%get_descriptor(), a_ctxt )

    Else

       ! Tranposed matrix - extract the transposed patch and then tranpose that to the right layout
       Call T2%create( nb, mb, A )
       Call pdgemr2d( nb, mb,  A%data, p, m,  A%matrix_map%get_descriptor(), &
                              T2%data, 1, 1, T2%matrix_map%get_descriptor(), a_ctxt )

       ! Use the addition routine to perform the transpose
       Call pdgeadd( 'T', mb, nb, 1.0_wp, T2%data, 1, 1, T2%matrix_map%get_descriptor(), &
                                  0.0_wp, T1%data, 1, 1, T1%matrix_map%get_descriptor() )
    End If

    B = T1

  End Function extract_real

  Function extract_complex( A, m, n, p, q ) Result( B )

    !! Extract a pach of one matrix into another matrix

    Use scalapack_interfaces, Only : pzgemr2d, pzgeadd

    Class( distributed_matrix ), Allocatable :: B
    
    Class( complex_distributed_matrix ), Intent( In    ) :: A
    Integer                         , Intent( In    ) :: m
    Integer                         , Intent( In    ) :: n
    Integer                         , Intent( In    ) :: p 
    Integer                         , Intent( In    ) :: q

    Type( complex_distributed_matrix ) :: T1, T2

    Integer :: mb
    Integer :: nb
    Integer :: a_ctxt

    mb = n - m + 1
    nb = q - p + 1

    Call T1%create( mb, nb, A )
    Call A%matrix_map%get_data( ctxt = a_ctxt )

    If( .Not. A%daggered ) Then

       ! Non-daggered matrix - just extract the patch
       Call pzgemr2d( mb, nb,  A%data, m, p,  A%matrix_map%get_descriptor(), &
                              T1%data, 1, 1, T1%matrix_map%get_descriptor(), a_ctxt )

    Else

       ! Daggered matrix - extract the daggered patch and then dagger that to the right layout
       Call T2%create( nb, mb, A )       
       Call pzgemr2d( nb, mb,  A%data, p, m,  A%matrix_map%get_descriptor(), &
                              T2%data, 1, 1, T2%matrix_map%get_descriptor(), a_ctxt )

       ! Use the addition routine to perform the transpose
       Call pzgeadd( 'C', mb, nb, ( 1.0_wp, 0.0_wp ), T2%data, 1, 1, T2%matrix_map%get_descriptor(), &
                                  ( 0.0_wp, 0.0_wp ), T1%data, 1, 1, T1%matrix_map%get_descriptor() )
    End If

    B = T1

  End Function extract_complex


  !##########################################################################
  ! Auxiliary routines
  
  Subroutine set_local_to_global( loc_to_glob, n, nb, myp, np, da )

    !! Generate an array mapping the local indices to global ones
    
    Integer, Dimension( : ), Allocatable, Intent(   Out ) :: loc_to_glob !! The indices
    Integer                             , Intent( In    ) :: n           !! The size of the dimension
    Integer                             , Intent( In    ) :: nb          !! the blocking factor
    Integer                             , Intent( In    ) :: myp         !! My processor coordinate
    Integer                             , Intent( In    ) :: np          !! The number of processors along this dimensions
    Integer                             , Intent( In    ) :: da          !! How big the mapping array need be
    
    Integer :: i_glob, i_loc, skip, start
    
    Allocate( loc_to_glob( 1:da ) )
    
    skip =  np * nb
    
    i_loc = 1
    start = myp * nb + 1
    Do While( start <= n )
       Do i_glob = start, Min( start + nb - 1, n )
          loc_to_glob( i_loc ) = i_glob
          i_loc = i_loc + 1
       End Do
       start = start + skip
    End Do
    
  End Subroutine set_local_to_global
  
  Subroutine set_global_to_local( glob_to_loc, n, nb, myp, np )
    
    !! Generate an array mapping the global indices to local ones
    !! If this process does not hold a given local index the entry is set to DISTRIBUTED_MATRIX_NOT_ME

    Integer, Dimension( : ), Allocatable, Intent(   Out ) :: glob_to_loc !! The indices
    Integer                             , Intent( In    ) :: n           !! The size of the dimension
    Integer                             , Intent( In    ) :: nb          !! the blocking factor
    Integer                             , Intent( In    ) :: myp         !! My processor coordinate
    Integer                             , Intent( In    ) :: np          !! The number of processors along this dimensions
    
    Integer :: i_glob, i_loc, skip, start
    
    Allocate( glob_to_loc( 1:n ) )
    
    glob_to_loc = DISTRIBUTED_MATRIX_NOT_ME
    
    skip =  np * nb
    
    i_loc = 1
    start = myp * nb + 1
    Do While( start <= n )
       Do i_glob = start, Min( start + nb - 1, n )
          glob_to_loc( i_glob ) = i_loc
          i_loc = i_loc + 1
       End Do
       start = start + skip
    End Do
    
  End Subroutine set_global_to_local

  Subroutine real_remap( A, is_A_dummy, parent_comm, B, is_B_dummy )

    !! Remap real data to another distribution

    Class( real_distributed_matrix ), Intent( In    ) :: A           !! Source Matrix
    Logical                         , Intent( In    ) :: is_A_dummy  !! If true the source is NOT on this process
    Integer                         , Intent( In    ) :: parent_comm !! A communicator that encompasses all the relevant processes
    Class(      distributed_matrix ), Intent( InOut ) :: B           !! Destination matrix
    Logical                         , Intent( In    ) :: is_B_dummy  !! If true the destination is NOT on this process

    Call B%real_remap( A, is_A_dummy, parent_comm, is_B_dummy )

  End Subroutine real_remap

  Subroutine complex_remap( A, is_A_dummy, parent_comm, B, is_B_dummy )

    !! Remap complex data to another distribution

    Class( complex_distributed_matrix ), Intent( In    ) :: A           !! Source Matrix
    Logical                            , Intent( In    ) :: is_A_dummy  !! If true the source is NOT on this process
    Integer                            , Intent( In    ) :: parent_comm !! A communicator that encompasses all the relevant processes
    Class(         distributed_matrix ), Intent( InOut ) :: B           !! Destination matrix
    Logical                            , Intent( In    ) :: is_B_dummy  !! If true the destination is NOT on this process

    Call B%complex_remap( A, is_A_dummy, parent_comm, is_B_dummy )

  End Subroutine complex_remap

  Subroutine complex_remap_real( A, is_A_dummy, parent_comm, B, is_B_dummy )

    !! Remap complex data to real data held in another distribution - ILLEGAL

    Class( complex_distributed_matrix ), Intent( In    ) :: A           !! Source Matrix
    Logical                            , Intent( In    ) :: is_A_dummy  !! If true the source is NOT on this process
    Integer                            , Intent( In    ) :: parent_comm !! A communicator that encompasses all the relevant processes
    Class(    real_distributed_matrix ), Intent( InOut ) :: B           !! Destination matrix
    Logical                            , Intent( In    ) :: is_B_dummy  !! If true the destination is NOT on this process

    Stop "Illegal combination of arguments in complex_remap_real"
    ! Silly lines to stop compiler warning about unused vars
    Write( *, * ) A%data, is_A_dummy, parent_comm, B%data, is_B_dummy
    
  End Subroutine complex_remap_real
  
  Subroutine real_remap_complex( A, is_A_dummy, parent_comm, B, is_B_dummy )

    !! Remap real data to complex data held in another distribution - ILLEGAL

    Class(    real_distributed_matrix ), Intent( In    ) :: A           !! Source Matrix
    Logical                            , Intent( In    ) :: is_A_dummy  !! If true the source is NOT on this process
    Integer                            , Intent( In    ) :: parent_comm !! A communicator that encompasses all the relevant processes
    Class( complex_distributed_matrix ), Intent( InOut ) :: B           !! Destination matrix
    Logical                            , Intent( In    ) :: is_B_dummy  !! If true the destination is NOT on this process

    Stop "Illegal combination of arguments in real_remap_complex"
    ! Silly lines to stop compiler warning about unused vars
    Write( *, * ) A%data, is_A_dummy, parent_comm, B%data, is_B_dummy
    
  End Subroutine real_remap_complex
  
  Subroutine real_remap_real( A, is_A_dummy, parent_comm, B, is_B_dummy )

    !! Remap a source real matrix in one distribution to another destination real
    !! matrix in another distribution
    ! A fairly horrible thing. The problems occur because the source and destination
    ! matrices may occupy different sets of processes. Thus we have the IS_?_DUMMY
    ! arguments which mean if set this process doesn't (currently) hold any data
    ! associated witht the matrix. If this is the case the ONLY use of the
    ! matrix is its type used to dispatch the approriate method - the routine should
    ! not attempt to read a matrix for which the corresponding IS_?_DUMMY argument is set.

    ! To keep some vague amount of sanity in all this one of A or B MUST be in this process
    ! This allows us to work out what type of thing is being redistributed.

    Use matrix_mapping_module, Only : matrix_mapping_comm_to_base
    Use Scalapack_interfaces , Only : pdgemr2d
    
    Class( real_distributed_matrix ), Intent( In    ) :: A           !! Source Matrix
    Logical                         , Intent( In    ) :: is_A_dummy  !! If true the source is NOT on this process
    Integer                         , Intent( In    ) :: parent_comm !! A communicator that encompasses all the relevant processes
    Class( real_distributed_matrix ), Intent( InOut ) :: B           !! Destination matrix
    Logical                         , Intent( In    ) :: is_B_dummy  !! If true the destination is NOT on this process

    Type( matrix_mapping ) :: mapping

    Real( wp ), Dimension( 1:1, 1:1 ) :: dum_A, dum_B

    Integer, Dimension( 1:9 ) :: desc_A, desc_B
    
    Integer :: m, n
    Integer :: m_A, n_A
    Integer :: m_B, n_B
    Integer :: parent_ctxt

    Logical :: p_A, p_B

    ! Are A or B "present", i.e. contain actual data. Hold over from when I tried to do
    ! this with both optional and allocatable arguments
    p_A = .Not. is_A_dummy
    p_B = .Not. is_B_dummy
    
    If( .Not. p_A .And. .Not. p_B ) Then
       Stop "In matrix_remap_data_real one of A or B must be supplied"
    End If
    
    ! Generate a context fron the parent_communicator
    Call matrix_mapping_comm_to_base( parent_comm, mapping )
    Call mapping%get_data( ctxt = parent_ctxt )

    ! Get sizes and descriptors for the matrices
    ! The redistrib routine use -1 to indicate no data on this process as that is what the Scalapack routine uses
    If( p_A ) Then
       Call A%matrix_map%get_data( m = m_A, n = n_A )
       desc_A = A%matrix_map%get_descriptor()
    Else
       m_A    = -1
       n_A    = -1
       desc_A = -1
    End If
    
    If( p_B ) Then
       Call B%matrix_map%get_data( m = m_B, n = n_B )
       desc_B = B%matrix_map%get_descriptor()
    Else
       m_B    = -1
       n_B    = -1
       desc_B = -1
    End If

    If( m_A /= -1 .And. n_A /= -1 .And. m_B /= -1 .And. n_B /= - 1 ) Then
       If( m_A /= m_B .Or. n_A /= n_B ) Then
          Stop "Inconsistent matrix sizes in matrix_remap_data_real"
       End If
    End If

    ! Work out the dimensions of the matrix we are redistributing
    If( p_A ) Then
       m = m_A
       n = n_A
    Else If( p_B ) Then
       m = m_B
       n = n_B
    Else
       Stop "In matrix_remap_data_real got to an impossible place!"
    End If

    ! Call the redistribution routine supplying dummy arrays as required
    If     (       p_A .And.       p_B ) Then
       Call pdgemr2d( m, n, A%data, 1, 1, desc_A, B%data, 1, 1, desc_B, parent_ctxt )

    Else If(       p_A .And. .Not. p_B ) Then
       Call pdgemr2d( m, n, A%data, 1, 1, desc_A, dum_B , 1, 1, desc_B, parent_ctxt )

    Else If( .Not. p_A .And.       p_B ) Then
       Call pdgemr2d( m, n, dum_A , 1, 1, desc_A, B%data, 1, 1, desc_B, parent_ctxt )

    Else If( .Not. p_A .And. .Not. p_B ) Then
       ! Shouldn't get here due to error check above
       Stop "In matrix_remap_data_real got to an impossible place!"

    End If

    ! Free the context we created
    Call mapping%free()
    
  End Subroutine real_remap_real
  
  Subroutine complex_remap_complex( A, is_A_dummy, parent_comm, B, is_B_dummy )

    !! Remap a source complex matrix in one distribution to another destination complex
    !! matrix in another distribution
    ! A fairly horrible thing. The problems occur because the source and destination
    ! matrices may occupy different sets of processes. Thus we have the IS_?_DUMMY
    ! arguments which mean if set this process doesn't (currently) hold any data
    ! associated witht the matrix. If this is the case the ONLY use of the
    ! matrix is its type used to dispatch the approriate method - the routine should
    ! not attempt to read a matrix for which the corresponding IS_?_DUMMY argument is set.

    ! To keep some vague amount of sanity in all this one of A or B MUST be in this process
    ! This allows us to work out what type of thing is being redistributed.

    Use matrix_mapping_module, Only : matrix_mapping_comm_to_base
    Use Scalapack_interfaces , Only : pzgemr2d

    Class( complex_distributed_matrix ), Intent( In    ) :: A           !! Source Matrix
    Logical                            , Intent( In    ) :: is_A_dummy  !! If true the source is NOT on this process
    Integer                            , Intent( In    ) :: parent_comm !! A communicator that encompasses all the relevant processes
    Class( complex_distributed_matrix ), Intent( InOut ) :: B           !! Destination matrix
    Logical                            , Intent( In    ) :: is_B_dummy  !! If true the destination is NOT on this process
    
    Type( matrix_mapping ) :: mapping

    Complex( wp ), Dimension( 1:1, 1:1 ) :: dum_A, dum_B

    Integer, Dimension( 1:9 ) :: desc_A, desc_B
    
    Integer :: m, n
    Integer :: m_A, n_A
    Integer :: m_B, n_B
    Integer :: parent_ctxt

    Logical :: p_A, p_B

    ! Are A or B "present", i.e. contain actual data. Hold over from when I tried to do
    ! this with both optional and allocatable arguments
    p_A = .Not. is_A_dummy
    p_B = .Not. is_B_dummy
    
    If( .Not. p_A .And. .Not. p_B ) Then
       Stop "In matrix_remap_data_complex one of A or B must be supplied"
    End If
    
    ! Generate a context fron the parent_communicator
    Call matrix_mapping_comm_to_base( parent_comm, mapping )
    Call mapping%get_data( ctxt = parent_ctxt )

    ! Get sizes and descriptors for the matrices
    ! The redistrib routine use -1 to indicate no data on this process as that is what the Scalapack routine uses
    If( p_A ) Then
       Call A%matrix_map%get_data( m = m_A, n = n_A )
       desc_A = A%matrix_map%get_descriptor()
    Else
       m_A    = -1
       n_A    = -1
       desc_A = -1
    End If
    
    If( p_B ) Then
       Call B%matrix_map%get_data( m = m_B, n = n_B )
       desc_B = B%matrix_map%get_descriptor()
    Else
       m_B    = -1
       n_B    = -1
       desc_B = -1
    End If

    If( m_A /= -1 .And. n_A /= -1 .And. m_B /= -1 .And. n_B /= - 1 ) Then
       If( m_A /= m_B .Or. n_A /= n_B ) Then
          Stop "Inconsistent matrix sizes in matrix_remap_data_complex"
       End If
    End If

    ! Work out the dimensions of the matrix we are redistributing
    If( p_A ) Then
       m = m_A
       n = n_A
    Else If( p_B ) Then
       m = m_B
       n = n_B
    Else
       Stop "In matrix_remap_data_real got to an impossible place!"
    End If

    ! Call the redistribution routine supplying dummy arrays as required
    If     (       p_A .And.       p_B ) Then
       Call pzgemr2d( m, n, A%data, 1, 1, desc_A, B%data, 1, 1, desc_B, parent_ctxt )

    Else If(       p_A .And. .Not. p_B ) Then
       Call pzgemr2d( m, n, A%data, 1, 1, desc_A, dum_B , 1, 1, desc_B, parent_ctxt )

    Else If( .Not. p_A .And.       p_B ) Then
       Call pzgemr2d( m, n, dum_A , 1, 1, desc_A, B%data, 1, 1, desc_B, parent_ctxt )

    Else If( .Not. p_A .And. .Not. p_B ) Then
       ! Shouldn't get here due to error check above
       Stop "In matrix_remap_data_complex got to an impossible place!"

    End If

    ! Free the context we created
    Call mapping%free()

  End Subroutine complex_remap_complex
  
End Module distributed_matrix_module
 
