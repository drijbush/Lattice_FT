Module ks_matrix_module

  !! Wrapper module for distributed matrix classes so we can implement an opaque array - 
  !! At this level we only consider 1 k/spin point 

  Use numbers_module           , Only : wp
  Use distributed_matrix_module, Only : distributed_matrix

  Implicit None

  Type, Public :: ks_matrix
     !! Wrapper type for distributed matices, i.e. either real or complex
     Class( distributed_matrix ), Allocatable, Private :: matrix
   Contains
     ! Public methods
     Procedure, Public :: create                 => ks_matrix_create                     !! Create a ks_matrix
     Generic  , Public :: Assignment( = )        => set_real_scalar                      !! Set all elements in the array to a real constant scalar
     Generic  , Public :: Operator( .Dagger. )   => dagger                               !! Dagger a ks_matrix
     Generic  , Public :: Operator( * )          => multiply                             !! Multiply 2 ks_matrix's
     Generic  , Public :: Operator( * )          => rscal_multiply                       !! Pre-multiply by a real scalar
     Generic  , Public :: Operator( * )          => multiply_rscal                       !! Post-multiply by a real scalar
     Generic  , Public :: Operator( + )          => plus                                 !! Unary plus operation
     Generic  , Public :: Operator( + )          => add                                  !! Add 2 ks_matrix's
     Generic  , Public :: Operator( + )          => add_diagonal                         !! Add a ks_matrix to a diagonal matrix
     Generic  , Public :: Operator( + )          => diagonal_add                         !! Add a ks_matrix to a diagonal matrix
     Generic  , Public :: Operator( - )          => minus                                !! Unary minus operation
     Generic  , Public :: Operator( - )          => subtract                             !! Subtract 2 ks_matrix's
     Generic  , Public :: Operator( - )          => subtract_diagonal                    !! Subtract a diagonal matrix from a ks_matrix 
     Generic  , Public :: Operator( - )          => diagonal_subtract                    !! Subtract a ks_matrix from a diagonal matrix
     Procedure, Public :: diag                   => ks_matrix_diag                       !! Diagonalise a ks_matrix
     Generic  , Public :: Operator( .Choleski. ) => choleski                             !! Choleski factor ks_matrix
     Generic  , Public :: Operator( .TrInv.    ) => tr_inv                               !! Invert a lower triangular ks_matrix
     Procedure, Public :: extract                => ks_matrix_extract                    !! Extract a patch of one matrix into another
     Procedure, Public :: size                   => ks_matrix_size                       !! Get the dimensions of the matrix
     Generic  , Public :: set_by_global          => set_global_real, set_global_complex  !! Set elements by global indices
     Generic  , Public :: get_by_global          => get_global_real, get_global_complex  !! Get elements using global indices
     Procedure, Public :: get_comm               => ks_matrix_communicator               !! Get the communicator containing the processes holding the matrix
     Procedure, Public :: global_to_local        => ks_matrix_global_to_local            !! Get an array for mapping global indices to local  ones
     Procedure, Public :: local_to_global        => ks_matrix_local_to_global            !! Get an array for mapping local  indices to global ones
     ! Private implementations
     Procedure,            Private :: set_real_scalar      => ks_matrix_set_real_scalar
     Procedure,            Private :: dagger               => ks_matrix_dagger
     Procedure,            Private :: multiply             => ks_matrix_mult
     Procedure, Pass( A ), Private :: rscal_multiply       => ks_matrix_rscal_mult
     Procedure,            Private :: multiply_rscal       => ks_matrix_mult_rscal
     Procedure,            Private :: plus                 => ks_matrix_plus
     Procedure,            Private :: add                  => ks_matrix_add
     Procedure,            Private :: add_diagonal         => ks_matrix_add_diagonal
     Procedure, Pass( A ), Private :: diagonal_add         => ks_matrix_diagonal_add
     Procedure,            Private :: minus                => ks_matrix_minus
     Procedure,            Private :: subtract             => ks_matrix_subtract
     Procedure,            Private :: subtract_diagonal    => ks_matrix_subtract_diagonal
     Procedure, Pass( A ), Private :: diagonal_subtract    => ks_matrix_diagonal_subtract
     Procedure,            Private :: choleski             => ks_matrix_choleski
     Procedure,            Private :: tr_inv               => ks_matrix_tr_inv
     Procedure,            Private :: set_global_real      => ks_matrix_set_global_real
     Procedure,            Private :: set_global_complex   => ks_matrix_set_global_complex
     Procedure,            Private :: get_global_real      => ks_matrix_get_global_real
     Procedure,            Private :: get_global_complex   => ks_matrix_get_global_complex
  End Type ks_matrix

  Public :: ks_matrix_init          !! Initialise the ks matrix system and optionally sets the default blocking factor
  Public :: ks_matrix_comm_to_base  !! Converts an MPI communicator into the data structures required to describe a ks matrix mapped onto it
  Public :: ks_matrix_remap_data    !! Remap the data held by A onto the distribution described by B
  Public :: ks_matrix_finalise      !! Finalise the matrix system

Contains

  !##########################################################################################################
  ! Non-type bound procedures

  Subroutine ks_matrix_init( nb )

    !! Initialise the matrix system and optionally sets the default blocking factor

    Use distributed_matrix_module, Only : distributed_matrix_init, &
         distributed_matrix_set_default_blocking

    Integer, Intent( In ), Optional :: nb !! Set a default blocking factor

    Call distributed_matrix_init
    
    If( Present( nb ) ) Then
       Call distributed_matrix_set_default_blocking( nb )
    End If

  End Subroutine ks_matrix_init

  Subroutine ks_matrix_comm_to_base( comm, base_matrix )

    !! Generate a base ks_matrix object from an MPI communicator
    
    Use distributed_matrix_module, Only : distributed_matrix_comm_to_base, &
         real_distributed_matrix

    Integer             , Intent( In    ) :: comm
    Type   ( ks_matrix ), Intent(   Out ) :: base_matrix

    Allocate( real_distributed_matrix:: base_matrix%matrix )

    Call distributed_matrix_comm_to_base( comm, base_matrix%matrix )
    
  End Subroutine ks_matrix_comm_to_base

  Subroutine ks_matrix_remap_data( A, parent_communicator, B )

    !! Remap the data held by A onto the distribution described by B
    !! Parent communicator contains both A and B

    ! Issues here because either in the source or remapped
    ! matrix a process may not actually hold any part of it
    ! and so the unallocated actual argument doesn't contain
    ! any information about what it is

    Type   ( ks_matrix ), Allocatable, Intent( In    ) :: A
    Integer             ,              Intent( In    ) :: parent_communicator
    Type   ( ks_matrix ), Allocatable, Intent( InOut ) :: B

    Type( ks_matrix ), Allocatable :: dummy_A, dummy_B

    Logical :: p_A, p_B

    p_A = Allocated( A )
    p_B = Allocated( B )

    ! One of A or B must be present to define the data type of what is being redistributed
    If( .Not. p_A .And. .Not. p_B ) Then
       Stop "Must specify one of the matrices in ks_matrix_remap_data"
    End If

    ! Both matrices MUST be of the same type. Thus if one is not specified
    ! create a dummy with the same type as one that is present
    If( .Not. p_A ) Then
       Allocate( dummy_A, Source = B )
    End If
    If( .Not. p_B ) Then
       Allocate( dummy_B, Source = A )
    End If

    ! Now can call remap data routines, carefully indicating which
    ! arguments are dummies because this process

    If     (       p_A .And.       p_B ) Then
       Call       A%matrix%remap( .False., parent_communicator,        B%matrix, .False. )

    Else If( .Not. p_A .And.       p_B ) Then
       Call dummy_A%matrix%remap( .True. , parent_communicator,        B%matrix, .False. )

    Else If(       p_A .And. .Not. p_B ) Then
       Call       A%matrix%remap( .False. , parent_communicator, dummy_B%matrix, .True.  )

    Else
       Stop "Must specify one of the matrices in ks_matrix_remap_data"
    End If

  End Subroutine ks_matrix_remap_data

  Subroutine ks_matrix_finalise

    !! Finalise the ks matrix system

    Use distributed_matrix_module, Only : distributed_matrix_finalise

    Call distributed_matrix_finalise
    
  End Subroutine ks_matrix_finalise

  !##########################################################################################################
  ! Type bound procedures

  Subroutine ks_matrix_create( A, is_complex, m, n, source_matrix )

    !! Create a distributed ks matrix

    Use distributed_matrix_module, Only : real_distributed_matrix, &
         complex_distributed_matrix

    Class  ( ks_matrix ), Intent(   Out ) :: A
    Logical             , Intent( In    ) :: is_complex
    Integer             , Intent( In    ) :: m
    Integer             , Intent( In    ) :: n
    Type   ( ks_matrix ), Intent( In    ) :: source_matrix

    If( is_complex ) Then
      Allocate( complex_distributed_matrix :: A%matrix ) 
    Else
      Allocate( real_distributed_matrix    :: A%matrix ) 
    End If

    Call A%matrix%create( m, n, source_matrix%matrix )
    
  End Subroutine ks_matrix_create

  Function ks_matrix_dagger( A ) Result( tA )

    !! Form the Hermitian conjugate of the matrix (Tranpose for real data)
    
    Type( ks_matrix ) :: tA

    Class( ks_matrix ), Intent( In ) :: A 

    tA%matrix = .Dagger. A%matrix

  End Function ks_matrix_dagger

  Function ks_matrix_rscal_mult( s, A ) Result( C )

    !! Multiply two matrices together
    
    Type( ks_matrix ) :: C
    
    Real( wp )        , Intent( In ) :: s
    Class( ks_matrix ), Intent( In ) :: A

    C%matrix = s * A%matrix

  End Function ks_matrix_rscal_mult

  Function ks_matrix_mult_rscal( A, s ) Result( C )

    !! Multiply two matrices together
    
    Type( ks_matrix ) :: C
    
    Class( ks_matrix ), Intent( In ) :: A
    Real( wp )        , Intent( In ) :: s

    C%matrix = s * A%matrix

  End Function ks_matrix_mult_rscal

  Function ks_matrix_mult( A, B ) Result( C )

    !! Multiply two matrices together
    
    Type( ks_matrix ) :: C

    Class( ks_matrix ), Intent( In ) :: A
    Type ( ks_matrix ), Intent( In ) :: B

    C%matrix = A%matrix * B%matrix

  End Function ks_matrix_mult

  Function ks_matrix_add( A, B ) Result( C )

    !! Add two matrices together
    
    Type( ks_matrix ) :: C

    Class( ks_matrix ), Intent( In ) :: A
    Type ( ks_matrix ), Intent( In ) :: B

    C%matrix = A%matrix + B%matrix

  End Function ks_matrix_add

  Function ks_matrix_add_diagonal( A, d ) Result( C )

    !! Add a matrix to a diagonal amtrix
    
    Type( ks_matrix ) :: C

    Class( ks_matrix )        , Intent( In ) :: A
    Real( wp ), Dimension( : ), Intent( In ) :: d

    C%matrix = A%matrix + d

  End Function ks_matrix_add_diagonal

  Function ks_matrix_diagonal_add( d, A ) Result( C )

    !! Add a matrix to a diagonal amtrix
    
    Type( ks_matrix ) :: C

    Real( wp ), Dimension( : ), Intent( In ) :: d
    Class( ks_matrix )        , Intent( In ) :: A

    C%matrix = d + A%matrix

  End Function ks_matrix_diagonal_add

  Function ks_matrix_subtract( A, B ) Result( C )

    !! Subtract two matrices together
    
    Type( ks_matrix ) :: C

    Class( ks_matrix ), Intent( In ) :: A
    Type ( ks_matrix ), Intent( In ) :: B

    C%matrix = A%matrix - B%matrix

  End Function ks_matrix_subtract

  Function ks_matrix_subtract_diagonal( A, d ) Result( C )

    !! Subtract a diagonal amtrix from a matrix
    
    Type( ks_matrix ) :: C

    Class( ks_matrix )        , Intent( In ) :: A
    Real( wp ), Dimension( : ), Intent( In ) :: d

    C%matrix = A%matrix - d

  End Function ks_matrix_subtract_diagonal

  Function ks_matrix_diagonal_subtract( d, A ) Result( C )

    !! Subtract a matrix from a diagonal amtrix
    
    Type( ks_matrix ) :: C

    Real( wp ), Dimension( : ), Intent( In ) :: d
    Class( ks_matrix )        , Intent( In ) :: A

    C%matrix = d - A%matrix

  End Function ks_matrix_diagonal_subtract

  Subroutine ks_matrix_diag( A, Q, E ) 

    !! Diagonalise a (assumed Hermitian) ks matrix
    
    Class( ks_matrix ),                      Intent( In    ) :: A
    Type ( ks_matrix ),                      Intent(   Out ) :: Q
    Real( wp ), Dimension( : ), Allocatable, Intent(   Out ) :: E

    ! Need to give a Type to Q as Intent( Out ) will have deallocatd all its components
    Q = A

    Call A%matrix%diag( Q%matrix, E )

  End Subroutine ks_matrix_diag

  Function ks_matrix_choleski( A ) Result( C )

    !! Choleski decompose A

    Type( ks_matrix ) :: C
    
    Class( ks_matrix ), Intent( In ) :: A

    C%matrix = .Choleski. A%matrix
    
  End Function ks_matrix_choleski

  Function ks_matrix_tr_inv( A ) Result( C )

    !! Choleski decompose A

    Type( ks_matrix ) :: C
    
    Class( ks_matrix ), Intent( In ) :: A

    C%matrix = .TrInv. A%matrix
    
  End Function ks_matrix_tr_inv

  Function ks_matrix_plus( A ) Result( C )

    !! Unary plus

    Type( ks_matrix ) :: C
    
    Class( ks_matrix ), Intent( In ) :: A

    C%matrix = + A%matrix
    
  End Function ks_matrix_plus

  Function ks_matrix_minus( A ) Result( C )

    !! Unary minus

    Type( ks_matrix ) :: C
    
    Class( ks_matrix ), Intent( In ) :: A

    C%matrix = - A%matrix
    
  End Function ks_matrix_minus

  Function ks_matrix_size( A, dim ) Result( n )

    !! Return the dimensions of the matrix

    Integer :: n
    
    Class( ks_matrix ), Intent( In ) :: A
    Integer           , Intent( In ) :: dim

    n = A%matrix%size( dim )
    
  End Function ks_matrix_size

  Function ks_matrix_extract( A, m, n, p, q ) Result( C )

    !! Choleski decompose A

    Type( ks_matrix ) :: C
    
    Class( ks_matrix ), Intent( In ) :: A
    Integer           , Intent( In ) :: m
    Integer           , Intent( In ) :: n
    Integer           , Intent( In ) :: p
    Integer           , Intent( In ) :: q

    C%matrix = A%matrix%extract( m, n, p, q )
    
  End Function ks_matrix_extract

  Subroutine ks_matrix_set_real_scalar( A, data )

    !! Set the whole of the matrix to a constant real scalar
    Class( ks_matrix ), Intent( InOut ) :: A
    Real( wp )        , Intent( In    ) :: data

    A%matrix = data
    
  End Subroutine ks_matrix_set_real_scalar

  Subroutine ks_matrix_set_global_real( A, m, n, p, q, data )

    !! Set the (m:n,p:q) patch of the matrix using global indices
    
    Class( ks_matrix )             , Intent( InOut ) :: A
    Integer                        , Intent( In    ) :: m
    Integer                        , Intent( In    ) :: n
    Integer                        , Intent( In    ) :: p
    Integer                        , Intent( In    ) :: q
    Real( wp ), Dimension( m:, p: ), Intent( In    ) :: data

    Call A%matrix%set_by_global( m, n, p, q, data )

  End Subroutine ks_matrix_set_global_real

  Subroutine ks_matrix_set_global_complex( A, m, n, p, q, data )
    
    !! Set the (m:n,p:q) patch of the matrix using global indices

    Class( ks_matrix )                , Intent( InOut ) :: A
    Integer                           , Intent( In    ) :: m
    Integer                           , Intent( In    ) :: n
    Integer                           , Intent( In    ) :: p
    Integer                           , Intent( In    ) :: q
    Complex( wp ), Dimension( m:, p: ), Intent( In    ) :: data

    Call A%matrix%set_by_global( m, n, p, q, data )

  End Subroutine ks_matrix_set_global_complex

  Subroutine ks_matrix_get_global_real( A, m, n, p, q, data )

    !! Get the (m:n,p:q) patch of the matrix using global indices
   
    Class( ks_matrix )             , Intent( In    ) :: A
    Integer                        , Intent( In    ) :: m
    Integer                        , Intent( In    ) :: n
    Integer                        , Intent( In    ) :: p
    Integer                        , Intent( In    ) :: q
    Real( wp ), Dimension( m:, p: ), Intent(   Out ) :: data

    Call A%matrix%get_by_global( m, n, p, q, data )

  End Subroutine ks_matrix_get_global_real

  Subroutine ks_matrix_get_global_complex( A, m, n, p, q, data )
    
    !! Get the (m:n,p:q) patch of the matrix using global indices

    Class( ks_matrix )                , Intent( In    ) :: A
    Integer                           , Intent( In    ) :: m
    Integer                           , Intent( In    ) :: n
    Integer                           , Intent( In    ) :: p
    Integer                           , Intent( In    ) :: q
    Complex( wp ), Dimension( m:, p: ), Intent(   Out ) :: data

    Call A%matrix%get_by_global( m, n, p, q, data )

  End Subroutine ks_matrix_get_global_complex

  Function ks_matrix_communicator( A ) Result( c )

    !! Get the communicator containing the processes holding the matrix

    Integer :: c

    Class( ks_matrix ), Intent( In ) :: A

    c = A%matrix%get_comm()
    
  End Function ks_matrix_communicator

  Function ks_matrix_global_to_local( A, what ) Result( gl_indexing )

    !! For the given matrix get an array for mapping global indices to local  ones
    !! Get row mapping arrays if WHAT='C' or 'c', row ones if it is equal to 'R' or 'r'

    Integer, Dimension( : ), Allocatable :: gl_indexing
    
    Class( ks_matrix   ), Intent( In ) :: A
    Character( Len = * ), Intent( In ) :: what

    gl_indexing = A%matrix%global_to_local( what )

  End Function ks_matrix_global_to_local
  
  Function ks_matrix_local_to_global( A, what ) Result( lg_indexing )

    !! For the given matrix get an array for mapping local indices to global  ones
    !! Get row mapping arrays if WHAT='C' or 'c', row ones if it is equal to 'R' or 'r'

    Integer, Dimension( : ), Allocatable :: lg_indexing
    
    Class( ks_matrix   ), Intent( In ) :: A
    Character( Len = * ), Intent( In ) :: what

    lg_indexing = A%matrix%local_to_global( what )

  End Function ks_matrix_local_to_global

End Module ks_matrix_module
