Module ks_array_module

  ! IRREPS NEED MORE THOUGHT
  !-------------------------
  ! Currently assume 1 - i.e. nosymada
  
  Use numbers_module  , Only : wp
  Use ks_matrix_module, Only : ks_matrix

  Implicit None

  Integer, Public, Parameter :: K_POINT_REAL = 0
  Integer, Public, Parameter :: K_POINT_COMPLEX = 1
  Integer, Public, Parameter :: K_POINT_NOT_EXIST = Huge( K_POINT_NOT_EXIST )

  Integer, Private, Parameter :: INVALID = -1
  Integer, Private, Parameter :: NOT_ME  = -2
  Integer, Private, Parameter :: NO_DATA = -3

  Type, Public :: ks_point_info
     !! A type to hold the basic info about a k point
     Integer                             , Public :: k_type
     Integer                             , Public :: spin
     Integer, Dimension( : ), Allocatable, Public :: k_indices
  End type ks_point_info

  Type, Private :: k_point_matrices
     !! A wrapper for data at a k point - will eventually be used to create arrays when we deal with irreps
     ! External label - e.g. irrep number
     Integer           :: label = INVALID
     Type( ks_matrix ) :: matrix
  End type k_point_matrices

  Type, Private :: k_point
     !! This wraps the info and data for a given ks point
     Type( ks_point_info     )                              :: info
     ! This splitting allows irreps within this k point
     Type( k_point_matrices  ), Dimension( : ), Allocatable :: data
     ! Want to hide eventually
     Integer                                                :: communicator = INVALID
  End type k_point

  Type, Public :: ks_array
     !! An array of ks matrices, operations on which are (almost always) independent and hence parallelisable
     Type( ks_point_info ), Dimension( : ), Allocatable, Private :: all_k_point_info                !! Info about all ks points held by the array
     ! this splitting allows multiple k points on this process
     Type( k_point       ), Dimension( : ), Allocatable, Private :: my_k_points                     !! Info and data for ks points held (in part) by this process
     Integer                                           , Private :: parent_communicator = INVALID   !! A parent communicator spanning all processes involved in the array
     Integer                                           , Private :: iterator_value      = INVALID   !! A value for the iterator - might need to rething as retainined on assignment of ks_arrays
   Contains
     ! Public Methods
     Procedure, Public :: create                  => ks_array_create                   !! Create a ks_array
     Procedure, Public :: split_ks                => ks_array_split_ks                 !! Split the distribution so k point //ism can be used
     Procedure, Public :: join_ks                 => ks_array_join_ks                  !! Rejoin a split_ks array into all k point mode
     Generic  , Public :: Operator( .Dagger. )    => dagger                            !! Dagger each element in the array
     Generic  , Public :: Operator( * )           => multiply                          !! Multiply each element of the array with the corresponding element in another array
     Generic  , Public :: Operator( * )           => rscal_multiply                    !! Pre-multiply by a real scalar
     Generic  , Public :: Operator( * )           => multiply_rscal                    !! Post-multiply by a real scalar
     Generic  , Public :: Operator( + )           => plus                              !! Unary plus operator
     Generic  , Public :: Operator( + )           => add                               !! Add each element of the array with the corresponding element in another array
     Generic  , Public :: Operator( + )           => add_diagonal                      !! Add each element of the array to a diagonal matrix
     Generic  , Public :: Operator( + )           => diagonal_add                      !! Add each element of the array to a diagonal matrix
     Generic  , Public :: Operator( - )           => minus                             !! Unary minus operator
     Generic  , Public :: Operator( - )           => subtract                          !! Subtract each element of the array with the corresponding element in another array
     Generic  , Public :: Operator( - )           => subtract_diagonal                 !! Subtract a diagonal matrix from a general one
     Generic  , Public :: Operator( - )           => diagonal_subtract                 !! Subtract a general matrix from a diagonal one
     Procedure, Public :: diag                    => ks_array_diag                     !! Diagonalise each matrix
     Generic  , Public :: Operator( .Choleski. )  => choleski                          !! Choleski decompose a matrix
     Generic  , Public :: Operator( .TrInv. )     => tr_inv                            !! Invert a lower traingular set of matrices
     Procedure, Public :: extract                 => ks_array_extract                  !! Extract a patch from the matrices and return a new ks_array holding it
     Generic  , Public :: set_by_global           => set_by_global_r, set_by_global_c  !! Set patches of an element
     Generic  , Public :: get_by_global           => get_by_global_r, get_by_global_c  !! Get patches of an element
     Procedure, Public :: global_to_local         => ks_array_g_to_l                   !! Get the global -> local  mapping arrays
     Procedure, Public :: local_to_global         => ks_array_l_to_g                   !! Get the local  -> global mapping arrays
     Procedure, Public :: size                    => ks_array_size                     !! Get the size of the matrix along a given dimension
     Procedure, Public :: print_info              => ks_array_print_info               !! print info about a KS_array
     Procedure, Public :: iterator_init           => ks_array_iterator_init            !! Initialise an iterator
     Procedure, Public :: iterator_reset          => ks_array_iterator_reset           !! Reset an iterator
     Procedure, Public :: iterator_next           => ks_array_iterator_next            !! Move to the next matrix in the ks_array
     Procedure, Public :: iterator_previous       => ks_array_iterator_previous        !! Move to the previous matrix in the ks_array
     ! Private implementations
     Procedure,            Private :: get_all_ks_index
     Procedure,            Private :: get_my_ks_index
     Procedure,            Private :: get_ks
     Procedure,            Private :: set_by_global_r      => ks_array_set_global_real
     Procedure,            Private :: set_by_global_c      => ks_array_set_global_complex
     Procedure,            Private :: get_by_global_r      => ks_array_get_global_real
     Procedure,            Private :: get_by_global_c      => ks_array_get_global_complex
     Procedure,            Private :: dagger               => ks_array_dagger
     Procedure,            Private :: multiply             => ks_array_mult
     Procedure, Pass( A ), Private :: rscal_multiply       => ks_array_rscal_mult
     Procedure,            Private :: multiply_rscal       => ks_array_mult_rscal
     Procedure,            Private :: plus                 => ks_array_plus
     Procedure,            Private :: add                  => ks_array_add
     Procedure,            Private :: add_diagonal         => ks_array_add_diagonal
     Procedure, Pass( A ), Private :: diagonal_add         => ks_array_diagonal_add
     Procedure,            Private :: minus                => ks_array_minus
     Procedure,            Private :: subtract             => ks_array_subtract
     Procedure,            Private :: subtract_diagonal    => ks_array_subtract_diagonal
     Procedure, Pass( A ), Private :: diagonal_subtract    => ks_array_diagonal_subtract
     Procedure,            Private :: choleski             => ks_array_choleski
     Procedure,            Private :: tr_inv               => ks_array_tr_inv
  End type ks_array
  
  Type, Public :: ks_eval_storage
     Type( ks_point_info )                   :: ks_point
     Real( wp ), Dimension( : ), Allocatable :: evals
  End type ks_eval_storage

  Public :: ks_array_init         !! Initalise the KS arrays
  Public :: ks_array_comm_to_base !! Turn an MPI communicator inot a base KS_array object
  Public :: ks_array_finalise     !! Finalise the KS array mechanism
  
  Private

Contains

  Subroutine ks_array_init( nb )

    !! Initalise the KS arrays
    
    Use ks_matrix_module, Only : ks_matrix_init
    
    Integer, Intent( In ), Optional :: nb !! Set a default blocking factor

    If( .Not. Present( nb ) ) Then
       Call ks_matrix_init
    Else
       Call ks_matrix_init( nb )
    End If

  End Subroutine ks_array_init

  Subroutine ks_array_comm_to_base( comm, n_spin, k_point_type, k_points, base_ks_array )

    !! Turn an MPI communicator inot a base KS_array object
    !! from which other KS_arrays can be created.
    !! Note each element of the array will be distributed across all processes in the communicator
    
    Use ks_matrix_module, Only : ks_matrix_comm_to_base
    
    Integer                   , Intent( In    ) :: comm             !! The communicator
    Integer                   , Intent( In    ) :: n_spin           !! The number of spins
    Integer, Dimension( :    ), Intent( In    ) :: k_point_type     !! the type of each k point
    Integer, Dimension( :, : ), Intent( In    ) :: k_points         !! The labels for each k point
    Type( ks_array )          , Intent(   Out ) :: base_ks_array    !! The resulting ks_array

    Type( ks_matrix ) :: base_matrix
    
    Integer :: n_k_points
    Integer :: s, k, ks
    
    n_k_points = Size( k_point_type)

    ! Set up the all k point data structure
    Allocate( base_ks_array%all_k_point_info( 1:n_spin * n_k_points ) )
    ks = 0
    Do s = 1, n_spin
       Do k = 1, n_k_points
          ks = ks + 1
          base_ks_array%all_k_point_info( ks )%k_type    = k_point_type( k )
          base_ks_array%all_k_point_info( ks )%spin      = s
          base_ks_array%all_k_point_info( ks )%k_indices = k_points( :, k )
       End Do
    End Do

    Call ks_matrix_comm_to_base( comm, base_matrix ) 

    ! Now my k points
    Allocate( base_ks_array%my_k_points( 1:n_spin * n_k_points ) )
    ks = 0
    Do s = 1, n_spin
       Do k = 1, n_k_points
          ks = ks + 1
          base_ks_array%my_k_points( ks )%info%k_type    = k_point_type( k )
          base_ks_array%my_k_points( ks )%info%spin      = s
          base_ks_array%my_k_points( ks )%info%k_indices = k_points( :, k )
          Allocate( base_ks_array%my_k_points( ks )%data( 1:1 ) )
          base_ks_array%my_k_points( ks )%data( 1 )%label = 1
          base_ks_array%my_k_points( ks )%data( 1 )%matrix = base_matrix
          base_ks_array%my_k_points( ks )%communicator = base_matrix%get_comm()
       End Do
    End Do

    base_ks_array%parent_communicator = base_matrix%get_comm()

  End Subroutine ks_array_comm_to_base

  Subroutine ks_array_finalise

    !! Finalise the KS array mechanism

    Use ks_matrix_module, Only : ks_matrix_finalise

    Call ks_matrix_finalise
    
  End Subroutine ks_array_finalise

  Subroutine ks_array_create( A, m, n, source )

    !! Create a KS array

    ! M and N should be arrays to allow different sizes at each ks point,
    ! or more likely irrep!!!!!!!!
    
    ! Create a matrix in all k point mode with no irreps parallelism
    
    ! Also want to think about what kind of object should be used as a source
    
    Class( ks_array ), Intent(   Out ) :: A
    Integer          , Intent( In    ) :: m
    Integer          , Intent( In    ) :: n
    Type ( ks_array ), Intent( In    ) :: source

    Integer :: n_all_ks, n_my_ks
    Integer :: ks
    
    ! Set up the all k point data structure
    n_all_ks = Size( source%all_k_point_info )
    Allocate( A%all_k_point_info( 1:n_all_ks ) )
    A%all_k_point_info = source%all_k_point_info

    ! Now my k points
    n_my_ks = Size( source%my_k_points )
    Allocate( A%my_k_points( 1:n_my_ks ) )
    Do ks = 1, n_my_ks
       A%my_k_points( ks )%info = source%my_k_points( ks )%info
       Allocate( A%my_k_points( ks )%data( 1:1 ) )
       A%my_k_points( ks )%data( 1 )%label = 1
       If( m /= NO_DATA .And. n /= NO_DATA ) Then
          Call A%my_k_points( ks )%data( 1 )%matrix%create( &
               A%my_k_points( ks )%info%k_type == K_POINT_COMPLEX, &
               m, n, source%my_k_points( ks )%data( 1 )%matrix )
       End If
       A%my_k_points( ks )%communicator = source%my_k_points( ks )%communicator 
    End Do

    A%parent_communicator = source%parent_communicator

    A%iterator_value = INVALID

  End Subroutine ks_array_create

  Subroutine ks_array_print_info( A, name, verbosity )

    !! Print information about an array, optionally with a title, and if the verbosity is high
    !! enough report (>99) also how it is distributed across the processes

    Use mpi, Only : MPI_Comm_rank, MPI_Comm_size, MPI_Allreduce, MPI_IN_PLACE, MPI_INTEGER, MPI_SUM

    Class( ks_array     ),           Intent( In ) :: A
    Character( Len = * ) , Optional, Intent( In ) :: name
    Integer              , Optional, Intent( In ) :: verbosity

    Integer, Dimension( :, : ), Allocatable :: dims
    
    Integer, Dimension( : ), Allocatable :: this_comm_ranks
    Integer, Dimension( : ), Allocatable :: packed_ranks
    Integer, Dimension( : ), Allocatable :: k
    
    Integer :: me, nproc, me_ks
    Integer :: error
    Integer :: n_ks
    Integer :: n_packed_ranks
    Integer :: ks, my_ks
    Integer :: s
    Integer :: i

    Logical :: in_range

    Call MPI_Comm_rank( A%parent_communicator, me   , error )
    Call MPI_Comm_size( A%parent_communicator, nproc, error )

    n_ks = Size( A%all_k_point_info )

    ! Get the size of each matrix - if the matrix is distributed in split mode
    ! we won't neccessarily know it for each k point so gather the data from the zero
    ! rank process resonsible at each k/spin
    Allocate( dims( 1:2, 1:n_ks ) )
    dims = 0
    Do ks = 1, n_ks
       ! Work out if I hold this k/spin point
       my_ks = A%get_my_ks_index( ks )
       If( my_ks /= NOT_ME ) Then
          k = A%all_k_point_info( ks )%k_indices
          s = A%all_k_point_info( ks )%spin
          ! I hold this k_point. If I am the zero ranked process store the data
          Call MPI_Comm_rank( A%my_k_points( my_ks )%communicator, me_ks, error )
          If( me_ks == 0 ) Then
             dims( 1, ks ) = A%size( k, s, 1 )
             dims( 2, ks ) = A%size( k, s, 2 )
          End If
       End If
    End Do
    ! Now replicate the data
    Call MPI_Allreduce( MPI_IN_PLACE, dims, Size( dims ), MPI_INTEGER, MPI_SUM, A%parent_communicator, error )
    
    If( me == 0 ) Then
       Write( *, * )
       If( Present( name ) ) Then
          Write( *, '( a, a )' ) 'Information on ks_array ', name
       Else
          Write( *, '( a )' ) 'Information on a ks_array'
       End If
       Write( *, '( a )' ) 'This ks array holds the following spins and k points'
       Write( *, '( a, t10, a, t30, a, t50, a )' ) 'Spin', 'Indices', 'Data Type', 'Dimensions'
       Do ks = 1, n_ks
          Write( *, '( i0, t10, "( ", i2, ", ", i2, ", ", i2, " )", t30, a, t45, i8, "x", i0 )' ) &
               A%all_k_point_info( ks )%spin, A%all_k_point_info( ks )%k_indices,                 &
               Merge( 'Real   ', 'Complex', A%all_k_point_info( ks )%k_type == K_POINT_REAL ),    &
               dims( :, ks )
       End Do
    End If

    ! Print Which procs hold which k point
    If( Present( verbosity ) ) Then
       If( verbosity > 99 ) Then
          If( me == 0 ) Then
             Write( *, '( a )' ) 'The ranks within the parent communicator map onto the k points as below:'
          End If
          ! For each k point in turn work out who owns it
          Do ks = 1, n_ks
             Allocate( this_comm_ranks( 0:nproc - 1 ) )
             this_comm_ranks = 0
             my_ks = A%get_my_ks_index( ks )
             If( my_ks /= NOT_ME ) Then
                this_comm_ranks( me ) = me + 1 ! Avoid problems when me == 0
             End If
             Call MPI_Allreduce( MPI_IN_PLACE, this_comm_ranks, Size( this_comm_ranks ), &
                  MPI_INTEGER, MPI_SUM, A%parent_communicator, error )
             If( me == 0 ) Then
                n_packed_ranks = Size( Pack( this_comm_ranks - 1, this_comm_ranks /= 0 ) )
                Allocate( packed_ranks( 0:n_packed_ranks  -1 ) )
                packed_ranks = Pack( this_comm_ranks - 1, this_comm_ranks /= 0 )
                Write( *, '( i0, t10, "( ", i2, ", ", i2, ", ", i2, " )", 3x )', Advance = 'No' )    &
                     A%all_k_point_info( ks )%spin, A%all_k_point_info( ks )%k_indices
                Write( *, '( i0 )', Advance = 'No' ) packed_ranks( 0 )
                If( Size( packed_ranks ) > 1 ) Then
                   in_range = .False.
                   Do i = Lbound( packed_ranks, Dim = 1 ) + 1, Ubound( packed_ranks, Dim = 1 ) - 1
                      If( packed_ranks( i ) - packed_ranks( i - 1 ) == 1 ) Then
                         If( .Not. in_range ) Then
                            Write( *, '( "-" )', Advance = 'No' )
                         End If
                         in_range = .True.
                      Else
                         If( in_range ) Then
                            Write( *, '( i0, ",", i0 )', Advance = 'No' ) &
                                 packed_ranks( i - 1 ), packed_ranks( i )
                            in_range = .False.
                         Else
                            Write( *, '( ",", i0, "," )', Advance = 'No' ) packed_ranks( i - 1 )
                         End If
                      End If
                   End Do
                   If( in_range ) Then
                      Write( *, '( i0 )' ) packed_ranks( i )
                   Else
                      Write( *, '( ",", i0 )' ) packed_ranks( i )
                   End If
                Else
                   Write( *, * )
                End If
                Deallocate( packed_ranks )
             End If
             Deallocate( this_comm_ranks )
          End Do
          If( me == 0 ) Then
             Write( *, * )
          End If
       End If
    End If
    
  End Subroutine ks_array_print_info

  Subroutine ks_array_split_ks( A, complex_weight, split_A, redistribute )

    !! Split a ks_array A so the resulting ks_array is ks point distributed

    Use mpi             , Only : MPI_Comm_size, MPI_Comm_rank, MPI_Comm_split, MPI_UNDEFINED
    Use ks_matrix_module, Only : ks_matrix_comm_to_base, ks_matrix_remap_data
    
    Class( ks_array     ), Intent( In    ) :: A                 !! The array to be split
    Real ( wp           ), Intent( In    ) :: complex_weight    !! The cost of complex points relative to real ones
    Type ( ks_array     ), Intent(   Out ) :: split_A           !! The k point distributed matrix
    Logical, Optional    , Intent( In    ) :: redistribute      !! By default the data associated with the matrix is redistributed

    Type( ks_matrix ), Allocatable :: this_ks_matrix
    Type( ks_matrix ), Allocatable :: split_ks_matrix

    Type( ks_matrix ) :: base_matrix

    Real( wp ), Dimension( : ), Allocatable :: weights

    Integer, Dimension( : ), Allocatable :: n_procs_ks
    Integer, Dimension( : ), Allocatable :: my_ks
    Integer, Dimension( : ), Allocatable :: this_k_indices

    Integer :: m, n
    Integer :: n_ks
    Integer :: n_procs_parent, me_parent, my_colour, k_comm, n_my_ks, n_left
    Integer :: this_k_type, this_s
    Integer :: top_rank
    Integer :: cost
    Integer :: error
    Integer :: ks, this_ks

    Logical :: loc_redist
    Logical :: is_split

    ! If A is already split across ks points no need to do very much except copy ...
    is_split = Size( A%all_k_point_info ) /= Size( A%my_k_points )
    If( is_split ) Then
       split_A = A
       ! .. and make sure the iterator in the new ks_array is reset
       Call split_A%iterator_reset()
       Return
    End If

    If( Present( redistribute ) ) Then
       loc_redist = redistribute
    Else
       loc_redist =.True.
    End If

    ! Set up stuff relevant to all k points
    split_A%all_k_point_info    = A%all_k_point_info
    split_A%parent_communicator = A%parent_communicator

    ! Now split the k points and return the new structure in split_A

    ! First split the communicator
    ! Generate an array for the weights
    n_ks = Size( split_A%all_k_point_info )
    Allocate( weights( 1:n_ks ) )
    Do ks = 1, n_ks
       weights( ks ) = Merge( 1.0_wp, complex_weight, split_A%all_k_point_info( ks )%k_type == K_POINT_REAL )
    End Do
    ! Makes a difference if somehow somebody has all k points complex
    weights = weights / Minval( weights )

    !THIS WHOLE SPLITTING STRATEGY PROBABLY NEEDS MORE THOUGHT
    Call MPI_Comm_size( split_A%parent_communicator, n_procs_parent, error )
    Call MPI_Comm_rank( split_A%parent_communicator,      me_parent, error )
    cost = Nint( Sum( weights ) )
    ! Scale up weights so fits inside the parent comm
    Allocate( n_procs_ks( 1:n_ks ) )
    n_procs_ks = Nint( weights * ( n_procs_parent / cost ) )

    ! Two possible cases

    k_split_strategy: If( cost <= n_procs_parent ) Then
       
       ! 1) There are sufficent processors for each k point to have its own, separate set if processors which work on it
       n_my_ks  = 1
       Allocate( my_ks( 1:n_my_ks ) ) 

       ! Decide which group ( if any ) I am in can probably write this more neatly but lets keep
       ! it explicit as I'm getting a little confused
       If( me_parent > Sum( n_procs_ks ) - 1 ) Then
          my_colour  = MPI_UNDEFINED
          my_ks( 1 ) = INVALID
          n_my_ks    = 0
       Else
          top_rank = 0
          Do ks = 1, n_ks
             top_rank = top_rank + n_procs_ks( ks )
             If( me_parent < top_rank ) Then
                ! Colour for comm spliting
                my_colour = ks
                ! As in this strategy we only have 1 k point store which one it is
                my_ks( 1 ) = ks
                Exit
             End If
          End Do
       End If

    Else

       ! 2) Not enough procs to make distibuting matrices wortwhile. Just use a simple
       ! round robin assignment and do operations in serial

       ! First work out how many ks points I will hold
       n_my_ks = n_ks / n_procs_parent
       n_left = n_ks - n_my_ks * n_procs_parent
       If( me_parent < n_left ) Then
          n_my_ks = n_my_ks + 1
       End If
       Allocate( my_ks( 1:n_my_ks ) )

       ! Now assign them
       this_ks = 0
       Do ks = me_parent + 1, n_ks, n_procs_parent
          this_ks = this_ks + 1
          my_ks( this_ks ) = ks
       End Do

       ! And colour my communicators, an mpi_comm_split is next
       my_colour = Merge( me_parent, MPI_UNDEFINED, me_parent < n_ks )

    End If k_split_strategy

    ! Can now split the communicator
    Call MPI_Comm_split( split_A%parent_communicator, my_colour, 0, k_comm, error )

    ! Now start setting up the k points held by this set of processes (if any!)
    Allocate( split_A%my_k_points( 1:n_my_ks ) )
    Allocate( this_k_indices( 1:Size( A%all_k_point_info( 1 )%k_indices ) ) )
    Do ks = 1, n_my_ks
       split_A%my_k_points( ks )%info = split_A%all_k_point_info( my_ks( ks ) )
       ! Irreps not split yet hence no split at this level
       Allocate( split_A%my_k_points( ks )%data( 1:1 ) )
       split_A%my_k_points( ks )%data( 1 )%label = 1
       this_k_type    = split_A%all_k_point_info( my_ks( ks ) )%k_type
       this_s         = split_A%all_k_point_info( my_ks( ks ) )%spin
       this_k_indices = split_A%all_k_point_info( my_ks( ks ) )%k_indices
       ! Now need to generate a source matrix from the communicator 
       Call ks_matrix_comm_to_base( k_comm, base_matrix )
       ! Need to get sizes for creation
       m = A%my_k_points( my_ks( ks ) )%data( 1 )%matrix%size( 1 )
       n = A%my_k_points( my_ks( ks ) )%data( 1 )%matrix%size( 2 )
       Call split_A%my_k_points( ks )%data( 1 )%matrix%create( this_k_type == K_POINT_COMPLEX, &
            m, n, base_matrix )
       split_A%my_k_points( ks )%communicator = base_matrix%get_comm()
    End Do

    ! If required redistribute the data from the all procs distribution into the split distribution
    ! Note all procs in parent comm must be involved as all hold data for the unsplit matrx
    If( loc_redist ) Then
       Allocate( this_ks_matrix )
       Do ks = 1, n_ks
          this_ks_matrix = A%my_k_points( ks )%data( 1 )%matrix
          ! Find out if I hold this in the split distribution
          ! Note indicate to the remap routine that we hod no data by having an unallocated array
          ! c.f. Null pointer in C type languages
          this_ks = split_A%get_my_ks_index( ks )
          If( this_ks /= NOT_ME ) Then
             Allocate( split_ks_matrix )
             split_ks_matrix = split_A%my_k_points( this_ks )%data( 1 )%matrix
          End If
          Call ks_matrix_remap_data( this_ks_matrix, A%parent_communicator, split_ks_matrix )
          If( this_ks /= NOT_ME ) Then
             split_A%my_k_points( this_ks )%data( 1 )%matrix = split_ks_matrix
             Deallocate( split_ks_matrix ) ! Important as an deallocated matrix indicates no data on this process in the remap routine
          End If
       End Do
       Deallocate( this_ks_matrix )
    End If

    ! NEED TO THINK ABOUT PRESERVING DAGGERS

    ! Make sure the iterator in the new ks_array is reset
    Call split_A%iterator_reset()

  End Subroutine ks_array_split_ks

  Subroutine ks_array_join_ks( A, joined_A, redistribute )

    !! Re-join a previously split ks_array A so that each matrix in it is distributed across all the processes

    Use mpi             , Only : MPI_Comm_rank, MPI_Allreduce, MPI_IN_PLACE, MPI_INTEGER, MPI_SUM
    Use ks_matrix_module, Only : ks_matrix, ks_matrix_comm_to_base, ks_matrix_remap_data
    
    Class( ks_array     ), Intent( In    ) :: A                 !! The new array
    Type ( ks_array     ), Intent(   Out ) :: joined_A          !! The new unsplit array
    Logical, Optional    , Intent( In    ) :: redistribute      !! By default the data associated with the matrix is redistributed

    Type( ks_matrix ), Allocatable :: this_ks_matrix
    Type( ks_matrix ), Allocatable :: split_ks_matrix

    Type( ks_matrix ) :: base
    
    Integer, Dimension( :, : ), Allocatable :: matrix_sizes

    Integer :: m, n
    Integer :: me_split
    Integer :: n_ks, this_ks, my_ks
    Integer :: ks
    Integer :: i
    Integer :: error
    
    Logical :: loc_redist
    Logical :: is_split

    ! If A is NOT already split across ks points no need to do very much except copy ...
    is_split = Size( A%all_k_point_info ) /= Size( A%my_k_points )
    If( .Not. is_split ) Then
       joined_A = A
       ! .. and make sure the iterator in the new ks_array is reset
       Call joined_A%iterator_reset()
       Return
    End If

    If( Present( redistribute ) ) Then
       loc_redist = redistribute
    Else
       loc_redist =.True.
    End If

    ! Have to get the sizes of the matrices. This is only known locally on the processes that
    ! actually hold data for the ks point. Thus we need to replicate those values.
    Allocate( matrix_sizes( 1:2, 1:Size( A%all_k_point_info ) ) )
    matrix_sizes = 0
    Do i = 1, Size( matrix_sizes, Dim = 2 )
       m = A%size( A%all_k_point_info( i )%k_indices, A%all_k_point_info( i )%spin, 1 )
       n = A%size( A%all_k_point_info( i )%k_indices, A%all_k_point_info( i )%spin, 2 )
       If( m /= NOT_ME ) Then
          ! This process holds data on this k point. If I am rank zero in this communicator fill in
          ! the size in the buffer
          my_ks = A%get_my_ks_index( i )
          Call MPI_Comm_rank( A%my_k_points( my_ks )%data( 1 )%matrix%get_comm(), me_split, error )
          If( me_split == 0 ) Then
             matrix_sizes( 1, i ) = m
             matrix_sizes( 2, i ) = n
          End If
       End If
    End Do
    Call MPI_Allreduce( MPI_IN_PLACE, matrix_sizes, Size( matrix_sizes ), MPI_Integer, MPI_SUM, A%parent_communicator, error )

    ! OK, now know the size of all of the matrices, can create a base object to act as a template
    Call ks_matrix_comm_to_base( A%parent_communicator, base )

    ! Now create the matrices in the joined array
    ! Once we can create arrays with matrices of different sizes this should be revisited as it can be simplified

    ! Set up stuff relevant to all ks points
    joined_A%all_k_point_info    = A%all_k_point_info
    joined_A%parent_communicator = A%parent_communicator

    ! As joined_A contains all the ks points setting up the list of my ks points is easy
    Allocate( joined_A%my_k_points( 1:Size( joined_A%all_k_point_info ) ) )
    joined_A%my_k_points( : )%info = joined_A%all_k_point_info

    n_ks = Size( joined_A%all_k_point_info )
    Do ks = 1, n_ks
       m = matrix_sizes( 1, ks )
       n = matrix_sizes( 2, ks )
       Allocate( joined_A%my_k_points( ks )%data( 1:1 ) )
       Call joined_A%my_k_points( ks )%data( 1 )%matrix%create( &
            joined_A%my_k_points( ks )%info%k_type == K_POINT_COMPLEX, &
            m, n, base )
       joined_A%my_k_points( ks )%communicator = joined_A%parent_communicator
    End Do

    If( loc_redist ) Then
       ! If required redistribute the data from the split matrix back to the unsplit one
       Do ks = 1, n_ks
          Allocate( this_ks_matrix )
          m = matrix_sizes( 1, ks )
          n = matrix_sizes( 2, ks )
          ! Give this_ks_matrix a datatype (i.e. real or complex )
          Call this_ks_matrix%create( joined_A%my_k_points( ks )%info%k_type == K_POINT_COMPLEX, m, n, &
               joined_A%my_k_points( ks )%data( 1 )%matrix )
          this_ks = A%get_my_ks_index( ks )
          If( this_ks /= NOT_ME ) Then
             Allocate( split_ks_matrix )
             split_ks_matrix = A%my_k_points( this_ks )%data( 1 )%matrix
          End If
          Call ks_matrix_remap_data( split_ks_matrix, A%parent_communicator, this_ks_matrix )
          If( this_ks /= NOT_ME ) Then
             Deallocate( split_ks_matrix ) ! Important as an deallocated matrix indicates no data on this process in the remap routine
          End If
          joined_A%my_k_points( ks )%data( 1 )%matrix = this_ks_matrix
          Deallocate( this_ks_matrix )
       End Do
    End If

    ! NEED TO THINK ABOUT PRESERVING DAGGERS
    
    ! Make sure the iterator in the new ks_array is reset
    Call joined_A%iterator_reset()

  End Subroutine ks_array_join_ks

  Function ks_array_dagger( A ) Result( tA )

    !! Form the Hermitian conjugate of the matrices (Tranpose for real data)

    Type( ks_array ) :: tA

    Class( ks_array ), Intent( In ) :: A

    Integer :: my_ks, my_irrep

    Call tA%create( NO_DATA, NO_DATA, A )
    
    Do my_ks = 1, Size( A%my_k_points )
       ! Irreps will need more thought - worrk currenly as burnt into as 1
       Do my_irrep = 1, Size( A%my_k_points( my_ks )%data )
          Associate( Aks  =>  A%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     tAks => tA%my_k_points( my_ks )%data( my_irrep )%matrix )
            tAks = .Dagger. Aks
          End Associate
       End Do
    End Do

  End Function ks_array_dagger

  Function ks_array_mult( A, B ) Result( C )

    !! Multiply the arays together element by element (i.e. matrix by matrix ) 

    Type( ks_array ) :: C

    Class( ks_array ), Intent( In ) :: A
    Type ( ks_array ), Intent( In ) :: B

    Integer :: my_ks, my_irrep

    Call C%create( NO_DATA, NO_DATA, A )
    
    Do my_ks = 1, Size( A%my_k_points )
       ! Irreps will need more thought - work currenly as burnt into as 1
       Do my_irrep = 1, Size( A%my_k_points( my_ks )%data )
          Associate( Aks => A%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     Bks => B%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     Cks => C%my_k_points( my_ks )%data( my_irrep )%matrix )
            Cks = Aks * Bks
          End Associate
       End Do
    End Do

  End Function ks_array_mult

  Function ks_array_rscal_mult( s, A ) Result( C )

    !! Multiply the arays together element by element (i.e. matrix by matrix ) 

    Type( ks_array ) :: C

    Real( wp )       , Intent( In ) :: s
    Class( ks_array ), Intent( In ) :: A

    Integer :: my_ks, my_irrep

    Call C%create( NO_DATA, NO_DATA, A )
    
    Do my_ks = 1, Size( A%my_k_points )
       ! Irreps will need more thought - work currenly as burnt into as 1
       Do my_irrep = 1, Size( A%my_k_points( my_ks )%data )
          Associate( Aks => A%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     Cks => C%my_k_points( my_ks )%data( my_irrep )%matrix )
            Cks = s * Aks
          End Associate
       End Do
    End Do

  End Function ks_array_rscal_mult

  Function ks_array_mult_rscal( A, s ) Result( C )

    !! Multiply the arays together element by element (i.e. matrix by matrix ) 

    Type( ks_array ) :: C

    Class( ks_array ), Intent( In ) :: A
    Real( wp )       , Intent( In ) :: s

    Integer :: my_ks, my_irrep

    Call C%create( NO_DATA, NO_DATA, A )
    
    Do my_ks = 1, Size( A%my_k_points )
       ! Irreps will need more thought - work currenly as burnt into as 1
       Do my_irrep = 1, Size( A%my_k_points( my_ks )%data )
          Associate( Aks => A%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     Cks => C%my_k_points( my_ks )%data( my_irrep )%matrix )
            Cks = s * Aks
          End Associate
       End Do
    End Do

  End Function ks_array_mult_rscal

  Function ks_array_add( A, B ) Result( C )

    !! Add the arays together element by element (i.e. matrix by matrix ) 

    Type( ks_array ) :: C

    Class( ks_array ), Intent( In ) :: A
    Type ( ks_array ), Intent( In ) :: B

    Integer :: my_ks, my_irrep

    Call C%create( NO_DATA, NO_DATA, A )
    
    Do my_ks = 1, Size( A%my_k_points )
       ! Irreps will need more thought - work currenly as burnt into as 1
       Do my_irrep = 1, Size( A%my_k_points( my_ks )%data )
          Associate( Aks => A%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     Bks => B%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     Cks => C%my_k_points( my_ks )%data( my_irrep )%matrix )
            Cks = Aks + Bks
          End Associate
       End Do
    End Do

  End Function ks_array_add

  Function ks_array_add_diagonal( A, d ) Result( C )

    !! Add a general matrix to a diagonal one

    Type( ks_array ) :: C

    Class( ks_array )         , Intent( In ) :: A
    Real( wp ), Dimension( : ), Intent( In ) :: d

    Integer :: my_ks, my_irrep

    Call C%create( NO_DATA, NO_DATA, A )
    
    Do my_ks = 1, Size( A%my_k_points )
       ! Irreps will need more thought - work currenly as burnt into as 1
       Do my_irrep = 1, Size( A%my_k_points( my_ks )%data )
          Associate( Aks => A%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     Cks => C%my_k_points( my_ks )%data( my_irrep )%matrix )
            Cks = Aks + d
          End Associate
       End Do
    End Do

  End Function ks_array_add_diagonal

  Function ks_array_diagonal_add( d, A ) Result( C )

    !! Add a general matrix to a diagonal one

    Type( ks_array ) :: C

    Real( wp ), Dimension( : ), Intent( In ) :: d
    Class( ks_array )         , Intent( In ) :: A
    
    Integer :: my_ks, my_irrep

    Call C%create( NO_DATA, NO_DATA, A )
    
    Do my_ks = 1, Size( A%my_k_points )
       ! Irreps will need more thought - work currenly as burnt into as 1
       Do my_irrep = 1, Size( A%my_k_points( my_ks )%data )
          Associate( Aks => A%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     Cks => C%my_k_points( my_ks )%data( my_irrep )%matrix )
            Cks = d + Aks
          End Associate
       End Do
    End Do

  End Function ks_array_diagonal_add

  Function ks_array_subtract( A, B ) Result( C )

    !! Subtract the arays together element by element (i.e. matrix by matrix ) 

    Type( ks_array ) :: C

    Class( ks_array ), Intent( In ) :: A
    Type ( ks_array ), Intent( In ) :: B

    Integer :: my_ks, my_irrep

    Call C%create( NO_DATA, NO_DATA, A )
    
    Do my_ks = 1, Size( A%my_k_points )
       ! Irreps will need more thought - work currenly as burnt into as 1
       Do my_irrep = 1, Size( A%my_k_points( my_ks )%data )
          Associate( Aks => A%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     Bks => B%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     Cks => C%my_k_points( my_ks )%data( my_irrep )%matrix )
            Cks = Aks - Bks
          End Associate
       End Do
    End Do

  End Function ks_array_subtract

  Function ks_array_subtract_diagonal( A, d ) Result( C )

    !! Subtract a diagonal matrix from a general matrix 

    Type( ks_array ) :: C

    Class( ks_array )         , Intent( In ) :: A
    Real( wp ), Dimension( : ), Intent( In ) :: d

    Integer :: my_ks, my_irrep

    Call C%create( NO_DATA, NO_DATA, A )
    
    Do my_ks = 1, Size( A%my_k_points )
       ! Irreps will need more thought - work currenly as burnt into as 1
       Do my_irrep = 1, Size( A%my_k_points( my_ks )%data )
          Associate( Aks => A%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     Cks => C%my_k_points( my_ks )%data( my_irrep )%matrix )
            Cks = Aks - d
          End Associate
       End Do
    End Do

  End Function ks_array_subtract_diagonal

  Function ks_array_diagonal_subtract( d, A ) Result( C )

    !! Subtract a general matrix from a diagonal one

    Type( ks_array ) :: C

    Real( wp ), Dimension( : ), Intent( In ) :: d
    Class( ks_array )         , Intent( In ) :: A
    
    Integer :: my_ks, my_irrep

    Call C%create( NO_DATA, NO_DATA, A )
    
    Do my_ks = 1, Size( A%my_k_points )
       ! Irreps will need more thought - work currenly as burnt into as 1
       Do my_irrep = 1, Size( A%my_k_points( my_ks )%data )
          Associate( Aks => A%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     Cks => C%my_k_points( my_ks )%data( my_irrep )%matrix )
            Cks = d - Aks
          End Associate
       End Do
    End Do

  End Function ks_array_diagonal_subtract

  Subroutine ks_array_diag( A, Q, E )

    !! Diagonalise each matrix

    Use mpi, Only : mpi_comm_rank, mpi_isend, mpi_recv, mpi_wait, mpi_bcast, &
         mpi_sizeof, mpi_type_match_size, MPI_INTEGER, MPI_ANY_SOURCE, MPI_STATUS_IGNORE, &
         MPI_TYPECLASS_REAL

    Class( ks_array        ),                 Intent( In    ) :: A
    Type ( ks_array        ),                 Intent(   Out ) :: Q
    Type ( ks_eval_storage ), Dimension( : ), Intent(   Out ) :: E

    Real( wp ) :: rdum

    Integer, Dimension( 1:2 ) :: buff_send, buff_recv
    
    Integer :: me, me_parent
    Integer :: ks_root, nb
    Integer :: error
    Integer :: request
    Integer :: rsize, handle
    Integer :: my_ks, ks
    Integer :: my_irrep

    Logical :: sending_data
    
    ! Make Q have the same set up as A
    Call Q%create( NO_DATA, NO_DATA, A )

    Do my_ks = 1, Size( A%my_k_points )
       ! Irreps will need more thought - worrk currenly as burnt into as 1
       Do my_irrep = 1, Size( A%my_k_points( my_ks )%data )
          ks = A%get_all_ks_index( my_ks )
          Associate( Aks => A%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     Qks => Q%my_k_points( my_ks )%data( my_irrep )%matrix )
            Call Aks%diag( Qks, E( ks )%evals )
          End Associate
       End Do
    End Do

    ! Replicate evals
    ! Again needs thought for ireps
    Call mpi_comm_rank( A%parent_communicator, me_parent, error )
    Do ks = 1, Size( A%all_k_point_info )
       E( ks )%ks_point = A%all_k_point_info( ks )
       my_ks = A%get_my_ks_index( ks )
       ! Work out who holds this set of evals and send how many there are
       ! the root node of the communicator holding them back to the root node of the parent communicator
       sending_data = .False.
       If( my_ks /= NOT_ME ) Then
          Call mpi_comm_rank( A%my_k_points( my_ks )%communicator, me, error )
          sending_data = me == 0
          If( sending_data ) then
             buff_send( 1 ) = me_parent
             buff_send( 2 ) = Size( E( ks )%evals )
             Call mpi_isend( buff_send, 2, MPI_INTEGER, 0, ks, A%parent_communicator, request, error )
          End If
       End If
       If( me_parent == 0 ) Then
          Call mpi_recv( buff_recv, 2, MPI_INTEGER, MPI_ANY_SOURCE, ks, A%parent_communicator, MPI_STATUS_IGNORE, error )
       End If
       If( sending_data ) Then
          Call mpi_wait( request, MPI_STATUS_IGNORE, error )
       End If
       ! Now on root of parent bcast back to all
       Call mpi_bcast( buff_recv, 2, MPI_INTEGER, 0, A%parent_communicator, error )
       ks_root = buff_recv( 1 )
       nb      = buff_recv( 2 )
       ! Now know how many evals we will recv - allocate memory if haven't done so already because I don't 'own' this k point
       If( .Not. Allocated( E( ks )%evals ) ) Then
          Allocate( E( ks )%evals( 1:nb ) )
       End If
       ! And finally bcast out the values from the root node for this set of evals
       Call mpi_sizeof( rdum, rsize, error )
       Call mpi_type_match_size( MPI_TYPECLASS_REAL, rsize, handle, error )
       Call mpi_bcast( E( ks )%evals, Size( E( ks )%evals ), handle, ks_root, A%parent_communicator, error )
    End Do
    
  End Subroutine ks_array_diag

  ! Choleski decompose each matrix
  
  Function ks_array_choleski( A ) Result( C )

    !! Choleski decompose the matrices

    Type( ks_array ) :: C

    Class( ks_array )         , Intent( In ) :: A
    
    Integer :: my_ks, my_irrep

    Call C%create( NO_DATA, NO_DATA, A )
    
    Do my_ks = 1, Size( A%my_k_points )
       ! Irreps will need more thought - work currenly as burnt into as 1
       Do my_irrep = 1, Size( A%my_k_points( my_ks )%data )
          Associate( Aks => A%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     Cks => C%my_k_points( my_ks )%data( my_irrep )%matrix )
            Cks = .Choleski. Aks
          End Associate
       End Do
    End Do

  End Function ks_array_choleski
  
  Function ks_array_tr_inv( A ) Result( C )

    !! Invert triangular matrix

    Type( ks_array ) :: C

    Class( ks_array )         , Intent( In ) :: A
    
    Integer :: my_ks, my_irrep

    Call C%create( NO_DATA, NO_DATA, A )
    
    Do my_ks = 1, Size( A%my_k_points )
       ! Irreps will need more thought - work currenly as burnt into as 1
       Do my_irrep = 1, Size( A%my_k_points( my_ks )%data )
          Associate( Aks => A%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     Cks => C%my_k_points( my_ks )%data( my_irrep )%matrix )
            Cks = .TrInv. Aks
          End Associate
       End Do
    End Do

  End Function ks_array_tr_inv
  
  Function ks_array_plus( A ) Result( C )

    !! Unary plus operation

    Type( ks_array ) :: C

    Class( ks_array )         , Intent( In ) :: A
    
    Integer :: my_ks, my_irrep

    Call C%create( NO_DATA, NO_DATA, A )
    
    Do my_ks = 1, Size( A%my_k_points )
       ! Irreps will need more thought - work currenly as burnt into as 1
       Do my_irrep = 1, Size( A%my_k_points( my_ks )%data )
          Associate( Aks => A%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     Cks => C%my_k_points( my_ks )%data( my_irrep )%matrix )
            Cks = + Aks
          End Associate
       End Do
    End Do

  End Function ks_array_plus
  
  Function ks_array_minus( A ) Result( C )

    !! Unary minus operation

    Type( ks_array ) :: C

    Class( ks_array )         , Intent( In ) :: A
    
    Integer :: my_ks, my_irrep

    Call C%create( NO_DATA, NO_DATA, A )
    
    Do my_ks = 1, Size( A%my_k_points )
       ! Irreps will need more thought - work currenly as burnt into as 1
       Do my_irrep = 1, Size( A%my_k_points( my_ks )%data )
          Associate( Aks => A%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     Cks => C%my_k_points( my_ks )%data( my_irrep )%matrix )
            Cks = - Aks
          End Associate
       End Do
    End Do

  End Function ks_array_minus
  
  Function ks_array_extract( A, r1, r2, c1, c2 ) Result( C )

    !! Extract a patch from each of the matrices and return a new matrix
    
    Type( ks_array ) :: C

    Class( ks_array ), Intent( In ) :: A
    ! Do we want the indices on the patches to be arrays so each ks point
    ! can extract a different patch??
    Integer          , Intent( In ) :: r1 
    Integer          , Intent( In ) :: r2
    Integer          , Intent( In ) :: c1 
    Integer          , Intent( In ) :: c2
    
    Integer :: my_ks, my_irrep

    Call C%create( NO_DATA, NO_DATA, A )
    
    Do my_ks = 1, Size( A%my_k_points )
       ! Irreps will need more thought - work currenly as burnt into as 1
       Do my_irrep = 1, Size( A%my_k_points( my_ks )%data )
          Associate( Aks => A%my_k_points( my_ks )%data( my_irrep )%matrix, &
                     Cks => C%my_k_points( my_ks )%data( my_irrep )%matrix )
            Cks = Aks%extract( r1, r2, c1, c2 )
          End Associate
       End Do
    End Do

  End Function ks_array_extract

  Pure Function get_all_ks_index( A, my_ks ) Result( ks )

    !! Given a local index into my list of ks points return the index in the full list
    
    Integer :: ks

    Class( ks_array ), Intent( In ) :: A
    Integer          , Intent( In ) :: my_ks

    Do ks = 1, Size( A%all_k_point_info )
       If( A%all_k_point_info( ks )%spin == A%my_k_points( my_ks )%info%spin ) Then
          If( All( A%all_k_point_info( ks )%k_indices == A%my_k_points( my_ks )%info%k_indices ) ) Then
             Exit
          End If
       End If
    End Do
    
  End Function get_all_ks_index

  Pure Function get_my_ks_index( A, ks ) Result( my_ks )

    !! Given an index into the full list of ks points return an index into my local list

    Integer :: my_ks

    Class( ks_array ), Intent( In ) :: A
    Integer          , Intent( In ) :: ks

    Integer :: k
    
    my_ks = NOT_ME

    Do k = 1, Size( A%my_k_points )
       If( A%all_k_point_info( ks )%spin == A%my_k_points( k )%info%spin ) Then
          If( All( A%all_k_point_info( ks )%k_indices == A%my_k_points( k )%info%k_indices ) ) Then
             my_ks = k
             Exit
          End If
       End If
    End Do

  End Function get_my_ks_index

  Pure Function get_ks( A, k, s ) Result( ks )

    !! Given the label for a k point and spin return the index in the full ks point list

    Integer :: ks

    Class( ks_array )        , Intent( In ) :: A
    Integer, Dimension( 1:3 ), Intent( In ) :: k
    Integer                  , Intent( In ) :: s

    Do ks = 1, Size( A%all_k_point_info )
       If( A%all_k_point_info( ks )%spin == s ) Then
          If( All( A%all_k_point_info( ks )%k_indices == k ) ) Then
             Exit
          End If
       End If
    End Do
    
  End Function get_ks

  Subroutine ks_array_set_global_real( A, k, s, m, n, p, q, data )

    !! For the matrix with spin label s and k-point label k set the (m:n,p:q) patch 
    ! Need to overload for irreps

    Class( ks_array )              , Intent( InOut ) :: A
    Integer    , Dimension( : )    , Intent( In    ) :: k
    Integer                        , Intent( In    ) :: s
    Integer                        , Intent( In    ) :: m
    Integer                        , Intent( In    ) :: n
    Integer                        , Intent( In    ) :: p
    Integer                        , Intent( In    ) :: q
    Real( wp ), Dimension( m:, p: ), Intent( In    ) :: data

    Integer :: ks, my_ks

    ks = A%get_ks( k, s )
    
    my_ks = A%get_my_ks_index( ks )

    If( my_ks /= NOT_ME ) Then
       Call A%my_k_points( my_ks )%data( 1 )%matrix%set_by_global( m, n, p, q, data )
    End If

  End Subroutine ks_array_set_global_real

  Subroutine ks_array_set_global_complex( A, k, s, m, n, p, q, data )

    !! For the matrix with spin label s and k-point label k set the (m:n,p:q) patch 
    ! Need to overload for irreps

    Class( ks_array )                 , Intent( InOut ) :: A
    Integer    , Dimension( : )       , Intent( In    ) :: k
    Integer                           , Intent( In    ) :: s
    Integer                           , Intent( In    ) :: m
    Integer                           , Intent( In    ) :: n
    Integer                           , Intent( In    ) :: p
    Integer                           , Intent( In    ) :: q
    Complex( wp ), Dimension( m:, p: ), Intent( In    ) :: data

    Integer :: ks, my_ks

    ks = A%get_ks( k, s )
    
    my_ks = A%get_my_ks_index( ks )

    If( my_ks /= NOT_ME ) Then
       Call A%my_k_points( my_ks )%data( 1 )%matrix%set_by_global( m, n, p, q, data )
    End If

  End Subroutine ks_array_set_global_complex  

  Subroutine ks_array_get_global_real( A, k, s, m, n, p, q, data )

    Use mpi, Only : MPI_Comm_rank, MPI_Isend, MPI_INTEGER, MPI_recv, MPI_ANY_SOURCE, &
         MPI_Wait, MPI_STATUS_IGNORE, MPI_Bcast, MPI_Sizeof, MPI_Type_match_size, &
         MPI_TYPECLASS_REAL

    !! For the matrix with spin label s and k-point label k get the (m:n,p:q) patch 
    ! Need to overload for irreps

    Class( ks_array )              , Intent( In    ) :: A
    Integer    , Dimension( : )    , Intent( In    ) :: k
    Integer                        , Intent( In    ) :: s
    Integer                        , Intent( In    ) :: m
    Integer                        , Intent( In    ) :: n
    Integer                        , Intent( In    ) :: p
    Integer                        , Intent( In    ) :: q
    Real( wp ), Dimension( m:, p: ), Intent(   Out ) :: data

    Real( wp ) :: rdum

    Integer :: me_parent, me
    Integer :: buff_send, buff_recv
    Integer :: ks, my_ks
    Integer :: request
    Integer :: ks_root
    Integer :: rsize, handle
    Integer :: error

    Logical :: sending_data
    
    ks = A%get_ks( k, s )
    
    my_ks = A%get_my_ks_index( ks )

    If( my_ks /= NOT_ME ) Then
       Call A%my_k_points( my_ks )%data( 1 )%matrix%get_by_global( m, n, p, q, data )
    End If

    ! Need to replicate data over parent communicator if split ks points
    
    ! Work out who holds this set of evals and send how many there are
    ! the root node of the communicator holding them back to the root node of the parent communicator
    Call mpi_comm_rank( A%parent_communicator, me_parent, error )
    sending_data = .False.
    If( my_ks /= NOT_ME ) Then
       Call mpi_comm_rank( A%my_k_points( my_ks )%communicator, me, error )
       sending_data = me == 0
       If( sending_data ) then
          buff_send = me_parent
          Call mpi_isend( buff_send, 1, MPI_INTEGER, 0, ks, A%parent_communicator, request, error )
       End If
    End If
    If( me_parent == 0 ) Then
       Call mpi_recv( buff_recv, 1, MPI_INTEGER, MPI_ANY_SOURCE, ks, A%parent_communicator, MPI_STATUS_IGNORE, error )
    End If
    If( sending_data ) Then
       Call mpi_wait( request, MPI_STATUS_IGNORE, error )
    End If
    ! Now on root of parent bcast back to all
    Call mpi_bcast( buff_recv, 1, MPI_INTEGER, 0, A%parent_communicator, error )
    ks_root = buff_recv 
    Call mpi_sizeof( rdum, rsize, error )
    Call mpi_type_match_size( MPI_TYPECLASS_REAL, rsize, handle, error )
    Call mpi_bcast( data, Size( data ), handle, ks_root, A%parent_communicator, error )    

  End Subroutine ks_array_get_global_real

  Subroutine ks_array_get_global_complex( A, k, s, m, n, p, q, data )

    !! For the matrix with spin label s and k-point label k set the (m:n,p:q) patch 

    Use mpi, Only : MPI_Comm_rank, MPI_Isend, MPI_INTEGER, MPI_recv, MPI_ANY_SOURCE, &
         MPI_Wait, MPI_STATUS_IGNORE, MPI_Bcast, MPI_Sizeof, MPI_Type_match_size, &
         MPI_TYPECLASS_COMPLEX

    ! Need to overload for irreps

    Class( ks_array )                 , Intent( In    ) :: A
    Integer    , Dimension( : )       , Intent( In    ) :: k
    Integer                           , Intent( In    ) :: s
    Integer                           , Intent( In    ) :: m
    Integer                           , Intent( In    ) :: n
    Integer                           , Intent( In    ) :: p
    Integer                           , Intent( In    ) :: q
    Complex( wp ), Dimension( m:, p: ), Intent(   Out ) :: data

    Complex( wp ) :: cdum

    Integer :: me_parent, me
    Integer :: buff_send, buff_recv
    Integer :: ks, my_ks
    Integer :: request
    Integer :: ks_root
    Integer :: rsize, handle
    Integer :: error

    Logical :: sending_data
    
    ks = A%get_ks( k, s )
    
    my_ks = A%get_my_ks_index( ks )

    If( my_ks /= NOT_ME ) Then
       Call A%my_k_points( my_ks )%data( 1 )%matrix%get_by_global( m, n, p, q, data )
    End If

    ! Need to replicate data over parent communicator if split ks points
    
    ! Work out who holds this set of evals and send how many there are
    ! the root node of the communicator holding them back to the root node of the parent communicator
    Call mpi_comm_rank( A%parent_communicator, me_parent, error )
    sending_data = .False.
    If( my_ks /= NOT_ME ) Then
       Call mpi_comm_rank( A%my_k_points( my_ks )%communicator, me, error )
       sending_data = me == 0
       If( sending_data ) then
          buff_send = me_parent
          Call mpi_isend( buff_send, 1, MPI_INTEGER, 0, ks, A%parent_communicator, request, error )
       End If
    End If
    If( me_parent == 0 ) Then
       Call mpi_recv( buff_recv, 1, MPI_INTEGER, MPI_ANY_SOURCE, ks, A%parent_communicator, MPI_STATUS_IGNORE, error )
    End If
    If( sending_data ) Then
       Call mpi_wait( request, MPI_STATUS_IGNORE, error )
    End If
    ! Now on root of parent bcast back to all
    Call mpi_bcast( buff_recv, 1, MPI_INTEGER, 0, A%parent_communicator, error )
    ks_root = buff_recv 
    Call mpi_sizeof( cdum, rsize, error )
    Call mpi_type_match_size( MPI_TYPECLASS_COMPLEX, rsize, handle, error )
    Call mpi_bcast( data, Size( data ), handle, ks_root, A%parent_communicator, error )    

  End Subroutine ks_array_get_global_complex

  Function ks_array_g_to_l( A, k, s, what ) Result( gl_indexing )

    !! For the matrix with spin label S and k-point label K get the global to local indexing array
    !! If WHAT == 'R' get the array for rows. If WHAT == 'C' get it for columns

    Integer, Dimension( : ), Allocatable :: gl_indexing

    Class( ks_array )          , Intent( In ) :: A
    Integer                    , Intent( In ) :: s
    Integer    , Dimension( : ), Intent( In ) :: k
    Character( Len = * )       , Intent( In ) :: what

    Integer :: ks, my_ks

    ks = A%get_ks( k, s )
    
    my_ks = A%get_my_ks_index( ks )

    If( my_ks /= NOT_ME ) Then
       gl_indexing = A%my_k_points( my_ks )%data( 1 )%matrix%global_to_local( what )
       Where( gl_indexing <= 0 )
          gl_indexing = NOT_ME
       End Where
    End If

  End Function ks_array_g_to_l

  Function ks_array_l_to_g( A, k, s, what ) Result( lg_indexing )

    !! For the matrix with spin label S and k-point label K get the local to global indexing array
    !! If WHAT == 'R' get the array for rows. If WHAT == 'C' get it for columns

    Integer, Dimension( : ), Allocatable :: lg_indexing

    Class( ks_array )          , Intent( In ) :: A
    Integer                    , Intent( In ) :: s
    Integer    , Dimension( : ), Intent( In ) :: k
    Character( Len = * )       , Intent( In ) :: what

    Integer :: ks, my_ks

    ks = A%get_ks( k, s )
    
    my_ks = A%get_my_ks_index( ks )

    If( my_ks /= NOT_ME ) Then
       lg_indexing = A%my_k_points( my_ks )%data( 1 )%matrix%local_to_global( what )
    End If

  End Function ks_array_l_to_g

  Function ks_array_size( A, k, s, dim ) Result( n )

    Integer :: n

    Class( ks_array )          , Intent( In )           :: A
    Integer                    , Intent( In )           :: s
    Integer    , Dimension( : ), Intent( In )           :: k
    Integer                    , Intent( In ), Optional :: dim

    Integer :: ks, my_ks

    ks = A%get_ks( k, s )
    
    my_ks = A%get_my_ks_index( ks )

    If( my_ks /= NOT_ME ) Then
       n = A%my_k_points( my_ks )%data( 1 )%matrix%size( dim )
    Else
       n = NOT_ME
    End If

  End Function ks_array_size

  Subroutine ks_array_iterator_init( A )

    !! Initialise an iterator

    Class( ks_array ), Intent( InOut ) :: A

    A%iterator_value = 0

  End Subroutine ks_array_iterator_init

  Subroutine ks_array_iterator_reset( A )

    !! Reset an iterator

    Class( ks_array ), Intent( InOut ) :: A

    A%iterator_value = INVALID

  End Subroutine ks_array_iterator_reset

  Function ks_array_iterator_next( A ) Result( ks )

    !! Move to the next matrix in the ks_array
    !! Return an unallocated info value if there is no next value, and set the iterator
    !! to ke beyond the end of the array

    Type( ks_point_info ), Allocatable :: ks
    
    Class( ks_array ), Intent( InOut ) :: A

    Allocate( ks )
    
    A%iterator_value = A%iterator_value + 1

    If( A%iterator_value <= Ubound( A%my_k_points, Dim = 1 ) ) Then
       ks = A%my_k_points( A%iterator_value )%info
    Else
       A%iterator_value = Ubound( A%my_k_points, Dim = 1 ) + 1
       ks%k_type = K_POINT_NOT_EXIST
    End If

  End Function ks_array_iterator_next

  Function ks_array_iterator_previous( A ) Result( ks )

    !! Move to the previous matrix in the ks_array
    !! Return an unallocated info value if there is no previous value, and set the iterator
    !! to ke beyond the end of the array

    Type( ks_point_info ), Allocatable :: ks
    
    Class( ks_array ), Intent( InOut ) :: A

    Allocate( ks )
    
    A%iterator_value = A%iterator_value - 1

    If( A%iterator_value >= Lbound( A%my_k_points, Dim = 1 ) ) Then
       ks = A%my_k_points( A%iterator_value )%info
    Else
       A%iterator_value = Lbound( A%my_k_points, Dim = 1 ) - 1
       ks%k_type = K_POINT_NOT_EXIST
    End If

  End Function ks_array_iterator_previous

End Module ks_array_module
