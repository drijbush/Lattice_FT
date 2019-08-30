Module matrix_mapping_module

  !! A module to describe how a matrix maps onto a set of processes
  !! This is the only level at which BLACS is accessed

  Use numbers_module     , Only : wp
  Use mpi                , Only : mpi_comm_size, mpi_comm_rank
  Use proc_mapping_module, Only : proc_mapping, proc_mapping_init, proc_mapping_comm_to_base, proc_mapping_finalise
  
  Implicit None

  Type, Public, Extends( proc_mapping )  :: matrix_mapping
     !! An instance of a matrix mapping
     Integer, Dimension( 1:9 ), Private :: descriptor                 !! The BLACS descriptor or the matrix
   Contains
     ! Public methods
     Procedure, Public  :: print          => print_matrix_mapping     !! Print information about a matrix mapping
     Procedure, Public  :: get_data       => get_matrix_mapping_data  !! Get information about the matrix mapping
     Procedure, Public  :: get_descriptor => get_matrix_descriptor    !! Returns the whole BLACS descriptor
     Generic  , Public  :: set            => set_matrix_mapping       !! Set a matrix mapping
     Procedure, Public  :: free           => free_matrix_mapping      !! Frees any data associated with a matrix mapping
     ! Private Implementations
     Procedure, Private :: set_matrix_mapping                         
  End type matrix_mapping

  Integer, Parameter, Private :: INVALID = -1

  Public :: matrix_mapping_init
  Public :: matrix_mapping_comm_to_base
  Public :: matrix_mapping_finalise
  
  Private

  ! Parameters to help internal use of BLACS
  ! What the elements of the descriptor mean
  Integer, Parameter, Private :: dtype_a = 1 ! descriptor type
  Integer, Parameter, Private :: ctxt_a  = 2 ! context
  Integer, Parameter, Private :: m_a     = 3 ! number of global rows
  Integer, Parameter, Private :: n_a     = 4 ! number of global columns
  Integer, Parameter, Private :: mb_a    = 5 ! row blocking factor
  Integer, Parameter, Private :: nb_a    = 6 ! column blocking factor
  Integer, Parameter, Private :: rsrc_a  = 7 ! first process row which holds a
  Integer, Parameter, Private :: csrc_a  = 8 ! first process col which golds a
  Integer, Parameter, Private :: lld_a   = 9 ! leading dimension of LOCAL a

Contains

  Subroutine matrix_mapping_init
    !! Initialise the matrix mapping system
    Call proc_mapping_init
  End Subroutine matrix_mapping_init

  Subroutine matrix_mapping_comm_to_base( comm, mapping )

    !! Converts a MPI communicator into a base mapping
    !! i.e. one which contains the details about the processes
    !! but nothing about the size of the matrix
    
    Integer               , Intent( In    ) :: comm
    Type( matrix_mapping ), Intent(   Out ) :: mapping

    Type( proc_mapping ) :: proc_mapping_base

    Integer :: context

    Call proc_mapping_comm_to_base( comm, proc_mapping_base )

    Call get_context_from_comm( comm, context )

    Call mapping%set( proc_mapping_base, context, &
         INVALID, INVALID, &
         INVALID, INVALID, &
         INVALID, INVALID, &
         INVALID )
    
  End Subroutine matrix_mapping_comm_to_base

  Subroutine matrix_mapping_finalise

    !! Finalise the matrix mapping system

!    Use blacs_interfaces, Only : blacs_exit
    
    Call proc_mapping_finalise

    ! Kill all blacs contexts I use, but not the MPI subsytem
    ! ( as indicated by the argument being non-zero )
    ! Unfortunately in practice this means restarting the BLACS seems impossible grrrr....
!    Call blacs_exit( 1 )
    
  End Subroutine matrix_mapping_finalise

  Subroutine print_matrix_mapping( map )

    !! Print information about a matrix mapping

    Class( matrix_mapping ), Intent( In ) :: map

    Integer :: rank
    Integer :: error
    
    Call mpi_comm_rank( map%get_comm(), rank, error )
    If( rank == 0 ) Then
       Write( *, '( a, 9( i0, 1x ) )' ) 'Descriptor for matrix mapping ', map%descriptor
       Call map%proc_mapping%print
    End If
    
  End Subroutine print_matrix_mapping

  Subroutine get_context_from_comm( comm, context )

    !! Generate a BLACS context from an MPI communicator
    
    Use blacs_interfaces, Only : blacs_gridinit

    Integer, Intent( In    ) :: comm
    Integer, Intent(   Out ) :: context

    Integer :: nproc, nprow, npcol
    Integer :: error
    
    Call mpi_comm_size( comm, nproc, error )
    Call factor( nproc, nprow, npcol )
    context = comm
    Call blacs_gridinit( context, 'C', nprow, npcol )
    
  End Subroutine get_context_from_comm
  
  Subroutine set_matrix_mapping( map, proc_map, context, m, n, mb, nb, rsrc, csrc, lld )

    !! Set a matrix mapping

    Class( matrix_mapping ), Intent(   Out ) :: map
    Type ( proc_mapping   ), Intent( In    ) :: proc_map
    Integer                , Intent( In    ) :: context
    Integer                , Intent( In    ) :: m
    Integer                , Intent( In    ) :: n
    Integer                , Intent( In    ) :: mb
    Integer                , Intent( In    ) :: nb
    Integer                , Intent( In    ) :: rsrc
    Integer                , Intent( In    ) :: csrc
    Integer                , Intent( In    ) :: lld

    map%proc_mapping = proc_map

    map%descriptor = INVALID
    
    map%descriptor( dtype_a ) = 1 ! Dense matrix

    map%descriptor( ctxt_a ) = context

    map%descriptor( m_a     ) = m ! Global rows
    map%descriptor( n_a     ) = n ! Global cols
    
    map%descriptor( mb_a    ) = mb ! Row block fac
    map%descriptor( nb_a    ) = nb ! Col block fac

    map%descriptor( rsrc_a  ) = rsrc ! First proc row
    map%descriptor( csrc_a  ) = csrc ! First proc col

    map%descriptor( lld_a   ) = lld ! Local leading dimension
    
  End Subroutine set_matrix_mapping

  Subroutine get_matrix_mapping_data( map, comm, nprow, npcol, myprow, mypcol,&
       ctxt, m, n, mb, nb, rsrc, csrc, lld ) 
    
    !! Get information about the matrix mapping

    Use blacs_interfaces, Only : blacs_gridinfo

    Class( matrix_mapping ), Intent( In    )           :: map     !! The map in questions
    Integer                , Intent(   Out ), Optional :: comm    !! The MPI communicator
    Integer                , Intent(   Out ), Optional :: nprow   !! Number of process rows
    Integer                , Intent(   Out ), Optional :: npcol   !! Number of process columns
    Integer                , Intent(   Out ), Optional :: myprow  !! My process row
    Integer                , Intent(   Out ), Optional :: mypcol  !! my process column 
    Integer                , Intent(   Out ), Optional :: ctxt    !! The BLACS context
    Integer                , Intent(   Out ), Optional :: m       !! The number of rows in the matrix
    Integer                , Intent(   Out ), Optional :: n       !! The number of columns in the matrix
    Integer                , Intent(   Out ), Optional :: mb      !! The row blocking factor
    Integer                , Intent(   Out ), Optional :: nb      !! the column blocking factor
    Integer                , Intent(   Out ), Optional :: rsrc    
    Integer                , Intent(   Out ), Optional :: csrc
    Integer                , Intent(   Out ), Optional :: lld     !! The leading dimension of the local storage

    Integer :: loc_nprow, loc_npcol
    Integer :: loc_myprow, loc_mypcol

    Call blacs_gridinfo( map%descriptor( ctxt_a ), loc_nprow, loc_npcol, loc_myprow, loc_mypcol )

    If( Present( comm ) ) Then
       comm = map%get_comm()
    End If

    If( Present( nprow ) ) Then
       nprow = loc_nprow
    End If

    If( Present( npcol ) ) Then
       npcol = loc_npcol
    End If

    If( Present( myprow ) ) Then
       myprow = loc_myprow
    End If

    If( Present( mypcol ) ) Then
       mypcol = loc_mypcol
    End If

    If( Present( ctxt ) ) Then
       ctxt = map%descriptor( ctxt_a )
    End If

    If( Present( m ) ) Then
       m = map%descriptor( m_a )
    End If

    If( Present( n ) ) Then
       n = map%descriptor( n_a )
    End If

    If( Present( mb ) ) Then
       mb = map%descriptor( mb_a )
    End If

    If( Present( nb ) ) Then
       nb = map%descriptor( nb_a )
    End If

    If( Present( rsrc ) ) Then
       rsrc = map%descriptor( rsrc_a )       
    End If

    If( Present( csrc ) ) Then
       csrc = map%descriptor( csrc_a )       
    End If

    If( Present( lld ) ) Then
       lld = map%descriptor( lld_a )       
    End If

  End Subroutine get_matrix_mapping_data

  Subroutine free_matrix_mapping( map )

    !! Frees any data associated with a matrix mapping

    Use blacs_interfaces, Only : blacs_gridexit

    Class( matrix_mapping ), Intent( InOut ) :: map

    Integer :: ctxt

    Call map%get_data( ctxt = ctxt )

    Call blacs_gridexit( ctxt )
    
  End Subroutine free_matrix_mapping

  Pure Function get_matrix_descriptor( map ) Result( descriptor )

    !! Returns the descriptor used in the mapping

    Integer, Dimension( 1:9 ) :: descriptor

    Class( matrix_mapping ), Intent( In ) :: map

    descriptor = map%descriptor
    
  End Function get_matrix_descriptor
    
  Subroutine factor( a, b, c )

    !! Factors a such that a = b * c

    Integer, Intent( In    ) :: a
    Integer, Intent(   Out ) :: b
    Integer, Intent(   Out ) :: c

    b = Nint( Sqrt( Real( a ) ) )

    Do
       c = a / b
       If( b * c == a ) Then
          Exit
       End If
       b = b - 1 ! Loop will terminate when b = 1
    End Do

  End Subroutine factor

End Module matrix_mapping_module
 
