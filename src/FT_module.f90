Module Numbers

  Implicit None
  
  !  A module to look after the data types used by Crystal.
  
  Integer, Parameter :: float = Selected_real_kind( 13, 100 )
  Integer, Parameter :: imag  = Kind((0.0_float,0.0_float))
  Integer, Parameter :: int32 = Kind( 1 )
  Integer, Parameter :: int64 = Selected_int_kind( 15 )
  Integer, Parameter :: short = Selected_int_kind( 8 )
  Integer, Parameter :: logic = Kind( .True. )
  
  !  Floating Cardinals 

End Module Numbers

Module MPP_Fourier_transform

  Implicit None

  Public :: MPP_Lattice_FT
  
  Private
  
  Integer, Private, Parameter :: NOT_ME = -2

Contains

  Subroutine MPP_Lattice_FT( n_ao, n_shells, size_shells, operator_G_space, operator_K_space )

    Use numbers        , Only : wp => float
    Use ks_array_module, Only : ks_array, ks_point_info, K_POINT_NOT_EXIST

    Integer                   , Intent( In    ) :: n_ao
    Integer                   , Intent( In    ) :: n_shells
    Integer   , Dimension( : ), Intent( In    ) :: size_shells
    Real( wp ), Dimension( : ), Intent( In    ) :: operator_G_space
    Type( ks_array )          , Intent( InOut ) :: operator_K_space

    Type( ks_point_info ) :: this_ks

    Integer, Dimension( :, : ), Allocatable :: global_to_local_row
    Integer, Dimension( :, : ), Allocatable :: global_to_local_col

    Integer, Dimension( :, : ), Allocatable :: shell_start_row
    Integer, Dimension( :, : ), Allocatable :: shell_start_col

    Integer, Dimension( :, : ), Allocatable :: shell_finish_row
    Integer, Dimension( :, : ), Allocatable :: shell_finish_col
    
    Integer :: max_size_shell
    Integer :: n_ks_points
    Integer :: shell_start, shell_finish
    Integer :: ks
    Integer :: i_shell
    
    Logical, Dimension( :, : ), Allocatable :: i_own_row
    Logical, Dimension( :, : ), Allocatable :: i_own_col

    ! Some basic information about the system
    max_size_shell = Maxval( size_shells( 1:n_shells ) )
    
    ! Set up ks point dependent stuff
    ! First count how many ks points are stored on this process
    Call operator_K_space%iterator_init()
    n_ks_points = 0
    this_ks = operator_K_space%iterator_next()
    Do While( this_ks%k_type /= K_POINT_NOT_EXIST )
       n_ks_points = n_ks_points + 1
       this_ks = operator_K_space%iterator_next()
    End Do
    Call operator_K_space%iterator_reset()

    ! Now get the mapping arrays for global indices to local indices
    Allocate( global_to_local_row( 1:n_ao, 1:n_ks_points ) )
    Allocate( global_to_local_col( 1:n_ao, 1:n_ks_points ) )
    ! And populate the array
    Call operator_K_space%iterator_init()
    ks = 0
    this_ks = operator_K_space%iterator_next()
    Do While( this_ks%k_type /= K_POINT_NOT_EXIST )
       ks = ks + 1
       global_to_local_row( :, ks ) = operator_K_space%global_to_local( this_ks%k_indices, this_ks%spin, 'R' )
       global_to_local_row( :, ks ) = operator_K_space%global_to_local( this_ks%k_indices, this_ks%spin, 'C' )
    End Do
    Call operator_K_space%iterator_reset()

    ! From this generate a list of shells this process holds at least part of
    Allocate( i_own_row( 1:n_shells, 1:n_ks_points ) )
    Allocate( i_own_col( 1:n_shells, 1:n_ks_points ) )
    ! And populate the array
    Call operator_K_space%iterator_init()
    ks = 0
    this_ks = operator_K_space%iterator_next()
    Do While( this_ks%k_type /= K_POINT_NOT_EXIST )
       ks = ks + 1
       shell_start = 1
       Do i_shell = 1, n_shells
          shell_finish = shell_start + size_shells( i_shell ) - 1
          i_own_row( i_shell, ks ) = Any( global_to_local_row( shell_start:shell_finish, ks ) > 0 )
          i_own_col( i_shell, ks ) = Any( global_to_local_col( shell_start:shell_finish, ks ) > 0 )
          shell_start = shell_finish + 1
       End Do
    End Do
    Call operator_K_space%iterator_reset()

    ! Now for each of the shells work out where my contribution starts, and where it finishes
    ! From this generate a list of shells this process holds at least part of
    Allocate( shell_start_row( 1:n_shells, 1:n_ks_points ) )
    Allocate( shell_start_col( 1:n_shells, 1:n_ks_points ) )
    shell_start_row = NOT_ME
    shell_start_col = NOT_ME
    Allocate( shell_finish_row( 1:n_shells, 1:n_ks_points ) )
    Allocate( shell_finish_col( 1:n_shells, 1:n_ks_points ) )
    shell_finish_row = NOT_ME
    shell_finish_col = NOT_ME
    ! And populate the array
    Call operator_K_space%iterator_init()
    ks = 0
    this_ks = operator_K_space%iterator_next()
    Do While( this_ks%k_type /= K_POINT_NOT_EXIST )
       ks = ks + 1
       shell_start = 1
       Do i_shell = 1, n_shells
          shell_finish = shell_start + size_shells( i_shell ) - 1
          If( i_own_row( i_shell, ks ) ) Then
             Call find_shell_start_and_finish( global_to_local_row( :, ks ), shell_start, shell_finish, &
                  shell_start_row ( i_shell, ks ), shell_finish_row( i_shell, ks ) )
          End If
          If( i_own_col( i_shell, ks ) ) Then
             Call find_shell_start_and_finish( global_to_local_col( :, ks ), shell_start, shell_finish, &
                  shell_start_col ( i_shell, ks ), shell_finish_col( i_shell, ks ) )
          End If
          shell_start = shell_finish + 1
       End Do
    End Do
    Call operator_K_space%iterator_reset()

    ! REMEMBER DO DIAG TERMS ONCE!! - or do we need to as only writing ... need to think
    ! ALSO NEED ARGUMENT FOR COMPLEX CONJUGATION FOR ONE TRAINGLE COMPARED TO OTHER
    ! Now can fourier transfrom and produce one triangle of the K space matrices

    ! For moment do all in one routine - thinks about splitting into upper and lower traingles later
    Call do_lattice_FT( n_ao, n_shells, size_shells, operator_G_space, &
         i_own_row, shell_start_row, shell_finish_row, &
         i_own_col, shell_start_col, shell_finish_col, &
         operator_K_space )
    ! Forget this idea for the moment and do it simply
!!$    ! And by transposing the ownership arrays we can do the other
!!$    Call do_lattice_FT( n_ao, n_shells, size_shells, operator_G_space, &
!!$         i_own_col, shell_start_col, shell_finish_col, &
!!$         i_own_row, shell_start_row, shell_finish_row, &
!!$         operator_K_space )

  End Subroutine MPP_Lattice_FT

  Subroutine find_shell_start_and_finish( global_to_local, shell_start, shell_finish, my_shell_start, my_shell_finish )

    Integer, Dimension( : ), Intent( In    ) :: global_to_local
    Integer                , Intent( In    ) :: shell_start
    Integer                , Intent( In    ) :: shell_finish
    Integer                , Intent(   Out ) :: my_shell_start
    Integer                , Intent(   Out ) :: my_shell_finish

    Integer :: first_index, last_index

    ! I own a bit of this shell, work out the first index I hold
    first_index = shell_start
    Do While( global_to_local( first_index ) <= 0 .And. first_index < shell_finish )
       first_index = first_index + 1
    End Do
    ! Now find last index we hold being careful not to over run the end
    last_index = first_index
    Do While( global_to_local(  last_index ) >  0 .And. last_index < shell_finish )
       last_index = last_index + 1
    End Do
    my_shell_start  = first_index
    my_shell_finish = last_index
       
  End Subroutine find_shell_start_and_finish

  Subroutine do_lattice_FT( n_ao, n_shells, size_shells, operator_G_space, &
       i_own_row, shell_start_row, shell_finish_row, &
       i_own_col, shell_start_col, shell_finish_col, &
       operator_K_space )

    Use numbers        , Only : wp => float
    Use ks_array_module, Only : ks_array, ks_point_info, K_POINT_NOT_EXIST

    Integer                      , Intent( In    ) :: n_ao
    Integer                      , Intent( In    ) :: n_shells
    Integer   , Dimension( :    ), Intent( In    ) :: size_shells
    Real( wp ), Dimension( :    ), Intent( In    ) :: operator_G_space
    Logical   , Dimension( :, : ), Intent( In    ) :: i_own_row
    Integer   , Dimension( :, : ), Intent( In    ) :: shell_start_row
    Integer   , Dimension( :, : ), Intent( In    ) :: shell_finish_row
    Logical   , Dimension( :, : ), Intent( In    ) :: i_own_col
    Integer   , Dimension( :, : ), Intent( In    ) :: shell_start_col
    Integer   , Dimension( :, : ), Intent( In    ) :: shell_finish_col
    Type( ks_array )             , Intent( InOut ) :: operator_K_space

    Type( ks_point_info ) :: this_ks

    Type FT_coeffs_container !!!!FINISH FROM HERE !!!!!!

    Integer :: ks, n_ks
    Integer :: n_G_vectors, start_G_vectors
    Integer :: la1, ll2, la2

    ! Calculate the FT coefficients at all K relevant to me
    n_ks = Size( i_own_row, Dim = 2 )
    Call operator_K_space%iterator_init()
    ks = 0
    this_ks = operator_K_space%iterator_next()
    Calc_FT_Coeffs_loop: Do While( this_ks%k_type /= K_POINT_NOT_EXIST )
       ks = ks + 1
       Call expu( this_ks%k_indices( 1 ), this_ks%k_indices( 2 ), this_ks%k_indices( 3 ) )
    End Do Calc_FT_Coeffs_loop

    ! Loop over all the shells
    Outer_shell_loop: Do la1 = 1, n_shells

       ! Work out if this process store any contribution to this shell
       Hold_outer_shell: If( Any( i_own_row( la1, : ) ) .Or. Any( i_own_col( la1, : ) ) ) Then

          ! Loop over all shells that couple with the outer shell
          Coupled_shells_loop: Do ll2 = icct( la1 ) + icc( la1 ), icct( la1 + 1 )

             ! Work out which inner shell this is
             la2 = ila12t( ll2 )

             ! Only need to do something if this process holds some data for this shell couple
             Hold_inner_shell: If( Any( i_own_row( la2, : ) ) .Or. Any( i_own_col( la2, : ) ) ) Then

                ! Only need to do something if there actually any any non-zero elements
                ! for this couple
                n_G_vectors = idimfc( jpoint( ll2 ) )
                Coupling_vectors_exist: If( n_G_vectors  /= 0 ) Then

                   ! Getting here means that there is data to FT held on this process

                   ! Find the list of G vectors relevant to this couple
                   start_G_vectors = iccs3( ll2 )
                   
                End If Coupling_vectors_exist
                
             End If Hold_inner_shell
             
          End Do Coupled_shells_loop
          
       End If Hold_outer_shell
       
    End Do

  End Subroutine do_lattice_FT
    
    
End Module MPP_Fourier_transform
