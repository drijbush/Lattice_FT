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

MODULE EXPO_MODULE
  USE NUMBERS, Only : float
  REAL(FLOAT),DIMENSION(:,:),ALLOCATABLE :: EX,EXVRS
END MODULE EXPO_MODULE

Module MPP_Fourier_transform

  Use numbers        , Only : wp => float
  
  Implicit None
    
  Public :: MPP_Lattice_FT
  
  Private
  
  Integer, Private, Parameter :: NOT_ME = -2

  Type FT_coeffs_container
     Private
     Real   ( wp ), Dimension( : ), Allocatable :: real_ft_coeffs
     Complex( wp ), Dimension( : ), Allocatable :: cmplx_ft_coeffs
  End type FT_coeffs_container


Contains

  Subroutine MPP_Lattice_FT( get_xg, n_ao, n_shells, size_shells, operator_G_space, &
       ft_coeffs, ila12t, idimfc, jpoint, ngshg, nqgshg, iccs3, icc, icct, &
       operator_K_space )

    Use numbers        , Only : wp => float
    Use ks_array_module, Only : ks_array, ks_point_info, K_POINT_NOT_EXIST

    Integer                   ,                    Intent( In    ) :: get_xg
    Integer                   ,                    Intent( In    ) :: n_ao
    Integer                   ,                    Intent( In    ) :: n_shells
    Integer                   , Dimension( :    ), Intent( In    ) :: size_shells
    Real( wp )                , Dimension( :    ), Intent( In    ) :: operator_G_space
    Integer                   , Dimension( :    ), Intent( In    ) :: ila12t
    Integer                   , Dimension( :    ), Intent( In    ) :: idimfc
    Integer                   , Dimension( :    ), Intent( In    ) :: jpoint
    Integer                   , Dimension( :    ), Intent( In    ) :: ngshg
    Integer                   , Dimension( :    ), Intent( In    ) :: nqgshg
    Integer                   , Dimension( :    ), Intent( In    ) :: iccs3
    Integer                   , Dimension( :    ), Intent( In    ) :: icc
    Integer                   , Dimension( :    ), Intent( In    ) :: icct
    Type( ks_array )          ,                    Intent( InOut ) :: operator_K_space

    Type( FT_coeffs_container ), Dimension( : ), Allocatable :: ft_coeffs
    
    Type( ks_point_info ) :: this_ks

    Integer, Dimension( :, : ), Allocatable :: global_to_local_row
    Integer, Dimension( :, : ), Allocatable :: global_to_local_col

    Integer, Dimension( :, : ), Allocatable :: shell_start_row
    Integer, Dimension( :, : ), Allocatable :: shell_start_col

    Integer, Dimension( :, : ), Allocatable :: shell_finish_row
    Integer, Dimension( :, : ), Allocatable :: shell_finish_col
    
    Integer :: max_size_shell
    Integer :: min_spin, max_spin
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
    min_spin =   Huge( min_spin )
    max_spin = - Huge( max_spin )
    Do While( this_ks%k_type /= K_POINT_NOT_EXIST )
       n_ks_points = n_ks_points + 1
       this_ks = operator_K_space%iterator_next()
       min_spin = Min( min_spin, this_ks%spin )
       max_spin = Max( max_spin, this_ks%spin )
    End Do
    Call operator_K_space%iterator_reset()

    ! Due to k point parallelism this process may not actually hold any data. We could at this point
    ! have a If( n_ks /= 0 ) BUT the code should work without that and it just adds complications. So
    ! let's miss it out - it will help find bugs.

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
       global_to_local_col( :, ks ) = operator_K_space%global_to_local( this_ks%k_indices, this_ks%spin, 'C' )
       this_ks = operator_K_space%iterator_next()
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
       this_ks = operator_K_space%iterator_next()
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
             shell_start_row ( i_shell, ks ) = shell_start
             shell_finish_row( i_shell, ks ) = shell_finish
          End If
          If( i_own_col( i_shell, ks ) ) Then
             shell_start_col ( i_shell, ks ) = shell_start
             shell_finish_col( i_shell, ks ) = shell_finish
          End If
          shell_start = shell_finish + 1
       End Do
       this_ks = operator_K_space%iterator_next()
    End Do
    Call operator_K_space%iterator_reset()

    ! Now calculate the FT coeffs (complex exponentials) for each k point
    Call calc_ft_coeffs( n_ks_points, operator_K_space, ft_coeffs )
    
    ! REMEMBER DO DIAG TERMS ONCE!! - or do we need to as only writing ... need to think
    ! ALSO NEED ARGUMENT FOR COMPLEX CONJUGATION FOR ONE TRAINGLE COMPARED TO OTHER
    ! Now can fourier transfrom and produce one triangle of the K space matrices

    ! For moment do all in one routine - thinks about splitting into upper and lower traingles later
    Call do_lattice_FT( get_xg, min_spin, max_spin, n_ao, n_shells, size_shells, operator_G_space, &
         i_own_row, shell_start_row, shell_finish_row, &
         i_own_col, shell_start_col, shell_finish_col, &
         ft_coeffs, ila12t, idimfc, jpoint, ngshg, nqgshg, iccs3, icc, icct, &
         operator_K_space )
    
    ! Forget this idea for the moment and do it simply
!!$    ! And by transposing the ownership arrays we can do the other
!!$    Call do_lattice_FT( n_ao, n_shells, size_shells, operator_G_space, &
!!$         i_own_col, shell_start_col, shell_finish_col, &
!!$         i_own_row, shell_start_row, shell_finish_row, &
!!$         operator_K_space )

  End Subroutine MPP_Lattice_FT

  Subroutine calc_ft_coeffs( n_ks, operator_K_space, ft_coeffs )

    ! Calculate the FT coeffs (complex exponentials) for each k point
    
    Use numbers        , Only : wp => float
    Use Expo_module    , Only : ex
    Use ks_array_module, Only : ks_array, ks_point_info, K_POINT_NOT_EXIST, K_POINT_COMPLEX

    Integer                                                 , Intent( In    ) :: n_ks
    Type( ks_array )                                        , Intent( InOut ) :: operator_K_space
    Type( FT_coeffs_container ), Dimension( : ), Allocatable, Intent(   Out ) :: ft_coeffs

    Interface
       Subroutine expt( j1, j2, j3 )
         Implicit None
         Integer, Intent( In ) :: j1, j2, j3
       End Subroutine expt
       Subroutine expu( j1, j2, j3 )
         Implicit None
         Integer, Intent( In ) :: j1, j2, j3
       End Subroutine expu
    End Interface
    
    Type( ks_point_info ) :: this_ks

    Integer :: ks
    
    Allocate( ft_coeffs( 1:n_ks ) )
    Call operator_K_space%iterator_init()
    ks = 0
    this_ks = operator_K_space%iterator_next()
    Calc_FT_Coeffs_loop: Do While( this_ks%k_type /= K_POINT_NOT_EXIST )
       ks = ks + 1
       Complex_or_real: If( this_ks%k_type == K_POINT_COMPLEX ) Then
          Call expu( this_ks%k_indices( 1 ), this_ks%k_indices( 2 ), this_ks%k_indices( 3 ) )
          ft_coeffs( ks )%cmplx_ft_coeffs = Cmplx( ex( 1, : ), ex( 2, : ), Kind = wp )
       Else
          Call expt( this_ks%k_indices( 1 ), this_ks%k_indices( 2 ), this_ks%k_indices( 3 ) )
          ft_coeffs( ks )%real_ft_coeffs = ex( 1, : )
       End If Complex_or_real
       this_ks = operator_K_space%iterator_next()
    End Do Calc_FT_Coeffs_loop
    Call operator_K_space%iterator_reset()

  End Subroutine calc_ft_coeffs

  Subroutine do_lattice_FT( get_xg, min_spin, max_spin, n_ao, n_shells, size_shells, operator_G_space, &
       i_own_row, shell_start_row, shell_finish_row, &
       i_own_col, shell_start_col, shell_finish_col, &
       ft_coeffs, ila12t, idimfc, jpoint, ngshg, nqgshg, iccs3, icc, icct, &
       operator_K_space )

    Use numbers        , Only : wp => float
    Use ks_array_module, Only : ks_array, ks_point_info, K_POINT_NOT_EXIST, K_POINT_COMPLEX

    Integer                                       , Intent( In    ) :: get_xg
    Integer                                       , Intent( In    ) :: min_spin
    Integer                                       , Intent( In    ) :: max_spin
    Integer                                       , Intent( In    ) :: n_ao
    Integer                                       , Intent( In    ) :: n_shells
    Integer                    , Dimension( :    ), Intent( In    ) :: size_shells
    Real( wp )                 , Dimension( :    ), Intent( In    ) :: operator_G_space
    Logical                    , Dimension( :, : ), Intent( In    ) :: i_own_row
    Integer                    , Dimension( :, : ), Intent( In    ) :: shell_start_row
    Integer                    , Dimension( :, : ), Intent( In    ) :: shell_finish_row
    Logical                    , Dimension( :, : ), Intent( In    ) :: i_own_col
    Integer                    , Dimension( :, : ), Intent( In    ) :: shell_start_col
    Integer                    , Dimension( :, : ), Intent( In    ) :: shell_finish_col
    Type( FT_coeffs_container ), Dimension( :    ), Intent( In    ) :: ft_coeffs
    Integer                    , Dimension( :    ), Intent( In    ) :: ila12t
    Integer                    , Dimension( :    ), Intent( In    ) :: idimfc
    Integer                    , Dimension( :    ), Intent( In    ) :: jpoint
    Integer                    , Dimension( :    ), Intent( In    ) :: ngshg
    Integer                    , Dimension( :    ), Intent( In    ) :: nqgshg
    Integer                    , Dimension( :    ), Intent( In    ) :: iccs3
    Integer                    , Dimension( :    ), Intent( In    ) :: icc
    Integer                    , Dimension( :    ), Intent( In    ) :: icct
    Type( ks_array )                              , Intent( InOut ) :: operator_K_space

    Complex( wp ), Dimension( :, : ), Allocatable :: fk_bit_2d_complex

    Real( wp ), Dimension( :, : ), Allocatable :: fg_red_bit
    
    Real( wp ), Dimension( :, : ), Allocatable :: fk_bit_2d_real

    Real( wp ), Dimension( : ), Allocatable :: fk_bit

    Real( wp ), Dimension( :, : ), Allocatable :: ex 
    
    Type( ks_point_info ) :: this_ks

    Integer :: max_size_shell
    Integer :: ks
    Integer :: n_G_vectors, start_G_vectors
    Integer :: spin
    Integer :: nbf_la1, nbf_la2, nbf12
    Integer :: la1, ll2, la2
    Integer :: ierr
    
    ! Some basic information about the system
    max_size_shell = Maxval( size_shells( 1:n_shells ) )

    ! Bit of memory to store the bit of the reducible G space represantation of the 
    ! matrix that we generate
    Allocate( fg_red_bit( 1:max_size_shell * max_size_shell, min_spin:max_spin ) )

    ! Similarly for the bit of the reciprocal space matrix. Factor of 2 for complex 
    Allocate( fk_bit( 1:2 * max_size_shell * max_size_shell ) )

    ! And the 2d representations of the bit of the reciprocal space matrix.
    Allocate( fk_bit_2d_complex( 1:max_size_shell, 1:max_size_shell ) )
    Allocate( fk_bit_2d_real   ( 1:max_size_shell, 1:max_size_shell ) )
    
    ! Initialise the k space representation of the operator
    operator_K_space = 0.0_wp
    
    ! Temporary storage for FT coeffs for a given shell couple at a given k point
    ! Avoid taking size of an unallocated array, which is not alloed by the standard
    ! The size of the first element is enough as all members should be the same size
    If( Allocated( ft_coeffs( 1 )%cmplx_ft_coeffs ) ) Then
       Allocate( ex( 1:Size( ft_coeffs( 1 )%cmplx_ft_coeffs ), 1:2 ) )
    Else
       Allocate( ex( 1:Size( ft_coeffs( 1 )%real_ft_coeffs ), 1:2  ) )
    End If

    ! Loop over all the shells
    Outer_shell_loop: Do la1 = 1, n_shells

       ! Number of basis functions in this shell
       nbf_la1 = size_shells( la1 )

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

                   ! Number of basis functions in this shell
                   nbf_la2 = size_shells( la2 )
                   nbf12   = nbf_la1 * nbf_la2
                   
                   ! Find the list of G vectors relevant to this couple
                   start_G_vectors = ngshg( iccs3( ll2 ) + 1 )
                   n_G_vectors     = ngshg( iccs3( ll2 ) + 1 + n_g_vectors ) - start_G_vectors

                   ! Calculate the bits of the reducible matrix we need
                   Real_space_spin_loop: Do spin = min_spin, max_spin
                      !!!! REALLY OUGHT TO FIX INTERFACE SO USES OPERATOR_G_SPACE and not through "common" !!!!!!
                      Call get_red_l1_l2( get_xg, spin, la2, la1, fg_red_bit( 1:nbf12, spin ), ierr )
                   End Do Real_space_spin_loop
                   
                   ! Loop over ks points and spins
                   Call operator_K_space%iterator_init()
                   ks = 0
                   this_ks = operator_K_space%iterator_next()
                   KS_point_loop: Do While( this_ks%k_type /= K_POINT_NOT_EXIST )
                      ks = ks + 1
                      spin = this_ks%spin
                      Select Case( this_ks%k_type )
                      Case( K_POINT_COMPLEX )
                         Call FT_complex
                      Case( K_POINT_REAL )
                         Call FT_real
                      End Select
                      this_ks = operator_K_space%iterator_next()
                   End Do KS_point_loop
                   Call operator_K_space%iterator_reset()
                
                End If Coupling_vectors_exist

             End If Hold_inner_shell
             
          End Do Coupled_shells_loop
          
       End If Hold_outer_shell
       
    End Do Outer_shell_loop

  Contains

    Subroutine FT_complex

      Integer :: this_G
      Integer :: shell_shell_block_size
      Integer :: ind_r, ind_c, ibf1, ibf2
      Integer :: i

      ! Get the FT coeffs that are needed at this k point
      Do i = 1, n_G_vectors
         this_G = nqgshg( start_G_vectors + i )
         ex( i, 1 ) = Real ( ft_coeffs( ks )%cmplx_ft_coeffs( this_G ), Kind = wp )
         ex( i, 2 ) = Aimag( ft_coeffs( ks )%cmplx_ft_coeffs( this_G )            )
      End Do

      ! Storage of reducible matrix differs for on and off diagonal blocks
      Diag_or_off_diag_block: If( la1 == la2 ) Then

         ! Diagonal Block
         shell_shell_block_size = ( ( nbf_la1 + 1 ) * ( nbf_la1 + 2 ) ) / 2
         ! Zero result matrix - * 2 for complex numbers
         fk_bit( 1:shell_shell_block_size * 2 ) = 0.0_wp
         ! Do the Lattice transform on this block
         Call mxmb( fg_red_bit( :, this_ks%spin ), 1, shell_shell_block_size, ex, 1, Size( ex, Dim = 1 ), &
              fk_bit, 1, shell_shell_block_size, shell_shell_block_size, n_G_vectors, 2 )
         ! Put into 2D matrix
         ind_r = 0                      ! index of real part in fk_bit
         ind_c = shell_shell_block_size ! index of complex part
         ! Fill in the upper triangle
         Do ibf2 = 1, nbf_la2
            Do ibf1 = 1, ibf2
               ind_r = ind_r + 1
               ind_c = ind_c + 1
               fk_bit_2d_complex( ibf1, ibf2 ) = Cmplx( fk_bit( ind_r ), fk_bit( ind_c ), Kind = wp )
            End Do
         End Do
         ! And the lower traing by Hermiticity
         Do ibf2 = 1, nbf_la2
            Do ibf1 = ibf2 + 1, nbf_la1
               fk_bit_2d_complex( ibf2, ibf1 ) = Conjg( fk_bit_2d_complex( ibf1, ibf2 ) )
            End Do
         End Do
         ! Now store
         Call operator_K_space%set_by_global( this_ks%k_indices, this_ks%spin, &
              shell_start_row( la1, ks ), shell_finish_row( la1, ks ), &
              shell_start_col( la1, ks ), shell_finish_col( la1, ks ), &
              fk_bit_2d_complex )

      Else

         ! Off diagonal block
         shell_shell_block_size = nbf_la1 * nbf_la2
         ! Zero result matrix - * 2 for complex numbers
         fk_bit( 1:shell_shell_block_size * 2 ) = 0.0_wp
         ! Do the Lattice transform on this block
         Call mxmb( fg_red_bit( 1:shell_shell_block_size, this_ks%spin ), 1, shell_shell_block_size, ex, 1, Size( ex, Dim = 1 ), &
              fk_bit, 1, shell_shell_block_size, shell_shell_block_size, n_G_vectors, 2 )
         ! Put into 2D matrix
         ind_r = 0                      ! index of real part in fk_bit
         ind_c = shell_shell_block_size ! index of complex part
         Do ibf2 = 1, nbf_la2
            Do ibf1 = 1, nbf_la1
               ind_r = ind_r + 1
               ind_c = ind_c + 1
               fk_bit_2d_complex( ibf1, ibf2 ) = Cmplx( fk_bit( ind_r ), fk_bit( ind_c ), Kind = wp )
            End Do
         End Do
         ! Now store
         Call operator_K_space%set_by_global( this_ks%k_indices, this_ks%spin, &
              shell_start_row( la1, ks ), shell_finish_row( la1, ks ), &
              shell_start_col( la2, ks ), shell_finish_col( la2, ks ), &
              fk_bit_2d_complex )
         Call operator_K_space%set_by_global( this_ks%k_indices, this_ks%spin, &
              shell_start_row( la2, ks ), shell_finish_row( la2, ks ), &
              shell_start_col( la1, ks ), shell_finish_col( la1, ks ), &
              Conjg( Transpose( fk_bit_2d_complex ) ) )
         
      End If Diag_or_off_diag_block
      
    End Subroutine FT_complex

    Subroutine FT_real

      Integer :: this_G
      Integer :: shell_shell_block_size
      Integer :: ind_r, ind_c, ibf1, ibf2
      Integer :: i

      ! Get the FT coeffs that are needed at this k point
      Do i = 1, n_G_vectors
         this_G = nqgshg( start_G_vectors + i )
         ex( i, 1 ) = ft_coeffs( ks )%real_ft_coeffs( this_G )
      End Do

      ! Storage of reducible matrix differs for on and off diagonal blocks
      Diag_or_off_diag_block: If( la1 == la2 ) Then

         ! Diagonal Block
         shell_shell_block_size = ( ( nbf_la1 + 1 ) * ( nbf_la1 + 2 ) ) / 2
         ! Zero result matrix
         fk_bit( 1:shell_shell_block_size ) = 0.0_wp
         ! Do the Lattice transform on this block
         Call mxmb( fg_red_bit( :, this_ks%spin ), 1, shell_shell_block_size, ex, 1, 1, &
              fk_bit, 1, 1, shell_shell_block_size, n_G_vectors, 1 )
         ! Put into 2D matrix
         ind_r = 0                      ! index of part in fk_bit
         ! Fill in the upper triangle
         Do ibf2 = 1, nbf_la2
            Do ibf1 = 1, ibf2
               ind_r = ind_r + 1
               fk_bit_2d_real( ibf1, ibf2 ) = fk_bit( ind_r )
            End Do
         End Do
         ! And the lower traing by Hermiticity
         Do ibf2 = 1, nbf_la2
            Do ibf1 = ibf2 + 1, nbf_la1
               fk_bit_2d_real( ibf2, ibf1 ) = fk_bit_2d_real( ibf1, ibf2 )
            End Do
         End Do
         ! Now store
         Call operator_K_space%set_by_global( this_ks%k_indices, this_ks%spin, &
              shell_start_row( la1, ks ), shell_finish_row( la1, ks ), &
              shell_start_col( la1, ks ), shell_finish_col( la1, ks ), &
              fk_bit_2d_real )

      Else

         ! Off diagonal block
         shell_shell_block_size = nbf_la1 * nbf_la2
         ! Zero result matrix
         fk_bit( 1:shell_shell_block_size ) = 0.0_wp
         ! Do the Lattice transform on this block
         Call mxmb( fg_red_bit( 1:shell_shell_block_size, this_ks%spin ), 1, shell_shell_block_size, ex, 1, 1, &
              fk_bit, 1, 1, shell_shell_block_size, n_G_vectors, 1 )
         ! Put into 2D matrix
         ind_r = 0                      ! index of part in fk_bit
         Do ibf2 = 1, nbf_la2
            Do ibf1 = 1, nbf_la1
               ind_r = ind_r + 1
               fk_bit_2d_real( ibf1, ibf2 ) = fk_bit( ind_r )
            End Do
         End Do
         ! Now store
         Call operator_K_space%set_by_global( this_ks%k_indices, this_ks%spin, &
              shell_start_row( la1, ks ), shell_finish_row( la1, ks ), &
              shell_start_col( la2, ks ), shell_finish_col( la2, ks ), &
              fk_bit_2d_real )
         Call operator_K_space%set_by_global( this_ks%k_indices, this_ks%spin, &
              shell_start_row( la2, ks ), shell_finish_row( la2, ks ), &
              shell_start_col( la1, ks ), shell_finish_col( la1, ks ), &
              Transpose( fk_bit_2d_real ) )
         
      End If Diag_or_off_diag_block

    End Subroutine FT_real

  End Subroutine do_lattice_FT
    
    
End Module MPP_Fourier_transform
