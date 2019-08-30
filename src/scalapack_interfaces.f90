Module scalapack_interfaces

  !! Interfaces to the linear algebra routines

  Use numbers_module, Only : wp

  Implicit None
  
  Interface

     !! Number of rows or columns distributed along this dimension
 
     Pure Function numroc( n, nb, iproc, isrcproc, nprocs  ) Result( loc_n )
       Implicit None
       Integer               :: loc_n
       Integer, Intent( In ) :: n         !! Length of Dimension
       Integer, Intent( In ) :: nb        !! Blocking Factor
       Integer, Intent( In ) :: iproc     !! Me
       Integer, Intent( In ) :: isrcproc  !! First proc holding the dim
       Integer, Intent( In ) :: nprocs    !! number of procs holding the dim
     End Function numroc

  End Interface

  Interface

     !! Matrix multiplies
  
     Subroutine pdgemm( ta, tb, m, n, k, alpha, A, ia, ja, descA, B, ib, jb, descB, beta, C, ic, jc, descC )
       Import :: wp
       Implicit None
       Character                   , Intent( In    ) :: ta
       Character                   , Intent( In    ) :: tb
       Integer                     , Intent( In    ) :: m
       Integer                     , Intent( In    ) :: n
       Integer                     , Intent( In    ) :: k
       Real( wp )                  , Intent( In    ) :: alpha
       Real( wp ), Dimension(  *  ), Intent( In    ) :: A
       Integer                     , Intent( In    ) :: ia
       Integer                     , Intent( In    ) :: ja
       Integer,    Dimension( 1:9 ), Intent( In    ) :: descA
       Real( wp ), Dimension(  *  ), Intent( In    ) :: B
       Integer                     , Intent( In    ) :: ib
       Integer                     , Intent( In    ) :: jb
       Integer,    Dimension( 1:9 ), Intent( In    ) :: descB
       Real( wp )                  , Intent( In    ) :: beta
       Real( wp ), Dimension(  *  ), Intent( InOut ) :: C
       Integer                     , Intent( In    ) :: ic
       Integer                     , Intent( In    ) :: jc
       Integer,    Dimension( 1:9 ), Intent( In    ) :: descC
     End Subroutine pdgemm

     Subroutine pzgemm( ta, tb, m, n, k, alpha, A, ia, ja, descA, B, ib, jb, descB, beta, C, ic, jc, descC )
       Import :: wp
       Implicit None
       Character                      , Intent( In    ) :: ta
       Character                      , Intent( In    ) :: tb
       Integer                        , Intent( In    ) :: m
       Integer                        , Intent( In    ) :: n
       Integer                        , Intent( In    ) :: k
       Complex( wp )                  , Intent( In    ) :: alpha
       Complex( wp ), Dimension(  *  ), Intent( In    ) :: A
       Integer                        , Intent( In    ) :: ia
       Integer                        , Intent( In    ) :: ja
       Integer,       Dimension( 1:9 ), Intent( In    ) :: descA
       Complex( wp ), Dimension(  *  ), Intent( In    ) :: B
       Integer                        , Intent( In    ) :: ib
       Integer                        , Intent( In    ) :: jb
       Integer,       Dimension( 1:9 ), Intent( In    ) :: descB
       Complex( wp )                  , Intent( In    ) :: beta
       Complex( wp ), Dimension(  *  ), Intent( InOut ) :: C
       Integer                        , Intent( In    ) :: ic
       Integer                        , Intent( In    ) :: jc
       Integer,       Dimension( 1:9 ), Intent( In    ) :: descC
     End Subroutine pzgemm

  End Interface

  Interface

     !! Matrix addition routines

     Subroutine pdgeadd( trans, m, n, alpha, a, ia, ja, desca, beta, c, ic, jc, descc )
       Import :: wp
       Implicit None
       Character                                     :: trans
       Integer                     , Intent( In    ) :: m
       Integer                     , Intent( In    ) :: n
       Real( wp )                  , Intent( In    ) :: alpha
       Real( wp ), Dimension( *   ), Intent( In    ) :: a
       Integer                     , Intent( In    ) :: ia
       Integer                     , Intent( In    ) :: ja
       Integer   , Dimension( 1:9 ), Intent( In    ) :: desca
       Real( wp )                  , Intent( In    ) :: beta
       Real( wp ), Dimension( *   ), Intent( InOut ) :: c
       Integer                     , Intent( In    ) :: ic
       Integer                     , Intent( In    ) :: jc
       Integer   , Dimension( 1:9 ), Intent( In    ) :: descc
     End Subroutine pdgeadd
     
     Subroutine pzgeadd( trans, m, n, alpha, a, ia, ja, desca, beta, c, ic, jc, descc )
       Import :: wp
       Implicit None
       Character                                        :: trans
       Integer                        , Intent( In    ) :: m
       Integer                        , Intent( In    ) :: n
       Complex( wp )                  , Intent( In    ) :: alpha
       Complex( wp ), Dimension( *   ), Intent( In    ) :: a
       Integer                        , Intent( In    ) :: ia
       Integer                        , Intent( In    ) :: ja
       Integer      , Dimension( 1:9 ), Intent( In    ) :: desca
       Complex( wp )                  , Intent( In    ) :: beta
       Complex( wp ), Dimension( *   ), Intent( InOut ) :: c
       Integer                        , Intent( In    ) :: ic
       Integer                        , Intent( In    ) :: jc
       Integer      , Dimension( 1:9 ), Intent( In    ) :: descc
     End Subroutine pzgeadd
     
  End Interface

  Interface

     !! Matrix diagonalisation routines

     Subroutine pdsyevd( jobz, uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, work, lwork, iwork, liwork, info )
       Import :: wp
       Implicit None
       Character                   , Intent( In    ) :: jobz
       Character                   , Intent( In    ) :: uplo
       Integer                     , Intent( In    ) :: n
       Real( wp ), Dimension( *   ), Intent( InOut ) :: a
       Integer                     , Intent( In    ) :: ia
       Integer                     , Intent( In    ) :: ja
       Integer   , Dimension( 1:9 ), Intent( In    ) :: desca
       Real( wp ), Dimension( *   ), Intent(   Out ) :: w
       Real( wp ), Dimension( *   ), Intent(   Out ) :: z
       Integer                     , Intent( In    ) :: iz
       Integer                     , Intent( In    ) :: jz
       Integer   , Dimension( 1:9 ), Intent( In    ) :: descz
       Real( wp ), Dimension( *   ), Intent(   Out ) :: work
       Integer                     , Intent( In    ) :: lwork
       Integer   , Dimension( *   ), Intent(   Out ) :: iwork
       Integer                     , Intent( In    ) :: liwork
       Integer                     , Intent(   Out ) :: info
     End Subroutine pdsyevd

     Subroutine pzheevd( jobz, uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, cwork, lcwork, rwork, lrwork, iwork, liwork, info )
       Import :: wp
       Implicit None
       Character                      , Intent( In    ) :: jobz
       Character                      , Intent( In    ) :: uplo
       Integer                        , Intent( In    ) :: n
       Complex( wp ), Dimension( *   ), Intent( InOut ) :: a
       Integer                        , Intent( In    ) :: ia
       Integer                        , Intent( In    ) :: ja
       Integer      , Dimension( 1:9 ), Intent( In    ) :: desca
       Real   ( wp ), Dimension( *   ), Intent(   Out ) :: w
       Complex( wp ), Dimension( *   ), Intent(   Out ) :: z
       Integer                        , Intent( In    ) :: iz
       Integer                        , Intent( In    ) :: jz
       Integer      , Dimension( 1:9 ), Intent( In    ) :: descz
       Complex( wp ), Dimension( *   ), Intent(   Out ) :: cwork
       Integer                        , Intent( In    ) :: lcwork
       Real( wp )   , Dimension( *   ), Intent(   Out ) :: rwork
       Integer                        , Intent( In    ) :: lrwork
       Integer      , Dimension( *   ), Intent(   Out ) :: iwork
       Integer                        , Intent( In    ) :: liwork
       Integer                        , Intent(   Out ) :: info
     End Subroutine pzheevd

  End Interface

  Interface

     !! Choleski decompostion routines

     Subroutine pdpotrf( uplo, n, a, ia, ja, desca, info )
       Import :: wp
       Character                   , Intent( In    ) :: uplo
       Integer                     , Intent( In    ) :: n
       Real( wp ), Dimension( *   ), Intent( InOut ) :: a
       Integer                     , Intent( In    ) :: ia
       Integer                     , Intent( In    ) :: ja
       Integer   , Dimension( 1:9 ), Intent( In    ) :: desca
       Integer                        , Intent(   Out ) :: info
     End Subroutine pdpotrf
     
     Subroutine pzpotrf( uplo, n, a, ia, ja, desca, info )
       Import :: wp
       Character                      , Intent( In    ) :: uplo
       Integer                        , Intent( In    ) :: n
       Complex( wp ), Dimension( *   ), Intent( InOut ) :: a
       Integer                        , Intent( In    ) :: ia
       Integer                        , Intent( In    ) :: ja
       Integer      , Dimension( 1:9 ), Intent( In    ) :: desca
       Integer                        , Intent(   Out ) :: info
     End Subroutine pzpotrf
     
  End Interface
  
  Interface

     !! Triangular invert routines

     Subroutine pdtrtri( uplo, diag, n, a, ia, ja, desca, info )
       Import :: wp
       Character                   , Intent( In    ) :: uplo
       Character                   , Intent( In    ) :: diag
       Integer                     , Intent( In    ) :: n
       Real( wp ), Dimension( *   ), Intent( InOut ) :: a
       Integer                     , Intent( In    ) :: ia
       Integer                     , Intent( In    ) :: ja
       Integer   , Dimension( 1:9 ), Intent( In    ) :: desca
       Integer                        , Intent(   Out ) :: info
     End Subroutine pdtrtri
     
     Subroutine pztrtri( uplo, diag, n, a, ia, ja, desca, info )
       Import :: wp
       Character                      , Intent( In    ) :: uplo
       Character                   , Intent( In    ) :: diag
       Integer                        , Intent( In    ) :: n
       Complex( wp ), Dimension( *   ), Intent( InOut ) :: a
       Integer                        , Intent( In    ) :: ia
       Integer                        , Intent( In    ) :: ja
       Integer      , Dimension( 1:9 ), Intent( In    ) :: desca
       Integer                        , Intent(   Out ) :: info
     End Subroutine pztrtri
     
  End Interface
  
  Interface

     !! Redistribution routines

     Subroutine pdgemr2d( m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt )
       Import wp
       Integer                     , Intent( In    ) :: m
       Integer                     , Intent( In    ) :: n
       Real( wp ), Dimension( *   ), Intent( In    ) :: a
       Integer                     , Intent( In    ) :: ia
       Integer                     , Intent( In    ) :: ja
       Integer   , Dimension( 1:9 ), Intent( In    ) :: desca
       Real( wp ), Dimension( *   ), Intent(   Out ) :: b
       Integer                     , Intent( In    ) :: ib
       Integer                     , Intent( In    ) :: jb
       Integer   , Dimension( 1:9 ), Intent( In    ) :: descb
       Integer                     , Intent( In    ) :: ictxt
     End Subroutine pdgemr2d
      
     Subroutine pzgemr2d( m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt )
       Import wp
       Integer                        , Intent( In    ) :: m
       Integer                        , Intent( In    ) :: n
       Complex( wp ), Dimension( *   ), Intent( In    ) :: a
       Integer                        , Intent( In    ) :: ia
       Integer                        , Intent( In    ) :: ja
       Integer      , Dimension( 1:9 ), Intent( In    ) :: desca
       Complex( wp ), Dimension( *   ), Intent(   Out ) :: b
       Integer                        , Intent( In    ) :: ib
       Integer                        , Intent( In    ) :: jb
       Integer      , Dimension( 1:9 ), Intent( In    ) :: descb
       Integer                        , Intent( In    ) :: ictxt
     End Subroutine pzgemr2d
      
  End Interface

End Module scalapack_interfaces
