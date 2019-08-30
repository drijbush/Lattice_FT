Module blacs_interfaces

  !! Interfaces for the BLACS routines we use
  
  Interface
     Subroutine blacs_gridinit( ctxt, order, nprow, npcol )
       Implicit None
       Integer             , Intent( InOut ) :: ctxt
       Character( Len = * ), Intent( In    ) :: order
       Integer             , Intent( In    ) :: nprow
       Integer             , Intent( In    ) :: npcol
     End Subroutine blacs_gridinit
     Subroutine blacs_gridinfo( ctxt, nprow, npcol, myprow, mypcol )
       Implicit None
       Integer             , Intent( In    ) :: ctxt
       Integer             , Intent(   Out ) :: nprow
       Integer             , Intent(   Out ) :: npcol
       Integer             , Intent(   Out ) :: myprow
       Integer             , Intent(   Out ) :: mypcol
     End Subroutine blacs_gridinfo
     Subroutine blacs_gridexit( ctxt )
       Implicit None
       Integer, Intent( InOut ) :: ctxt
     End Subroutine blacs_gridexit
     Subroutine blacs_exit( cont )
       Implicit None
       Integer, Intent( In ) :: cont
     End Subroutine blacs_exit
  End Interface

End Module blacs_interfaces
