Module numbers_module

  !! Little module to set the precision of real numbers that we shall use
  
  Implicit None 
  
  Integer, Public, Parameter :: wp = Selected_real_kind( 14, 70 )  !! Kind for real numbers

  Private
  
End Module numbers_module

