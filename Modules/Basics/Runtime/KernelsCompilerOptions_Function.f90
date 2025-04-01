!-- KernelsCompilerOptions return the compiler options used compiling 
!   for kernels, commands, and functions.

module KernelsCompilerOptions_Function
  
  use iso_fortran_env

  implicit none
  private
  
  public :: &
    KernelsCompilerOptions
    
contains


  function KernelsCompilerOptions ( ) result ( KCO )
    
    character ( len = : ), allocatable :: &
      KCO
    
    KCO = compiler_options ( )
    
  end function KernelsCompilerOptions


end module KernelsCompilerOptions_Function
