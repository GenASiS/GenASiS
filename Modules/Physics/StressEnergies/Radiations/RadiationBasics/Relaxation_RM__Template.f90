module Relaxation_RM__Template

  use Basics
  use Mathematics

  implicit none
  private
  
  type, public, abstract :: Relaxation_RM_Template
    integer ( KDI ) :: &
      IGNORABILITY = 0
    character ( LDL ) :: &
      Name = '', &
      Type = ''
    type ( LinearEquations_LAPACK_Form ), dimension ( : ), allocatable :: &
      LinearEquations
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure, public, pass :: &
      FinalizeTemplate
  end type Relaxation_RM_Template

contains

  subroutine InitializeTemplate ( R, Name )

    class ( Relaxation_RM_Template ), intent ( inout ) :: &
      R
    character ( * ), intent ( in ) :: &
      Name

    if ( R % Type == '' ) &
      R % Type  =  'a Relaxation_RM'

    R % Name  =  Name

    R % IGNORABILITY = CONSOLE % INFO_1 
    call Show ( 'Initializing ' // trim ( R % Type ), R % IGNORABILITY )
    call Show ( R % Name, 'Name', R % IGNORABILITY )

  end subroutine InitializeTemplate


  subroutine FinalizeTemplate ( R )

    class ( Relaxation_RM_Template ), intent ( inout ) :: &
      R

    if ( allocated ( R % LinearEquations ) ) &
      deallocate ( R % LinearEquations )

    call Show ( 'Finalizing ' // trim ( R % Type ), R % IGNORABILITY )
    call Show ( R % Name, 'Name', R % IGNORABILITY )

  end subroutine FinalizeTemplate


end module Relaxation_RM__Template
