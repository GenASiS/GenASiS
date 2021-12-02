module Coarsening_C__Form

  !-- Coarsening_Central_Form

  use Basics
  use Manifolds
  use FieldSets
  use Geometry_F__Form

  implicit none
  private

  type, public, extends ( FieldSetForm ) :: Coarsening_C_Form
    integer ( KDI ) :: &
      COARSENING_2 = 0, &
      COARSENING_3 = 0
  contains
    procedure, private, pass :: &
      Initialize_C
    generic, public :: &
      Initialize => Initialize_C
    final :: &
      Finalize
  end type Coarsening_C_Form

    private :: &
      SetCoarseningPolar

contains


  subroutine Initialize_C ( C, G )

    class ( Coarsening_C_Form ), intent ( inout ), target :: &
      C
    class ( Geometry_F_Form ), intent ( in ), target :: &
      G

    if ( C % Type == '' ) &
      C % Type  =  'a Coarsening_C' 
    
    C % COARSENING_2  =  1
    C % COARSENING_3  =  2

    call C % FieldSetForm % Initialize &
           ( G % Atlas, &
             FieldOption = [ 'Coarsening_2', 'Coarsening_3' ], &
             NameOption = 'Coarsening', &
             DeviceMemoryOption = G % DeviceMemory, &
             PinnedMemoryOption = G % PinnedMemory, &
             DevicesCommunicateOption = G % DevicesCommunicate, &
             nFieldsOption = 2 )

    select type ( A  =>  G % Atlas )
    class is ( Atlas_SCG_CC_Form )
      call SetCoarseningPolar &
             ( C   = A % Chart_GS_CC, &
               R   = G % Storage_GS % Value ( :, G % CENTER_U_1 ), &
               C_2 = C % Storage_GS % Value ( :, C % COARSENING_2 ) )
    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Coarsening_C_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Initialize_C', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

  end subroutine Initialize_C


  subroutine SetCoarseningPolar ( C, R, C_2 )
    
    class ( Chart_GS_C_Form ), intent ( in ) :: &
      C
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      R
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      C_2

    integer ( KDI ) :: &
      iV
    real ( KDR ) :: &
      dTheta

    dTheta  =  CONSTANT % PI  /  C % nCellsPolar

    do iV = 1, size ( C_2 )
      if ( .not. C % ProperCell ( iV ) ) &
        cycle
      C_2 ( iV )  =  1.0_KDR
      Coarsen_2: do
        if ( C_2 ( iV )  *  R ( iV )  *  dTheta  >  C % MinWidth ) &
          exit Coarsen_2
        C_2 ( iV )  =  2.0_KDR  *  C_2 ( iV )
      end do Coarsen_2
    end do !-- iV

  end subroutine SetCoarseningPolar


  impure elemental subroutine Finalize ( C )

    type ( Coarsening_C_Form ), intent ( inout ) :: &
      C

  end subroutine Finalize

  
end module Coarsening_C__Form
