module Coarsening_C_F__Form

  !-- Coarsening_Central_Fluid_Form

  use Basics
  use Mathematics
  use Fluid_D__Form

  implicit none
  private

  type, public, extends ( Coarsening_C_Form ) :: Coarsening_C_F_Form
    integer ( KDI ) :: &
      nCellsZero
    class ( Fluid_D_Form ), pointer :: &
      Fluid => null ( )
  contains
    procedure, private, pass :: &
      Initialize_C_F
    generic, public :: &
      Initialize => Initialize_C_F
    procedure, public, pass ( C ) :: &
      Compute
    final :: &
      Finalize
  end type Coarsening_C_F_Form


contains


  subroutine Initialize_C_F ( C, F, G, nCellsZero )

    class ( Coarsening_C_F_Form ), intent ( inout ) :: &
      C
    class ( Fluid_D_Form ), intent ( in ), target :: &
      F
    class ( Geometry_F_Form ), intent ( in ) :: &
      G
    integer ( KDI ), intent ( in ) :: &
      nCellsZero

    if ( C % Type == '' ) &
      C % Type  =  'a Coarsening_C_F' 
    
    call C % Coarsening_C_Form % Initialize ( G )

    C % Fluid       =>  F
    C % nCellsZero  =   nCellsZero

  end subroutine Initialize_C_F


  subroutine Compute ( FS, C )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS
    class ( Coarsening_C_F_Form ), intent ( in ) :: &
      C

    integer ( KDI ) :: &
      iMomentum_2, &
      iMomentum_3
    real ( KDR ), dimension ( :, :, :, : ), pointer :: &
      FS_4D

    select type ( A  =>  FS % Atlas )
      class is ( Atlas_SCG_CC_Form )
    associate &
      ( F        =>  C % Fluid, &
        C_GS_CC  =>  A % Chart_GS_CC )

    call Search &
           ( F % iaBalanced, F % MOMENTUM_DENSITY_D_2, iMomentum_2 )
    call Search &
           ( F % iaBalanced, F % MOMENTUM_DENSITY_D_3, iMomentum_3 )

    call C % Coarsening_C_Form % Compute ( FS )

    call C_GS_CC % SetFieldPointer &
           ( FS % Storage_GS % Value, FS_4D )

    if ( C_GS_CC % iaBrick ( 1 )  ==  1 ) then
      FS_4D ( 1, :, :, iMomentum_2 )  =  0.0_KDR
      FS_4D ( 1, :, :, iMomentum_3 )  =  0.0_KDR
    end if

    end associate !-- F, etc.

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Coarsening_C_F_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

  end subroutine Compute


  impure elemental subroutine Finalize ( C )

    type ( Coarsening_C_F_Form ), intent ( inout ) :: &
      C

    nullify ( C % Fluid )

  end subroutine Finalize

  
end module Coarsening_C_F__Form
