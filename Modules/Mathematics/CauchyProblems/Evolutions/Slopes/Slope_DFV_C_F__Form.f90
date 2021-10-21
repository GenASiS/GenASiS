module Slope_DFV_C_F__Form

  !-- Slope_DivergenceFiniteVolume_Connection_Flat__Form

  use Basics
  use Manifolds
  use Fields
  use Slope_H__Form

  implicit none
  private

  type, public, extends ( Slope_H_Form ) :: Slope_DFV_C_F_Form
    type ( FieldSetForm ), allocatable :: &
      Stress
    class ( CurrentSetForm ), pointer :: &
      CurrentSet => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_C_F
    generic, public :: &
      Initialize => InitializeAllocate_C_F
    final :: &
      Finalize
  end type Slope_DFV_C_F_Form


contains


  subroutine InitializeAllocate_C_F ( S, CS, SuffixOption, IgnorabilityOption )

    class ( Slope_DFV_C_F_Form ), intent ( inout ) :: &
      S
    class ( CurrentSetForm ), intent ( in ), target :: &
      CS
    character ( * ), intent ( in ), optional :: &
      SuffixOption    
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Slope_DFV_C_F' 
    
    if ( S % TimerName  ==  '' ) &
      S % TimerName  =  'S_DFV_C_F_' // trim ( CS % Name )

    Name  =  'S_DFV_C_F_' // trim ( CS % Name )
    if ( present ( SuffixOption ) ) &
      Name  =  trim ( Name ) // '_' // trim ( SuffixOption )

    S % CurrentSet  =>  CS

    call S % Slope_H_Form % Initialize &
           ( CS % Atlas, &
             FieldOption = CS % Balanced, &
             NameOption = Name, &
             DeviceMemoryOption = CS % DeviceMemory, &
             PinnedMemoryOption = CS % PinnedMemory, &
             DevicesCommunicateOption = CS % DevicesCommunicate, &
             nFieldsOption = CS % nBalanced, &
             IgnorabilityOption = IgnorabilityOption )

    allocate ( S % Stress )
    associate ( SS  =>  S % Stress )
    call SS % Initialize &
           ( CS % Atlas, &
             FieldOption = [ 'Stress_22', 'Stress_33' ], &
             NameOption = 'Stress', &
             DeviceMemoryOption = CS % DeviceMemory, &
             nFieldsOption = 2, &
             IgnorabilityOption = CS % IGNORABILITY + 1 )
    end associate !-- SSS

  end subroutine InitializeAllocate_C_F


  impure elemental subroutine Finalize ( S )

    type ( Slope_DFV_C_F_Form ), intent ( inout ) :: &
      S

    nullify ( S % CurrentSet )
    
    if ( allocated ( S % Stress ) ) &
      deallocate ( S % Stress )

  end subroutine Finalize



end module Slope_DFV_C_F__Form
