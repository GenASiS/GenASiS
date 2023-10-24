module Interactions_C__Form

  !-- Interactions_Constant_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( Interactions_BM_Form ) :: Interactions_C_Form
    real ( KDR ) :: &
      OpacityAbsorption = 0.0_KDR
  contains
    procedure, private, pass :: &
      InitializeAllocate_I
    procedure, public, pass :: &
      SetOpacityAbsorption
    final :: &
      Finalize
  end type Interactions_C_Form


contains


  subroutine InitializeAllocate_I &
               ( I, F, Units_R, FieldOption, NameOption, UnitOption, &
                 nFieldsOption, IgnorabilityOption )

    class ( Interactions_C_Form ), intent ( inout ) :: &
      I
    class ( Fluid_D_Form ), intent ( in ), target :: &
      F
    class ( Units_R_Form ), dimension ( : ), intent ( in ) :: &
      Units_R
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( QuantityForm ), dimension ( :, : ), intent ( in ), optional :: &
      UnitOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption, &
      IgnorabilityOption

    if ( I % Type  ==  '' ) &
      I % Type  =  'an Interactions_C' 
    
    call I % Interactions_BM_Form % Initialize &
           ( F, Units_R, &
             FieldOption = FieldOption, &
             NameOption = NameOption, &
             UnitOption = UnitOption, &
             nFieldsOption = nFieldsOption, &
             IgnorabilityOption = IgnorabilityOption )

  end subroutine InitializeAllocate_I


  subroutine SetOpacityAbsorption ( I, Kappa_A )

    class ( Interactions_C_Form ), intent ( inout ) :: &
      I
    real ( KDR ), intent ( in ) :: &
      Kappa_A

    I % OpacityAbsorption  =  Kappa_A

    call Show ( 'Setting BaryonDensityMin of an Interactions_C', &
                I % IGNORABILITY + 1 )
    call Show ( I % Name, 'Name', &
                I % IGNORABILITY + 1 )
    call Show ( I % OpacityAbsorption, 'OpacityAbsorption', &
                I % IGNORABILITY + 1 )

  end subroutine SetOpacityAbsorption


  impure elemental subroutine Finalize ( I )

    type ( Interactions_C_Form ), intent ( inout ) :: &
      I

  end subroutine Finalize


end module Interactions_C__Form
