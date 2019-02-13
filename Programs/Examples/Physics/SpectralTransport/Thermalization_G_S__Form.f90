module Thermalization_G_S__Form

  !-- Thermalization_Generic_Grey_Form

  use GenASiS
  use Interactions_G_S_BSLL_ASC_CSLD__Form

  implicit none
  private

  type, public, extends ( UniverseTemplate ) :: Thermalization_G_S_Form
  contains
    procedure, public, pass :: &
      Initialize
  end type Thermalization_G_S_Form

    real ( KDR ), private :: &
      TemperatureMin, &
      TemperatureMax, &
      OpacityAbsorption

contains


  subroutine Initialize ( T, Name )

    class ( Thermalization_G_S_Form ), intent ( inout ) :: &
      T
    character ( * ), intent ( in )  :: &
      Name


    if ( T % Type == '' ) &
      T % Type = 'a Thermalization_G_S'

!    Thermalization => T

    call T % InitializeTemplate ( Name )

    TemperatureMin  =   0.1_KDR
    TemperatureMax  =  10.0_KDR
    call PROGRAM_HEADER % GetParameter ( TemperatureMin, 'TemperatureMin' )
    call PROGRAM_HEADER % GetParameter ( TemperatureMax, 'TemperatureMax' )

    OpacityAbsorption = 1.0_KDR
    call PROGRAM_HEADER % GetParameter &
           ( OpacityAbsorption, 'OpacityAbsorption' )


    !-- Integrator

    allocate ( SpectralRadiationBoxForm :: T % Integrator )
    select type ( SRB => T % Integrator )
    type is ( SpectralRadiationBoxForm )
    call SRB % Initialize &
           ( RadiationName = [ 'Radiation' ], RadiationType = [ 'GENERIC' ], &
             Name = Name, ApplyStreamingOption = .false. )
!     SRB % SetReference => SetReference

    select type ( PS => SRB % PositionSpace )
    class is ( Atlas_SC_Form )

    select type ( MS => SRB % MomentumSpace )
    class is ( Bundle_SLL_ASC_CSLD_Form )

    select type ( RMB => SRB % Current_BSLL_ASC_CSLD_1D ( 1 ) % Element )
    class is ( RadiationMoments_BSLL_ASC_CSLD_Form )


    !-- Interactions

    allocate ( Interactions_G_S_BSLL_ASC_CSLD_Form :: &
                 SRB % Interactions_BSLL_ASC_CSLD )
    select type ( IB => SRB % Interactions_BSLL_ASC_CSLD )
    class is ( Interactions_G_S_BSLL_ASC_CSLD_Form )
    call IB % Initialize ( MS, 'GENERIC_SPECTRAL' )
    call IB % Set ( OpacityAbsorption = OpacityAbsorption )
    call RMB % SetInteractions ( IB )
    end select !-- IB


    !-- Cleanup

    end select !-- RMB
    end select !-- MS
    end select !-- PS
    end select !-- SRB

  end subroutine Initialize


end module Thermalization_G_S__Form
