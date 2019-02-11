module Thermalization_G_G__Form

  !-- Thermalization_Generic_Grey_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( UniverseTemplate ) :: Thermalization_G_G_Form
  contains
    procedure, public, pass :: &
      Initialize
  end type Thermalization_G_G_Form

    real ( KDR ), private :: &
      TemperatureMin, &
      TemperatureMax, &
      OpacityAbsorption

contains


  subroutine Initialize ( T, Name )

    class ( Thermalization_G_G_Form ), intent ( inout ) :: &
      T
    character ( * ), intent ( in )  :: &
      Name

    if ( T % Type == '' ) &
      T % Type = 'a Thermalization_G_G'

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

    allocate ( GreyRadiationBoxForm :: T % Integrator )
    select type ( GRB => T % Integrator )
    type is ( GreyRadiationBoxForm )
    call GRB % Initialize &
           ( RadiationName = [ 'Radiation' ], RadiationType = [ 'GENERIC' ], &
             Name = Name )
!    GRB % SetReference => SetReference

    select type ( PS => GRB % PositionSpace )
    class is ( Atlas_SC_Form )

    select type ( RMA => GRB % Current_ASC_1D ( 1 ) % Element )
    class is ( RadiationMoments_ASC_Form )


    !-- Cleanup

    end select !-- RMA
    end select !-- PS
    end select !-- GRB

  end subroutine Initialize


end module Thermalization_G_G__Form
