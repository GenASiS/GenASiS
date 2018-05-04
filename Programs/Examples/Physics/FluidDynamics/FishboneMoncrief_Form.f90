module FishboneMoncrief_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( UniverseTemplate ) :: FishboneMoncriefForm
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type FishboneMoncriefForm

contains


  subroutine Initialize ( FM, Name )

    class ( FishboneMoncriefForm ), intent ( inout ) :: &
      FM
    character ( * ), intent ( in ) :: &
      Name

    if ( FM % Type == '' ) &
      FM % Type = 'a FishboneMoncrief'

    call FM % InitializeTemplate ( Name )


    !-- Integrator

    allocate ( FluidCentralExcisionForm :: FM % Integrator )
    select type ( FCE => FM % Integrator )
    type is ( FluidCentralExcisionForm )
    call FCE % Initialize &
           ( Name, FluidType = 'IDEAL', &
             GeometryType = 'NEWTONIAN', &
             RadialRatioOption = 3.0_KDR )
    end select !-- FCE

  end subroutine Initialize


  impure elemental subroutine Finalize ( FM )
    
    type ( FishboneMoncriefForm ), intent ( inout ) :: &
      FM

    call FM % FinalizeTemplate ( )

  end subroutine Finalize


end module FishboneMoncrief_Form
