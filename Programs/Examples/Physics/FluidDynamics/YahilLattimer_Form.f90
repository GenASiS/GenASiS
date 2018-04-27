module YahilLattimer_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( UniverseTemplate ) :: YahilLattimerForm
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type YahilLattimerForm

    private :: &
      SetFluid!, &
!      SetReference

!      private :: &
!        SetFluidKernel

contains


  subroutine Initialize ( YL, Name )

    class ( YahilLattimerForm ), intent ( inout ) :: &
      YL
    character ( * ), intent ( in ) :: &
      Name

    class ( Fluid_P_I_Form ), pointer :: &
      F

    if ( YL % Type == '' ) &
      YL % Type = 'a YahilLattimer'

    call YL % InitializeTemplate ( Name )


    !-- Integrator

    allocate ( FluidCentralCoreForm :: YL % Integrator )
    select type ( FCC => YL % Integrator )
    type is ( FluidCentralCoreForm )
    call FCC % Initialize &
           ( Name, FluidType = 'IDEAL', &
             GeometryType = 'NEWTONIAN' )
!    FB % SetReference => SetReference

   select type ( PS => FCC % PositionSpace )
   class is ( Atlas_SC_Form )

   select type ( FA => FCC % Current_ASC )
   class is ( Fluid_ASC_Form )


    !-- Initial Conditions

!-- FIXME: Set default parameters, read parameters

    F => FA % Fluid_P_I ( )
    call SetFluid ( YL, F, Time = 0.0_KDR )


    !-- Cleanup

    end select !-- FA
    end select !-- PS
    end select !-- FCC
    nullify ( F )

  end subroutine Initialize


  impure elemental subroutine Finalize ( YL )
    
    type ( YahilLattimerForm ), intent ( inout ) :: &
      YL

    call YL % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetFluid ( YL, F, Time )

    class ( YahilLattimerForm ), intent ( inout ) :: &
      YL
    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      Time

  end subroutine SetFluid


end module YahilLattimer_Form
