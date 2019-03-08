module PhotonMoments_CSL__Form

  !-- PhotonMoments_ChartSingleLevel__Form

  use Basics
  use Mathematics
  use RadiationBasics
  use PhotonMoments_G__Form
  use PhotonMoments_S__Form

  implicit none
  private

  type, public, extends ( RadiationMoments_CSL_Form ) :: PhotonMoments_CSL_Form
  contains
    procedure, private, pass :: &
      SetField
  end type PhotonMoments_CSL_Form

contains


  subroutine SetField ( FC )

    class ( PhotonMoments_CSL_Form ), intent ( inout ) :: &
      FC

    allocate ( FC % FieldOutput )

    select case ( trim ( FC % RadiationMomentsType ) )
    case ( 'PHOTONS_GREY' )
      allocate ( PhotonMoments_G_Form :: FC % Field )
      select type ( PM => FC % Field )
      type is ( PhotonMoments_G_Form )
        call PM % Initialize &
               ( FC % RiemannSolverType, FC % UseLimiter, &
                 FC % Velocity_U_Unit, FC % MomentumDensity_U_Unit, &
                 FC % MomentumDensity_D_Unit, FC % EnergyDensityUnit, &
                 FC % TemperatureUnit, FC % LimiterParameter, FC % nValues, &
                 NameOption = FC % NameShort )
        call PM % SetPrimitiveConserved ( )
        call PM % SetOutput ( FC % FieldOutput )
      end select !-- PM
    case ( 'PHOTONS_SPECTRAL' )
      allocate ( PhotonMoments_S_Form :: FC % Field )
      select type ( PM => FC % Field )
      type is ( PhotonMoments_S_Form )
        call PM % Initialize &
               ( FC % RiemannSolverType, FC % UseLimiter, &
                 FC % Velocity_U_Unit, FC % MomentumDensity_U_Unit, &
                 FC % MomentumDensity_D_Unit, FC % EnergyDensityUnit, &
                 FC % LimiterParameter, FC % nValues, &
                 NameOption = FC % NameShort )
        call PM % SetPrimitiveConserved ( )
        call PM % SetOutput ( FC % FieldOutput )
      end select !-- PM
    case default
      call Show ( 'RadiationMomentsType not recognized', CONSOLE % ERROR )
      call Show ( FC % RadiationMomentsType, 'RadiationMomentsType', &
                  CONSOLE % ERROR )
      call Show ( 'PhotonMoments_CSL__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetField', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- RadiationMomentsType

  end subroutine SetField


end module PhotonMoments_CSL__Form
