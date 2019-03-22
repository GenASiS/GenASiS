module PhotonMoments_CSL__Form

  !-- PhotonMoments_ChartSingleLevel__Form

  use Basics
  use Mathematics
  use StressEnergyBasics
  use RadiationBasics
  use PhotonMoments_G__Form
  use PhotonMoments_S__Form

  implicit none
  private

  type, public, extends ( RadiationMoments_CSL_Form ) :: PhotonMoments_CSL_Form
  contains
    procedure, public, pass :: &
      Initialize
    procedure, private, pass :: &
      SetField
  end type PhotonMoments_CSL_Form

contains


  subroutine Initialize &
               ( RMC, C, NameShort, RadiationMomentsType, RiemannSolverType, &
                 ReconstructedType, UseLimiter, Units, LimiterParameter, &
                 nValues, IgnorabilityOption )

    class ( PhotonMoments_CSL_Form ), intent ( inout ) :: &
      RMC
    class ( ChartTemplate ), intent ( in ) :: &
      C
    character ( * ), intent ( in ) :: &
      NameShort, &
      RadiationMomentsType, &
      RiemannSolverType, &
      ReconstructedType
    logical ( KDL ), intent ( in ) :: &
      UseLimiter
    class ( StressEnergyUnitsForm ), intent ( in ), target :: &
      Units
    real ( KDR ), intent ( in ) :: &
      LimiterParameter
    integer ( KDI ), intent ( in ) :: &
      nValues
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( RMC % Type == '' ) &
      RMC % Type = 'a PhotonMoments_CSL'

    call RMC % RadiationMoments_CSL_Form % Initialize &
           ( C, NameShort, RadiationMomentsType, RiemannSolverType, &
             ReconstructedType, UseLimiter, Units, LimiterParameter, &
             nValues, IgnorabilityOption )

  end subroutine Initialize


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
               ( FC % RadiationMomentsType, FC % RiemannSolverType, &
                 FC % ReconstructedType, FC % UseLimiter, FC % Units, &
                 FC % LimiterParameter, FC % nValues, &
                 NameOption = FC % NameShort )
        call PM % SetPrimitiveConserved ( )
        call PM % SetReconstructed ( )
        call PM % SetOutput ( FC % FieldOutput )
      end select !-- PM
    case ( 'PHOTONS_SPECTRAL' )
      allocate ( PhotonMoments_S_Form :: FC % Field )
      select type ( PM => FC % Field )
      type is ( PhotonMoments_S_Form )
        call PM % Initialize &
               ( FC % RadiationMomentsType, FC % RiemannSolverType, &
                 FC % ReconstructedType, FC % UseLimiter, FC % Units, &
                 FC % LimiterParameter, FC % nValues, &
                 NameOption = FC % NameShort )
        call PM % SetPrimitiveConserved ( )
        call PM % SetReconstructed ( )
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
