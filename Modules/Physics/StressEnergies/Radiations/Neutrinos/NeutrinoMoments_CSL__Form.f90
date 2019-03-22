module NeutrinoMoments_CSL__Form

  !-- NeutrinoMoments_ChartSingleLevel__Form

  use Basics
  use Mathematics
  use StressEnergyBasics
  use RadiationBasics
  use NeutrinoMoments_G__Form
  use NeutrinoMoments_S__Form

  implicit none
  private

  type, public, extends ( RadiationMoments_CSL_Form ) :: &
    NeutrinoMoments_CSL_Form
  contains
    procedure, public, pass :: &
      Initialize
    procedure, private, pass :: &
      SetField
  end type NeutrinoMoments_CSL_Form

contains


  subroutine Initialize &
               ( RMC, C, NameShort, RadiationMomentsType, RiemannSolverType, &
                 ReconstructedType, UseLimiter, Units, LimiterParameter, &
                 nValues, IgnorabilityOption )

    class ( NeutrinoMoments_CSL_Form ), intent ( inout ) :: &
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
      RMC % Type = 'a NeutrinoMoments_CSL'

    call RMC % RadiationMoments_CSL_Form % Initialize &
           ( C, NameShort, RadiationMomentsType, RiemannSolverType, &
             ReconstructedType, UseLimiter, Units, LimiterParameter, &
             nValues, IgnorabilityOption )

  end subroutine Initialize


  subroutine SetField ( FC )

    class ( NeutrinoMoments_CSL_Form ), intent ( inout ) :: &
      FC

    allocate ( FC % FieldOutput )

    select case ( trim ( FC % RadiationMomentsType ) )
    case ( 'NEUTRINOS_ELECTRON_GREY', 'ANTINEUTRINOS_ELECTRON_GREY', &
           'NEUTRINOS_X_GREY' )
      allocate ( NeutrinoMoments_G_Form :: FC % Field )
      select type ( NM => FC % Field )
      type is ( NeutrinoMoments_G_Form )
        call NM % Initialize &
               ( FC % RadiationMomentsType, FC % RiemannSolverType, &
                 FC % ReconstructedType, FC % UseLimiter, FC % Units, &
                 FC % LimiterParameter, FC % nValues, &
                 NameOption = FC % NameShort )
        call NM % SetPrimitiveConserved ( )
        call NM % SetReconstructed ( )
        call NM % SetOutput ( FC % FieldOutput )
      end select !-- NM
    case ( 'NEUTRINOS_ELECTRON_SPECTRAL', 'ANTINEUTRINOS_ELECTRON_SPECTRAL', &
           'NEUTRINOS_X_SPECTRAL' )
      allocate ( NeutrinoMoments_S_Form :: FC % Field )
      select type ( NM => FC % Field )
      type is ( NeutrinoMoments_S_Form )
        call NM % Initialize &
               ( FC % RadiationMomentsType, FC % RiemannSolverType, &
                 FC % ReconstructedType, FC % UseLimiter, FC % Units, &
                 FC % LimiterParameter, FC % nValues, &
                 NameOption = FC % NameShort )
        call NM % SetPrimitiveConserved ( )
        call NM % SetReconstructed ( )
        call NM % SetOutput ( FC % FieldOutput )
      end select !-- NM
    case default
      call Show ( 'RadiationMomentsType not recognized', CONSOLE % ERROR )
      call Show ( FC % RadiationMomentsType, 'RadiationMomentsType', &
                  CONSOLE % ERROR )
      call Show ( 'NeutrinoMoments_CSL__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetField', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- RadiationMomentsType

  end subroutine SetField


end module NeutrinoMoments_CSL__Form
