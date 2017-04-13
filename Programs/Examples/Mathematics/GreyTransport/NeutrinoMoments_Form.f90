module NeutrinoMoments_Form

  use Basics
  use RadiationMoments_Form

  implicit none
  private

  type, public, extends ( RadiationMomentsForm ) :: NeutrinoMomentsForm
    real ( KDR ) :: &
      DegeneracyParameter_Infty = 0.0_KDR
  contains
    procedure, public, pass :: &
      InitializeAllocate_NM
    procedure, public, pass :: &
      SetDegeneracyParameter_Infty
    procedure, public, pass ( RM ) :: &
      ComputeSpectralParameters      
  end type NeutrinoMomentsForm

contains


  subroutine InitializeAllocate_NM &
               ( NM, NeutrinoType, RiemannSolverType, UseLimiter, &
                 Velocity_U_Unit, MomentumDensity_U_Unit, &
                 MomentumDensity_D_Unit, EnergyDensityUnit, TemperatureUnit, &
                 LimiterParameter, nValues, VariableOption, VectorOption, &
                 NameOption, ClearOption, UnitOption, VectorIndicesOption )

    class ( NeutrinoMomentsForm ), intent ( inout ) :: &
      NM
    character ( * ), intent ( in ) :: &
      NeutrinoType, &
      RiemannSolverType
    logical ( KDL ), intent ( in ) :: &
      UseLimiter
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      Velocity_U_Unit, &
      MomentumDensity_U_Unit, &
      MomentumDensity_D_Unit
    type ( MeasuredValueForm ), intent ( in ) :: &
      EnergyDensityUnit, &
      TemperatureUnit
    real ( KDR ), intent ( in ) :: &
      LimiterParameter
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ClearOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption

    NM % Type = NeutrinoType

    call NM % RadiationMomentsForm % Initialize &
           ( RiemannSolverType, UseLimiter, Velocity_U_Unit, &
             MomentumDensity_U_Unit, MomentumDensity_D_Unit, &
             EnergyDensityUnit, TemperatureUnit, LimiterParameter, &
             nValues, VariableOption, VectorOption, NameOption, &
             ClearOption, UnitOption, VectorIndicesOption )

  end subroutine InitializeAllocate_NM


  subroutine SetDegeneracyParameter_Infty ( F, DegeneracyParameter_Infty )

    class ( NeutrinoMomentsForm ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      DegeneracyParameter_Infty

    F % DegeneracyParameter_Infty = DegeneracyParameter_Infty

  end subroutine SetDegeneracyParameter_Infty


  subroutine ComputeSpectralParameters ( T, Eta, RM, J, FF )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      T, &
      Eta
    class ( NeutrinoMomentsForm ), intent ( in ) :: &
      RM
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      J, &
      FF

    call ComputeSpectralParametersKernel &
           ( T, Eta, J, FF, RM % Value ( :, RM % DEGENERACY_PARAMETER_EQ ), &
             RM % DegeneracyParameter_Infty )

  end subroutine ComputeSpectralParameters


  subroutine ComputeSpectralParametersKernel &
               ( T, Eta, J, FF, Eta_EQ, Eta_Infty )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      T, &
      Eta
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      J, &
      FF, &
      Eta_EQ
    real ( KDR ), intent ( in ) :: &
      Eta_Infty

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      Factor, &
      Fermi_3, &
      fdeta, fdeta2, &
      fdtheta, fdtheta2, &
      fdetadtheta

    nValues = size ( T )

    associate &
      ( k      => CONSTANT % BOLTZMANN, &
        c      => CONSTANT % SPEED_OF_LIGHT, &
        hBar   => CONSTANT % PLANCK_REDUCED, &
        FourPi => 4.0_KDR * CONSTANT % PI, &
        TwoPi  => 2.0_KDR * CONSTANT % PI )

    Factor  =  FourPi  *  k ** 4  /  ( TwoPi * hBar * c ) ** 3

    end associate !-- k, etc.

    !$OMP parallel do private ( iV )
    do iV = 1, nValues

      Eta ( iV )  =  Eta_EQ ( iV )  *  ( 1.0_KDR  -  FF ( iV ) ) &
                     +   Eta_Infty  *  FF ( iV )

      Fermi_3 = 1.0_KDR
      !-- FIXME: call Compute Fermi_3
      call DFERMI &
           ( 3.0_KDR, Eta ( iV ), 0.0_KDR, Fermi_3, &
            fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta)

      T  ( iV )  =  ( J ( iV )  /  ( Factor * Fermi_3 ) ) ** ( 0.25_KDR )

    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeSpectralParametersKernel


end module NeutrinoMoments_Form
