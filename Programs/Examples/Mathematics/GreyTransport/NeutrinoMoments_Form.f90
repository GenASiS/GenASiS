module NeutrinoMoments_Form

  use Basics
  use RadiationMoments_Form

  implicit none
  private

  type, public, extends ( RadiationMomentsForm ) :: NeutrinoMomentsForm
    real ( KDR ) :: &
      DegeneracyInfinity
  contains
    procedure, public, pass :: &
      InitializeAllocate_NM
    procedure, public, pass ( RM ) :: &
      ComputeSpectralParameters      
  end type NeutrinoMomentsForm

contains


  subroutine InitializeAllocate_NM &
               ( NM, NeutrinoType, RiemannSolverType, UseLimiter, &
                 Velocity_U_Unit, MomentumDensity_U_Unit, &
                 MomentumDensity_D_Unit, EnergyDensityUnit, TemperatureUnit, &
                 DegeneracyInfinity, LimiterParameter, nValues, &
                 VariableOption, VectorOption, NameOption, ClearOption, &
                 UnitOption, VectorIndicesOption )

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
      DegeneracyInfinity, &
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

    NM % DegeneracyInfinity = DegeneracyInfinity

    call NM % RadiationMomentsForm % Initialize &
           ( RiemannSolverType, UseLimiter, Velocity_U_Unit, &
             MomentumDensity_U_Unit, MomentumDensity_D_Unit, &
             EnergyDensityUnit, TemperatureUnit, LimiterParameter, &
             nValues, VariableOption, VectorOption, NameOption, &
             ClearOption, UnitOption, VectorIndicesOption )

  end subroutine InitializeAllocate_NM


  subroutine ComputeSpectralParameters &
               ( T, Eta, E_Ave, J_EQ, RM, J, FF, T_EQ, Eta_EQ )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      T, &
      Eta, &
      E_Ave, &
      J_EQ
    class ( NeutrinoMomentsForm ), intent ( in ) :: &
      RM
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      J, &
      FF, &
      T_EQ, &
      Eta_EQ

    call ComputeSpectralParametersKernel &
           ( T, Eta, E_Ave, J_EQ, J, FF, T_EQ, Eta_EQ, &
             RM % DegeneracyInfinity )

  end subroutine ComputeSpectralParameters


  subroutine ComputeSpectralParametersKernel &
               ( T, Eta, E_Ave, J_EQ, J, FF, T_EQ, Eta_EQ, Eta_Infty )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      T, &
      Eta, &
      E_Ave, &
      J_EQ
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      J, &
      FF, &
      T_EQ, &
      Eta_EQ
    real ( KDR ), intent ( in ) :: &
      Eta_Infty

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      k, &
      Factor, &
      Fermi_2, Fermi_3, &
      fdeta, fdeta2, &
      fdtheta, fdtheta2, &
      fdetadtheta

    nValues = size ( T )

    k = CONSTANT % BOLTZMANN

    associate &
      ( c      => CONSTANT % SPEED_OF_LIGHT, &
        hBar   => CONSTANT % PLANCK_REDUCED, &
        FourPi => 4.0_KDR * CONSTANT % PI, &
        TwoPi  => 2.0_KDR * CONSTANT % PI )

    Factor  =  FourPi  *  k ** 4  /  ( TwoPi * hBar * c ) ** 3

    end associate !-- k, etc.

    !$OMP parallel do private ( iV )
    do iV = 1, nValues

      Eta ( iV )  =  Eta_EQ ( iV )  *  ( 1.0_KDR  -  FF ( iV ) ) &
                     +   Eta_Infty  *  FF ( iV )

      call DFERMI ( 2.0_KDR, Eta ( iV ), 0.0_KDR, Fermi_2, &
                    fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
      call DFERMI ( 3.0_KDR, Eta ( iV ), 0.0_KDR, Fermi_3, &
                    fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )

      T  ( iV )     =  ( J ( iV )  /  ( Factor * Fermi_3 ) ) ** ( 0.25_KDR )

      E_Ave ( iV )  =  ( k * T ( iV ) ) * Fermi_3 / Fermi_2

      J_EQ ( iV )   =  Factor  *  T_EQ ( iV ) ** 4  *  Fermi_3

    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeSpectralParametersKernel


end module NeutrinoMoments_Form
