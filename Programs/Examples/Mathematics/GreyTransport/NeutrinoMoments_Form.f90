module NeutrinoMoments_Form

  use Basics
  use RadiationMoments_Form

  implicit none
  private

  type, public, extends ( RadiationMomentsForm ) :: NeutrinoMomentsForm
  contains
    procedure, public, pass :: &
      InitializeAllocate_NM
    procedure, public, pass ( RM ) :: &
      ComputeSpectralParameters      
  end type NeutrinoMomentsForm

    private :: &
      ComputeSpectralParametersKernel

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


  subroutine ComputeSpectralParameters &
               ( T, Eta, E_Ave, F_Ave, J_EQ, RM, J, N, T_EQ, Eta_EQ )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      T, &
      Eta, &
      E_Ave, &
      F_Ave, &
      J_EQ
    class ( NeutrinoMomentsForm ), intent ( in ) :: &
      RM
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      J, &
      N, &
      T_EQ, &
      Eta_EQ

    call ComputeSpectralParametersKernel &
           ( T, Eta, E_Ave, F_Ave, J_EQ, J, N, T_EQ, Eta_EQ )

  end subroutine ComputeSpectralParameters


  subroutine ComputeSpectralParametersKernel &
               ( T, Eta, E_Ave, F_Ave, J_EQ, J, N, T_EQ, Eta_EQ )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      T, &
      Eta, &
      E_Ave, &
      F_Ave, &
      J_EQ
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      J, &
      N, &
      T_EQ, &
      Eta_EQ

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      k, &
      Factor_Eta, &
      Factor_T, &
      Fermi_3, &
      Fermi_3_EQ, &
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

    Factor_Eta =  FourPi ** ( 1.0_KDR / 3.0_KDR ) / ( TwoPi * hBar * c )
    Factor_T   =  FourPi  *  k ** 4  /  ( TwoPi * hBar * c ) ** 3

    end associate !-- k, etc.

    !$OMP parallel do private ( iV )
    do iV = 1, nValues

      if ( N ( iV ) > 0.0_KDR ) then
        Eta ( iV )  &
          =  - 3 * log ( Factor_Eta * ( J ( iV ) / 6.0_KDR ) &
                         * ( 2.0_KDR / N ( iV ) ) ** ( 4.0_KDR / 3.0_KDR ) )
      end if

      call DFERMI ( 3.0_KDR, Eta ( iV ), 0.0_KDR, Fermi_3, &
                    fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
      call DFERMI ( 3.0_KDR, Eta_EQ ( iV ), 0.0_KDR, Fermi_3_EQ, &
                    fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )

      T  ( iV )     =  ( J ( iV )  /  ( Factor_T * Fermi_3 ) ) ** ( 0.25_KDR )

      E_Ave ( iV )  =  J ( iV )  /  max ( N ( iV ), tiny ( 0.0_KDR ) )

      J_EQ ( iV )   =  Factor_T  *  T_EQ ( iV ) ** 4  *  Fermi_3_EQ

    end do !-- iV
    !$OMP end parallel do

!call Show ( J, '>>> J' )
!call Show ( N, '>>> N' )
!call Show ( E_Ave, '>>> E_Ave' )

    ! !$OMP parallel do private ( iV )
    ! do iV = 1, nValues
    !   if ( T ( iV ) > 0.0_KDR ) then
    !     F_Ave ( iV )  &
    !       =  1.0_KDR &
    !          / ( exp ( E_Ave ( iV ) / ( k * T ( iV ) )  -  Eta ( iV ) ) &
    !              + 1.0_KDR )
    !   else
    !     F_Ave ( iV ) = 0.0_KDR
    !   end if
    ! end do !-- iV
    ! !$OMP end parallel do

  end subroutine ComputeSpectralParametersKernel


end module NeutrinoMoments_Form
