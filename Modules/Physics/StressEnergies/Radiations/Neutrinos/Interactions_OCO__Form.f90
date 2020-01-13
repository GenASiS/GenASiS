module Interactions_OCO__Form

  !-- Interactions_OConnorOtt__Form

  use NULIBTABLE, only: &
    NULIBTABLE_NUMBER_GROUPS, &
    NULIBTABLE_NUMBER_EASVARIABLES, &
    NULIBTABLE_LOGRHO_MIN
  use Basics
  use Mathematics
  use StressEnergyBasics
  use Fluids
  use RadiationBasics
  use NeutrinoMoments_S__Form

  implicit none
  private

  type, public, extends ( InteractionsTemplate ) :: Interactions_OCO_Form
    real ( KDR ), dimension ( :, : ), private, allocatable :: &
      EmissionAbsorptionScattering  
      !-- 1st dimension: Energy
      !-- 2nd dimension: Variable 1-3
      !     1: Emissivity * bin widths (ergs/cm^3/s/srad)
      !     2: Absorption opacity (cm^-1)  !-- Includes "stimulated absorption"
      !     3: Scattering opacity (cm^-1)
    logical ( KDL ), private :: &
      Include_NES, &
      IncludePairs
  contains
    procedure, private, pass :: &
      InitializeAllocate_I
    procedure, private, pass :: &
      Set_S
    generic, public :: &
      Set => Set_S
    procedure, public, pass :: &
      Compute
    procedure, public, pass :: &
      ComputeTimeScale
    final :: &
      Finalize
    procedure, private, pass :: &
      ComputeTimeScaleKernel_S
  end type Interactions_OCO_Form

    real ( KDR ), private, protected :: &
      MassDensity_CGS, &
      MeV, &
      cm

contains


  subroutine InitializeAllocate_I &
               ( I, MomentsType, Units, nValues, VariableOption, NameOption, &
                 ClearOption, UnitOption )

    class ( Interactions_OCO_Form ), intent ( inout ) :: &
      I
    character ( * ), intent ( in ) :: &
      MomentsType
    class ( StressEnergyUnitsForm ), intent ( in ) :: &
      Units
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ClearOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption

    if ( I % Type == '' ) &
      I % Type = 'an Interactions_OCO'

    call I % InitializeTemplate &
           ( MomentsType, Units, nValues, VariableOption, NameOption, &
             ClearOption, UnitOption )

  end subroutine InitializeAllocate_I


  subroutine Set_S ( I, Fluid, Energy, d3_Energy, Include_NES, IncludePairs, &
                     iBaseCell )

    class ( Interactions_OCO_Form ), intent ( inout ) :: &
      I
    class ( Fluid_P_Template ), intent ( in ), target :: &
      Fluid
    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      Energy, &
      d3_Energy
    logical ( KDL ), intent ( in ) :: &
      Include_NES, &
      IncludePairs
    integer ( KDI ), intent ( in ) :: &
      iBaseCell

    I % iBaseCell     =   iBaseCell
    I % Energy        =>  Energy
    I % d3_Energy     =>  d3_Energy
    I % Fluid         =>  Fluid

    I % Include_NES   =  Include_NES
    I % IncludePairs  =  IncludePairs

    allocate ( I % EmissionAbsorptionScattering &
                 ( NULIBTABLE_NUMBER_GROUPS, &
                   NULIBTABLE_NUMBER_EASVARIABLES ) )

    MassDensity_CGS  =  UNIT % MASS_DENSITY_CGS
    MeV              =  UNIT % MEGA_ELECTRON_VOLT
    cm               =  UNIT % CENTIMETER

  end subroutine Set_S


  subroutine Compute ( I, R )

    class ( Interactions_OCO_Form ), intent ( inout ) :: &
      I
    class ( CurrentTemplate ), intent ( in ) :: &
      R

    associate &
      (   F => I % Fluid, &
        iBC => I % iBaseCell )
    select type ( R )
    class is ( RadiationMomentsForm )

    select case ( trim ( I % MomentsType ) )
    case ( 'GREY' )
      ! select type ( R )
      ! class is ( PhotonMoments_G_Form )
      !   call I % ComputeKernel_G &
      !          ( R % Value ( :, R % TEMPERATURE_PARAMETER ), &
      !            F % Value ( :, F % BARYON_MASS ), &
      !            F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
      !            F % Value ( :, F % TEMPERATURE ), &
      !            I % Value ( :, I % EMISSIVITY_J ), &
      !            I % Value ( :, I % OPACITY_J ), &
      !            I % Value ( :, I % OPACITY_H ), &
      !            I % Value ( :, I % EQUILIBRIUM_J ) )
      ! end select !-- R
      call Show ( 'GREY not implemented', CONSOLE % ERROR )
      call Show ( 'Interactions_OCO_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    case ( 'SPECTRAL' )

      ! call SetPlanckSpectrum &
      !        ( I % Energy, &
      !          F % Value ( iBC, F % TEMPERATURE ), &
      !          I % Value ( :, I % EQUILIBRIUM_J ) )
      ! call I % ComputeKernel_S &
      !        ( I % Value ( :, I % EQUILIBRIUM_J ), &
      !          F % Value ( iBC, F % BARYON_MASS ), &
      !          F % Value ( iBC, F % COMOVING_BARYON_DENSITY ), &
      !          F % Value ( iBC, F % TEMPERATURE ), &
      !          I % Value ( :, I % EMISSIVITY_J ), &
      !          I % Value ( :, I % OPACITY_J ), &
      !          I % Value ( :, I % OPACITY_H ) )

      select case ( trim ( R % RadiationType ) )
      case ( 'NEUTRINOS_E' )

      case ( 'NEUTRINOS_E_BAR' )

      case ( 'NEUTRINOS_X' )

      case default
        call Show ( 'Radiation Type not recognized', CONSOLE % ERROR )
        call Show ( R % RadiationType, 'RadiationType', CONSOLE % ERROR )
        call Show ( 'Interactions_OCO_Form', 'module', CONSOLE % ERROR )
        call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- R % RadiationType

    end select !-- MomentsType

    end select !-- R
    end associate !-- F, etc.

  end subroutine Compute


  subroutine ComputeTimeScale ( I, R )

    class ( Interactions_OCO_Form ), intent ( inout ) :: &
      I
    class ( CurrentTemplate ), intent ( in ) :: &
      R

    integer ( KDI ) :: &
      iSpecies_OCO

    associate ( iBC => I % iBaseCell )
    select type ( F => I % Fluid )
    class is ( Fluid_P_HN_Form )
    select type ( SF => F % Sources )
    class is ( Sources_F_Form )
    select type ( R )
    class is ( RadiationMomentsForm )

    select case ( trim ( I % MomentsType ) )
    case ( 'GREY' )
      ! select type ( R )
      ! class is ( PhotonMoments_G_Form )
      !   call I % ComputeTimeScaleKernel_G &
      !          (  R % Value ( :,  R % TEMPERATURE_PARAMETER ), &
      !             F % Value ( :,  F % BARYON_MASS ), &
      !             F % Value ( :,  F % COMOVING_BARYON_DENSITY ), &
      !             F % Value ( :,  F % INTERNAL_ENERGY ), &
      !             F % Value ( :,  F % TEMPERATURE ), &
      !             R % Value ( :,  R % COMOVING_ENERGY ), &
      !            SF % Value ( :, SF % RADIATION_TIME ) )
      ! end select !-- R
      call Show ( 'GREY not implemented', CONSOLE % ERROR )
      call Show ( 'Interactions_OCO_Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeTimeScale', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    case ( 'SPECTRAL' )

      select type ( R )
      class is ( NeutrinoMoments_S_Form )

      !   call I % ComputeTimeScaleKernel_S &
      !          (  I % Value ( :, I % EQUILIBRIUM_J ), &
      !             R % Value ( :, R % COMOVING_ENERGY ), &
      !             I % d3_Energy, &
      !             F % Value ( iBC,  F % BARYON_MASS ), &
      !             F % Value ( iBC,  F % COMOVING_BARYON_DENSITY ), &
      !             F % Value ( iBC,  F % INTERNAL_ENERGY ), &
      !             F % Value ( iBC,  F % TEMPERATURE ), &
      !            SF % Value ( iBC, SF % RADIATION_TIME ) )

      select case ( trim ( R % RadiationType ) )
      case ( 'NEUTRINOS_E' )

        iSpecies_OCO  =  1  

        call SetFermiDiracSpectrum &
               ( I % Energy, &
                 F % Value ( iBC, F % TEMPERATURE ), &
                 F % Value ( iBC, F % CHEMICAL_POTENTIAL_E ) &
                   -  F % Value ( iBC, F % CHEMICAL_POTENTIAL_N_P ), &
                 I % Value ( :, I % EQUILIBRIUM_J ) )

      case ( 'NEUTRINOS_E_BAR' )

        iSpecies_OCO  =  2  

        call SetFermiDiracSpectrum &
               ( I % Energy, &
                 F % Value ( iBC, F % TEMPERATURE ), &
                 F % Value ( iBC, F % CHEMICAL_POTENTIAL_N_P ) &
                   -  F % Value ( iBC, F % CHEMICAL_POTENTIAL_E ), &
                 I % Value ( :, I % EQUILIBRIUM_J ) )

      case default
        call Show ( 'Radiation Type not recognized', CONSOLE % ERROR )
        call Show ( R % RadiationType, 'RadiationType', CONSOLE % ERROR )
        call Show ( 'Interactions_OCO_Form', 'module', CONSOLE % ERROR )
        call Show ( 'ComputeTimeScale', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- R % RadiationType

      call I % ComputeTimeScaleKernel_S &
             (  I % Value ( :, I % EQUILIBRIUM_J ), &
                R % Value ( :, R % COMOVING_ENERGY ), &
                I % d3_Energy, &
                F % Value ( iBC,  F % BARYON_MASS ), &
                F % Value ( iBC,  F % COMOVING_BARYON_DENSITY ), &
                F % Value ( iBC,  F % INTERNAL_ENERGY ), &
                F % Value ( iBC,  F % TEMPERATURE ), &
                F % Value ( iBC,  F % ELECTRON_FRACTION ), &
                iS  =  iSpecies_OCO, &
                RT  =  SF % Value ( iBC, SF % RADIATION_TIME ) ) 
             
      end select !-- R

    end select !-- MomentsType

    end select !-- R
    end select !-- SF
    end select !-- F
    end associate !-- iBC

  end subroutine ComputeTimeScale


  impure elemental subroutine Finalize ( I )

    type ( Interactions_OCO_Form ), intent ( inout ) :: &
      I

    if ( allocated ( I % EmissionAbsorptionScattering ) ) &
      deallocate ( I % EmissionAbsorptionScattering )

  end subroutine Finalize


  subroutine ComputeTimeScaleKernel_S ( I, J_EQ, J, dV, M, N, U, T, Y, iS, RT )

    class ( Interactions_OCO_Form ), intent ( in ) :: &
      I
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      J_EQ, &
      J, &
      dV
    real ( KDR ), intent ( in ) :: &
      M, &
      N, &
      U, &
      T, &
      Y
    integer ( KDI ), intent ( in ) :: &
      iS
    real ( KDR ), intent ( out ) :: &
      RT

    real ( KDR ) :: &
      Rho_CGS, &
      T_MeV, &
      Q, &
      SqrtTiny
    real ( KDR ), dimension ( size ( J ) ) :: &
      Kappa_Star  !-- includes stimulated absoprtion, Burrows et al. (2006)

    Rho_CGS  =  M * N / MassDensity_CGS
    T_MeV    =  T / MeV

    SqrtTiny = sqrt ( tiny ( 0.0_KDR ) )

    if ( Rho_CGS  >  10.0_KDR ** NULIBTABLE_LOGRHO_MIN ) then

      call NULIBTABLE_SINGLE_SPECIES_RANGE_ENERGY &
             ( Rho_CGS, T_MeV, Y, iS, I % EmissionAbsorptionScattering, &
               NULIBTABLE_NUMBER_GROUPS, NULIBTABLE_NUMBER_EASVARIABLES )

      Kappa_Star  =  I % EmissionAbsorptionScattering ( :, 2 )  *  cm ** (-1)

      Q   =  abs ( sum ( Kappa_Star * ( J_EQ - J ) * dV ) )
      RT  =  U / max ( Q, SqrtTiny )

    else      
      RT  =  huge ( 1.0_KDR )
    end if

  end subroutine ComputeTimeScaleKernel_S


end module Interactions_OCO__Form
