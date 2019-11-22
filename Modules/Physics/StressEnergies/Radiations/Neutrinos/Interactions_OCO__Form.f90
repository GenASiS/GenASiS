module Interactions_OCO__Form

  !-- Interactions_OConnorOtt__Form

  use Basics
  use Mathematics
  use StressEnergyBasics
  use Fluids
  use RadiationBasics
  use NeutrinoMoments_S__Form

  implicit none
  private

  type, public, extends ( InteractionsTemplate ) :: Interactions_OCO_Form
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
  end type Interactions_OCO_Form

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

        call SetFermiDiracSpectrum &
               ( I % Energy, &
                 F % Value ( iBC, F % TEMPERATURE ), &
                 F % Value ( iBC, F % CHEMICAL_POTENTIAL_E ) &
                   -  F % Value ( iBC, F % CHEMICAL_POTENTIAL_N_P ), &
                 I % Value ( :, I % EQUILIBRIUM_J ) )

      case ( 'NEUTRINOS_E_BAR' )

        call SetFermiDiracSpectrum &
               ( I % Energy, &
                 F % Value ( iBC, F % TEMPERATURE ), &
                 F % Value ( iBC, F % CHEMICAL_POTENTIAL_N_P ) &
                   -  F % Value ( iBC, F % CHEMICAL_POTENTIAL_E ), &
                 I % Value ( :, I % EQUILIBRIUM_J ) )

      case ( 'NEUTRINOS_X' )

      case default
        call Show ( 'Radiation Type not recognized', CONSOLE % ERROR )
        call Show ( R % RadiationType, 'RadiationType', CONSOLE % ERROR )
        call Show ( 'Interactions_OCO_Form', 'module', CONSOLE % ERROR )
        call Show ( 'ComputeTimeScale', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- R % RadiationType

      end select !-- R

    end select !-- MomentsType

    end select !-- R
    end select !-- SF
    end select !-- F
    end associate !-- iBC

  end subroutine ComputeTimeScale


end module Interactions_OCO__Form
