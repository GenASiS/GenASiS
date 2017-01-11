module Fluid_ASC__Form

  !-- Fluid_AtlasSingleChart__Form

  use Basics
  use Mathematics
  use Fluid_D__Form
  use Fluid_P_P__Form
  use Fluid_CSL__Form
  use Tally_F_D__Form
  use Tally_F_P__Form
  
  implicit none
  private
  
  type, public, extends ( Current_ASC_Template ) :: Fluid_ASC_Form
    type ( MeasuredValueForm ) :: &
      MassDensityUnit, &
      EnergyDensityUnit, &
      NumberDensityUnit, &
      TemperatureUnit
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      VelocityUnit
    character ( LDF ) :: &
      FluidType = ''
  contains
    procedure, public, pass :: &
      Initialize
    procedure, private, pass :: &
      Fluid_D_CSL
    generic, public :: &
      Fluid_D => Fluid_D_CSL
    procedure, private, pass :: &
      Fluid_P_P_CSL
    generic, public :: &
      Fluid_P_P => Fluid_P_P_CSL
    final :: &
      Finalize
    procedure, public, pass :: &
      SetField
  end type Fluid_ASC_Form

contains


  subroutine Initialize &
               ( FA, A, FluidType, NameOutputOption, VelocityUnitOption, &
                 MassDensityUnitOption, EnergyDensityUnitOption, &
                 NumberDensityUnitOption, TemperatureUnitOption, &
                 MassUnitOption, EnergyUnitOption, MomentumUnitOption, &
                 AngularMomentumUnitOption )

    class ( Fluid_ASC_Form ), intent ( inout ) :: &
      FA
    class ( Atlas_SC_Form ), intent ( in ) :: &
      A
    character ( * ), intent ( in ) :: &
      FluidType
    character ( * ), intent ( in ), optional :: &
      NameOutputOption
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ), optional :: &
      VelocityUnitOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      MassDensityUnitOption, &
      EnergyDensityUnitOption, &
      NumberDensityUnitOption, &
      TemperatureUnitOption, &
      MassUnitOption, &
      EnergyUnitOption, &
      MomentumUnitOption, &
      AngularMomentumUnitOption

    integer ( KDI ) :: &
      iB  !-- iBoundary

    FA % Type = 'a Fluid_ASC'
    FA % FluidType = FluidType

    if ( present ( MassDensityUnitOption ) ) &
      FA % MassDensityUnit = MassDensityUnitOption
    if ( present ( EnergyDensityUnitOption ) ) &
      FA % EnergyDensityUnit = EnergyDensityUnitOption
    if ( present ( NumberDensityUnitOption ) ) &
      FA % NumberDensityUnit = NumberDensityUnitOption
    if ( present ( TemperatureUnitOption ) ) &
      FA % TemperatureUnit = TemperatureUnitOption
    if ( present ( VelocityUnitOption ) ) &
      FA % VelocityUnit = VelocityUnitOption

    call FA % InitializeTemplate_ASC &
           ( A, NameOutputOption = NameOutputOption )

    if ( .not. allocated ( FA % TallyInterior ) ) then
      select case ( trim ( FluidType ) )
      case ( 'DUST' )
        allocate ( Tally_F_D_Form :: FA % TallyInterior )
        allocate ( Tally_F_D_Form :: FA % TallyTotal )
        allocate ( Tally_F_D_Form :: FA % TallyChange )
        allocate ( FA % TallyBoundaryLocal  ( A % nBoundaries ) )
        allocate ( FA % TallyBoundaryGlobal ( A % nBoundaries ) )
        do iB = 1, A % nBoundaries 
          allocate &
            ( Tally_F_D_Form :: FA % TallyBoundaryLocal  ( iB ) % Element )
          allocate &
            ( Tally_F_D_Form :: FA % TallyBoundaryGlobal ( iB ) % Element )
        end do !-- iB
      case ( 'POLYTROPIC' )
        allocate ( Tally_F_P_Form :: FA % TallyInterior )
        allocate ( Tally_F_P_Form :: FA % TallyTotal )
        allocate ( Tally_F_P_Form :: FA % TallyChange )
        allocate ( FA % TallyBoundaryLocal  ( A % nBoundaries ) )
        allocate ( FA % TallyBoundaryGlobal ( A % nBoundaries ) )
        do iB = 1, A % nBoundaries
          allocate &
            ( Tally_F_P_Form :: FA % TallyBoundaryLocal  ( iB ) % Element )
          allocate &
            ( Tally_F_P_Form :: FA % TallyBoundaryGlobal ( iB ) % Element )
        end do !-- iB
      case default
        call Show ( 'FluidType not recognized', CONSOLE % ERROR )
        call Show ( 'Fluid_ASC__Form', 'module', CONSOLE % ERROR )
        call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- FluidType   
    end if !-- allocated TallyInterior
      
    select type ( TI => FA % TallyInterior )
    class is ( Tally_F_D_Form )
      call TI % Initialize &
             ( A, MassUnitOption = MassUnitOption, &
               EnergyUnitOption = EnergyUnitOption, &
               MomentumUnitOption = MomentumUnitOption, &
               AngularMomentumUnitOption = AngularMomentumUnitOption )
    end select !-- TI

    select type ( TT => FA % TallyTotal )
    class is ( Tally_F_D_Form )
      call TT % Initialize &
             ( A, MassUnitOption = MassUnitOption, &
               EnergyUnitOption = EnergyUnitOption, &
               MomentumUnitOption = MomentumUnitOption, &
               AngularMomentumUnitOption = AngularMomentumUnitOption )
    end select !-- TT

    select type ( TC => FA % TallyChange )
    class is ( Tally_F_D_Form )
      call TC % Initialize &
             ( A, MassUnitOption = MassUnitOption, &
               EnergyUnitOption = EnergyUnitOption, &
               MomentumUnitOption = MomentumUnitOption, &
               AngularMomentumUnitOption = AngularMomentumUnitOption )
    end select !-- TC

    do iB = 1, A % nBoundaries
      select type ( TB => FA % TallyBoundaryLocal ( iB ) % Element )
      class is ( Tally_F_D_Form )
        call TB % Initialize &
               ( A, MassUnitOption = MassUnitOption, &
                 EnergyUnitOption = EnergyUnitOption, &
                 MomentumUnitOption = MomentumUnitOption, &
                 AngularMomentumUnitOption = AngularMomentumUnitOption )
      end select !-- TB
      select type ( TB => FA % TallyBoundaryGlobal ( iB ) % Element )
      class is ( Tally_F_D_Form )
        call TB % Initialize &
               ( A, MassUnitOption = MassUnitOption, &
                 EnergyUnitOption = EnergyUnitOption, &
                 MomentumUnitOption = MomentumUnitOption, &
                 AngularMomentumUnitOption = AngularMomentumUnitOption )
      end select !-- TB
    end do !-- iB

  end subroutine Initialize


  function Fluid_D_CSL ( FA ) result ( F )

    class ( Fluid_ASC_Form ), intent ( in ) :: &
      FA
    class ( Fluid_D_Form ), pointer :: &
      F

    select type ( FC => FA % Chart )
    class is ( Fluid_CSL_Form )
      F => FC % Fluid_D ( )
    class default
      call Show ( 'Fluid Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Fluid_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Fluid_D_CSL', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FC

  end function Fluid_D_CSL


  function Fluid_P_P_CSL ( FA ) result ( F )

    class ( Fluid_ASC_Form ), intent ( in ) :: &
      FA
    class ( Fluid_P_P_Form ), pointer :: &
      F

    select type ( FC => FA % Chart )
    class is ( Fluid_CSL_Form )
      F => FC % Fluid_P_P ( )
    class default
      call Show ( 'Fluid Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Fluid_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Fluid_P_P_CSL', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FC

  end function Fluid_P_P_CSL


  impure elemental subroutine Finalize ( FA )

    type ( Fluid_ASC_Form ), intent ( inout ) :: &
      FA

    call FA % FinalizeTemplate_ASC_C ( )

  end subroutine Finalize


  subroutine SetField ( FA, NameOutputOption )

    class ( Fluid_ASC_Form ), intent ( inout ) :: &
      FA
    character ( * ), intent ( in ), optional :: &
      NameOutputOption

    select type ( A => FA % Atlas )
    class is ( Atlas_SC_Template )

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
    associate ( nValues => C % nProperCells + C % nGhostCells )

    allocate ( Fluid_CSL_Form :: FA % Chart )

    select type ( FC => FA % Chart )
    class is ( Fluid_CSL_Form )
      call FC % Initialize &
             ( C, FA % FluidType, FA % VelocityUnit, FA % MassDensityUnit, &
               FA % EnergyDensityUnit, FA % NumberDensityUnit, &
               FA % TemperatureUnit, nValues, &
               NameOutputOption = NameOutputOption )
    end select !-- FC

    call A % AddField ( FA )

    end associate !-- nValues
    end select !-- C
    end select !-- A

  end subroutine SetField


end module Fluid_ASC__Form
