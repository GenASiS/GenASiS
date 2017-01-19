module RadiationMoments_BSLL_ASC_CSLD__Form

  use Basics
  use Manifolds
  use RadiationMoments_Form
  use RadiationMoments_CSL__Form
  use RadiationMoments_ASC__Form

  implicit none
  private

  type, public, extends ( Field_BSLL_ASC_CSLD_Template ) :: &
    RadiationMoments_BSLL_ASC_CSLD_Form
      type ( MeasuredValueForm ) :: &
        EnergyDensityUnit, &
        EnergyUnit, &
        MomentumUnit, &
        AngularMomentumUnit
      type ( MeasuredValueForm ), dimension ( 3 ) :: &
        VelocityUnit
      character ( LDF ) :: &
        RadiationType = ''
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      SetField
  end type RadiationMoments_BSLL_ASC_CSLD_Form

contains


  subroutine Initialize &
               ( RMB, B, RadiationType, NameOutputOption, VelocityUnitOption, &
                 EnergyDensityUnitOption, EnergyUnitOption, &
                 MomentumUnitOption, AngularMomentumUnitOption )

    class ( RadiationMoments_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      RMB
    class ( Bundle_SLL_ASC_CSLD_Form ), intent ( in ), target :: &
      B
    character ( * ), intent ( in )  :: &
      RadiationType
    character ( * ), intent ( in ), optional :: &
      NameOutputOption
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ), optional :: &
      VelocityUnitOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      EnergyDensityUnitOption, &
      EnergyUnitOption, &
      MomentumUnitOption, &
      AngularMomentumUnitOption

    ! class ( GeometryFlatForm ), pointer :: &
    !   GF

    if ( RMB % Type == '' ) &
      RMB % Type = 'a RadiationMoments_BSLL_ASC_CSLD'
    RMB % RadiationType = RadiationType

    if ( present ( EnergyDensityUnitOption ) ) &
      RMB % EnergyDensityUnit = EnergyDensityUnitOption
    if ( present ( VelocityUnitOption ) ) &
      RMB % VelocityUnit = VelocityUnitOption
    if ( present ( EnergyUnitOption ) ) &
      RMB % EnergyUnit = EnergyUnitOption
    if ( present ( MomentumUnitOption ) ) &
      RMB % MomentumUnit = MomentumUnitOption
    if ( present ( AngularMomentumUnitOption ) ) &
      RMB % AngularMomentumUnit = AngularMomentumUnitOption

    ! select type ( B )
    ! class is ( Bundle_SLL_ASC_CSLD_Form )

    ! RMB % nFibers = B % nFibers

    ! RMB % nBaseValues &
    !   = B % Base_CSLD % nProperCells  +  B % Base_CSLD % nGhostCells
    ! RMB % ChartBase => B % Base_CSLD

    ! allocate ( RMB % iaBaseCell ( size ( B % iaBaseCell ) ) )
    ! RMB % iaBaseCell = B % iaBaseCell

    ! GF => B % GeometryFiber ( )
    ! associate ( Energy => GF % Value ( :, GF % CENTER ( 1 ) ) )
    ! RMB % nEnergyValues = size ( Energy )
    ! allocate ( RMB % Energy ( RMB % nEnergyValues ) )
    ! RMB % Energy = Energy
    ! end associate !-- Energy

    ! RMB % EnergyUnit = GF % Unit ( GF % CENTER ( 1 ) )

    ! associate ( AF => B % FiberMaster )
    ! select type ( CF => AF % Chart )
    ! class is ( Chart_SLL_Form )
    !   RMB % ChartFiber => CF
    ! end select !--C
    ! end associate !-- AF

    ! class default
    !   call Show ( 'Bundle type not recognized', CONSOLE % ERROR )
    !   call Show ( 'RadiationMoments_BSLL_ASC_CSLD__Form', 'module', &
    !               CONSOLE % ERROR )
    !   call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
    !   call PROGRAM_HEADER % Abort ( )
    ! end select !-- B

    call RMB % InitializeTemplate_BSLL &
           ( B, NameOutputOption = NameOutputOption )

    ! nullify ( GF )

  end subroutine Initialize


  subroutine SetField ( FB, NameOutputOption )

    class ( RadiationMoments_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      FB
    character ( * ), intent ( in ), optional :: &
      NameOutputOption

    integer ( KDI ) :: &
      iF  !-- iFiber

    select type ( B => FB % Bundle )
    class is ( Bundle_SLL_ASC_CSLD_Form )

    allocate ( FB % Fiber )
    associate ( FBF => FB % Fiber )
    call FBF % Initialize ( B % nFibers )

    do iF = 1, B % nFibers
      allocate ( RadiationMoments_ASC_Form :: FBF % Atlas ( iF ) % Element )
      select type ( RMA => FBF % Atlas ( iF ) % Element )
      class is ( RadiationMoments_ASC_Form )

      select type ( AF => B % Fiber % Atlas ( iF ) % Element )
      class is ( Atlas_SC_Form )

      call RMA % Initialize &
             ( AF, FB % RadiationType, NameOutputOption = NameOutputOption, &
               VelocityUnitOption = FB % VelocityUnit, &
               EnergyDensityUnitOption = FB % EnergyDensityUnit, &
               EnergyUnitOption = FB % EnergyUnit, &
               MomentumUnitOption = FB % MomentumUnit, &
               AngularMomentumUnitOption = FB % AngularMomentumUnit )

      end select !-- AF
      end select !-- RMA
    end do !-- iF
    
    end associate !-- FBF
    end select !-- B

  end subroutine SetField


end module RadiationMoments_BSLL_ASC_CSLD__Form
