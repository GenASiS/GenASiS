module RadiationMoments_BSLL_ASC_CSLD__Form

  use Basics
  use Mathematics
  use RadiationMoments_Form
  use RadiationMoments_ASC__Form

  implicit none
  private

  type, public, extends ( Current_BSLL_ASC_CSLD_Template ) :: &
    RadiationMoments_BSLL_ASC_CSLD_Form
      integer ( KDI ) :: &
        nEnergyValues = 0
      real ( KDR ), dimension ( : ), allocatable :: &
        Energy
      type ( MeasuredValueForm ) :: &
        EnergyDensityUnit, &
        EnergyUnit, &
        MomentumUnit, &
        AngularMomentumUnit
      type ( MeasuredValueForm ), dimension ( 3 ) :: &
        Velocity_U_Unit, &
        MomentumDensity_U_Unit, &
        MomentumDensity_D_Unit
      character ( LDF ) :: &
        RadiationType = ''
      class ( RadiationMoments_ASC_Form ), allocatable :: &
        EnergyIntegral
      class ( Field_BSLL_ASC_CSLD_Template ), pointer :: &
        Interactions_BSLL_ASC_CSLD => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      RadiationMomentsFiber
    procedure, public, pass :: &
      ComputeTally
    procedure, public, pass :: &
      SetInteractions
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type RadiationMoments_BSLL_ASC_CSLD_Form

    private :: &
      ComputeEnergyIntegral

contains


  subroutine Initialize &
               ( RMB, B, RadiationType, NameShortOption, &
                 Velocity_U_UnitOption, MomentumDensity_U_UnitOption, &
                 MomentumDensity_D_UnitOption, EnergyDensityUnitOption, &
                 EnergyUnitOption, MomentumUnitOption, &
                 AngularMomentumUnitOption )

    class ( RadiationMoments_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      RMB
    class ( Bundle_SLL_ASC_CSLD_Form ), intent ( in ), target :: &
      B
    character ( * ), intent ( in )  :: &
      RadiationType
    character ( * ), intent ( in ), optional :: &
      NameShortOption
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ), optional :: &
      Velocity_U_UnitOption, &
      MomentumDensity_U_UnitOption, &
      MomentumDensity_D_UnitOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      EnergyDensityUnitOption, &
      EnergyUnitOption, &
      MomentumUnitOption, &
      AngularMomentumUnitOption

    character ( LDL ) :: &
      NameShort
    class ( GeometryFlatForm ), pointer :: &
      GF

    if ( RMB % Type == '' ) &
      RMB % Type = 'a RadiationMoments_BSLL_ASC_CSLD'
    RMB % RadiationType = RadiationType

    if ( present ( EnergyDensityUnitOption ) ) &
      RMB % EnergyDensityUnit = EnergyDensityUnitOption
    if ( present ( Velocity_U_UnitOption ) ) &
      RMB % Velocity_U_Unit = Velocity_U_UnitOption
    if ( present ( MomentumDensity_U_UnitOption ) ) &
      RMB % MomentumDensity_U_Unit = MomentumDensity_U_UnitOption
    if ( present ( MomentumDensity_D_UnitOption ) ) &
      RMB % MomentumDensity_D_Unit = MomentumDensity_D_UnitOption

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

    GF => B % GeometryFiber ( )
    associate ( Energy => GF % Value ( :, GF % CENTER ( 1 ) ) )
    RMB % nEnergyValues = size ( Energy )
    allocate ( RMB % Energy ( RMB % nEnergyValues ) )
    RMB % Energy = Energy
    end associate !-- Energy

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

    NameShort = 'RadiationMoments'
    if ( present ( NameShortOption ) ) &
      NameShort = NameShortOption

    call RMB % InitializeTemplate_BSLL_ASC_CSLD ( B, NameShort )

    nullify ( GF )

  end subroutine Initialize


  subroutine ComputeTally ( CB, ComputeChangeOption )

    class ( RadiationMoments_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      CB
    logical ( KDL ), intent ( in ), optional :: &
      ComputeChangeOption

!    call CB % ComputeTallySections ( ComputeChangeOption )

    call ComputeEnergyIntegral ( CB )
    call CB % EnergyIntegral % ComputeTally ( ComputeChangeOption )

  end subroutine ComputeTally


  function RadiationMomentsFiber ( RMB, iFiber ) result ( RMF )

    class ( RadiationMoments_BSLL_ASC_CSLD_Form ), intent ( in ) :: &
      RMB
    integer ( KDI ), intent ( in ) :: &
      iFiber
    class ( RadiationMomentsForm ), pointer :: &
      RMF

    select type ( RMA => RMB % Fiber % Atlas ( iFiber ) % Element )
    class is ( Field_ASC_Template )
      select type ( RMC => RMA % Chart )
      class is ( Field_CSL_Template )   
        select type ( RM => RMC % Field )
        class is ( RadiationMomentsForm )
          RMF => RM
        end select !-- RM
      end select !-- RMC
    end select !-- RMA

  end function RadiationMomentsFiber


  subroutine SetInteractions ( RMB, IB )

    class ( RadiationMoments_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      RMB
    class ( Field_BSLL_ASC_CSLD_Template ), intent ( in ), target :: &
      IB

    integer ( KDI ) :: &
      iF  !-- iFiber

    RMB % Interactions_BSLL_ASC_CSLD => IB

    select type ( MS => RMB % Bundle )
    class is ( Bundle_SLL_ASC_CSLD_Form )

    do iF = 1, MS % nFibers
      select type ( RMA => RMB % Fiber % Atlas ( iF ) % Element )
      class is ( RadiationMoments_ASC_Form )

      select type ( IA => IB % Fiber % Atlas ( iF ) % Element )
      class is ( Field_ASC_Template )
      call RMA % SetInteractions ( IA )
      end select !-- I

      end select !-- RMA
    end do !-- iF

    end select !-- MS

  end subroutine SetInteractions


  impure elemental subroutine Finalize ( RMB )

    type ( RadiationMoments_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      RMB

    if ( allocated ( RMB % Energy ) ) &
      deallocate ( RMB % Energy )

    call RMB % FinalizeTemplate_BSLL_ASC_CSLD ( )

  end subroutine Finalize


  subroutine SetField ( FB )

    class ( RadiationMoments_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      FB

    integer ( KDI ) :: &
      iF, &  !-- iFiber
      iE     !-- iEnergy
    character ( 1 + 2 ) :: &
      EnergyNumber

    associate ( B => FB % Bundle_SLL_ASC_CSLD )

    !-- Fibers

    allocate ( FB % Fiber )
    associate ( FBF => FB % Fiber )
    call FBF % Initialize ( FB % nFibers )

    do iF = 1, FB % nFibers
      allocate ( RadiationMoments_ASC_Form :: FBF % Atlas ( iF ) % Element )
      select type ( RMA => FBF % Atlas ( iF ) % Element )
      class is ( RadiationMoments_ASC_Form )

      select type ( AF => B % Fiber % Atlas ( iF ) % Element )
      class is ( Atlas_SC_Form )

      call RMA % Initialize &
             ( AF, FB % RadiationType, NameShortOption = FB % NameShort, &
               Velocity_U_UnitOption = FB % Velocity_U_Unit, &
               MomentumDensity_U_UnitOption = FB % MomentumDensity_U_Unit, &
               MomentumDensity_D_UnitOption = FB % MomentumDensity_D_Unit, &
               EnergyDensityUnitOption = FB % EnergyDensityUnit, &
               EnergyUnitOption = FB % EnergyUnit, &
               MomentumUnitOption = FB % MomentumUnit, &
               AngularMomentumUnitOption = FB % AngularMomentumUnit )

      end select !-- AF
      end select !-- RMA
    end do !-- iF

    end associate !-- FBF

    !-- Sections

    allocate ( FB % Section )
    associate ( FBS => FB % Section )
    call FBS % Initialize ( FB % nSections )

    do iE = 1, FB % nEnergyValues
      write ( EnergyNumber, fmt = '(a1,i2.2)' ) '_', iE
      allocate ( RadiationMoments_ASC_Form :: FBS % Atlas ( iE ) % Element )
      select type ( RMA => FBS % Atlas ( iE ) % Element )
      class is ( RadiationMoments_ASC_Form )
        call RMA % Initialize &
               ( B % Base_ASC, FB % RadiationType, &
                 NameShortOption = trim ( FB % NameShort ) // EnergyNumber, &
                 Velocity_U_UnitOption = FB % Velocity_U_Unit, &
                 MomentumDensity_U_UnitOption = FB % MomentumDensity_U_Unit, &
                 MomentumDensity_D_UnitOption = FB % MomentumDensity_D_Unit, &
                 EnergyDensityUnitOption = FB % EnergyDensityUnit, &
                 EnergyUnitOption = FB % EnergyUnit, &
                 MomentumUnitOption = FB % MomentumUnit, &
                 AngularMomentumUnitOption = FB % AngularMomentumUnit, &
                 IgnorabilityOption = CONSOLE % INFO_5 )
      end select !-- RMA
    end do !-- iE

    end associate !-- FBS

    !-- EnergyIntegral

    allocate ( FB % EnergyIntegral )
    associate ( EI => FB % EnergyIntegral )
      call EI % Initialize &
             ( B % Base_ASC, FB % RadiationType, &
               NameShortOption = trim ( FB % NameShort ) // '_Integral' )
    end associate !-- EI

    end associate !-- B

  end subroutine SetField


  subroutine ComputeEnergyIntegral ( RMB )

    type ( RadiationMoments_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      RMB

    integer ( KDI ) :: &
      iC, &  !-- iConserved
      iF     !-- iFiber
    real ( KDR ), dimension ( : ), allocatable :: &
      Integral
    type ( Real_1D_Form ), dimension ( : ), allocatable :: &
      Integrand
    class ( GeometryFlatForm ), pointer :: &
      G
    type ( VolumeIntegralForm ) :: &
      VI
    class ( RadiationMomentsForm ), pointer :: &
      RMEI, &
      RMF

    associate ( MS => RMB % Bundle_SLL_ASC_CSLD )

    RMEI => RMB % EnergyIntegral % RadiationMoments ( )
    G    => MS % Base_CSLD % Geometry ( )

    allocate ( Integral  ( RMEI % N_CONSERVED ) )
    allocate ( Integrand ( RMEI % N_CONSERVED ) )
    do iC = 1, RMEI % N_CONSERVED
      call Integrand ( iC ) % Initialize ( RMB % nEnergyValues )
    end do !-- iC

    do iF = 1, MS % nFibers
      associate &
        ( iaC => RMEI % iaConserved, &
          iBC => MS % iaBaseCell ( iF ), &
          CF  => MS % Fiber_CSLL )

      RMF => RMB % RadiationMomentsFiber ( iF )
      do iC = 1, RMF % N_CONSERVED
        Integrand ( iC ) % Value = RMF % Value ( :, iaC ( iC ) )
      end do !-- iC

      call VI % Compute ( CF, Integrand, Integral )

      do iC = 1, RMEI % N_CONSERVED
        RMEI % Value ( iBC, iaC ( iC ) ) = Integral ( iC ) 
      end do !-- iC
      
      end associate !-- iaC, etc.
    end do !-- iF

    call RMEI % ComputeFromConserved ( G )

    end associate !-- MS
    nullify ( G, RMEI, RMF )
    
  end subroutine ComputeEnergyIntegral


end module RadiationMoments_BSLL_ASC_CSLD__Form
