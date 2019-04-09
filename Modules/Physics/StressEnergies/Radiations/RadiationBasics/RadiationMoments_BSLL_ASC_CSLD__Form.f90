module RadiationMoments_BSLL_ASC_CSLD__Form

  use Basics
  use Mathematics
  use StressEnergyBasics
  use RadiationMoments_Form
  use RadiationMoments_ASC__Form

  implicit none
  private

  type, public, extends ( Current_BSLL_ASC_CSLD_Template ) :: &
    RadiationMoments_BSLL_ASC_CSLD_Form
      integer ( KDI ) :: &
        nEnergyValues = 0
      real ( KDR ) :: &
        LimiterParameter
      real ( KDR ), dimension ( : ), allocatable :: &
        Energy
      logical ( KDL ) :: &
        UseLimiter
      type ( StressEnergyUnitsForm ) :: &
        UnitsSpectral
      class ( StressEnergyUnitsForm ), pointer :: &
        Units => null ( )
      character ( LDF ) :: &
        RadiationType = ''
      class ( FieldAtlasTemplate ), allocatable :: &
        EnergyIntegral
      class ( Field_BSLL_ASC_CSLD_Template ), pointer :: &
        Interactions_BSLL_ASC_CSLD => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      RadiationMoments
    procedure, public, pass :: &
      ComputeTally
    procedure, public, pass :: &
      SetInteractions
    final :: &
      Finalize
    procedure, private, pass :: &
      SetType
    procedure, private, pass :: &
      SetField
    procedure, private, pass :: &
      AllocateField
    procedure, public, pass :: &
      ComputeEnergyIntegral
  end type RadiationMoments_BSLL_ASC_CSLD_Form

contains


  subroutine Initialize &
               ( RMB, B, RadiationType, Units, NameShortOption, &
                 UseLimiterOption, LimiterParameterOption )

    class ( RadiationMoments_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      RMB
    class ( Bundle_SLL_ASC_CSLD_Form ), intent ( in ), target :: &
      B
    character ( * ), intent ( in )  :: &
      RadiationType
    class ( StressEnergyUnitsForm ), intent ( in ), target :: &
      Units
    character ( * ), intent ( in ), optional :: &
      NameShortOption
    logical ( KDL ), intent ( in ), optional :: &
      UseLimiterOption
    real ( KDR ), intent ( in ), optional :: &
      LimiterParameterOption

    type ( MeasuredValueForm ) :: &
      ParticleEnergyUnit
    character ( LDL ) :: &
      NameShort
    class ( GeometryFlatForm ), pointer :: &
      GF

    call RMB % SetType ( )

    RMB % RadiationType = RadiationType

    RMB % Units => Units

    associate ( AF => B % FiberMaster )
    select type ( CF => AF % Chart )
    class is ( Chart_SLL_Form )
      ParticleEnergyUnit  =  CF % CoordinateUnit ( 1 )
    end select !--CF
    end associate !-- AF
    RMB % UnitsSpectral  =  Units
    RMB % UnitsSpectral % NumberDensity &
      =  Units % NumberDensity  *  ParticleEnergyUnit ** (-3)
    RMB % UnitsSpectral % EnergyDensity &
      =  Units % EnergyDensity  *  ParticleEnergyUnit ** (-3)
    RMB % UnitsSpectral % MomentumDensity_U &
      =  Units % MomentumDensity_U  *  ParticleEnergyUnit ** (-3)
    RMB % UnitsSpectral % MomentumDensity_D &
      =  Units % MomentumDensity_D  *  ParticleEnergyUnit ** (-3)

    RMB % UseLimiter = .false.
    if ( present ( UseLimiterOption ) ) &
      RMB % UseLimiter = UseLimiterOption
    call PROGRAM_HEADER % GetParameter &
           ( RMB % UseLimiter, 'UseLimiter' )

    RMB % LimiterParameter = 1.4_KDR
    if ( present ( LimiterParameterOption ) ) &
      RMB % LimiterParameter = LimiterParameterOption
    call PROGRAM_HEADER % GetParameter &
           ( RMB % LimiterParameter, 'LimiterParameter' )

    GF => B % GeometryFiber ( )
    associate ( Energy => GF % Value ( :, GF % CENTER_U ( 1 ) ) )
    RMB % nEnergyValues = size ( Energy )
    allocate ( RMB % Energy ( RMB % nEnergyValues ) )
    RMB % Energy = Energy
    end associate !-- Energy

    NameShort = 'RadiationMoments'
    if ( present ( NameShortOption ) ) &
      NameShort = NameShortOption

    call RMB % InitializeTemplate_BSLL_ASC_CSLD ( B, NameShort )

    call Show ( RMB % RadiationType, 'RadiationType', &
                RMB % IGNORABILITY )

    nullify ( GF )

  end subroutine Initialize


  function RadiationMoments ( RMB, iFiber ) result ( RMF )

    class ( RadiationMoments_BSLL_ASC_CSLD_Form ), intent ( in ), target :: &
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

  end function RadiationMoments


  subroutine ComputeTally ( CB, ComputeChangeOption, IgnorabilityOption )

    class ( RadiationMoments_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      CB
    logical ( KDL ), intent ( in ), optional :: &
      ComputeChangeOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

!    call CB % ComputeTallySections &
!           ( ComputeChangeOption = ComputeChangeOption, &
!             IgnorabilityOption  = IgnorabilityOption )

    call CB % ComputeEnergyIntegral ( )
    select type ( RMA => CB % EnergyIntegral )
    class is ( RadiationMoments_ASC_Form )
      call RMA % ComputeTally &
             ( ComputeChangeOption = ComputeChangeOption, &
               IgnorabilityOption  = IgnorabilityOption )
    end select !-- RMA

  end subroutine ComputeTally


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

    nullify ( RMB % Units )

    call RMB % FinalizeTemplate_BSLL_ASC_CSLD ( )

  end subroutine Finalize


  subroutine SetType ( RMB )

    class ( RadiationMoments_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      RMB

    if ( RMB % Type == '' ) &
      RMB % Type = 'a RadiationMoments_BSLL_ASC_CSLD'

  end subroutine SetType


  subroutine SetField ( FB )

    class ( RadiationMoments_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      FB

    logical ( KDL ) :: &
      SuppressWrite
    integer ( KDI ) :: &
      iF, &  !-- iFiber
      iE, &  !-- iEnergy
      Ignorability
    character ( 1 + 2 ) :: &
      EnergyNumber

    associate ( B => FB % Bundle_SLL_ASC_CSLD )

    !-- Fibers

    allocate ( FB % Fiber )
    associate ( FBF => FB % Fiber )
    call FBF % Initialize ( FB % nFibers )

    do iF = 1, FB % nFibers
      call FB % AllocateField ( FBF % Atlas ( iF ) % Element )
      select type ( RMA => FBF % Atlas ( iF ) % Element )
      class is ( RadiationMoments_ASC_Form )

      select type ( AF => B % Fiber % Atlas ( iF ) % Element )
      class is ( Atlas_SC_Form )

      call RMA % Initialize &
             ( AF, RadiationType = FB % RadiationType, &
               MomentsType = 'SPECTRAL', Units = FB % UnitsSpectral, &
               NameShortOption = FB % NameShort, &
               SuppressWriteSourcesOption = .false. )

      end select !-- AF
      end select !-- RMA
    end do !-- iF

    end associate !-- FBF

    !-- Sections

    allocate ( FB % Section )
    associate ( FBS => FB % Section )
    call FBS % Initialize ( FB % nSections )

    do iE = 1, FB % nEnergyValues
      SuppressWrite  =  ( mod ( iE - 1, B % sSectionsWrite ) /= 0 )
      write ( EnergyNumber, fmt = '(a1,i2.2)' ) '_', iE
      call FB % AllocateField ( FBS % Atlas ( iE ) % Element )
      select type ( RMA => FBS % Atlas ( iE ) % Element )
      class is ( RadiationMoments_ASC_Form )

      Ignorability = CONSOLE % INFO_5
      if ( iE == 1 ) then
        Ignorability = B % Base_ASC % Ignorability
        call Show ( 'Representative position space section', Ignorability )
      end if
   
      call RMA % Initialize &
             ( B % Base_ASC, RadiationType = FB % RadiationType, &
               MomentsType = 'SPECTRAL', Units = FB % UnitsSpectral, &
               NameShortOption = trim ( FB % NameShort ) // EnergyNumber, &
               UseLimiterOption = FB % UseLimiter, &
               SuppressWriteOption = SuppressWrite, &
               SuppressWriteSourcesOption = .true., &
               LimiterParameterOption = FB % LimiterParameter, &
               IgnorabilityOption = Ignorability )

      end select !-- RMA
    end do !-- iE

    end associate !-- FBS

    !-- EnergyIntegral

    call FB % AllocateField ( FB % EnergyIntegral )
    select type ( EI => FB % EnergyIntegral )
    class is ( RadiationMoments_ASC_Form )
      call EI % Initialize &
             ( B % Base_ASC, RadiationType = FB % RadiationType, &
               MomentsType = 'GREY', Units = FB % Units, &
               NameShortOption = trim ( FB % NameShort ) // '_Integral', &
               IgnorabilityOption = CONSOLE % INFO_5 )
    end select !-- EI

    end associate !-- B

  end subroutine SetField


  subroutine AllocateField ( RMB, RMA )

    class ( RadiationMoments_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      RMB
    class ( FieldAtlasTemplate ), intent ( out ), allocatable :: &
      RMA

    allocate ( RadiationMoments_ASC_Form :: RMA )

  end subroutine AllocateField


  subroutine ComputeEnergyIntegral ( RMB )

    class ( RadiationMoments_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
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

    select type ( RMA => RMB % EnergyIntegral )
    class is ( RadiationMoments_ASC_Form )
      RMEI => RMA % RadiationMoments ( )
    end select !-- RMA

    G => MS % Base_CSLD % Geometry ( )

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

      RMF => RMB % RadiationMoments ( iF )
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
