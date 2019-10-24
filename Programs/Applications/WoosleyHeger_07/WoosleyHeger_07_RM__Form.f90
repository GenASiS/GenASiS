module WoosleyHeger_07_RM__Form

  !-- WoosleyHeger_07_RadiationMoments_Form

  use GenASiS
  use WoosleyHeger_07__Template

  implicit none
  private

  type, public, extends ( WoosleyHeger_07_Template ) :: WoosleyHeger_07_RM_Form
  contains
    procedure, private, pass :: &
      Initialize_WH
    generic, public :: &
      Initialize => Initialize_WH
    final :: &
      Finalize
  end type WoosleyHeger_07_RM_Form

    private :: &
      InitializeRadiationCentralCore, &
      SetProblem

      private :: &
        InitializeInteractions

    character ( LDF ), private :: &
      InteractionsType

contains


  subroutine Initialize_WH ( WH, MomentsType, Name )

    class ( WoosleyHeger_07_RM_Form ), intent ( inout ), target :: &
      WH
    character ( * ), intent ( in )  :: &
      MomentsType, &
      Name

    if ( WH % Type == '' ) &
      WH % Type = 'a WoosleyHeger_07_RM'

    call InitializeRadiationCentralCore ( WH, MomentsType, Name )
    call SetProblem ( WH, MomentsType )

  end subroutine Initialize_WH


  subroutine Finalize ( WH )

    type ( WoosleyHeger_07_RM_Form ), intent ( inout ) :: &
      WH

  end subroutine Finalize


  subroutine InitializeRadiationCentralCore ( WH, MomentsType, Name )

    class ( WoosleyHeger_07_RM_Form ), intent ( inout ), target :: &
      WH
    character ( * ), intent ( in ) :: &
      MomentsType, &
      Name

    real ( KDR ) :: &
      MaxEnergy, &
      MinWidthEnergy, &
      EnergyScale
    character ( LDL ) :: &
      GeometryType
    character ( LDL ), dimension ( : ), allocatable :: &
      RadiationName, &
      RadiationType

    GeometryType = 'NEWTONIAN'
    call PROGRAM_HEADER % GetParameter ( GeometryType, 'GeometryType' )

    allocate ( RadiationName ( 2 ) )
    allocate ( RadiationType ( 2 ) )
    RadiationName  =  [ 'Nu_E   ', 'NuBar_E' ]
    RadiationType  =  [ 'NEUTRINOS_E    ', 'NEUTRINOS_E_BAR' ]

    MinWidthEnergy  =  2.0_KDR    *  UNIT % MEGA_ELECTRON_VOLT  !-- Geometric
    MaxEnergy       =  310.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT

    EnergyScale  =  10.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT  !-- Compactified

    call WH % Initialize &
           ( RadiationName = RadiationName, RadiationType = RadiationType, &
             MomentsType = MomentsType, FluidType = 'HEAVY_NUCLEUS', &
             GeometryType = GeometryType, Name = Name, &
             ShockThresholdOption = 1.0_KDR, &
             MinWidthEnergyOption = MinWidthEnergy, &
             MaxEnergyOption = MaxEnergy, &
             EnergyScaleOption = EnergyScale, &
             nWriteOption = 30 )

  end subroutine InitializeRadiationCentralCore


  subroutine SetProblem ( WH, MomentsType )

    class ( WoosleyHeger_07_RM_Form ), intent ( inout ) :: &
      WH
    character ( * ), intent ( in )  :: &
      MomentsType

    call WH % SetFluid ( )

  end subroutine SetProblem


  subroutine InitializeInteractions ( WH, MomentsType )

    class ( WoosleyHeger_07_RM_Form ), intent ( inout ) :: &
      WH
    character ( * ), intent ( in )  :: &
      MomentsType

  end subroutine InitializeInteractions


end module WoosleyHeger_07_RM__Form
