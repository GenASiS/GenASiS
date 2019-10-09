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
      InitializeRadiationCentralCore

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
    call WH % SetFluid ( )

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

    call WH % Initialize &
           ( RadiationName = RadiationName, RadiationType = RadiationType, &
             MomentsType = MomentsType, FluidType = 'HEAVY_NUCLEUS', &
             GeometryType = GeometryType, Name = Name, &
             ShockThresholdOption = 1.0_KDR, nWriteOption = 30 )

  end subroutine InitializeRadiationCentralCore


end module WoosleyHeger_07_RM__Form
