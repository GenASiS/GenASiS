module WoosleyHeger_07_NM__Form

  !-- WoosleyHeger_07_NeutrinoMoments_Form

  use GenASiS
  use WoosleyHeger_07__Form

  implicit none
  private

  type, public, extends ( WoosleyHeger_07_Form ) :: WoosleyHeger_07_NM_Form
  contains
    procedure, public, pass :: &
      Initialize_NM
    final :: &
      Finalize
  end type WoosleyHeger_07_NM_Form

    private :: &
      InitializeUniverse, &
      SetInitial


contains


  subroutine Initialize_NM ( U, FormalismType, Name )

    class ( WoosleyHeger_07_NM_Form ), intent ( inout ), target :: &
      U
    character ( * ), intent ( in ) :: &
      FormalismType, &
      Name

    if ( U % Type == '' ) &
      U % Type = 'a WoosleyHeger_07_NM'

    call InitializeUniverse ( U, FormalismType, Name )

  end subroutine Initialize_NM


  subroutine Finalize ( WH )

    type ( WoosleyHeger_07_NM_Form ), intent ( inout ) :: &
      WH

  end subroutine Finalize


  subroutine InitializeUniverse ( WH, FormalismType, Name )

    class ( WoosleyHeger_07_NM_Form ), intent ( inout ), target :: &
      WH
    character ( * ), intent ( in )  :: &
      FormalismType, &
      Name

    real ( KDR ) :: &
      FinishTime
    logical ( KDL ) :: &
      Neutrinos_EB, &
      Neutrinos_HL
    character ( LDL ), dimension ( : ), allocatable :: &
      RadiationName, &
      RadiationType

    FinishTime  =  0.7_KDR  *  UNIT % SECOND

    Neutrinos_EB  =  .true.
    Neutrinos_HL  =  .true.
    call PROGRAM_HEADER % GetParameter &
           ( Neutrinos_EB, 'Neutrinos_EB' )
    call PROGRAM_HEADER % GetParameter &
           ( Neutrinos_HL, 'Neutrinos_HL' )

    if ( Neutrinos_EB .and. Neutrinos_HL ) then
      allocate ( Interactions_NM_G_Form :: WH % Interactions_NM_G ( 3 ) )
      allocate ( RadiationName ( 3 ) )
      allocate ( RadiationType ( 3 ) )
      RadiationName  =  [ 'Neutrinos_E ',   &
                          'Neutrinos_EB', &
                          'Neutrinos_HL' ]
      RadiationType  =  [ 'NEUTRINOS_E ',   &
                          'NEUTRINOS_EB', &
                          'NEUTRINOS_HL' ]
    else if ( Neutrinos_EB ) then
      allocate ( Interactions_NM_G_Form :: WH % Interactions_NM_G ( 2 ) )
      allocate ( RadiationName ( 2 ) )
      allocate ( RadiationType ( 2 ) )
      RadiationName  =  [ 'Neutrinos_E ',   &
                          'Neutrinos_EB' ]
      RadiationType  =  [ 'NEUTRINOS_E ',   &
                          'NEUTRINOS_EB' ]
    else
      allocate ( Interactions_NM_G_Form :: WH % Interactions_NM_G ( 1 ) )
      allocate ( RadiationName ( 1 ) )
      allocate ( RadiationType ( 1 ) )
      RadiationName  =  [ 'Neutrinos_E' ]
      RadiationType  =  [ 'NEUTRINOS_E ' ]
    end if

    !-- Initialization

    call WH % Initialize &
           ( RadiationName = RadiationName, &
             RadiationType = RadiationType, &
             FormalismType = FormalismType, &
             FluidType = 'HEAVY_NUCLEUS', &
             GravitationType = 'NEWTON_SG', &
             Name = Name, &
             UnitsTypeOption = 'ASTROPHYSICS', &
             FinishTimeOption = FinishTime, &
             nCellsPolarOption = 128, &
             nWriteOption = 10 )

    WH % Integrator % SetInitial  =>  SetInitial
    WH % Integrator % System      =>  WH

  end subroutine InitializeUniverse


  subroutine SetInitial ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    select type ( WH  =>  I % System )
      class is ( WoosleyHeger_07_Form )

    call WH % SetFluid ( )

    end select !-- WH

  end subroutine SetInitial


end module WoosleyHeger_07_NM__Form
