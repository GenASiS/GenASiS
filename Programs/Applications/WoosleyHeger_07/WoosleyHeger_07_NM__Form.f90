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

    FinishTime  =  0.7_KDR  *  UNIT % SECOND
    call PROGRAM_HEADER % GetParameter ( FinishTime, 'FinishTime' )

    allocate ( Interactions_NM_G_Form :: WH % Interactions_NM_G ( 3 ) )

    call WH % Initialize &
           ( RadiationName = [ 'Neutrinos_E ',   &
                               'Neutrinos_EB', &
                               'Neutrinos_HL' ], &
             RadiationType = [ 'NEUTRINOS_E ',   &
                               'NEUTRINOS_EB', &
                               'NEUTRINOS_HL' ], &
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
