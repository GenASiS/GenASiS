module WoosleyHeger_07_A__Form

  !-- WoosleyHeger_07_Adiabatic__Form

  use GenASiS
  use WoosleyHeger_07__Form

  implicit none
  private

  type, public, extends ( WoosleyHeger_07_Form ) :: WoosleyHeger_07_A_Form
    type ( FieldSet_BM_Form ), allocatable :: &
      Pressure
    type ( GradientForm ), allocatable :: &
      GradientPressure
 contains
    procedure, private, pass :: &
      Initialize_H
    final :: &
      Finalize
  end type WoosleyHeger_07_A_Form

    private :: &
      InitializeUniverse, &
      InitializeDiagnostics, &
      SetInitial, &
      SetReference


contains


  subroutine Initialize_H ( U, Name, CommunicatorOption, UnitsTypeOption )

    class ( WoosleyHeger_07_A_Form ), intent ( inout ), target :: &
      U
    character ( * ), intent ( in ) :: &
      Name
    type ( CommunicatorForm ), intent ( in ), target, optional :: &
      CommunicatorOption
    character ( * ), intent ( in ), optional :: &
      UnitsTypeOption

    if ( U % Type == '' ) &
      U % Type = 'a WoosleyHeger_07_A'

    call InitializeUniverse ( U, Name )
    call InitializeDiagnostics ( U )

  end subroutine Initialize_H


  subroutine Finalize ( WH )

    type ( WoosleyHeger_07_A_Form ), intent ( inout ) :: &
      WH

    if ( allocated ( WH % GradientPressure ) ) &
      deallocate ( WH % GradientPressure )
    if ( allocated ( WH % Pressure ) ) &
      deallocate ( WH % Pressure )

  end subroutine Finalize


  subroutine InitializeUniverse ( WH, Name )

    class ( WoosleyHeger_07_A_Form ), intent ( inout ), target :: &
      WH
    character ( * ), intent ( in )  :: &
      Name

    real ( KDR ) :: &
      FinishTime

    FinishTime  =  1.5_KDR  *  UNIT % SECOND

    call WH % Initialize &
           ( RadiationName = [ 'None' ], &
             RadiationType = [ 'NONE' ], &
             FormalismType = 'NONE', &
             FluidType = 'HEAVY_NUCLEUS', &
             GravitationType = 'NEWTON_SG', &
             Name = Name, &
             UnitsTypeOption = 'ASTROPHYSICS', &
             FinishTimeOption = FinishTime, &
             nCellsPolarOption = 128, &
             nWriteOption = 30 )

    WH % Integrator % SetInitial    =>  SetInitial
    WH % Integrator % SetReference  =>  SetReference
    WH % Integrator % System        =>  WH

  end subroutine InitializeUniverse


  subroutine InitializeDiagnostics ( WH )

    class ( WoosleyHeger_07_A_Form ), intent ( inout ) :: &
      WH

      allocate &
        ( WH % Pressure, &
        WH % GradientPressure )
      select type ( I  =>  WH % Integrator )
        class is ( Integrator_CS_Form )
      select type ( F  =>  I % CurrentSet_X )
        class is ( Fluid_P_HN_Form )
      associate &
        (  P  =>  WH % Pressure, &
          GP  =>  WH % GradientPressure, &
           G  =>  I % Geometry_X, &
           S  =>  I % Checkpoint_X )
      call  P % Initialize ( F, NameOption = 'Pressure', &
                             iaSelected = [ F % PRESSURE ] )
      call GP % Initialize ( G, P )
      call  S % AddFieldSet ( GP )
      end associate !-- P, etc.
      end select !-- F
      end select !-- I

  end subroutine InitializeDiagnostics


  subroutine SetInitial ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    select type ( WH  =>  I % System )
      class is ( WoosleyHeger_07_Form )

    call WH % SetFluid ( )

    end select !-- WH

  end subroutine SetInitial


  subroutine SetReference ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    select type ( WH  =>  I % System )
      class is ( WoosleyHeger_07_A_Form )
    if ( allocated ( WH % GradientPressure ) ) then
      associate ( GP   =>  WH % GradientPressure )
      associate ( GPV  =>  GP % Storage_GS % Value )
      call GP % Compute ( iD = 1 )
      GPV  =  -1.0_KDR  *  GPV
      end associate !-- GPV
      end associate !-- GP
    end if
    end select !-- WH

  end subroutine SetReference


end module WoosleyHeger_07_A__Form
