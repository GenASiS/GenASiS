module Integrator_CS_1D_C_CS__Form

  !-- Integrator_CurrentSet_1D_Collected_CurrentSet__Form

  use Basics
  use Fields
  use Steps
  use Series_CS_1D_C_CS__Form
  use Integrator_H__Form
  use Integrator_CS_1D_CS__Form

  implicit none
  private

  type, public, extends ( Integrator_CS_1D_CS_Form ) :: &
    Integrator_CS_1D_C_CS_Form
      class ( CurrentSetForm ), dimension ( : ), allocatable :: &
        CurrentSet_X_1D
      class ( EigenspeedSet_F_Form ), dimension ( :, : ), allocatable :: &
        EigenspeedSet_X_1D
  contains
    procedure, private, pass :: &  !-- 1
      Initialize_CS_1D_C_CS
    generic, public :: &
      Initialize => Initialize_CS_1D_C_CS
    final :: &
      Finalize
    procedure, public, pass :: &  !-- 2
      ShowParameters
    procedure, public, pass :: &   !-- 2
      ShowFields
    procedure, public, pass :: &   !-- 2
      PrepareEvolution
    procedure, public, pass :: &   !-- 3
      UpdateHost => UpdateHost_CS_1D
    procedure, private, pass :: &
      ComputeTally_1D
  end type Integrator_CS_1D_C_CS_Form

    private :: &
      InitializeSeries_CS_1D_C_CS


contains


  subroutine Initialize_CS_1D_C_CS &
               ( I, nCurrentSets_1D, CommunicatorOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, Unit_T_Option, T_FinishOption, &
                 nWriteOption )

    class ( Integrator_CS_1D_C_CS_Form ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ) :: &
      nCurrentSets_1D
    type ( CommunicatorForm ), intent ( in ), target, optional :: &
      CommunicatorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      DeviceMemoryOption, &
      PinnedMemoryOption, &
      DevicesCommunicateOption
    type ( QuantityForm ), intent ( in ), optional :: &
      Unit_T_Option
    real ( KDR ), intent ( in ), optional :: &
      T_FinishOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    integer ( KDI ) :: &
      iCS, &
      iD
    character ( 1 ) :: &
      Label_CS, &
      Label_D

    if ( I % Type == '' ) &
      I % Type = 'an Integrator_CS_1D_C_CS'

    call I % Integrator_CS_1D_CS_Form % Initialize &
           ( CommunicatorOption = CommunicatorOption, &
             NameOption = NameOption, &
             DeviceMemoryOption = DeviceMemoryOption, &
             PinnedMemoryOption = PinnedMemoryOption, &
             DevicesCommunicateOption = DevicesCommunicateOption, &
             Unit_T_Option = Unit_T_Option, & 
             T_FinishOption = T_FinishOption, &
             nWriteOption = nWriteOption )

    I % nCurrentSets  =  nCurrentSets_1D

    !-- CurrentSet, if necessary

    if ( .not. allocated ( I % CurrentSet_X_1D ) ) then
      allocate ( I % CurrentSet_X_1D ( nCurrentSets_1D ) )
      do iCS  =  1,  nCurrentSets_1D
      associate &
        ( CS  =>  I % CurrentSet_X_1D ( iCS ), &
           G  =>  I % Geometry_X )
      write ( Label_CS, fmt = '(i1.1)' ) iCS
      call CS % Initialize ( G, NameOption = 'CurrentSet_X_1D_' // Label_CS )
      end associate !-- CS, etc.
      end do !-- iCS 
    end if

    !-- EigenspeedSet

    allocate ( I % EigenspeedSet_X_1D ( 3, nCurrentSets_1D ) )
    do iCS  =  1, nCurrentSets_1D
      do iD  =  1, 3
        associate &
          ( ES  =>  I % EigenspeedSet_X_1D ( iD, iCS ), &
            CS  =>  I % CurrentSet_X_1D ( iCS ) )
        write ( Label_CS, fmt = '(i1.1)' ) iCS
        write ( Label_D,  fmt = '(i1.1)' ) iD
        call ES % Initialize &
               ( CS, CS, SuffixOption = Label_CS // '_' // Label_D ) 
        end associate !-- ES, etc.
      end do !-- iD
    end do !-- iCS

    !-- Stream

    do iCS  =  1, nCurrentSets_1D

      associate &
        ( CS_X_1D  =>  I % CurrentSet_X_1D ( iCS ), &
           S_X     =>  I % Checkpoint_X )

      call CS_X_1D % SetStream ( S_X )

      ! if ( allocated ( I % Step_1D ) ) then
      !   associate ( S  =>  I % Step_1D )
      !   call S % SetStream ( S_X )
      !   end associate !-- S
      ! end if

      end associate !-- CS_X, etc.

    end do !-- iCS

    !-- Series

    I % InitializeSeries  =>  InitializeSeries_CS_1D_C_CS
    
  end subroutine Initialize_CS_1D_C_CS


  impure elemental subroutine Finalize ( I )

    type ( Integrator_CS_1D_C_CS_Form ), intent ( inout ) :: &
      I

    if ( allocated ( I % EigenspeedSet_X_1D ) ) &
      deallocate ( I % EigenspeedSet_X_1D )
    if ( allocated ( I % CurrentSet_X_1D ) ) &
      deallocate ( I % CurrentSet_X_1D )

  end subroutine Finalize


  subroutine ShowParameters ( I )

    class ( Integrator_CS_1D_C_CS_Form ), intent ( in ) :: &
      I

    integer ( KDI ) :: &
      iCS

    call I % Integrator_CS_1D_CS_Form % ShowParameters ( )

  end subroutine ShowParameters


  subroutine ShowFields ( I )

    class ( Integrator_CS_1D_C_CS_Form ), intent ( in ) :: &
      I

    integer ( KDI ) :: &
      iCS, &
      iD

    call I % Integrator_CS_Form % ShowFields ( )

    do iCS  =  1,  I % nCurrentSets

      call I % CurrentSet_X_1D ( iCS ) % Show ( )

      do iD  =  1, 3
        call I % EigenspeedSet_X_1D ( iD, iCS ) % Show ( )
      end do !-- iD

    end do !-- iCS

  end subroutine ShowFields


  subroutine PrepareEvolution ( I )

    class ( Integrator_CS_1D_C_CS_Form ), intent ( inout ) :: &
      I

    integer ( KDI ) :: &
      iCS

    call I % Integrator_CS_Form % PrepareEvolution ( )

    if ( .not. allocated ( I % CurrentSet_X_1D ) ) &
      return

    do iCS  =  1,  I % nCurrentSets
      associate ( CS  =>  I % CurrentSet_X_1D ( iCS ) )
      call CS % UpdateDevice ( )
      call CS % ExchangeGhostData ( )
      call CS % ComputeFromInitial ( )
      call CS % ApplyBoundaryConditions ( )
      call CS % UpdateHost ( )
      end associate !-- CS
    end do !-- iCS

  end subroutine PrepareEvolution
  
  
  subroutine UpdateHost_CS_1D ( I )

    class ( Integrator_CS_1D_C_CS_Form ), intent ( inout ) :: &
      I

    integer ( KDI ) :: &
      iCS

    call I % Integrator_CS_Form % UpdateHost ( )

    do iCS  =  1,  I % nCurrentSets
      associate ( CS  =>  I % CurrentSet_X_1D ( iCS ) )
      call CS % UpdateHost ( )
      end associate !-- CS
    end do

  end subroutine UpdateHost_CS_1D


  subroutine ComputeTally_1D ( I, ChangeOption, IgnorabilityOption )

    class ( Integrator_CS_1D_C_CS_Form ), intent ( inout ) :: &
      I
    logical ( KDL ), intent ( in ), optional :: &
      ChangeOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    integer ( KDI ) :: &
      iCS

    if ( .not. allocated ( I % CurrentSet_X_1D ) ) &
      return

    do iCS  =  1,  I % nCurrentSets
      associate ( CS => I % CurrentSet_X_1D ( iCS ) )
      call CS % ComputeTally &
             ( ChangeOption = ChangeOption, &
               IgnorabilityOption = IgnorabilityOption )
      end associate !-- CS
    end do

  end subroutine ComputeTally_1D


  subroutine InitializeSeries_CS_1D_C_CS ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    allocate ( Series_CS_1D_C_CS_Form :: I % Series )

    select type ( I )
      class is ( Integrator_CS_1D_C_CS_Form )
    select type ( S  =>  I % Series )
      class is ( Series_CS_1D_C_CS_Form )
    call S % Initialize &
      ( I % CurrentSet_X_1D, I % CurrentSet_X, I % GridImageStream, &
        I % dT_Label, I % Unit_T, I % dT_Candidate, I % T, &
        I % Communicator % Rank, I % nWrite, I % iCycle )
    end select !-- S
    end select !-- I

  end subroutine InitializeSeries_CS_1D_C_CS


end module Integrator_CS_1D_C_CS__Form
