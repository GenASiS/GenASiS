module ImposedAdvection_VL__Form

  !-- ImposedAdvection_VL__Form

  use Basics
  use Mathematics
  use CurrentSet_VLC__Form

  implicit none
  private

  type, public :: ImposedAdvection_VL_Form
    real ( KDR ), private :: &
      SpeedInitial, &
      LengthInitial, &
      DensityInitial, &
      DimensionFactor, &
      T_Finish
    character ( LDL ), private :: &
      CoordinateSystem, &
      AdvectionType
    type ( Integrator_CS_Form ), allocatable :: &
      Integrator
    type ( CurrentSet_VLC_Form ), allocatable :: &
      Reference, &
      Difference
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Evolve
    final :: &
      Finalize
  end type ImposedAdvection_VL_Form

    private :: &
      SetParameters, &
      SetInitial, &
      ShowSystem, &
      SetReference

      private :: &
        SetCurrentSet_VLC, &
        ComputeError


contains


  subroutine Initialize ( IA, CoordinateSystem, AdvectionType )

    class ( ImposedAdvection_VL_Form ), intent ( inout ), target :: &
      IA
    character ( * ), intent ( in ) :: &
      CoordinateSystem, &
      AdvectionType
      
    logical ( KDL ) :: &
      DeviceMemory, &   
      PinnedMemory, &   
      DevicesCommunicate

    call Show ( 'Initializing an ImposedAdvection_VL' )
    
    DeviceMemory  =  OffloadEnabled ( )  .and.  NumberOfDevices ( ) >= 1
    call PROGRAM_HEADER % GetParameter ( DeviceMemory, 'DeviceMemory' )

    PinnedMemory        =  DeviceMemory    
    DevicesCommunicate  =  DeviceMemory    
    call PROGRAM_HEADER % GetParameter &   
           ( PinnedMemory, 'PinnedMemory' )
    call PROGRAM_HEADER % GetParameter &
           ( DevicesCommunicate, 'DevicesCommunicate' )

    IA % CoordinateSystem  =  CoordinateSystem
    IA % AdvectionType     =  AdvectionType

    allocate ( IA % Integrator )
    associate ( I  =>  IA % Integrator )

    select case ( trim ( CoordinateSystem ) )
    case ( 'RECTANGULAR' )
      allocate ( Atlas_SCG_Form :: I % X )
      select type ( A  =>  I % X )
        class is ( Atlas_SCG_Form )
      call A % Initialize &
             ( CommunicatorOption = PROGRAM_HEADER % Communicator, &
               NameOption = 'X' )
      end select !-- A
    case ( 'SPHERICAL' )
      allocate ( Atlas_SCG_CC_Form :: I % X )
      select type ( A  =>  I % X )
        class is ( Atlas_SCG_CC_Form )
      call A % Initialize &
             ( RadiusMax = 10.0_KDR, &
               RadiusCore = 10.0_KDR / 8.0_KDR, &
               CommunicatorOption = PROGRAM_HEADER % Communicator, &
               NameOption = 'X' )
      end select !-- A
    end select !-- CoordinateSystem

    call SetParameters ( IA )

    allocate ( I % Geometry_X )
    associate ( G  =>  I % Geometry_X )
    call G % Initialize &
           ( I % X, &
             DeviceMemoryOption = DeviceMemory, &
             PinnedMemoryOption = PinnedMemory, &
             DevicesCommunicateOption = DevicesCommunicate )
    end associate !-- G

    allocate ( CurrentSet_VLC_Form :: I % CurrentSet_X )
    select type ( CS  =>  I % CurrentSet_X )
      class is ( CurrentSet_VLC_Form )
    call CS % Initialize &
           ( I % Geometry_X, IA % SpeedInitial, IA % LengthInitial )
    call CS % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'OUTFLOW   ' ], iC = 1, iD = 1 )
    select case ( trim ( CoordinateSystem ) )
    case ( 'RECTANGULAR' )
      call CS % SetBoundaryConditionsFace &
             ( [ 'PERIODIC', 'PERIODIC' ], iC = 1, iD = 2 )
      call CS % SetBoundaryConditionsFace &
             ( [ 'PERIODIC', 'PERIODIC' ], iC = 1, iD = 3 )
    case ( 'SPHERICAL' )
      call CS % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'REFLECTING' ], iC = 1, iD = 2 )
      call CS % SetBoundaryConditionsFace &
             ( [ 'PERIODIC', 'PERIODIC' ], iC = 1, iD = 3 )
    end select !-- CoordinateSystem
    end select !-- CS

    call I % Initialize ( T_FinishOption = IA % T_Finish )
    I % System        =>  IA
    I % SetInitial    =>  SetInitial
    I % ShowSystem    =>  ShowSystem
    I % SetReference  =>  SetReference

    allocate ( IA % Reference )
    allocate ( IA % Difference )
    associate &
      ( CS_R  =>  IA % Reference, &
        CS_D  =>  IA % Difference, &
         G    =>  I % Geometry_X, &
         S    =>  I % Checkpoint_X )
    call CS_R % Initialize ( G, NameOption = 'Reference' )
    call CS_D % Initialize ( G, NameOption = 'Difference' )
    call CS_R % SetStream ( S )
    call CS_D % SetStream ( S )
    end associate !-- CS_R, etc.

    end associate !-- I

  end subroutine Initialize


  subroutine Evolve ( IA )

    class ( ImposedAdvection_VL_Form ), intent ( inout ) :: &
      IA

    associate ( I  =>  IA % Integrator )
    call I % Evolve ( )
    end associate !-- I
   
  end subroutine Evolve


  impure elemental subroutine Finalize ( IA )

    type ( ImposedAdvection_VL_Form ), intent ( inout ) :: &
      IA

    if ( allocated ( IA % Difference ) ) &
      deallocate ( IA % Difference )
    if ( allocated ( IA % Reference ) ) &
      deallocate ( IA % Reference )
    if ( allocated ( IA % Integrator ) ) &
      deallocate ( IA % Integrator )

    call Show ( 'Finalizing an ImposedAdvection_VL' )

  end subroutine Finalize


  subroutine SetParameters ( IA )

    class ( ImposedAdvection_VL_Form ), intent ( inout ) :: &
      IA

    associate &
      ( I  =>  IA % Integrator )
    select type ( A  =>  I % X )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS )
    associate &
      ( V_0  =>  IA % SpeedInitial, &
        L_0  =>  IA % LengthInitial, &
        N_0  =>  IA % DensityInitial, &
        D    =>  IA % DimensionFactor, &
        T_F  =>  IA % T_Finish )

    select case ( trim ( IA % CoordinateSystem ) )
    case ( 'RECTANGULAR' )
      D  =  1.0_KDR
    case ( 'SPHERICAL' )
      D  =  3.0_KDR
    end select

    select case ( trim ( IA % AdvectionType ) )
    case ( 'CONTRACTION' )

      V_0  =  -1.0_KDR
      L_0  =   C % MaxCoordinate ( 1 )  *  0.8_KDR
      call PROGRAM_HEADER % GetParameter ( V_0, 'SpeedInitial'  )
      call PROGRAM_HEADER % GetParameter ( L_0, 'LengthInitial' )

      T_F  =  log ( 0.25_KDR ** D )  *  L_0 / ( D * V_0 )

    case ( 'EXPANSION' )

      V_0  =  +1.0_KDR
      L_0  =   C % MaxCoordinate ( 1 )  *  0.2_KDR
      call PROGRAM_HEADER % GetParameter ( V_0, 'SpeedInitial'  )
      call PROGRAM_HEADER % GetParameter ( L_0, 'LengthInitial' )

      T_F  =  log ( 4.0_KDR ** D )  *  L_0 / ( D * V_0 )

   end select !-- AdvectionType

    N_0  =  1.0_KDR
    call PROGRAM_HEADER % GetParameter ( N_0, 'DensityInitial' )

    end associate !-- V_0, etc.
    end associate !-- C
    end select !-- A
    end associate !-- I

  end subroutine SetParameters


  subroutine SetInitial ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    select type ( IA  =>  I % System )
      class is ( ImposedAdvection_VL_Form )
    select type ( I )
      class is ( Integrator_CS_Form )
    select type ( CS  =>  I % CurrentSet_X )
      class is ( CurrentSet_VLC_Form )

    call SetCurrentSet_VLC ( CS, IA, T = 0.0_KDR )

    end select !-- CS
    end select !-- I
    end select !-- IA

  end subroutine SetInitial


  subroutine ShowSystem ( I )

    class ( Integrator_H_Form ), intent ( in ) :: &
      I

    select type ( IA  =>  I % System )
      class is ( ImposedAdvection_VL_Form )
    call Show ( 'ImposedAdvection_VL Parameters' )
    call Show ( IA % CoordinateSystem, 'CoordinateSystem' )
    call Show ( IA % AdvectionType,    'AdvectionType' )
    call Show ( IA % SpeedInitial, 'SpeedInitial' )
    call Show ( IA % LengthInitial, 'LengthInitial' )
    call Show ( IA % DensityInitial, 'DensityInitial' )
    end select !-- IA

  end subroutine ShowSystem


  subroutine SetReference ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    select type ( IA  =>  I % System )
      class is ( ImposedAdvection_VL_Form )
    select type ( I )
      class is ( Integrator_CS_Form )
    select type ( CS  =>  I % CurrentSet_X )
      class is ( CurrentSet_VLC_Form )
    associate &
      ( CS_R  =>  IA % Reference, &
        CS_D  =>  IA % Difference )
    
    call SetCurrentSet_VLC ( CS_R, IA, I % T )
    
    call CS_D % MultiplyAdd ( CS, CS_R, -1.0_KDR, UseDeviceOption = .false. )

    call ComputeError ( CS_D, CS_R )

    end associate !-- CS_R, etc.
    end select !-- CS
    end select !-- I
    end select !-- IA

  end subroutine SetReference


  subroutine SetCurrentSet_VLC ( CS, IA, T )

    class ( CurrentSet_VLC_Form ), intent ( inout ) :: &
      CS
    class ( ImposedAdvection_VL_Form ), intent ( in ) :: &
      IA
    real ( KDR ), intent ( in ) :: &
      T

    real ( KDR ) :: &
      L_T, &
      N_T, &
      E_T, &
      Dim   

    associate &
      ( G  =>  CS % Geometry )
    associate &
      (  GV  =>   G % Storage_GS % Value, &
        CSV  =>  CS % Storage_GS % Value )
    associate &
      ( X_1  =>   GV ( :,  G % CENTER_U_1 ), &
        N    =>  CSV ( :, CS % DENSITY_CS ) )
    associate &
      ( V_0  =>  IA % SpeedInitial, &
        L_0  =>  IA % LengthInitial, &
        N_0  =>  IA % DensityInitial, &
        D    =>  IA % DimensionFactor )

    call CS % SetVelocityLinear ( V_0, L_0 )

    E_T  =  exp ( - D * V_0 / L_0  *  T )

    N_T  =  N_0  *  E_T
    L_T  =  L_0  *  E_T ** ( - 1.0_KDR / D )

    where ( X_1  <=  L_T )
      N  =  N_T
    elsewhere
      N  =  0.0_KDR
    end where

    end associate !-- V_0, etc.
    end associate !-- X_1, etc.
    end associate !-- GV, etc.
    end associate !-- G

  end subroutine SetCurrentSet_VLC


  subroutine ComputeError ( FS_D, FS_R )

    class ( FieldSetForm ), intent ( in ) :: &
      FS_D, &
      FS_R

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iF     !-- iField
    real ( KDR ), dimension ( FS_R % nFields ) :: &
      Norm_D, &
      Norm_R
    type ( CollectiveOperation_R_Form ) :: &
      CO

    call Show ( 'Computing error' )
    call Show ( FS_R % Name, 'FieldSet' )

    associate ( nF  =>  FS_R % nFields )

    select type ( A  =>  FS_R % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS )

    associate ( Cr  =>  C % Communicator )
    call CO % Initialize &
           ( Cr, nOutgoing = [ 2 * nF ], nIncoming = [ 2 * nF ] )
    end associate !-- Cr
    
    do iS  =  1, nF
      iF  =  FS_R % iaSelected ( iS )
      associate &
        ( F_D  =>  FS_D % Storage_GS % Value ( :, iF ), &
          F_R  =>  FS_R % Storage_GS % Value ( :, iF ) )
      
      !-- proper cells only
      CO % Outgoing % Value ( iS )  &
        =  sum ( pack ( abs ( F_D ), mask = C % ProperCell ) )
      CO % Outgoing % Value ( nF + iS )  &
        =  sum ( pack ( abs ( F_R ), mask = C % ProperCell ) )

!      !-- with ghost cells
!      CO % Outgoing % Value ( iF )  &
!        =  sum ( abs ( F  -  F_R ) )

      end associate !-- F_R, etc.
    end do !-- iS
    
    call CO % Reduce ( REDUCTION % SUM )
    Norm_D  =        CO % Incoming % Value (      1 :      nF )
    Norm_R  =  max ( CO % Incoming % Value ( nF + 1 : nF + nF ), &
                     sqrt ( tiny ( 0.0_KDR ) ) )
                         
    call Show ( Norm_D  /  Norm_R , 'L1 Error' )

    end associate !-- C
    end select !-- A
    end associate !-- nF

  end subroutine ComputeError


end module ImposedAdvection_VL__Form
