module Measures_F_CC__Form

  !-- Measures_Fluid_CentralCore__Form

  use Basics
  use Mathematics
  use Gravitations
  use Fluids

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_MEASURES_D     =  8, &
      N_MEASURES_P     = 11, &
      N_MEASURES_P_HN  = 12, &
      N_MEASURES_F_MAX = 12

  type, public :: Measures_F_CC_Form
    integer ( KDI ) :: &
      N_MEASURES_D     = N_MEASURES_D, &
      N_MEASURES_P     = N_MEASURES_P, &
      N_MEASURES_P_HN  = N_MEASURES_P_HN, &
      N_MEASURES_F_MAX = N_MEASURES_F_MAX
    integer ( KDI ) :: &
      nMeasures_F
    integer ( KDI ) :: &
      nMeasures = 0
    real ( KDR ) :: &
      VelocityMax, &
      Radius_V_Max, &
      Baryons_V_Max, &
      Mass_V_Max, &
      BaryonDensity_V_Max, &
      MassDensity_V_Max, &
      BaryonDensity_C, &
      MassDensity_C, &
      Temperature_C, &
      EntropyPerBaryon_C, &
      RadiusShock, &
      ElectronFraction_C
    real ( KDR ), dimension ( : ), allocatable :: &
      Value
    real ( KDR ), dimension ( : ), allocatable:: &
      Radius
    type ( QuantityForm ), dimension ( : ), allocatable :: &
      Unit
    character ( LDL ), dimension ( : ), allocatable :: &
      Name
    class ( Atlas_H_Form ), pointer :: &
      Atlas_SA => null ( )
    class ( FieldSet_BM_Form ), pointer :: &
      Geometry_SA => null ( ), &
      Fluid_SA    => null ( )
    class ( Units_F_Form ), pointer :: &
      Units_F => null ( )
  contains
    procedure, private, pass :: &
      Initialize_F_CC
    generic, public :: &
      Initialize => Initialize_F_CC
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Measures_F_CC_Form

    real ( KDR ), dimension ( : ), allocatable, private :: &
       R, &
      dV, &
       N, &
       V, &
       T, &
       S, &
       SS, &
       Y
    type ( CollectiveOperation_R_Form ), allocatable, private :: &
      CO


contains


  subroutine Initialize_F_CC &
               ( M, F_SA, G_SA, A_SA, Units_F, nMeasuresAddOption )

    class ( Measures_F_CC_Form ), intent ( inout ) :: &
      M
    class ( FieldSet_BM_Form ), intent ( in ), target :: &
      F_SA, &
      G_SA
    class ( Atlas_H_Form ), intent ( in ), target :: &
      A_SA
    class ( Units_F_Form ), intent ( in ), target :: &
      Units_F
    integer ( KDI ), intent ( in ), optional :: &
      nMeasuresAddOption

    M % Atlas_SA     =>  A_SA
    M % Fluid_SA     =>  F_SA
    M % Geometry_SA  =>  G_SA
    M % Units_F      =>  Units_F

    select type ( F_SA )
    class is ( Fluid_P_HN_Form )
      M % nMeasures_F  =  N_MEASURES_P_HN
    class is ( Fluid_P_Form )
      M % nMeasures_F  =  N_MEASURES_P
    class is ( Fluid_D_Form )
      M % nMeasures_F  =  N_MEASURES_D
    class default
      call Show ( 'Fluid type not found', CONSOLE % ERROR )
      call Show ( 'Measures_F_CC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- F_SA

    M % nMeasures  =  M % nMeasures_F
    if ( present ( nMeasuresAddOption ) ) &
      M % nMeasures  =  M % nMeasures  +  nMeasuresAddOption

    allocate ( M % Name ( M % nMeasures ) )   
    allocate ( M % Unit ( M % nMeasures ) )
    allocate ( M % Value ( M % nMeasures ) )
   
    associate ( UF  =>  M % Units_F )

    M % Name ( 1 : M % N_MEASURES_D ) &
      = [ 'VelocityMax        ', &
          'Radius_V_Max       ', &
          'Baryons_V_Max      ', &
          'Mass_V_Max         ', &
          'BaryonDensity_V_Max', &
          'MassDensity_V_Max  ', &
          'BaryonDensity_C    ', &
          'MassDensity_C      ' ]
    M % Unit ( 1 : M % N_MEASURES_D ) &
      =  [ UF % Velocity_U ( 1 ), &     !-- VelocityMax
           UF % Coordinate_PS ( 1 ), &  !-- Radius_V_Max
           UF % Number, &               !-- Baryons_V_Max
           UF % Mass, &                 !-- Mass_V_Max
           UF % NumberDensity, &        !-- BaryonDensity_V_Max
           UF % MassDensity, &          !-- MassDensity_V_Max
           UF % NumberDensity, &        !-- BaryonDensity_C
           UF % MassDensity ]           !-- MassDensity_C

    if ( M % nMeasures_F  >  M % N_MEASURES_D ) then
      M % Name ( M % N_MEASURES_D + 1 : M % N_MEASURES_P ) &
        = [ 'Temperature_C      ', &
            'EntropyPerBaryon_C ', &
            'RadiusShock        ' ]
      M % Unit ( M % N_MEASURES_D + 1 : M % N_MEASURES_P ) &
        = [ UF % Temperature, &           !-- Temperature_C
            UF % EnergyDensity &
                /  UF % NumberDensity  &
                /  UF % Temperature,  &   !-- EntropyPerBaryon_C
            UF % Coordinate_PS ( 1 ) ]    !-- RadiusShock
    end if !-- N_MEASURES_P

    if ( M % nMeasures_F  >  M % N_MEASURES_P ) then
      M % Name ( M % N_MEASURES_P + 1 : M % N_MEASURES_P_HN ) &
        = [ 'ElectronFraction_C ' ]
      M % Unit ( M % N_MEASURES_P + 1 : M % N_MEASURES_P_HN ) &
        = [ UNIT % IDENTITY ]            !-- ElectronFraction_C
    end if !-- N_MEASURES_P_HN

    end associate !-- UF

  end subroutine Initialize_F_CC


  subroutine Compute ( M )

    class ( Measures_F_CC_Form ), intent ( inout ) :: &
      M

    integer ( KDI ) :: &
      iB, &  !-- iBrick
      iR, &  !-- iRadius
      iC, &  !-- iCell
      iM, &  !-- iMeasure
      oC, &  !-- oCell
      oI, &  !-- oIncoming
      nF
    real ( KDR ) :: &
      MachNumber, &
      EntropyShock, &
      SqrtTiny
    real ( KDR ), dimension ( : ), pointer :: &
      T_P, &
      S_P, &
      SS_P, &
      Y_P
    real ( KDR ), dimension ( :, : ), pointer :: &
      Outgoing_2D, &
      Incoming_2D
    logical ( KDL ) :: &
      Supersonic

    select type ( F_SA  =>  M % Fluid_SA )
      class is ( Fluid_D_Form )
    select type ( G_SA  =>  M % Geometry_SA )
      class is ( Gravitation_G_Form )
    select type ( A_SA  =>  M % Atlas_SA )
      class is ( Atlas_SCG_Form )
    associate &
      ( C_SA    =>  A_SA % Chart_GS, &
        G_SA_V  =>  G_SA % Storage_GS % Value, &
        F_SA_V  =>  F_SA % Storage_GS % Value, &
        M_B     =>  F_SA % BaryonMass, &
        N_Min   =>  F_SA % BaryonDensityMin )
    associate &
      ( nGL   =>  C_SA % nGhostLayers ( 1 ), &
        nC    =>  C_SA % nCells ( 1 ), &
        nCB   =>  C_SA % nCellsBrick ( 1 ), &
        nCBG  =>  C_SA % nCellsBrickGlobal ( 1 ), &
        nB    =>  C_SA % nBricks ( 1 ), &
         R_P  =>  G_SA_V ( :, G_SA % CENTER_U_1 ), &
        dV_P  =>  G_SA_V ( :, G_SA % VOLUME ), &
         N_P  =>  F_SA_V ( :, F_SA % BARYON_DENSITY_C ), &
         V_P  =>  F_SA_V ( :, F_SA % VELOCITY_U_1 ) )

    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

    !-- Gather spherically averaged fluid fields
    !   (assume decomposition in spherical shells)

    nF    =   4
     T_P  =>  null ( )
     S_P  =>  null ( )
    SS_P  =>  null ( )
     Y_P  =>  null ( )

    select type ( F_SA )
    class is ( Fluid_P_Form )
       nF   =   nF + 3
       T_P  =>  F_SA_V ( :, F_SA % TEMPERATURE )
       S_P  =>  F_SA_V ( :, F_SA % ENTROPY_PER_BARYON )
      SS_P  =>  F_SA_V ( :, F_SA % SOUND_SPEED )
    end select !-- F_SA

    select type ( F_SA )
    class is ( Fluid_P_HN_Form )
      nF   =   nF + 1
      Y_P  =>  F_SA_V ( :, F_SA % ELECTRON_FRACTION )
    end select !-- F_SA

    if ( .not. allocated ( CO ) ) then
      allocate ( CO )
      call CO % Initialize &
             ( C_SA % Communicator, &
               nOutgoing = [ nF * nCB ], nIncoming = nF * nCBG % Value )
    end if

    Outgoing_2D ( 1 : nCB, 1 : nF )  =>  CO % Outgoing % Value
    Outgoing_2D ( 1 : nCB, 1 )  =   R_P ( nGL + 1 : nGL + nCB )
    Outgoing_2D ( 1 : nCB, 2 )  =  dV_P ( nGL + 1 : nGL + nCB )
    Outgoing_2D ( 1 : nCB, 3 )  =   N_P ( nGL + 1 : nGL + nCB )
    Outgoing_2D ( 1 : nCB, 4 )  =   V_P ( nGL + 1 : nGL + nCB )
    if ( associated ( T_P ) ) &
      Outgoing_2D ( 1 : nCB, 5 )  =  T_P ( nGL + 1 : nGL + nCB )
    if ( associated ( S_P ) ) &
      Outgoing_2D ( 1 : nCB, 6 )  =  S_P ( nGL + 1 : nGL + nCB )
    if ( associated ( SS_P ) ) &
      Outgoing_2D ( 1 : nCB, 7 )  =  SS_P ( nGL + 1 : nGL + nCB )
    if ( associated ( Y_P ) ) &
      Outgoing_2D ( 1 : nCB, 8 )  =  Y_P ( nGL + 1 : nGL + nCB )

    call CO % Gather_V ( )

    if ( .not. allocated ( R ) ) &
      allocate ( R ( nC ), dV ( nC ), N ( nC ), V ( nC ) )
    if ( associated ( T_P ) .and. .not. allocated ( T ) ) &
      allocate ( T ( nC ) )
    if ( associated ( S_P ) .and. .not. allocated ( S ) ) &
      allocate ( S ( nC ) )
    if ( associated ( SS_P ) .and. .not. allocated ( SS ) ) &
      allocate ( SS ( nC ) )
    if ( associated ( Y_P ) .and. .not. allocated ( Y ) ) &
      allocate ( Y ( nC ) )
    
    do iB  =  1,  nB
      if ( iB == 1 ) then
        oC = 0
      else
        oC  = oC + nCBG % Value ( iB - 1 )
      end if
      oI  =  oC * nF
      associate ( nCBG_V => nCBG % Value ( iB ) )
      Incoming_2D ( 1 : nCBG_V, 1 : nF )  &
        =>  CO % Incoming % Value ( oI + 1 : oI + nCBG_V * nF )
       R ( oC + 1 : oC + nCBG_V )  =  Incoming_2D ( : , 1 )
      dV ( oC + 1 : oC + nCBG_V )  =  Incoming_2D ( : , 2 )
       N ( oC + 1 : oC + nCBG_V )  =  Incoming_2D ( : , 3 )
       V ( oC + 1 : oC + nCBG_V )  =  Incoming_2D ( : , 4 )
      if ( allocated ( T ) ) &
         T ( oC + 1 : oC + nCBG_V )  =  Incoming_2D ( : , 5 )
      if ( allocated ( S ) ) &
         S ( oC + 1 : oC + nCBG_V )  =  Incoming_2D ( : , 6 )
      if ( allocated ( SS ) ) &
        SS ( oC + 1 : oC + nCBG_V )  =  Incoming_2D ( : , 7 )
      if ( allocated ( Y ) ) &
         Y ( oC + 1 : oC + nCBG_V )  =  Incoming_2D ( : , 8 )
      end associate   !-- nCBG_V
    end do !-- iB

    if ( .not. allocated ( M % Radius ) ) &
      allocate ( M % Radius ( nC ) )
    M % Radius  =  R

    associate &
      (   V_Max    =>  M % VelocityMax, &
          R_V_Max  =>  M % Radius_V_Max, &
          B_V_Max  =>  M % Baryons_V_Max, &
          M_V_Max  =>  M % Mass_V_Max, &
          N_V_Max  =>  M % BaryonDensity_V_Max, &
        Rho_V_Max  =>  M % MassDensity_V_Max, &
          N_C      =>  M % BaryonDensity_C, &
        Rho_C      =>  M % MassDensity_C, &
          T_C      =>  M % Temperature_C, &
          S_C      =>  M % EntropyPerBaryon_C, &
          R_S      =>  M % RadiusShock, &
          Y_C      =>  M % ElectronFraction_C )

    !-- VelocityMax

    V_Max  =  maxval ( abs ( V ) )

    iR  =  nC
    do iC  =  nC, 1, -1
      if ( N ( iC )  >  1.01_KDR  *  N_Min ) then
        if ( abs ( V ( iC ) )  ==  V_Max ) then
          iR  =  iC
          exit
        end if
      end if
    end do !-- iC

      R_V_Max  =  R ( iR )
      B_V_Max  =  sum ( N ( : iR )  *  dV ( : iR ) )
      M_V_Max  =  M_B * B_V_Max
      N_V_Max  =  B_V_Max  /  sum ( dV ( : iR ) )
    Rho_V_Max  =  M_B * N_V_Max

    !-- Center

      N_C  =  N ( 1 )
    Rho_C  =  M_B  *  N_C
    if ( allocated ( S ) )  S_C  =  S ( 1 )
    if ( allocated ( T ) )  T_C  =  T ( 1 )
    if ( allocated ( Y ) )  Y_C  =  Y ( 1 )

    !-- Shock

    associate ( UF  =>  M % Units_F )
    EntropyShock  =  3.0_KDR  *  UF % EnergyDensity &
                                 /  UF % NumberDensity  &
                                 /  UF % Temperature 
    end associate !-- UF


    iR  =  1
    Supersonic  =  .false.

    do iC  =  nC, 1, -1
      MachNumber  =  abs ( V ( iC ) ) / max ( SS ( iC ), SqrtTiny )
      if ( MachNumber  >  1.0_KDR ) then
        Supersonic = .true.
      end if
      if ( Supersonic .and. MachNumber  <  1.0_KDR &
                      .and.   S ( iC )  >  EntropyShock ) &
      then
        iR  =  iC
        exit
      end if
    end do !-- iC

    R_S  =  R ( iR )

    !-- Record

    M % Value (  1 )  =      V_Max
    M % Value (  2 )  =    R_V_Max
    M % Value (  3 )  =    B_V_Max
    M % Value (  4 )  =    M_V_Max
    M % Value (  5 )  =    N_V_Max
    M % Value (  6 )  =  Rho_V_Max
    M % Value (  7 )  =    N_C
    M % Value (  8 )  =  Rho_C
    if ( M % nMeasures_F  >  M % N_MEASURES_D ) then
      M % Value (  9 )  =  T_C
      M % Value ( 10 )  =  S_C
      M % Value ( 11 )  =  R_S
    end if !-- N_MEASURES_P
    if ( M % nMeasures_F  >  M % N_MEASURES_P ) then
      M % Value ( 12 )  =  Y_C
    end if !-- N_MEASURES_P_HN

    !-- Display

    call Show ( 'Fluid_CentralCore Measures' )
    do iM  =  1,  M % nMeasures_F
      call Show ( M % Value ( iM ), M % Unit ( iM ), M % Name ( iM ) )
    end do !-- iM

    !-- Cleanup

    end associate !-- V_Max, etc.
    end associate !-- nGL, etc.
    end associate !-- C_SA, etc.
    end select !-- A_SA
    end select !-- G_SA
    end select !-- F_SA

  end subroutine Compute


  impure elemental subroutine Finalize ( M )

    type ( Measures_F_CC_Form ), intent ( inout ) :: &
      M

    if ( allocated ( M % Name ) ) &
      deallocate ( M % Name )
    if ( allocated ( M % Unit ) ) &
      deallocate ( M % Unit )
    if ( allocated ( M % Value ) ) &
      deallocate ( M % Value )
    if ( allocated ( M % Radius ) ) &
      deallocate ( M % Radius )

    if ( allocated ( CO ) ) &
      deallocate ( CO )    

    nullify ( M % Units_F )
    nullify ( M % Fluid_SA )
    nullify ( M % Geometry_SA )
    nullify ( M % Atlas_SA )

  end subroutine Finalize


end module Measures_F_CC__Form
