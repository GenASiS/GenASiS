module Measures_F_CC__Form

  use Basics
  use Mathematics
  use Gravitations
  use Fluids

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_MEASURES_D    =  8, &
      N_MEASURES_P    = 10, &
      N_MEASURES_P_HN = 11, &
      N_MEASURES      = 11

  type, public :: Measures_F_CC_Form
    integer ( KDI ) :: &
      nMeasures = N_MEASURES
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
      ElectronFraction_C
    real ( KDR ), dimension ( N_MEASURES ) :: &
      Measure
    type ( MeasuredValueForm ), dimension ( N_MEASURES ) :: &
      Unit
    character ( LDL ), dimension ( N_MEASURES ) :: &
      Name
    class ( Atlas_H_Form ), pointer :: &
      Atlas_SA => null ( )
    class ( FieldSetForm ), pointer :: &
      Geometry_SA => null ( ), &
      Fluid_SA    => null ( )
    class ( Units_F_Form ), pointer :: &
      Units_F => null ( )
  contains
    procedure, public, pass :: &
      Initialize
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
       Y
    type ( CollectiveOperation_R_Form ), allocatable, private :: &
      CO


contains


  subroutine Initialize ( M, F_SA, G_SA, A_SA, Units_F )

    class ( Measures_F_CC_Form ), intent ( inout ) :: &
      M
    class ( FieldSetForm ), intent ( in ), target :: &
      F_SA, &
      G_SA
    class ( Atlas_H_Form ), intent ( in ), target :: &
      A_SA
    class ( Units_F_Form ), intent ( in ), target :: &
      Units_F

    M % Atlas_SA     =>  A_SA
    M % Fluid_SA     =>  F_SA
    M % Geometry_SA  =>  G_SA
    M % Units_F      =>  Units_F

    M % Name = [ 'VelocityMax        ', &
                 'Radius_V_Max       ', &
                 'Baryons_V_Max      ', &
                 'Mass_V_Max         ', &
                 'BaryonDensity_V_Max', &
                 'MassDensity_V_Max  ', &
                 'BaryonDensity_C    ', &
                 'MassDensity_C      ', &
                 'Temperature_C      ', &
                 'EntropyPerBaryon_C ', &
                 'ElectronFraction_C ' ]

    associate ( UF  =>  M % Units_F )
    M % Unit  =  [ UF % Velocity_U ( 1 ), &     !-- VelocityMax
                   UF % Coordinate_PS ( 1 ), &  !-- Radius_V_Max
                   UF % Number, &               !-- Baryons_V_Max
                   UF % Mass, &                 !-- Mass_V_Max
                   UF % NumberDensity, &        !-- BaryonDensity_V_Max
                   UF % MassDensity, &          !-- MassDensity_V_Max
                   UF % NumberDensity, &        !-- BaryonDensity_C
                   UF % MassDensity, &          !-- MassDensity_C
                   UF % Temperature, &          !-- Temperature_C
                   UF % EnergyDensity &
                     /  UF % NumberDensity  &
                     /  UF % Temperature, &     !-- EntropyPerBaryon_C
                   UNIT % IDENTITY ]            !-- ElectronFraction_C
    end associate !-- UF

    select type ( F_SA )
    class is ( Fluid_P_HN_Form )
      M % nMeasures  =  N_MEASURES_P_HN
    class is ( Fluid_P_Form )
      M % nMeasures  =  N_MEASURES_P
    class is ( Fluid_D_Form )
      M % nMeasures  =  N_MEASURES_D
    class default
      call Show ( 'Fluid type not found', CONSOLE % ERROR )
      call Show ( 'Measures_F_CC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- F_SA

  end subroutine Initialize


  subroutine Compute ( M )

    class ( Measures_F_CC_Form ), intent ( inout ) :: &
      M

    integer ( KDI ) :: &
      iP, &  !-- iProcess
      iR, &  !-- iRadius
      iC, &  !-- iCell
      iM, &  !-- iMeasure
      oC, &  !-- oCell
      oI, &  !-- oIncoming
      nF
    real ( KDR ), dimension ( : ), pointer :: &
      T_P, &
      S_P, &
      Y_P
    real ( KDR ), dimension ( :, : ), pointer :: &
      Outgoing_2D, &
      Incoming_2D

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
        nP    =>  C_SA % Communicator % Size, &
         R_P  =>  G_SA_V ( :, G_SA % CENTER_U_1 ), &
        dV_P  =>  G_SA_V ( :, G_SA % VOLUME ), &
         N_P  =>  F_SA_V ( :, F_SA % BARYON_DENSITY_C ), &
         V_P  =>  F_SA_V ( :, F_SA % VELOCITY_U_1 ) )

    nF   =   4
    T_P  =>  null ( )
    S_P  =>  null ( )
    Y_P  =>  null ( )

    select type ( F_SA )
    class is ( Fluid_P_Form )
      nF   =   nF + 2
      T_P  =>  F_SA_V ( :, F_SA % TEMPERATURE )
      S_P  =>  F_SA_V ( :, F_SA % ENTROPY_PER_BARYON )
    end select !-- F_SA

    select type ( F_SA )
    class is ( Fluid_P_HN_Form )
      nF   =   nF + 1
      Y_P  =>  F_SA_V ( :, F_SA % ELECTRON_FRACTION )
    end select !-- F_SA

    !-- Gather spherically averaged density and velocity
    !   (assume decomposition in spherical shells)

    if ( .not. allocated ( CO ) ) then
      allocate ( CO )
      call CO % Initialize &
             ( C_SA % Communicator, &
               nOutgoing = [ nF * nCB ], nIncoming = [ nF * nC ] )
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
    if ( associated ( Y_P ) ) &
      Outgoing_2D ( 1 : nCB, 7 )  =  Y_P ( nGL + 1 : nGL + nCB )

    call CO % Gather ( )

    if ( .not. allocated ( R ) ) &
      allocate ( R ( nC ), dV ( nC ), N ( nC ), V ( nC ) )
    if ( associated ( T_P ) .and. .not. allocated ( T ) ) &
      allocate ( T ( nC ) )
    if ( associated ( S_P ) .and. .not. allocated ( S ) ) &
      allocate ( S ( nC ) )
    if ( associated ( Y_P ) .and. .not. allocated ( Y ) ) &
      allocate ( Y ( nC ) )
    do iP  =  0,  nP - 1
      oC  =  iP * nCB
      oI  =  oC * nF
      Incoming_2D ( 1 : nCB, 1 : nF )  &
        =>  CO % Incoming % Value ( oI + 1 : oI + nCB * nF ) 
       R ( oC + 1 : oC + nCB )  =  Incoming_2D ( 1 : nCB, 1 )
      dV ( oC + 1 : oC + nCB )  =  Incoming_2D ( 1 : nCB, 2 )
       N ( oC + 1 : oC + nCB )  =  Incoming_2D ( 1 : nCB, 3 )
       V ( oC + 1 : oC + nCB )  =  Incoming_2D ( 1 : nCB, 4 )
      if ( allocated ( T ) ) &
        T ( oC + 1 : oC + nCB )  =  Incoming_2D ( 1 : nCB, 5 )
      if ( allocated ( S ) ) &
        S ( oC + 1 : oC + nCB )  =  Incoming_2D ( 1 : nCB, 6 )
      if ( allocated ( Y ) ) &
        Y ( oC + 1 : oC + nCB )  =  Incoming_2D ( 1 : nCB, 7 )
    end do !-- iP

    associate &
      (  UF        =>  M % Units_F, &
          V_Max    =>  M % VelocityMax, &
          R_V_Max  =>  M % Radius_V_Max, &
          B_V_Max  =>  M % Baryons_V_Max, &
          M_V_Max  =>  M % Mass_V_Max, &
          N_V_Max  =>  M % BaryonDensity_V_Max, &
        Rho_V_Max  =>  M % MassDensity_V_Max, &
          N_C      =>  M % BaryonDensity_C, &
        Rho_C      =>  M % MassDensity_C, &
          T_C      =>  M % Temperature_C, &
          S_C      =>  M % EntropyPerBaryon_C, &
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

    !-- Record

    M % Measure (  1 )  =      V_Max
    M % Measure (  2 )  =    R_V_Max
    M % Measure (  3 )  =    B_V_Max
    M % Measure (  4 )  =    M_V_Max
    M % Measure (  5 )  =    N_V_Max
    M % Measure (  6 )  =  Rho_V_Max
    M % Measure (  7 )  =    N_C
    M % Measure (  8 )  =  Rho_C
    M % Measure (  9 )  =    T_C
    M % Measure ( 10 )  =    S_C
    M % Measure ( 11 )  =    Y_C

    !-- Display

    call Show ( 'Fluid_CentralCore Measures' )
    do iM  =  1,  M % nMeasures
      call Show ( M % Measure ( iM ), M % Unit ( iM ), M % Name ( iM ) )
    end do !-- iM

    !-- Cleanup

    end associate !-- UF, etc.
    end associate !-- nGL, etc.
    end associate !-- C_SA
    end select !-- A_SA
    end select !-- G_SA
    end select !-- F_SA

  end subroutine Compute


  impure elemental subroutine Finalize ( M )

    type ( Measures_F_CC_Form ), intent ( inout ) :: &
      M

    if ( allocated ( CO ) ) deallocate ( CO )

    nullify ( M % Units_F )
    nullify ( M % Fluid_SA )
    nullify ( M % Geometry_SA )
    nullify ( M % Atlas_SA )

  end subroutine Finalize


end module Measures_F_CC__Form
