module Step_RK_NM_G_1D_C__Form

  !-- Step_RungeKutta_NeutrinoMoments_Grey_1D_Collected__Form

  use Basics
  use Mathematics
  use Gravitations
  use Fluids
  use NeutrinoMoments_G__Form
  use Interactions_NM_G__Form

  implicit none
  private

  type, public, extends ( Step_RK_CS_1D_C_CS_Form ) :: Step_RK_NM_G_1D_C_Form
    real ( KDR ), dimension ( : ), allocatable :: &
      Residual_J_Eq_E,  Residual_N_Eq_E, &
      Residual_J_Eq_EB, Residual_N_Eq_EB
  contains
    procedure, private, pass :: &
      Initialize_CS_1D_C_CS
    final :: &
      Finalize
    procedure, public, pass :: &
      SolveUpdateImplicit
  end type Step_RK_NM_G_1D_C_Form

    private :: &
      SetBalancedIndices, &
      SetStoragePointers_F, &
      SetStoragePointers_R, &
      SetFieldPointers_FS_B, &
      SetFieldPointers_F, &
      SetFieldPointers_R, &
      SetFieldPointers_I, &
      SolveKernel

contains


  subroutine Initialize_CS_1D_C_CS &
               ( S, CS_1D, CS, NameOption, ImplicitExplicitOption, &
                 nStagesOption )

    class ( Step_RK_NM_G_1D_C_Form ), intent ( inout ) :: &
      S
    class ( CurrentSetForm ), dimension ( : ), intent ( in ), target :: &
      CS_1D
    class ( CurrentSetForm ), intent ( in ), target :: &
      CS
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ImplicitExplicitOption
    integer ( KDI ), intent ( in ), optional :: &
      nStagesOption

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Step_RK_NM_G_1D_C'

    call S % Step_RK_CS_1D_C_CS_Form % Initialize &
           ( CS_1D, CS, NameOption, ImplicitExplicitOption, nStagesOption )

    associate ( mII  =>  S % MaxImplicitIterations )

    allocate ( S % Residual_J_Eq_E ( mII ) )
    allocate ( S % Residual_N_Eq_E ( mII ) )

    allocate ( S % Residual_J_Eq_EB ( mII ) )
    allocate ( S % Residual_N_Eq_EB ( mII ) )

    end associate !-- nNM, etc.

  end subroutine Initialize_CS_1D_C_CS


  impure elemental subroutine Finalize ( S )

    type ( Step_RK_NM_G_1D_C_Form ), intent ( inout ) :: &
      S

    if ( allocated ( S % Residual_N_Eq_EB ) ) &
      deallocate ( S % Residual_N_Eq_EB )
    if ( allocated ( S % Residual_J_Eq_EB ) ) &
      deallocate ( S % Residual_J_Eq_EB )

    if ( allocated ( S % Residual_N_Eq_E ) ) &
      deallocate ( S % Residual_N_Eq_E )
    if ( allocated ( S % Residual_J_Eq_E ) ) &
      deallocate ( S % Residual_J_Eq_E )

  end subroutine Finalize


  subroutine SolveUpdateImplicit  ( S, T, dT, iS )

    class ( Step_RK_NM_G_1D_C_Form ), intent ( inout ), target :: &
      S
    real ( KDR ), intent ( in ) :: &
       T, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iS

    integer ( KDI ) :: &
      iR, &  !-- iRadiation
      iC     !-- iChart
    integer ( KDI ) :: &
      iEnergy_R, iNumber_R, &
      iEnergy_F, iNumber_F
    integer ( KDI ), dimension ( 3 ) :: &
      iMomentum_R, iMomentum_F
    !-- Field pointers
    real ( KDR ), dimension ( : ), pointer :: &
      Error, nIterations, Omega, Residual
    real ( KDR ), dimension ( : ), pointer :: &
      KK_F_E,   KK_F_S_1,  KK_F_S_2,  KK_F_S_3,  KK_F_D, &
      KK_E_E ,  KK_E_S_1,  KK_E_S_2,  KK_E_S_3,  KK_E_D, & 
      KK_EB_E,  KK_EB_S_1, KK_EB_S_2, KK_EB_S_3, KK_EB_D 
    real ( KDR ), dimension ( : ), pointer :: &
      E_F_0,  S_F_1_0,  S_F_2_0,  S_F_3_0,  D_F_0, &
      E_E_0,  S_E_1_0,  S_E_2_0,  S_E_3_0,  D_E_0, &
      E_EB_0, S_EB_1_0, S_EB_2_0, S_EB_3_0, D_EB_0
    real ( KDR ), dimension ( : ), pointer :: &
      M_DD_11, M_DD_22, M_DD_33
    real ( KDR ), dimension ( : ), pointer :: &
      E_F, S_F_1, S_F_2, S_F_3, D_F
    real ( KDR ), dimension ( : ), pointer :: &
      J_E, H_E_1, H_E_2, H_E_3, N_E, &
      E_E, S_E_1, S_E_2, S_E_3, D_E, &
      J_Eq_E, N_Eq_E
    real ( KDR ), dimension ( : ), pointer :: &
      J_EB, H_EB_1, H_EB_2, H_EB_3, N_EB, &
      E_EB, S_EB_1, S_EB_2, S_EB_3, D_EB, &
      J_Eq_EB, N_Eq_EB
    real ( KDR ), dimension ( : ), pointer :: &
      Xi_J_E,  Xi_H_E,  Xi_N_E,  Chi_J_E,  Chi_H_E,  Chi_N_E, &
      Xi_J_EB, Xi_H_EB, Xi_N_EB, Chi_J_EB, Chi_H_EB, Chi_N_EB
    !-- Storage % Value pointers
    real ( KDR ), dimension ( :, : ), pointer :: &
      ID_V
    real ( KDR ), dimension ( :, : ), pointer :: &
      KK_F_V, KK_E_V, KK_EB_V
    real ( KDR ), dimension ( :, : ), pointer :: &
      Y_I_F_V, Y_I_E_V, Y_I_EB_V
    real ( KDR ), dimension ( :, : ), pointer :: &
      G_V
    real ( KDR ), dimension ( :, : ), pointer :: &
      F_V
    real ( KDR ), dimension ( :, : ), pointer :: &
      R_E_V, R_EB_V, &
      I_E_V, I_EB_V
    !-- FieldSet pointers
    class ( ImplicitDiagnosticsForm ), pointer :: &
      ID
    class ( FieldSet_BM_Form ), pointer :: &
      KK_F, KK_E, KK_EB
    class ( FieldSet_BM_Form ), pointer :: &
      Y_I_F, Y_I_E, Y_I_EB
    class ( Gravitation_N_SG_Form ), pointer :: &
      G_N
    class ( Fluid_P_HN_Form ), pointer :: &
      F_HN
    class ( NeutrinoMoments_G_Form ), pointer :: &
      R_E, R_EB
    class ( Interactions_NM_G_Form ), pointer :: &
      I_E, I_EB

integer ( KDI ) :: &
  iV

    call Show ( 'SolveUpdateImplicit', CONSOLE % INFO_5 )
    call Show ( S % Name, 'Step', CONSOLE % INFO_5 )

    associate &
      (  S_R  =>  S % Step_CS_1D ( : ), &
         S_F  =>  S % Step_CS, &
         Res_J_Eq_E   =>  S % Residual_J_Eq_E, &
         Res_N_Eq_E   =>  S % Residual_N_Eq_E, &
         Res_J_Eq_EB  =>  S % Residual_J_Eq_EB, &
         Res_N_Eq_EB  =>  S % Residual_N_Eq_EB, &
         AA   =>  S % AA ( iS ) % Value ( iS ), &
         Tol  =>  S % ImplicitTolerance, &
        mII   =>  S % MaxImplicitIterations, &
        mRI   =>  S % MaxRelaxationIterations, & 
        nR    =>  S % nCurrentSets_1D )

    !-- FieldSet pointers

    ID  =>  S % ImplicitDiagnostics ( iS )

    select type ( F  =>  S_F % CurrentSet )
    class is ( Fluid_P_HN_Form )
      F_HN   =>  F 
      Y_I_F  =>  S_F % Intermediate
       KK_F  =>  S_F % SlopeStageImplicit ( iS ) % Element
    end select !-- F

    select type ( G  =>  F_HN % Geometry )
    class is ( Gravitation_N_SG_Form )
      G_N  =>  G
    end select !-- G

    do iR  =  1,  nR
      select type ( R  =>  S_R ( iR ) % CurrentSet )
        class is ( NeutrinoMoments_G_Form )
      select type ( I  =>  R % Interactions )
        class is ( Interactions_NM_G_Form )
      select case ( trim ( R % RadiationType ) )
      case ( 'NEUTRINOS_E' )
          I_E  =>  I
          R_E  =>  R
        Y_I_E  =>  S_R ( iR ) % Intermediate
         KK_E  =>  S_R ( iR ) % SlopeStageImplicit ( iS ) % Element
      case ( 'NEUTRINOS_E_BAR' )
          I_EB  =>  I
          R_EB  =>  R
        Y_I_EB  =>  S_R ( iR ) % Intermediate
         KK_EB  =>  S_R ( iR ) % SlopeStageImplicit ( iS ) % Element
      end select !-- RadiationType
      end select !-- I
      end select !-- R
    end do !-- iR

    call SetBalancedIndices &
           ( R_E, F_HN, iMomentum_R, iMomentum_F, iEnergy_R, iEnergy_F, &
             iNumber_R, iNumber_F )

    do iC  =  1,  S % Atlas % nCharts
      select type ( C  =>  S % Atlas % Chart ( iC ) % Element )
        class is ( Chart_GS_Form )

      !-- Storage % Value pointers

      ID_V  =>   ID % Storage ( iC ) % Value

       G_V  =>  G_N % Storage ( iC ) % Value

      call SetStoragePointers_F &
             ( F_HN, Y_I_F,   KK_F,  iC, &
               F_V,  Y_I_F_V, KK_F_V )
      call SetStoragePointers_R &
             ( I_E,   R_E,   Y_I_E,   KK_E,  iC, &
               I_E_V, R_E_V, Y_I_E_V, KK_E_V )
      call SetStoragePointers_R &
             ( I_EB,   R_EB,   Y_I_EB,   KK_EB,  iC, &
               I_EB_V, R_EB_V, Y_I_EB_V, KK_EB_V )

      !-- Field pointers

            Error  =>  ID_V ( :, ID % ERROR )
      nIterations  =>  ID_V ( :, ID % N_ITERATIONS )
            Omega  =>  ID_V ( :, ID % RELAXATION )
         Residual  =>  ID_V ( :, ID % RESIDUAL )

      M_DD_11  =>  G_V ( :, G_N % METRIC_F_DD_11 )
      M_DD_22  =>  G_V ( :, G_N % METRIC_F_DD_22 )
      M_DD_33  =>  G_V ( :, G_N % METRIC_F_DD_33 )

      call SetFieldPointers_FS_B &
             ( KK_F_V, iMomentum_F, iEnergy_F, iNumber_F, &
               KK_F_S_1, KK_F_S_2, KK_F_S_3, KK_F_E, KK_F_D )
      call SetFieldPointers_FS_B &
             ( KK_E_V, iMomentum_R, iEnergy_R, iNumber_R, &
               KK_E_S_1, KK_E_S_2, KK_E_S_3, KK_E_E, KK_E_D )
      call SetFieldPointers_FS_B &
             ( KK_EB_V, iMomentum_R, iEnergy_R, iNumber_R, &
               KK_EB_S_1, KK_EB_S_2, KK_EB_S_3, KK_EB_E, KK_EB_D )

      call SetFieldPointers_FS_B &
             ( Y_I_F_V, iMomentum_F, iEnergy_F, iNumber_F, &
               S_F_1_0, S_F_2_0, S_F_3_0, E_F_0, D_F_0 )
      call SetFieldPointers_FS_B &
             ( Y_I_E_V, iMomentum_R, iEnergy_R, iNumber_R, &
               S_E_1_0, S_E_2_0, S_E_3_0, E_E_0, D_E_0 )
      call SetFieldPointers_FS_B &
             ( Y_I_EB_V, iMomentum_R, iEnergy_R, iNumber_R, &
               S_EB_1_0, S_EB_2_0, S_EB_3_0, E_EB_0, D_EB_0 )

      call SetFieldPointers_F &
             ( F_HN, F_V, E_F, S_F_1, S_F_2, S_F_3, D_F )

      call SetFieldPointers_R &
             ( R_E, R_E_V, &
               J_E, H_E_1, H_E_2, H_E_3, N_E, &
               E_E, S_E_1, S_E_2, S_E_3, D_E, J_Eq_E, N_Eq_E )
      call SetFieldPointers_R &
             ( R_EB, R_EB_V, &
               J_EB, H_EB_1, H_EB_2, H_EB_3, N_EB, &
               E_EB, S_EB_1, S_EB_2, S_EB_3, D_EB, J_Eq_EB, N_Eq_EB )

      call SetFieldPointers_I &
             ( I_E, I_E_V, &
               Xi_J_E, Xi_H_E, Xi_N_E, Chi_J_E, Chi_H_E, Chi_N_E )
      call SetFieldPointers_I &
             ( I_EB, I_EB_V, &
               Xi_J_EB, Xi_H_EB, Xi_N_EB, Chi_J_EB, Chi_H_EB, Chi_N_EB )

if ( iS == 2 ) then
  associate &
    ( Y_E  =>  S_R ( 1 ) % Solution, &
      K_E  =>  S_R ( 1 ) % SlopeStageExplicit ( iS - 1 ) % Element )
  associate &
    ( Y_E_V  =>  Y_E % Storage ( iC ) % Value, &
      K_E_V  =>  K_E % Storage ( iC ) % Value )
  associate &
    ( R    =>  G_V ( :, G_N % CENTER_U_1 ), &
      J_N  =>  Y_E_V ( :, iEnergy_R ), &
      H_N  =>  Y_E_V ( :, iMomentum_R ( 1 ) ), &
      K_E_E   =>  K_E_V ( :, iEnergy_R ), &
      K_E_S_1 =>  K_E_V ( :, iMomentum_R ( 1 ) ) )
  call Show ( '>>> Stage' )
  call Show ( iS, '>>> iS' )
  call Show ( dT, '>>> dT' )
  iV  =  40
  call Show ( '>>> Cell' )
  call Show ( iV, '>>> iV' )
  call Show ( R ( iV ), UNIT % KILOMETER, '>>> R' )
  call Show ( '>>> Before implicit solve' )
    call Show ( '>>> Energy' )
    call Show ( J_N ( iV ), '>>> J_(1)' )
    call Show ( dT * K_E_E ( iV ), '>>> dT * K_E_E_(1)' )
    call Show ( E_E_0 ( iV ), '>>> J_(1+)' )
    call Show ( J_N ( iV ) + dT * K_E_E ( iV ), '>>> J_(1+) check' )
    call Show ( K_E_E ( iV ), '>>> K_E_E' )
    call Show ( - 2. * H_N ( iV ) / R ( iV ), '>>> - 2H/R' )
    call Show ( - ( H_N ( iV + 1 ) - H_N ( iV - 1 ) ) &
                  / ( R ( iV + 1 ) - R ( iV - 1 ) ), &
                '>>> - dH/dR' )
    call Show ( - 2. * H_N ( iV ) / R ( iV ) &
                - ( H_N ( iV + 1 ) - H_N ( iV - 1 ) ) &
                  / ( R ( iV + 1 ) - R ( iV - 1 ) ), &
                '>>> - ( 2H/R + dH/dR )' )
    call Show ( '>>> Momentum' )
    call Show ( H_N ( iV ), '>>> H_(1)' )
    call Show ( dT * K_E_S_1 ( iV ), '>>> dT * K_E_S_1_(1)' )
    call Show ( S_E_1_0 ( iV ), '>>> H_(1+)' )
    call Show ( H_N ( iV ) + dT * K_E_S_1 ( iV ), '>>> H_(1+) check' )
    call Show ( K_E_S_1 ( iV ), '>>> K_E_S_1' )
    call Show ( - ( J_N ( iV + 1 ) - J_N ( iV - 1 ) ) &
                  / ( 3. * ( R ( iV + 1 ) - R ( iV - 1 ) ) ), &
                '>>> -(1/3) dJ/dR' )
  end associate !-- R, etc.
  end associate !-- Y_E_V
  end associate !-- Y_E
end if

      call SolveKernel &
             ( I_E, I_EB, R_E, R_EB, F_HN, &
               Xi_J_E, Xi_H_E, Xi_N_E, Chi_J_E, Chi_H_E, Chi_N_E, &
               Xi_J_EB, Xi_H_EB, Xi_N_EB, Chi_J_EB, Chi_H_EB, Chi_N_EB, &
               J_E, H_E_1, H_E_2, H_E_3, N_E, &
               E_E, S_E_1, S_E_2, S_E_3, D_E, J_Eq_E, N_Eq_E, &
               J_EB, H_EB_1, H_EB_2, H_EB_3, N_EB, &
               E_EB, S_EB_1, S_EB_2, S_EB_3, D_EB, J_Eq_EB, N_Eq_EB, &
               E_F, S_F_1, S_F_2, S_F_3, D_F, &
               Error, nIterations, Omega, Residual, &
               C % ProperCell, &
               E_E_0,  S_E_1_0,  S_E_2_0,  S_E_3_0,  D_E_0, &
               E_EB_0, S_EB_1_0, S_EB_2_0, S_EB_3_0, D_EB_0, &
               E_F_0,  S_F_1_0,  S_F_2_0,  S_F_3_0,  D_F_0, &
               M_DD_11, M_DD_22, M_DD_33, &
               AA, Tol, dT, mRI, mII, iC, &
               KK_E_E,  KK_E_S_1,  KK_E_S_2,  KK_E_S_3,  KK_E_D, & 
               KK_EB_E, KK_EB_S_1, KK_EB_S_2, KK_EB_S_3, KK_EB_D, & 
               KK_F_E,  KK_F_S_1,  KK_F_S_2,  KK_F_S_3,  KK_F_D, &
               Res_J_Eq_E,  Res_N_Eq_E, &
               Res_J_Eq_EB, Res_N_Eq_EB )

if ( iS == 2 ) then
  associate &
    ( K_E  =>  S_R ( 1 ) % SlopeStageExplicit ( iS - 1 ) % Element )
  associate &
    ( K_E_V  =>  K_E % Storage ( iC ) % Value )
  associate &
    ( K_E_S_1 =>  K_E_V ( :, iMomentum_R ( 1 ) ) )
  call Show ( '>>> After implicit solve' )
    call Show ( '>>> Energy' )
    call Show ( dT * KK_E_E ( iV ), '>>> dT * KK_E_E_(2)' )
    call Show ( E_E_0 ( iV ) +  dT * KK_E_E ( iV ), '>>> J_(2)' )
    call Show ( '>>> Momentum' )
    call Show ( dT * KK_E_S_1 ( iV ), '>>> dT * KK_E_S_1_(2)' )
    call Show ( S_E_1_0 ( iV ) + dT * KK_E_S_1 ( iV ), '>>> H_(2)' )
    call Show ( K_E_S_1 ( iV ) / Chi_H_E ( iV ), 'K_E_S_1 / Chi_H_E' )
    call Show ( Chi_H_E ( iV ), '>>> Chi_H_E' )
  end associate !-- K_E
  end associate !-- K_E_V
  end associate !-- K_E_S_1
end if


      class default
        call Show ( 'Chart type not recognized', CONSOLE % ERROR )
        call Show ( 'Step_RK_NM_G_1D__Form', 'module', CONSOLE % ERROR )
        call Show ( 'SolveUpdateImplicit', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- C
    end do !-- iC

    end associate !-- S_R, etc.

  end subroutine SolveUpdateImplicit


  subroutine SetBalancedIndices &
               ( R, F, iMomentum_R, iMomentum_F, iEnergy_R, iEnergy_F, &
                 iNumber_R, iNumber_F )

    class ( NeutrinoMoments_G_Form ), intent ( in ) :: &
      R
    class ( Fluid_P_HN_Form ), intent ( in ) :: &
      F
    integer ( KDI ), dimension ( 3 ), intent ( out ) :: &
      iMomentum_R, &
      iMomentum_F
    integer ( KDI ), intent ( out ) :: &
      iEnergy_R, &
      iEnergy_F, &
      iNumber_R, &
      iNumber_F

    call Search &
           ( R % iaBalanced, R % ENERGY_DENSITY_B, iEnergy_R )
    call Search &
           ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_1, iMomentum_R ( 1 ) )
    call Search &
           ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_2, iMomentum_R ( 2 ) )
    call Search &
           ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_3, iMomentum_R ( 3 ) )
    call Search &
           ( R % iaBalanced, R % NUMBER_DENSITY_B, iNumber_R )

    call Search &
           ( F % iaBalanced, F % ENERGY_DENSITY_B, iEnergy_F )
    call Search &
           ( F % iaBalanced, F % MOMENTUM_DENSITY_D_1, iMomentum_F ( 1 ) )
    call Search &
           ( F % iaBalanced, F % MOMENTUM_DENSITY_D_2, iMomentum_F ( 2 ) )
    call Search &
           ( F % iaBalanced, F % MOMENTUM_DENSITY_D_3, iMomentum_F ( 3 ) )
    call Search &
           ( F % iaBalanced, F % ELECTRON_DENSITY_B, iNumber_F )

  end subroutine SetBalancedIndices


  subroutine SetStoragePointers_F ( F, Y_I, KK, iC, F_V, Y_I_V, KK_V )

    class ( FieldSet_BM_Form ), intent ( in ), target :: &
      F, Y_I, KK
    integer ( KDI ), intent ( in ) :: &
      iC
    real ( KDR ), dimension ( :, : ), pointer, intent ( out ) :: &
      F_V, Y_I_V, KK_V

      F_V  =>    F % Storage ( iC ) % Value
    Y_I_V  =>  Y_I % Storage ( iC ) % Value
     KK_V  =>   KK % Storage ( iC ) % Value

  end subroutine SetStoragePointers_F


  subroutine SetStoragePointers_R &
               ( I, R, Y_I, KK, iC, I_V, R_V, Y_I_V, KK_V )

    class ( FieldSet_BM_Form ), intent ( in ), target :: &
      I, R, Y_I, KK
    integer ( KDI ), intent ( in ) :: &
      iC
    real ( KDR ), dimension ( :, : ), pointer, intent ( out ) :: &
      I_V, R_V, Y_I_V, KK_V

      I_V  =>    I % Storage ( iC ) % Value
      R_V  =>    R % Storage ( iC ) % Value
    Y_I_V  =>  Y_I % Storage ( iC ) % Value
     KK_V  =>   KK % Storage ( iC ) % Value

  end subroutine SetStoragePointers_R


  subroutine SetFieldPointers_FS_B &  !-- FieldSet_Balanced
               ( FS_V, iMomentum, iEnergy, iNumber, &
                 FS_S_1, FS_S_2, FS_S_3, FS_E, FS_D )
    
    real ( KDR ), dimension ( :, : ), intent ( in ), target :: &
      FS_V
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      iMomentum
    integer ( KDI ), intent ( in ) :: &
      iEnergy, &
      iNumber
    real ( KDR ), dimension ( : ), intent ( out ), pointer :: &
      FS_S_1, FS_S_2, FS_S_3, &
      FS_E, &
      FS_D

    FS_S_1  =>  FS_V ( :, iMomentum ( 1 ) ) 
    FS_S_2  =>  FS_V ( :, iMomentum ( 2 ) ) 
    FS_S_3  =>  FS_V ( :, iMomentum ( 3 ) ) 
    
    FS_E    =>  FS_V ( :, iEnergy )
    FS_D    =>  FS_V ( :, iNumber )

  end subroutine SetFieldPointers_FS_B


  subroutine SetFieldPointers_F &
               ( F, F_V, E, S_1, S_2, S_3, D )

    class ( Fluid_P_HN_Form ), intent ( in ) :: &
      F
    real ( KDR ), dimension ( :, : ), intent ( in ), target :: &
      F_V
    real ( KDR ), dimension ( : ), intent ( out ), pointer :: &
      E, S_1, S_2, S_3, D

      E     =>  F_V ( :, F % ENERGY_DENSITY_B )
      S_1   =>  F_V ( :, F % MOMENTUM_DENSITY_D_1 )
      S_2   =>  F_V ( :, F % MOMENTUM_DENSITY_D_2 )
      S_3   =>  F_V ( :, F % MOMENTUM_DENSITY_D_3 )
      D     =>  F_V ( :, F % ELECTRON_DENSITY_B )

  end subroutine SetFieldPointers_F


  subroutine SetFieldPointers_R &
               ( R, R_V, J, H_1, H_2, H_3, N, E, S_1, S_2, S_3, D, J_Eq, N_Eq )

    class ( NeutrinoMoments_G_Form ), intent ( in ) :: &
      R
    real ( KDR ), dimension ( :, : ), intent ( in ), target :: &
      R_V
    real ( KDR ), dimension ( : ), intent ( out ), pointer :: &
      J, H_1, H_2, H_3, N, &
      E, S_1, S_2, S_3, D, &
      J_Eq, N_Eq

      J     =>  R_V ( :, R % ENERGY_DENSITY_C )
      H_1   =>  R_V ( :, R % MOMENTUM_DENSITY_C_U_1 )
      H_2   =>  R_V ( :, R % MOMENTUM_DENSITY_C_U_2 )
      H_3   =>  R_V ( :, R % MOMENTUM_DENSITY_C_U_3 )
      N     =>  R_V ( :, R % NUMBER_DENSITY_C )
      E     =>  R_V ( :, R % ENERGY_DENSITY_B )
      S_1   =>  R_V ( :, R % MOMENTUM_DENSITY_B_D_1 )
      S_2   =>  R_V ( :, R % MOMENTUM_DENSITY_B_D_2 )
      S_3   =>  R_V ( :, R % MOMENTUM_DENSITY_B_D_3 )
      D     =>  R_V ( :, R % NUMBER_DENSITY_B )
      J_Eq  =>  R_V ( :, R % ENERGY_DENSITY_C_EQ )
      N_Eq  =>  R_V ( :, R % NUMBER_DENSITY_C_EQ )

  end subroutine SetFieldPointers_R


  subroutine SetFieldPointers_I &
               ( I, I_V, Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N )

    class ( Interactions_NM_G_Form ), intent ( in ) :: &
      I
    real ( KDR ), dimension ( :, : ), intent ( in ), target :: &
      I_V
    real ( KDR ), dimension ( : ), intent ( out ), pointer :: &
       Xi_J,  Xi_H,  Xi_N, &
      Chi_J, Chi_H, Chi_N

     Xi_J  =>  I_V ( :, I % EMISSIVITY_J )
     Xi_H  =>  I_V ( :, I % EMISSIVITY_H )
     Xi_N  =>  I_V ( :, I % EMISSIVITY_N )
    Chi_J  =>  I_V ( :, I % OPACITY_J )
    Chi_H  =>  I_V ( :, I % OPACITY_H )
    Chi_N  =>  I_V ( :, I % OPACITY_N )
 
  end subroutine SetFieldPointers_I


  ! subroutine SolveKernel &
  !              ( I_E, I_EB, R_E, R_EB, F_HN, &
  !                Xi_J_E, Xi_H_E, Xi_N_E, Chi_J_E, Chi_H_E, Chi_N_E, &
  !                Xi_J_EB, Xi_H_EB, Xi_N_EB, Chi_J_EB, Chi_H_EB, Chi_N_EB, &
  !                J_E, H_E_1, H_E_2, H_E_3, N_E, &
  !                E_E, S_E_1, S_E_2, S_E_3, D_E, J_Eq_E, N_Eq_E, &
  !                J_EB, H_EB_1, H_EB_2, H_EB_3, N_EB, &
  !                E_EB, S_EB_1, S_EB_2, S_EB_3, D_EB, J_Eq_EB, N_Eq_EB, &
  !                E_F, S_F_1, S_F_2, S_F_3, D_F, &
  !                Error, nIterations, Omega, Residual, &
  !                ProperCell, &
  !                E_E_0,  S_E_1_0,  S_E_2_0,  S_E_3_0,  D_E_0, &
  !                E_EB_0, S_EB_1_0, S_EB_2_0, S_EB_3_0, D_EB_0, &
  !                E_F_0,  S_F_1_0,  S_F_2_0,  S_F_3_0,  D_F_0, &
  !                M_DD_11, M_DD_22, M_DD_33, &
  !                AA, Tol, dT, mRI, mII, iC, &
  !                KK_E_E,  KK_E_S_1,  KK_E_S_2,  KK_E_S_3,  KK_E_D, & 
  !                KK_EB_E, KK_EB_S_1, KK_EB_S_2, KK_EB_S_3, KK_EB_D, & 
  !                KK_F_E,  KK_F_S_1,  KK_F_S_2,  KK_F_S_3,  KK_F_D, &
  !                Res_J_Eq_E,  Res_N_Eq_E, &
  !                Res_J_Eq_EB, Res_N_Eq_EB )

  !   class ( Interactions_NM_G_Form ), intent ( inout ) :: &
  !     I_E, I_EB
  !   class ( NeutrinoMoments_G_Form ), intent ( inout ) :: &
  !     R_E, R_EB
  !   class ( Fluid_P_HN_Form ), intent ( inout ) :: &
  !     F_HN
  !   real ( KDR ), dimension ( : ), intent ( inout ) :: &
  !     Xi_J_E,  Xi_H_E,  Xi_N_E,  Chi_J_E,  Chi_H_E,  Chi_N_E, &
  !     Xi_J_EB, Xi_H_EB, Xi_N_EB, Chi_J_EB, Chi_H_EB, Chi_N_EB
  !   real ( KDR ), dimension ( : ), intent ( inout ) :: &
  !     J_E, H_E_1, H_E_2, H_E_3, N_E, &
  !     E_E, S_E_1, S_E_2, S_E_3, D_E, J_Eq_E, N_Eq_E
  !   real ( KDR ), dimension ( : ), intent ( inout ) :: &
  !     J_EB, H_EB_1, H_EB_2, H_EB_3, N_EB, &
  !     E_EB, S_EB_1, S_EB_2, S_EB_3, D_EB, J_Eq_EB, N_Eq_EB
  !   real ( KDR ), dimension ( : ), intent ( inout ) :: &
  !     E_F, S_F_1, S_F_2, S_F_3, D_F
  !   real ( KDR ), dimension ( : ), intent ( inout ) :: &
  !     Error, nIterations, Omega, Residual
  !   logical ( KDL ), dimension ( : ), intent ( in ) :: &
  !     ProperCell
  !   real ( KDR ), dimension ( : ), intent ( in ) :: &
  !     E_E_0,  S_E_1_0,  S_E_2_0,  S_E_3_0,  D_E_0, &
  !     E_EB_0, S_EB_1_0, S_EB_2_0, S_EB_3_0, D_EB_0, &
  !     E_F_0,  S_F_1_0,  S_F_2_0,  S_F_3_0,  D_F_0
  !   real ( KDR ), dimension ( : ), intent ( in ) :: &
  !     M_DD_11, M_DD_22, M_DD_33
  !   real ( KDR ), intent ( in ) :: &
  !     AA, Tol, dT
  !   integer ( KDI ), intent ( in ) :: &
  !     mRI, &
  !     mII, &
  !     iC
  !   real ( KDR ), dimension ( : ), intent ( out ) :: &
  !     KK_E_E,  KK_E_S_1,  KK_E_S_2,  KK_E_S_3,  KK_E_D, & 
  !     KK_EB_E, KK_EB_S_1, KK_EB_S_2, KK_EB_S_3, KK_EB_D, & 
  !     KK_F_E,  KK_F_S_1,  KK_F_S_2,  KK_F_S_3,  KK_F_D
  !   real ( KDR ), dimension ( : ), intent ( out ) :: &
  !     Res_J_Eq_E,  Res_N_Eq_E, &
  !     Res_J_Eq_EB, Res_N_Eq_EB

  !   integer ( KDI ) :: &
  !     iV, &  !-- iValue
  !     iR, &  !-- iRelaxation
  !     iI, &  !-- iIteration
  !     nV
  !   real ( KDR ) :: &
  !     J_Eq_E_0,  N_Eq_E_0,  &  !-- upon entry
  !     J_Eq_EB_0, N_Eq_EB_0
  !   real ( KDR ) :: &
  !     J_Eq_E_P,  N_Eq_E_P,  &  !-- previous iteration
  !     J_Eq_EB_P, N_Eq_EB_P
  !   real ( KDR ) :: &
  !     E_E_P,  E_E_N,  D_E_P,  D_E_N, &  !-- previous, new
  !     E_EB_P, E_EB_N, D_EB_P, D_EB_N
  !   real ( KDR ) :: &
  !     dOmega, &
  !     SqrtTiny

  !   SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

  !   nV  =  size ( ProperCell )

  !   do iV = 1, nV
  !     if ( ProperCell ( iV ) ) then      

  !       !-- Iterate radiation and fluid energy and number to convergence

  !       J_Eq_E_0   =  J_Eq_E  ( iV )
  !       N_Eq_E_0   =  N_Eq_E  ( iV )

  !       J_Eq_EB_0  =  J_Eq_EB ( iV )
  !       N_Eq_EB_0  =  N_Eq_EB ( iV )

  !       dOmega  =  1.0_KDR  /  mRI
  !       if ( Omega ( iV )  ==  0.0_KDR ) then
  !         Omega ( iV )  =  1.0_KDR
  !       else
  !         Omega ( iV )  =  min ( Omega ( iV )  +  dOmega,  1.0_KDR )
  !       end if

  !       iR  =  0
  !       Relaxation: do 

  !         iR  =  iR + 1

  !         if ( iR  >  1 ) &
  !           Omega ( iV )  =  Omega ( iV )  -  dOmega

  !         if ( Omega ( iV )  <  0.99 * dOmega ) then

  !           !-- Solve fails with vanishing relaxation parameter Omegax

  !           !-- Reset

  !           E_E ( iV )  =  E_E_0 ( iV )
  !           D_E ( iV )  =  D_E_0 ( iV )
  !           call R_E % ComputeFromBalanced ( iC, iV )

  !           E_EB ( iV )  =  E_EB_0 ( iV )
  !           D_EB ( iV )  =  D_EB_0 ( iV )
  !           call R_EB % ComputeFromBalanced ( iC, iV )

  !           E_F ( iV )  =  E_F_0 ( iV )
  !           D_F ( iV )  =  D_F_0 ( iV )
  !           call F_HN % ComputeFromBalanced ( iC, iV )

  !           J_Eq_E ( iV )  =  J_Eq_E_0
  !           N_Eq_E ( iV )  =  N_Eq_E_0

  !           J_Eq_EB ( iV )  =  J_Eq_E_0
  !           N_Eq_EB ( iV )  =  N_Eq_E_0

  !           !-- Nullify updates
          
  !           KK_E_E   ( iV )  =  0.0_KDR
  !           KK_E_S_1 ( iV )  =  0.0_KDR
  !           KK_E_S_2 ( iV )  =  0.0_KDR
  !           KK_E_S_3 ( iV )  =  0.0_KDR
  !           KK_E_D   ( iV )  =  0.0_KDR

  !           KK_EB_E   ( iV )  =  0.0_KDR
  !           KK_EB_S_1 ( iV )  =  0.0_KDR
  !           KK_EB_S_2 ( iV )  =  0.0_KDR
  !           KK_EB_S_3 ( iV )  =  0.0_KDR
  !           KK_EB_D   ( iV )  =  0.0_KDR

  !           KK_F_E   ( iV )  =  0.0_KDR
  !           KK_F_S_1 ( iV )  =  0.0_KDR
  !           KK_F_S_2 ( iV )  =  0.0_KDR
  !           KK_F_S_3 ( iV )  =  0.0_KDR
  !           KK_F_D   ( iV )  =  0.0_KDR

  !           !-- Abort solve

  !           Error ( iV )  =  1.0_KDR
  !           exit Relaxation

  !         end if

  !         iI  =  0
  !         Implicit: do 

  !           iI  =  iI + 1

  !           !-- Compute interactions

  !           call I_E  % Compute ( iC, iV )
  !           call I_EB % Compute ( iC, iV )

  !           !-- For vanishing radiation initial conditions

  !           if ( J_Eq_E_0  ==  0.0_KDR )  &
  !             J_Eq_E_0  =  J_Eq_E ( iV )
  !           if ( N_Eq_E_0  ==  0.0_KDR )  &
  !             N_Eq_E_0  =  N_Eq_E ( iV )

  !           if ( J_Eq_EB_0  ==  0.0_KDR )  &
  !             J_Eq_EB_0  =  J_Eq_EB ( iV )
  !           if ( N_Eq_EB_0  ==  0.0_KDR )  &
  !             N_Eq_EB_0  =  N_Eq_EB ( iV )

  !           !-- Compute radiation energy and number updates

  !           KK_E_E ( iV )  &
  !             =  ( Xi_J_E ( iV )  -  Chi_J_E ( iV )  *  J_E ( iV ) ) &
  !                /  ( 1.0_KDR  +  Chi_J_E ( iV ) * dT )
  !           KK_E_D ( iV )  &
  !             =  ( Xi_N_E ( iV )  -  Chi_N_E ( iV )  *  N_E ( iV ) ) &
  !                /  ( 1.0_KDR  +  Chi_N_E ( iV ) * dT )

  !           KK_EB_E ( iV )  &
  !             =  ( Xi_J_EB ( iV )  -  Chi_J_EB ( iV )  *  J_EB ( iV ) ) &
  !                /  ( 1.0_KDR  +  Chi_J_EB ( iV ) * dT )
  !           KK_EB_D ( iV )  &
  !             =  ( Xi_N_EB ( iV )  -  Chi_N_EB ( iV )  *  N_EB ( iV ) ) &
  !                /  ( 1.0_KDR  +  Chi_N_EB ( iV ) * dT )

  !           !-- Previous radiation values

  !           E_E_P  =  E_E ( iV )
  !           D_E_P  =  D_E ( iV )

  !           E_EB_P  =  E_EB ( iV )
  !           D_EB_P  =  D_EB ( iV )

  !           !-- New radiation values

  !           E_E_N  =  E_E_0 ( iV )  +  dT * AA * KK_E_E ( iV )  
  !           D_E_N  =  D_E_0 ( iV )  +  dT * AA * KK_E_D ( iV )

  !           E_EB_N  =  E_EB_0 ( iV )  +  dT * AA * KK_EB_E ( iV )  
  !           D_EB_N  =  D_EB_0 ( iV )  +  dT * AA * KK_EB_D ( iV )

  !           !-- New radiation values with underrelaxation

  !           E_E ( iV )  =  ( 1.0_KDR  -  Omega ( iV ) )  *  E_E_P  &
  !                                     +  Omega ( iV )    *  E_E_N
  !           D_E ( iV )  =  ( 1.0_KDR  -  Omega ( iV ) )  *  D_E_P  &
  !                                     +  Omega ( iV )    *  D_E_N

  !           E_EB ( iV )  =  ( 1.0_KDR  -  Omega ( iV ) )  *  E_EB_P  &
  !                                      +  Omega ( iV )    *  E_EB_N
  !           D_EB ( iV )  =  ( 1.0_KDR  -  Omega ( iV ) )  *  D_EB_P  &
  !                                      +  Omega ( iV )    *  D_EB_N

  !           !-- Adjusted radiation updates

  !           KK_E_E ( iV )  =  ( E_E ( iV )  -  E_E_0 ( iV ) )  &
  !                             /  ( dT * AA )
  !           KK_E_D ( iV )  =  ( D_E ( iV )  -  D_E_0 ( iV ) )  &
  !                             /  ( dT * AA )

  !           KK_EB_E ( iV )  =  ( E_EB ( iV )  -  E_EB_0 ( iV ) )  &
  !                             /  ( dT * AA )
  !           KK_EB_D ( iV )  =  ( D_EB ( iV )  -  D_EB_0 ( iV ) )  &
  !                             /  ( dT * AA )

  !           !-- Fluid updates

  !           KK_F_E ( iV )  =  - KK_E_E ( iV )  -  KK_EB_E ( iV )
  !           KK_F_D ( iV )  =  - KK_E_D ( iV )  +  KK_EB_D ( iV )

  !           !-- New fluid values

  !           E_F ( iV )  =  E_F_0 ( iV )  +  dT * AA * KK_F_E ( iV )
  !           D_F ( iV )  =  D_F_0 ( iV )  +  dT * AA * KK_F_D ( iV )

  !           !-- If negative radiation density, abort for this value of Omega

  !           if (      E_E  ( iV )  <  0.0_KDR .or. D_E  ( iV )  <  0.0_KDR  &
  !                .or. E_EB ( iV )  <  0.0_KDR .or. D_EB ( iV )  <  0.0_KDR )  &
  !            then
  !             exit Implicit
  !           end if

  !           !-- Exit test

  !           if ( iI  >  1 ) then

  !             associate &
  !               ( dJ_Eq_E   =>  Res_J_Eq_E  ( iI ), &
  !                 dN_Eq_E   =>  Res_N_Eq_E  ( iI ), &
  !                 dJ_Eq_EB  =>  Res_J_Eq_EB ( iI ), &
  !                 dN_Eq_EB  =>  Res_N_Eq_EB ( iI ) )

  !             dJ_Eq_E  =  abs ( J_Eq_E ( iV )  -  J_Eq_E_P )  &
  !                         /  max ( abs ( J_Eq_E_0 ), SqrtTiny )
  !             dN_Eq_E  =  abs ( N_Eq_E ( iV )  -  N_Eq_E_P )  &
  !                         /  max ( abs ( N_Eq_E_0 ), SqrtTiny )

  !             dJ_Eq_EB  =  abs ( J_Eq_EB ( iV )  -  J_Eq_EB_P )  &
  !                          /  max ( abs ( J_Eq_EB_0 ), SqrtTiny )
  !             dN_Eq_EB  =  abs ( N_Eq_EB ( iV )  -  N_Eq_EB_P )  &
  !                          /  max ( abs ( N_Eq_EB_0 ), SqrtTiny )

  !             nIterations ( iV )  =  iI
  !                Residual ( iV )  =  max ( dJ_Eq_E,  dN_Eq_E, &
  !                                          dJ_Eq_EB, dN_Eq_EB )

  !             if (       dJ_Eq_E   <  Tol .and. dN_Eq_E   <  Tol  &
  !                  .and. dJ_Eq_EB  <  Tol .and. dN_Eq_EB  <  Tol )  &
  !             then
  !               Error ( iV )  =  0.0_KDR
  !               exit Relaxation
  !             else if ( iI  ==  mII ) then
  !               exit Implicit
  !             end if
  !             end associate !-- dJ_Eq_E, etc.

  !           end if !-- iI > 1 (exit test)

  !           !-- Prepare for next iteration

  !           J_Eq_E_P  =  J_Eq_E ( iV )
  !           N_Eq_E_P  =  N_Eq_E ( iV )

  !           J_Eq_EB_P  =  J_Eq_EB ( iV )
  !           N_Eq_EB_P  =  N_Eq_EB ( iV )

  !           call R_E  % ComputeFromBalanced ( iC, iV )
  !           call R_EB % ComputeFromBalanced ( iC, iV )
  !           call F_HN % ComputeFromBalanced ( iC, iV )

  !         end do Implicit

  !         !-- Reset, try again with lower relaxation parameter Omega

  !         E_E ( iV )  =  E_E_0 ( iV )
  !         D_E ( iV )  =  D_E_0 ( iV )
  !         call R_E % ComputeFromBalanced ( iC, iV )

  !         E_EB ( iV )  =  E_EB_0 ( iV )
  !         D_EB ( iV )  =  D_EB_0 ( iV )
  !         call R_EB % ComputeFromBalanced ( iC, iV )

  !         E_F ( iV )  =  E_F_0 ( iV )
  !         D_F ( iV )  =  D_F_0 ( iV )
  !         call F_HN % ComputeFromBalanced ( iC, iV )

  !         J_Eq_E ( iV )  =  J_Eq_E_0
  !         N_Eq_E ( iV )  =  N_Eq_E_0

  !         J_Eq_EB ( iV )  =  J_Eq_E_0
  !         N_Eq_EB ( iV )  =  N_Eq_E_0

  !       end do Relaxation

  !       !--- Momentum update

  !       if ( Error ( iV )  ==  0.0_KDR ) then

  !         KK_E_S_1 ( iV )  &
  !           =  ( Xi_H_E ( iV )  &
  !                   -  Chi_H_E ( iV )  *  M_DD_11 ( iV ) * H_E_1 ( iV ) ) &
  !              /  ( 1.0_KDR  +  Chi_H_E ( iV ) * dT )
  !         KK_E_S_2 ( iV )  &
  !           =  ( Xi_H_E ( iV )  &
  !                   -  Chi_H_E ( iV )  *  M_DD_22 ( iV ) * H_E_2 ( iV ) ) &
  !              /  ( 1.0_KDR  +  Chi_H_E ( iV ) * dT )
  !         KK_E_S_3 ( iV )  &
  !           =  ( Xi_H_E ( iV )  &
  !                   -  Chi_H_E ( iV )  *  M_DD_33 ( iV ) * H_E_3 ( iV ) ) &
  !              /  ( 1.0_KDR  +  Chi_H_E ( iV ) * dT )

  !         KK_EB_S_1 ( iV )  &
  !           =  ( Xi_H_EB ( iV )  &
  !                   -  Chi_H_EB ( iV )  *  M_DD_11 ( iV ) * H_EB_1 ( iV ) ) &
  !              /  ( 1.0_KDR  +  Chi_H_EB ( iV ) * dT )
  !         KK_EB_S_2 ( iV )  &
  !           =  ( Xi_H_EB ( iV )  &
  !                   -  Chi_H_EB ( iV )  *  M_DD_22 ( iV ) * H_EB_2 ( iV ) ) &
  !              /  ( 1.0_KDR  +  Chi_H_EB ( iV ) * dT )
  !         KK_EB_S_3 ( iV )  &
  !           =  ( Xi_H_EB ( iV )  &
  !                   -  Chi_H_EB ( iV )  *  M_DD_33 ( iV ) * H_EB_3 ( iV ) ) &
  !              /  ( 1.0_KDR  +  Chi_H_EB ( iV ) * dT )

  !         KK_F_S_1 ( iV )  =  - KK_E_S_1 ( iV )  -  KK_EB_S_1 ( iV )
  !         KK_F_S_2 ( iV )  =  - KK_E_S_2 ( iV )  -  KK_EB_S_2 ( iV )
  !         KK_F_S_3 ( iV )  =  - KK_E_S_3 ( iV )  -  KK_EB_S_3 ( iV )

  !       end if !-- Error = 0

  !     else !-- .not. ProperCell

  !       KK_E_E   ( iV )  =  0.0_KDR
  !       KK_E_S_1 ( iV )  =  0.0_KDR
  !       KK_E_S_2 ( iV )  =  0.0_KDR
  !       KK_E_S_3 ( iV )  =  0.0_KDR
  !       KK_E_D   ( iV )  =  0.0_KDR

  !       KK_EB_E   ( iV )  =  0.0_KDR
  !       KK_EB_S_1 ( iV )  =  0.0_KDR
  !       KK_EB_S_2 ( iV )  =  0.0_KDR
  !       KK_EB_S_3 ( iV )  =  0.0_KDR
  !       KK_EB_D   ( iV )  =  0.0_KDR

  !       KK_F_E   ( iV )  =  0.0_KDR
  !       KK_F_S_1 ( iV )  =  0.0_KDR
  !       KK_F_S_2 ( iV )  =  0.0_KDR
  !       KK_F_S_3 ( iV )  =  0.0_KDR
  !       KK_F_D   ( iV )  =  0.0_KDR

  !     end if !-- ProperCell
  !   end do !-- iV

  ! end subroutine SolveKernel


  subroutine SolveKernel &
               ( I_E, I_EB, R_E, R_EB, F_HN, &
                 Xi_J_E, Xi_H_E, Xi_N_E, Chi_J_E, Chi_H_E, Chi_N_E, &
                 Xi_J_EB, Xi_H_EB, Xi_N_EB, Chi_J_EB, Chi_H_EB, Chi_N_EB, &
                 J_E, H_E_1, H_E_2, H_E_3, N_E, &
                 E_E, S_E_1, S_E_2, S_E_3, D_E, J_Eq_E, N_Eq_E, &
                 J_EB, H_EB_1, H_EB_2, H_EB_3, N_EB, &
                 E_EB, S_EB_1, S_EB_2, S_EB_3, D_EB, J_Eq_EB, N_Eq_EB, &
                 E_F, S_F_1, S_F_2, S_F_3, D_F, &
                 Error, nIterations, Omega, Residual, &
                 ProperCell, &
                 E_E_0,  S_E_1_0,  S_E_2_0,  S_E_3_0,  D_E_0, &
                 E_EB_0, S_EB_1_0, S_EB_2_0, S_EB_3_0, D_EB_0, &
                 E_F_0,  S_F_1_0,  S_F_2_0,  S_F_3_0,  D_F_0, &
                 M_DD_11, M_DD_22, M_DD_33, &
                 AA, Tol, dT, mRI, mII, iC, &
                 KK_E_E,  KK_E_S_1,  KK_E_S_2,  KK_E_S_3,  KK_E_D, & 
                 KK_EB_E, KK_EB_S_1, KK_EB_S_2, KK_EB_S_3, KK_EB_D, & 
                 KK_F_E,  KK_F_S_1,  KK_F_S_2,  KK_F_S_3,  KK_F_D, &
                 Res_J_Eq_E,  Res_N_Eq_E, &
                 Res_J_Eq_EB, Res_N_Eq_EB )

    class ( Interactions_NM_G_Form ), intent ( inout ) :: &
      I_E, I_EB
    class ( NeutrinoMoments_G_Form ), intent ( inout ) :: &
      R_E, R_EB
    class ( Fluid_P_HN_Form ), intent ( inout ) :: &
      F_HN
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Xi_J_E,  Xi_H_E,  Xi_N_E,  Chi_J_E,  Chi_H_E,  Chi_N_E, &
      Xi_J_EB, Xi_H_EB, Xi_N_EB, Chi_J_EB, Chi_H_EB, Chi_N_EB
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      J_E, H_E_1, H_E_2, H_E_3, N_E, &
      E_E, S_E_1, S_E_2, S_E_3, D_E, J_Eq_E, N_Eq_E
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      J_EB, H_EB_1, H_EB_2, H_EB_3, N_EB, &
      E_EB, S_EB_1, S_EB_2, S_EB_3, D_EB, J_Eq_EB, N_Eq_EB
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      E_F, S_F_1, S_F_2, S_F_3, D_F
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Error, nIterations, Omega, Residual
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      ProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      E_E_0,  S_E_1_0,  S_E_2_0,  S_E_3_0,  D_E_0, &
      E_EB_0, S_EB_1_0, S_EB_2_0, S_EB_3_0, D_EB_0, &
      E_F_0,  S_F_1_0,  S_F_2_0,  S_F_3_0,  D_F_0
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M_DD_11, M_DD_22, M_DD_33
    real ( KDR ), intent ( in ) :: &
      AA, Tol, dT
    integer ( KDI ), intent ( in ) :: &
      mRI, &
      mII, &
      iC
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      KK_E_E,  KK_E_S_1,  KK_E_S_2,  KK_E_S_3,  KK_E_D, & 
      KK_EB_E, KK_EB_S_1, KK_EB_S_2, KK_EB_S_3, KK_EB_D, & 
      KK_F_E,  KK_F_S_1,  KK_F_S_2,  KK_F_S_3,  KK_F_D
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      Res_J_Eq_E,  Res_N_Eq_E, &
      Res_J_Eq_EB, Res_N_Eq_EB

    integer ( KDI ) :: &
      iV, &  !-- iValue
      iR, &  !-- iRelaxation
      iI, &  !-- iIteration
      nV
    real ( KDR ) :: &
      J_Eq_E_0,  N_Eq_E_0,  &  !-- upon entry
      J_Eq_EB_0, N_Eq_EB_0
    real ( KDR ) :: &
      J_Eq_E_P,  N_Eq_E_P,  &  !-- previous iteration
      J_Eq_EB_P, N_Eq_EB_P
    real ( KDR ) :: &
      E_E_P,  E_E_N,  D_E_P,  D_E_N, &  !-- previous, new
      E_EB_P, E_EB_N, D_EB_P, D_EB_N
    real ( KDR ) :: &
      dOmega, &
      SqrtTiny

    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

    nV  =  size ( ProperCell )

    do iV = 1, nV
      if ( ProperCell ( iV ) ) then      

        !-- Iterate radiation and fluid energy and number to convergence

        J_Eq_E_0   =  J_Eq_E  ( iV )
        N_Eq_E_0   =  N_Eq_E  ( iV )

        J_Eq_EB_0  =  J_Eq_EB ( iV )
        N_Eq_EB_0  =  N_Eq_EB ( iV )

        ! dOmega  =  1.0_KDR  /  mRI
        ! if ( Omega ( iV )  ==  0.0_KDR ) then
          Omega ( iV )  =  1.0_KDR
        ! else
        !   Omega ( iV )  =  min ( Omega ( iV )  +  dOmega,  1.0_KDR )
        ! end if

        ! iR  =  0
        ! Relaxation: do 

        !   iR  =  iR + 1

        !   if ( iR  >  1 ) &
        !     Omega ( iV )  =  Omega ( iV )  -  dOmega

        !   if ( Omega ( iV )  <  0.99 * dOmega ) then

        !     !-- Solve fails with vanishing relaxation parameter Omegax

        !     !-- Reset

        !     E_E ( iV )  =  E_E_0 ( iV )
        !     D_E ( iV )  =  D_E_0 ( iV )
        !     call R_E % ComputeFromBalanced ( iC, iV )

        !     E_EB ( iV )  =  E_EB_0 ( iV )
        !     D_EB ( iV )  =  D_EB_0 ( iV )
        !     call R_EB % ComputeFromBalanced ( iC, iV )

        !     E_F ( iV )  =  E_F_0 ( iV )
        !     D_F ( iV )  =  D_F_0 ( iV )
        !     call F_HN % ComputeFromBalanced ( iC, iV )

        !     J_Eq_E ( iV )  =  J_Eq_E_0
        !     N_Eq_E ( iV )  =  N_Eq_E_0

        !     J_Eq_EB ( iV )  =  J_Eq_E_0
        !     N_Eq_EB ( iV )  =  N_Eq_E_0

        !     !-- Nullify updates
          
        !     KK_E_E   ( iV )  =  0.0_KDR
        !     KK_E_S_1 ( iV )  =  0.0_KDR
        !     KK_E_S_2 ( iV )  =  0.0_KDR
        !     KK_E_S_3 ( iV )  =  0.0_KDR
        !     KK_E_D   ( iV )  =  0.0_KDR

        !     KK_EB_E   ( iV )  =  0.0_KDR
        !     KK_EB_S_1 ( iV )  =  0.0_KDR
        !     KK_EB_S_2 ( iV )  =  0.0_KDR
        !     KK_EB_S_3 ( iV )  =  0.0_KDR
        !     KK_EB_D   ( iV )  =  0.0_KDR

        !     KK_F_E   ( iV )  =  0.0_KDR
        !     KK_F_S_1 ( iV )  =  0.0_KDR
        !     KK_F_S_2 ( iV )  =  0.0_KDR
        !     KK_F_S_3 ( iV )  =  0.0_KDR
        !     KK_F_D   ( iV )  =  0.0_KDR

        !     !-- Abort solve

        !     Error ( iV )  =  1.0_KDR
        !     exit Relaxation

        !   end if

          iI  =  0
          Implicit: do 

            iI  =  iI + 1

            !-- Compute interactions

            call I_E  % Compute ( iC, iV )
            call I_EB % Compute ( iC, iV )

            !-- For vanishing radiation initial conditions

            if ( J_Eq_E_0  ==  0.0_KDR )  &
              J_Eq_E_0  =  J_Eq_E ( iV )
            if ( N_Eq_E_0  ==  0.0_KDR )  &
              N_Eq_E_0  =  N_Eq_E ( iV )

            if ( J_Eq_EB_0  ==  0.0_KDR )  &
              J_Eq_EB_0  =  J_Eq_EB ( iV )
            if ( N_Eq_EB_0  ==  0.0_KDR )  &
              N_Eq_EB_0  =  N_Eq_EB ( iV )

            !-- Compute radiation energy and number updates

            KK_E_E ( iV )  &
              =  ( Xi_J_E ( iV )  -  Chi_J_E ( iV )  *  E_E_0 ( iV ) ) &
                 /  ( 1.0_KDR  +  Chi_J_E ( iV ) * dT )
            KK_E_D ( iV )  &
              =  ( Xi_N_E ( iV )  -  Chi_N_E ( iV )  *  D_E_0 ( iV ) ) &
                 /  ( 1.0_KDR  +  Chi_N_E ( iV ) * dT )

            KK_EB_E ( iV )  &
              =  ( Xi_J_EB ( iV )  -  Chi_J_EB ( iV )  *  E_EB_0 ( iV ) ) &
                 /  ( 1.0_KDR  +  Chi_J_EB ( iV ) * dT )
            KK_EB_D ( iV )  &
              =  ( Xi_N_EB ( iV )  -  Chi_N_EB ( iV )  *  D_EB_0 ( iV ) ) &
                 /  ( 1.0_KDR  +  Chi_N_EB ( iV ) * dT )

            !-- Previous radiation values

            E_E_P  =  E_E ( iV )
            D_E_P  =  D_E ( iV )

            E_EB_P  =  E_EB ( iV )
            D_EB_P  =  D_EB ( iV )

            !-- New radiation values

            E_E_N  =  E_E_0 ( iV )  +  dT * AA * KK_E_E ( iV )  
            D_E_N  =  D_E_0 ( iV )  +  dT * AA * KK_E_D ( iV )

            E_EB_N  =  E_EB_0 ( iV )  +  dT * AA * KK_EB_E ( iV )  
            D_EB_N  =  D_EB_0 ( iV )  +  dT * AA * KK_EB_D ( iV )

            !-- New radiation values with underrelaxation

            E_E ( iV )  =  ( 1.0_KDR  -  Omega ( iV ) )  *  E_E_P  &
                                      +  Omega ( iV )    *  E_E_N
            D_E ( iV )  =  ( 1.0_KDR  -  Omega ( iV ) )  *  D_E_P  &
                                      +  Omega ( iV )    *  D_E_N

            E_EB ( iV )  =  ( 1.0_KDR  -  Omega ( iV ) )  *  E_EB_P  &
                                       +  Omega ( iV )    *  E_EB_N
            D_EB ( iV )  =  ( 1.0_KDR  -  Omega ( iV ) )  *  D_EB_P  &
                                       +  Omega ( iV )    *  D_EB_N

            !-- Adjusted radiation updates

            KK_E_E ( iV )  =  ( E_E ( iV )  -  E_E_0 ( iV ) )  &
                              /  ( dT * AA )
            KK_E_D ( iV )  =  ( D_E ( iV )  -  D_E_0 ( iV ) )  &
                              /  ( dT * AA )

            KK_EB_E ( iV )  =  ( E_EB ( iV )  -  E_EB_0 ( iV ) )  &
                              /  ( dT * AA )
            KK_EB_D ( iV )  =  ( D_EB ( iV )  -  D_EB_0 ( iV ) )  &
                              /  ( dT * AA )

            !-- Fluid updates

            KK_F_E ( iV )  =  - KK_E_E ( iV )  -  KK_EB_E ( iV )
            KK_F_D ( iV )  =  - KK_E_D ( iV )  +  KK_EB_D ( iV )

            !-- New fluid values

            E_F ( iV )  =  E_F_0 ( iV )  +  dT * AA * KK_F_E ( iV )
            D_F ( iV )  =  D_F_0 ( iV )  +  dT * AA * KK_F_D ( iV )

            !-- If negative radiation density, abort for this value of Omega

            if (      E_E  ( iV )  <  0.0_KDR .or. D_E  ( iV )  <  0.0_KDR  &
                 .or. E_EB ( iV )  <  0.0_KDR .or. D_EB ( iV )  <  0.0_KDR )  &
             then
              exit Implicit
            end if

            !-- Exit test

            if ( iI  >  1 ) then

              associate &
                ( dJ_Eq_E   =>  Res_J_Eq_E  ( iI ), &
                  dN_Eq_E   =>  Res_N_Eq_E  ( iI ), &
                  dJ_Eq_EB  =>  Res_J_Eq_EB ( iI ), &
                  dN_Eq_EB  =>  Res_N_Eq_EB ( iI ) )

              dJ_Eq_E  =  abs ( J_Eq_E ( iV )  -  J_Eq_E_P )  &
                          /  max ( abs ( J_Eq_E_0 ), SqrtTiny )
              dN_Eq_E  =  abs ( N_Eq_E ( iV )  -  N_Eq_E_P )  &
                          /  max ( abs ( N_Eq_E_0 ), SqrtTiny )

              dJ_Eq_EB  =  abs ( J_Eq_EB ( iV )  -  J_Eq_EB_P )  &
                           /  max ( abs ( J_Eq_EB_0 ), SqrtTiny )
              dN_Eq_EB  =  abs ( N_Eq_EB ( iV )  -  N_Eq_EB_P )  &
                           /  max ( abs ( N_Eq_EB_0 ), SqrtTiny )

              nIterations ( iV )  =  iI
                 Residual ( iV )  =  max ( dJ_Eq_E,  dN_Eq_E, &
                                           dJ_Eq_EB, dN_Eq_EB )

              if (       dJ_Eq_E   <  Tol .and. dN_Eq_E   <  Tol  &
                   .and. dJ_Eq_EB  <  Tol .and. dN_Eq_EB  <  Tol )  &
              then
                Error ( iV )  =  0.0_KDR
!                exit Relaxation
                exit Implicit
              else if ( iI  ==  mII ) then
                Error ( iV )  =  1.0_KDR
                exit Implicit
              end if
              end associate !-- dJ_Eq_E, etc.

            end if !-- iI > 1 (exit test)

            !-- Prepare for next iteration

            J_Eq_E_P  =  J_Eq_E ( iV )
            N_Eq_E_P  =  N_Eq_E ( iV )

            J_Eq_EB_P  =  J_Eq_EB ( iV )
            N_Eq_EB_P  =  N_Eq_EB ( iV )

!            call R_E  % ComputeFromBalanced ( iC, iV )
!            call R_EB % ComputeFromBalanced ( iC, iV )
            call F_HN % ComputeFromBalanced ( iC, iV )

          end do Implicit

        !   !-- Reset, try again with lower relaxation parameter Omega

        !   E_E ( iV )  =  E_E_0 ( iV )
        !   D_E ( iV )  =  D_E_0 ( iV )
        !   call R_E % ComputeFromBalanced ( iC, iV )

        !   E_EB ( iV )  =  E_EB_0 ( iV )
        !   D_EB ( iV )  =  D_EB_0 ( iV )
        !   call R_EB % ComputeFromBalanced ( iC, iV )

        !   E_F ( iV )  =  E_F_0 ( iV )
        !   D_F ( iV )  =  D_F_0 ( iV )
        !   call F_HN % ComputeFromBalanced ( iC, iV )

        !   J_Eq_E ( iV )  =  J_Eq_E_0
        !   N_Eq_E ( iV )  =  N_Eq_E_0

        !   J_Eq_EB ( iV )  =  J_Eq_E_0
        !   N_Eq_EB ( iV )  =  N_Eq_E_0

        ! end do Relaxation

        !--- Momentum update

!        if ( Error ( iV )  ==  0.0_KDR ) then

          KK_E_S_1 ( iV )  &
            =  ( Xi_H_E ( iV )  -  Chi_H_E ( iV )  *  S_E_1_0 ( iV ) ) &
               /  ( 1.0_KDR  +  Chi_H_E ( iV ) * dT )
          KK_E_S_2 ( iV )  &
            =  ( Xi_H_E ( iV )  -  Chi_H_E ( iV )  *  S_E_2_0 ( iV ) ) &
               /  ( 1.0_KDR  +  Chi_H_E ( iV ) * dT )
          KK_E_S_3 ( iV )  &
            =  ( Xi_H_E ( iV )  -  Chi_H_E ( iV )  *  S_E_3_0 ( iV ) ) &
               /  ( 1.0_KDR  +  Chi_H_E ( iV ) * dT )

          KK_EB_S_1 ( iV )  &
            =  ( Xi_H_EB ( iV )  -  Chi_H_EB ( iV )  *  S_EB_1_0 ( iV ) ) &
               /  ( 1.0_KDR  +  Chi_H_EB ( iV ) * dT )
          KK_EB_S_2 ( iV )  &
            =  ( Xi_H_EB ( iV )  -  Chi_H_EB ( iV )  *  S_EB_2_0 ( iV ) ) &
               /  ( 1.0_KDR  +  Chi_H_EB ( iV ) * dT )
          KK_EB_S_3 ( iV )  &
            =  ( Xi_H_EB ( iV )  -  Chi_H_EB ( iV )  *  S_EB_3_0 ( iV ) ) &
               /  ( 1.0_KDR  +  Chi_H_EB ( iV ) * dT )

          KK_F_S_1 ( iV )  =  - KK_E_S_1 ( iV )  -  KK_EB_S_1 ( iV )
          KK_F_S_2 ( iV )  =  - KK_E_S_2 ( iV )  -  KK_EB_S_2 ( iV )
          KK_F_S_3 ( iV )  =  - KK_E_S_3 ( iV )  -  KK_EB_S_3 ( iV )

!        end if !-- Error = 0

      else !-- .not. ProperCell

        KK_E_E   ( iV )  =  0.0_KDR
        KK_E_S_1 ( iV )  =  0.0_KDR
        KK_E_S_2 ( iV )  =  0.0_KDR
        KK_E_S_3 ( iV )  =  0.0_KDR
        KK_E_D   ( iV )  =  0.0_KDR

        KK_EB_E   ( iV )  =  0.0_KDR
        KK_EB_S_1 ( iV )  =  0.0_KDR
        KK_EB_S_2 ( iV )  =  0.0_KDR
        KK_EB_S_3 ( iV )  =  0.0_KDR
        KK_EB_D   ( iV )  =  0.0_KDR

        KK_F_E   ( iV )  =  0.0_KDR
        KK_F_S_1 ( iV )  =  0.0_KDR
        KK_F_S_2 ( iV )  =  0.0_KDR
        KK_F_S_3 ( iV )  =  0.0_KDR
        KK_F_D   ( iV )  =  0.0_KDR

      end if !-- ProperCell
    end do !-- iV

  end subroutine SolveKernel


end module Step_RK_NM_G_1D_C__Form
