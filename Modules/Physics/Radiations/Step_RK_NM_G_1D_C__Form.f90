module Step_RK_NM_G_1D_C__Form

  !-- Step_RungeKutta_NeutrinoMoments_Grey_1D_Collected__Form

  use Basics
  use Mathematics
  use Gravitations
  use Fluids
  use NeutrinoMoments_G__Form
  use Interactions_NM_G__Form
  use ImplicitDiagnostics_NM_G__Form

  implicit none
  private

  type, public, extends ( Step_RK_CS_1D_C_CS_Form ) :: Step_RK_NM_G_1D_C_Form
    real ( KDR ), dimension ( :, : ), allocatable :: &
      Residual_J_Eq, &
      Residual_N_Eq
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
      SetFieldPointers_R, &
      SetFieldPointers_KK, &
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

    integer ( KDI ) :: &
      iNM, &
      iS

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Step_RK_NM_G_1D_C'

    call S % Step_RK_CS_1D_C_CS_Form % Initialize &
           ( CS_1D, CS, NameOption, ImplicitExplicitOption, nStagesOption )

    associate &
      ( nNM  =>  S % nCurrentSets_1D, &
        nS   =>  S % nStages, &
        mII  =>  S % MaxImplicitIterations )

    allocate ( ImplicitDiagnostics_NM_G_Form  &
               :: S % ImplicitDiagnostics ( 2 : nS, 1 : nNM ) )
    do iNM  =  1, nNM
      do iS  =  2, nS
        associate &
          ( NM  =>  CS_1D ( iNM ) )
        select type ( ID  =>  S % ImplicitDiagnostics ( iS, iNM ) )
          class is ( ImplicitDiagnostics_NM_G_Form )
        call ID % Initialize &
               ( S % Atlas, NM % Name, iS, &
                 DeviceMemoryOption = NM % DeviceMemory, &
                 PinnedMemoryOption = NM % PinnedMemory, &
                 DevicesCommunicateOption = NM % DevicesCommunicate )
        end select !-- ID
        end associate !-- NM
      end do !-- iS
    end do !-- iNM

    allocate ( S % Residual_J_Eq ( mII, nNM ) )
    allocate ( S % Residual_N_Eq ( mII, nNM ) )

    end associate !-- nNM, etc.

  end subroutine Initialize_CS_1D_C_CS


  impure elemental subroutine Finalize ( S )

    type ( Step_RK_NM_G_1D_C_Form ), intent ( inout ) :: &
      S

    if ( allocated ( S % Residual_N_Eq ) ) &
      deallocate ( S % Residual_N_Eq )
    if ( allocated ( S % Residual_J_Eq ) ) &
      deallocate ( S % Residual_J_Eq )

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
      iEnergy_R, &
      iEnergy_F, &
      iNumber_R, &
      iNumber_F!, &
!       ErrorRank
! integer ( KDI ) :: &
!   iShow
    integer ( KDI ), dimension ( 3 ) :: &
      iMomentum_R, &
      iMomentum_F
    !-- Field pointers
    real ( KDR ), dimension ( : ), pointer :: &
      KK_E_E,  KK_E_S_1,  KK_E_S_2,  KK_E_S_3,  KK_E_N, & 
      KK_EB_E, KK_EB_S_1, KK_EB_S_2, KK_EB_S_3, KK_EB_N, & 
      KK_F_E,  KK_F_S_1,  KK_F_S_2,  KK_F_S_3,  KK_F_N
    real ( KDR ), dimension ( : ), pointer :: &
      J_Eq_E,  N_Eq_E, &
      J_Eq_EB, N_Eq_EB
    !-- Storage % Value pointers
    real ( KDR ), dimension ( :, : ), pointer :: &
      Y_I_F_V => null ( ), Y_I_E_V => null ( ), Y_I_EB_V => null ( )
    real ( KDR ), dimension ( :, : ), pointer :: &
      KK_F_V => null ( ), KK_E_V => null ( ), KK_EB_V => null ( )
    real ( KDR ), dimension ( :, : ), pointer :: &
      F_V => null ( )
    real ( KDR ), dimension ( :, : ), pointer :: &
      R_E_V => null ( ), R_EB_V => null ( ), &
      I_E_V => null ( ), I_EB_V => null ( )
    !-- FieldSet pointers
    class ( FieldSet_BM_Form ), pointer :: &
      Y_I_F => null ( ), Y_I_E => null ( ), Y_I_EB => null ( )
    class ( FieldSet_BM_Form ), pointer :: &
      KK_F => null ( ), KK_E => null ( ), KK_EB => null ( )
    class ( Fluid_P_HN_Form ), pointer :: &
      F_HN => null ( )
    class ( NeutrinoMoments_G_Form ), pointer :: &
      R_E => null ( ), R_EB => null ( )
    class ( Interactions_NM_G_Form ), pointer :: &
      I_E => null ( ), I_EB => null ( )

    call Show ( 'SolveUpdateImplicit', CONSOLE % INFO_5 )
    call Show ( S % Name, 'Step', CONSOLE % INFO_5 )

    associate &
      (  S_R  =>  S % Step_CS_1D ( : ), &
         S_F  =>  S % Step_CS, &
        nR    =>  S % nCurrentSets_1D )

    select type ( F  =>  S_F % CurrentSet )
    class is ( Fluid_P_HN_Form )
      F_HN   =>  F 
      Y_I_F  =>  S_F % Intermediate
       KK_F  =>  S_F % SlopeStageImplicit ( iS ) % Element
    end select !-- F

    do iR  =  1,  nR
      select type ( R  =>  S_R ( iR ) % CurrentSet )
        class is ( NeutrinoMoments_G_Form )
      select type ( I  =>  R % Interactions )
        class is ( Interactions_NM_G_Form )
      select type ( ID =>  S % ImplicitDiagnostics ( iS, iR ) )
        class is ( ImplicitDiagnostics_NM_G_Form )
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
      end select !-- ID
      end select !-- I
      end select !-- R
    end do !-- iR

    call SetBalancedIndices &
           ( R_E, F_HN, iMomentum_R, iMomentum_F, iEnergy_R, iEnergy_F, &
             iNumber_R, iNumber_F )

    do iC  =  1,  S % Atlas % nCharts
      select type ( C  =>  S % Atlas % Chart ( iC ) % Element )
        class is ( Chart_GS_Form )

      call SetStoragePointers_F &
             ( F_HN, Y_I_F,   KK_F,  iC, &
               F_V,  Y_I_F_V, KK_F_V )
      call SetStoragePointers_R &
             ( I_E,   R_E,   Y_I_E,   KK_E,  iC, &
               I_E_V, R_E_V, Y_I_E_V, KK_E_V )
      call SetStoragePointers_R &
             ( I_EB,   R_EB,   Y_I_EB,   KK_EB,  iC, &
               I_EB_V, R_EB_V, Y_I_EB_V, KK_EB_V )

      call SetFieldPointers_KK &
               ( KK_F_V, iMomentum_F, iEnergy_F, iNumber_F, &
                 KK_F_S_1, KK_F_S_2, KK_F_S_3, KK_F_E, KK_F_N )
      call SetFieldPointers_KK &
               ( KK_E_V, iMomentum_R, iEnergy_R, iNumber_R, &
                 KK_E_S_1, KK_E_S_2, KK_E_S_3, KK_E_E, KK_E_N )
      call SetFieldPointers_KK &
               ( KK_EB_V, iMomentum_R, iEnergy_R, iNumber_R, &
                 KK_EB_S_1, KK_EB_S_2, KK_EB_S_3, KK_EB_E, KK_EB_N )

      call SetFieldPointers_R ( R_E,  R_E_V,  J_Eq_E,  N_Eq_E  )
      call SetFieldPointers_R ( R_EB, R_EB_V, J_Eq_EB, N_Eq_EB )

      call SolveKernel &
             ( C % ProperCell, &
               J_Eq_E, N_Eq_E, J_Eq_EB, N_Eq_EB, &
               KK_E_E,  KK_E_S_1,  KK_E_S_2,  KK_E_S_3,  KK_E_N, & 
               KK_EB_E, KK_EB_S_1, KK_EB_S_2, KK_EB_S_3, KK_EB_N, & 
               KK_F_E,  KK_F_S_1,  KK_F_S_2,  KK_F_S_3,  KK_F_N )

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

    class ( FieldSet_BM_Form ), intent ( in ) :: &
      F, Y_I, KK
    integer ( KDI ), intent ( in ) :: &
      iC
    real ( KDR ), dimension ( :, : ), pointer, intent ( out ) :: &
      F_V, Y_I_V, KK_V

      F_V  =>    F % Storage ( iC ) % Value
    Y_I_V  =>  Y_I % Storage ( iC ) % Value
     KK_V  =>   KK % Storage ( iC ) % Value

  end subroutine SetStoragePointers_F


  subroutine SetStoragePointers_R ( I, R, Y_I, KK, iC, I_V, R_V, Y_I_V, KK_V )

    class ( FieldSet_BM_Form ), intent ( in ) :: &
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


  subroutine SetFieldPointers_KK &
               ( KK_V, iMomentum, iEnergy, iNumber, &
                 KK_S_1, KK_S_2, KK_S_3, KK_E, KK_N )
    
    real ( KDR ), dimension ( :, : ), intent ( in ), target :: &
      KK_V
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      iMomentum
    integer ( KDI ), intent ( in ) :: &
      iEnergy, &
      iNumber
    real ( KDR ), dimension ( : ), intent ( out ), pointer :: &
      KK_S_1, KK_S_2, KK_S_3, &
      KK_E, &
      KK_N

    KK_S_1  =>  KK_V ( :, iMomentum ( 1 ) ) 
    KK_S_2  =>  KK_V ( :, iMomentum ( 2 ) ) 
    KK_S_3  =>  KK_V ( :, iMomentum ( 3 ) ) 
    
    KK_E    =>  KK_V ( :, iEnergy )
    KK_N    =>  KK_V ( :, iNumber )

  end subroutine SetFieldPointers_KK


  subroutine SetFieldPointers_R ( R, R_V, J_Eq, N_Eq )

    class ( NeutrinoMoments_G_Form ), intent ( in ) :: &
      R
    real ( KDR ), dimension ( :, : ), intent ( in ), target :: &
      R_V
    real ( KDR ), dimension ( : ), intent ( out ), pointer :: &
      J_Eq, N_Eq

      J_Eq  =>  R_V ( :, R % ENERGY_DENSITY_C_EQ )
      N_Eq  =>  R_V ( :, R % NUMBER_DENSITY_C_EQ )

  end subroutine SetFieldPointers_R


  subroutine SolveKernel &
               ( ProperCell, &
                 J_Eq_E, N_Eq_E, J_Eq_EB, N_Eq_EB, &
                 KK_E_E,  KK_E_S_1,  KK_E_S_2,  KK_E_S_3,  KK_E_N, & 
                 KK_EB_E, KK_EB_S_1, KK_EB_S_2, KK_EB_S_3, KK_EB_N, & 
                 KK_F_E,  KK_F_S_1,  KK_F_S_2,  KK_F_S_3,  KK_F_N )

    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      ProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      J_Eq_E,  N_Eq_E, &
      J_Eq_EB, N_Eq_EB
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      KK_E_E,  KK_E_S_1,  KK_E_S_2,  KK_E_S_3,  KK_E_N, & 
      KK_EB_E, KK_EB_S_1, KK_EB_S_2, KK_EB_S_3, KK_EB_N, & 
      KK_F_E,  KK_F_S_1,  KK_F_S_2,  KK_F_S_3,  KK_F_N

    integer ( KDI ) :: &
      iV, &  !-- iValue
      iR, &  !-- iRelaxation
      iI, &  !-- iIteration
      nV
    real ( KDR ) :: &
      J_Eq_E_0,  N_Eq_E_0,  &  !-- upon entry
      J_Eq_EB_0, N_Eq_EB_0
!      J_Eq_P, &  !-- previous iteration
!      N_Eq_P, &
!      E_R_P, E_R_N, &  !-- previous, new
!      N_R_P, N_R_N, &
    real ( KDR ) :: &
!      dOmega, &
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

      else

        KK_E_E   ( iV )  =  0.0_KDR
        KK_E_S_1 ( iV )  =  0.0_KDR
        KK_E_S_2 ( iV )  =  0.0_KDR
        KK_E_S_3 ( iV )  =  0.0_KDR
        KK_E_N   ( iV )  =  0.0_KDR

        KK_EB_E   ( iV )  =  0.0_KDR
        KK_EB_S_1 ( iV )  =  0.0_KDR
        KK_EB_S_2 ( iV )  =  0.0_KDR
        KK_EB_S_3 ( iV )  =  0.0_KDR
        KK_EB_N   ( iV )  =  0.0_KDR

        KK_F_E   ( iV )  =  0.0_KDR
        KK_F_S_1 ( iV )  =  0.0_KDR
        KK_F_S_2 ( iV )  =  0.0_KDR
        KK_F_S_3 ( iV )  =  0.0_KDR
        KK_F_N   ( iV )  =  0.0_KDR

      end if !-- ProperCell
    end do !-- iV

  end subroutine SolveKernel


end module Step_RK_NM_G_1D_C__Form
