module Step_RK_NM_G__Form

  !-- Step_RungeKutta_NeutrinoMoments_Grey_Form

  use Basics
  use Mathematics
  use Gravitations
  use Fluids
  use NeutrinoMoments_G__Form
  use Interactions_NM_G__Form
  use ImplicitDiagnostics_NM_G__Form

  implicit none
  private

  type, public, extends ( Step_RK_CS_CS_Form ) :: Step_RK_NM_G_Form
    ! real ( KDR ), dimension ( : ), allocatable :: &
    !   Residual_R_E, &
    !   Residual_R_N, &
    !   Residual_F_E, &
    !   Residual_F_N
    real ( KDR ), dimension ( : ), allocatable :: &
      Residual_J_Eq, &
      Residual_N_Eq
  contains
    procedure, private, pass :: &
      Initialize_CS_CS
    final :: &
      Finalize
    procedure, public, pass :: &
      SolveUpdateImplicit
  end type Step_RK_NM_G_Form


contains


  subroutine Initialize_CS_CS &
               ( S, CS_1, CS_2, NameOption, ImplicitExplicitOption, &
                 nStagesOption )

    class ( Step_RK_NM_G_Form ), intent ( inout ) :: &
      S
    class ( CurrentSetForm ), intent ( in ), target :: &
      CS_1, &
      CS_2
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ImplicitExplicitOption
    integer ( KDI ), intent ( in ), optional :: &
      nStagesOption

    integer ( KDI ) :: &
      iS

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Step_RK_NM_G'

    call S % Step_RK_CS_CS_Form % Initialize &
           ( CS_1, CS_2, NameOption, ImplicitExplicitOption, nStagesOption )

    associate ( nS  =>  S % nStages )
    allocate ( S % ImplicitDiagnostics ( 2 : nS ) )
    do iS  =  2, nS
      allocate ( ImplicitDiagnostics_NM_G_Form &
                 :: S % ImplicitDiagnostics ( iS ) % Element )
      select type ( ID  =>  S % ImplicitDiagnostics ( iS ) % Element )
        class is ( ImplicitDiagnostics_NM_G_Form )
      call ID % Initialize &
             ( S % Atlas, CS_1 % Name, iS, &
               DeviceMemoryOption = CS_1 % DeviceMemory, &
               PinnedMemoryOption = CS_1 % PinnedMemory, &
               DevicesCommunicateOption = CS_1 % DevicesCommunicate )
      end select !-- ID
    end do !-- iS
    end associate !-- nS

    ! allocate ( S % Residual_R_E ( S % MaxImplicitIterations ) )
    ! allocate ( S % Residual_R_N ( S % MaxImplicitIterations ) )
    ! allocate ( S % Residual_F_E ( S % MaxImplicitIterations ) )
    ! allocate ( S % Residual_F_N ( S % MaxImplicitIterations ) )
    ! S % Residual_R_E  =  0.0_KDR
    ! S % Residual_R_N  =  0.0_KDR
    ! S % Residual_F_E  =  0.0_KDR
    ! S % Residual_F_N  =  0.0_KDR
    allocate ( S % Residual_J_Eq ( S % MaxImplicitIterations ) )
    allocate ( S % Residual_N_Eq ( S % MaxImplicitIterations ) )

  end subroutine Initialize_CS_CS


  impure elemental subroutine Finalize ( S )

    type ( Step_RK_NM_G_Form ), intent ( inout ) :: &
      S

    ! if ( allocated ( S % Residual_F_N ) ) &
    !   deallocate ( S % Residual_F_N )
    ! if ( allocated ( S % Residual_F_E ) ) &
    !   deallocate ( S % Residual_F_E )
    ! if ( allocated ( S % Residual_R_N ) ) &
    !   deallocate ( S % Residual_R_N )
    ! if ( allocated ( S % Residual_R_E ) ) &
    !   deallocate ( S % Residual_R_E )
    if ( allocated ( S % Residual_N_Eq ) ) &
      deallocate ( S % Residual_N_Eq )
    if ( allocated ( S % Residual_J_Eq ) ) &
      deallocate ( S % Residual_J_Eq )

  end subroutine Finalize


!   subroutine SolveUpdateImplicit  ( S, T, dT, iS )

!     class ( Step_RK_NM_G_Form ), intent ( inout ) :: &
!       S
!     real ( KDR ), intent ( in ) :: &
!        T, &
!       dT
!     integer ( KDI ), intent ( in ) :: &
!       iS

!     integer ( KDI ) :: &
!       iC, &
!       iV, &
!       iI, &
!       nV
! integer ( KDI ) :: &
!   iShow
!     integer ( KDI ) :: &
!       iEnergy_R, &
!       iEnergy_F, &
!       iNumber_R, &
!       iNumber_F, &
!       ErrorRank
!     integer ( KDI ), dimension ( 3 ) :: &
!       iMomentum_R, &
!       iMomentum_F
!     real ( KDR ) :: &
!       E_R_P, &  !-- Energy_Radiation_Previous
!       N_R_P, &  !-- Number_Radiation_Previous
!       E_F_P, &  !-- Energy_Fluid_Previous
!       N_F_P, &  !-- Number_Fluid_Previous
!       SqrtTiny
!     logical ( KDL ) :: &
!       Edit_E, &
!       Edit_N

!     call Show ( 'SolveUpdateImplicit', CONSOLE % INFO_5 )
!     call Show ( S % Name, 'Step', CONSOLE % INFO_5 )

!     SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

!     associate &
!       ( S_R  =>  S % Step_CS_1, &
!         S_F  =>  S % Step_CS_2 )
!     select type ( R  =>  S_R % CurrentSet )
!       class is ( NeutrinoMoments_G_Form )
!     select type ( F  =>  S_F % CurrentSet )
!       class is ( Fluid_P_HN_Form )
!     select type ( I  =>  R % Interactions )
!       class is ( Interactions_NM_G_Form )
!     associate &
!       ( Y_I_R     =>  S_R % Intermediate, &
!         Y_I_F     =>  S_F % Intermediate, &
!          KK_R     =>  S_R % SlopeStageImplicit ( iS ) % Element, &
!          KK_F     =>  S_F % SlopeStageImplicit ( iS ) % Element, &
!          AA       =>  S_R % AA ( iS ) % Value ( iS ), &
!          Tol      =>  S   % ImplicitTolerance, &
!          Max_I    =>  S   % MaxImplicitIterations, &
!          Res_R_E  =>  S   % Residual_R_E, &
!          Res_R_N  =>  S   % Residual_R_N, &
!          Res_F_E  =>  S   % Residual_F_E, &
!          Res_F_N  =>  S   % Residual_F_N )
!     select type ( ID  =>  S % ImplicitDiagnostics ( iS ) % Element )
!       class is ( ImplicitDiagnostics_NM_G_Form )

!     call Search &
!            ( R % iaBalanced, R % ENERGY_DENSITY_B, iEnergy_R )
!     call Search &
!            ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_1, iMomentum_R ( 1 ) )
!     call Search &
!            ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_2, iMomentum_R ( 2 ) )
!     call Search &
!            ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_3, iMomentum_R ( 3 ) )
!     call Search &
!            ( R % iaBalanced, R % NUMBER_DENSITY_B, iNumber_R )

!     call Search &
!            ( F % iaBalanced, F % ENERGY_DENSITY_B, iEnergy_F )
!     call Search &
!            ( F % iaBalanced, F % MOMENTUM_DENSITY_D_1, iMomentum_F ( 1 ) )
!     call Search &
!            ( F % iaBalanced, F % MOMENTUM_DENSITY_D_2, iMomentum_F ( 2 ) )
!     call Search &
!            ( F % iaBalanced, F % MOMENTUM_DENSITY_D_3, iMomentum_F ( 3 ) )
!     call Search &
!            ( F % iaBalanced, F % ELECTRON_DENSITY_B, iNumber_F )

!     do iC  =  1,  S % Atlas % nCharts
!       select type ( C  =>  S % Atlas % Chart ( iC ) % Element )
!         class is ( Chart_GS_Form )
!       associate &
!         (     I_V  =>      I % Storage ( iC ) % Value, &
!               R_V  =>      R % Storage ( iC ) % Value, &
!               F_V  =>      F % Storage ( iC ) % Value, &
!           Y_I_R_V  =>  Y_I_R % Storage ( iC ) % Value, &
!           Y_I_F_V  =>  Y_I_F % Storage ( iC ) % Value, &
!            KK_R_V  =>   KK_R % Storage ( iC ) % Value, &
!            KK_F_V  =>   KK_F % Storage ( iC ) % Value, &
!              ID_V  =>     ID % Storage ( iC ) % Value )
!       associate &
!         (   Xi_J      =>  I_V ( :, I % EMISSIVITY_J ), &
!             Xi_H      =>  I_V ( :, I % EMISSIVITY_H ), &
!             Xi_N      =>  I_V ( :, I % EMISSIVITY_N ), &
!            Chi_J      =>  I_V ( :, I % OPACITY_J ), &
!            Chi_H      =>  I_V ( :, I % OPACITY_H ), &
!            Chi_N      =>  I_V ( :, I % OPACITY_N ), &
!              J        =>  R_V ( :, R % ENERGY_DENSITY_C ), &
!              H_1      =>  R_V ( :, R % MOMENTUM_DENSITY_C_U_1 ), &
!              H_2      =>  R_V ( :, R % MOMENTUM_DENSITY_C_U_2 ), &
!              H_3      =>  R_V ( :, R % MOMENTUM_DENSITY_C_U_3 ), &
!              N        =>  R_V ( :, R % NUMBER_DENSITY_C ), &
!              E_R      =>  R_V ( :, R % ENERGY_DENSITY_B ), &
!              S_R_1    =>  R_V ( :, R % MOMENTUM_DENSITY_B_D_1 ), &
!              S_R_2    =>  R_V ( :, R % MOMENTUM_DENSITY_B_D_2 ), &
!              S_R_3    =>  R_V ( :, R % MOMENTUM_DENSITY_B_D_3 ), &
!              N_R      =>  R_V ( :, R % NUMBER_DENSITY_B ), &
!              E_F      =>  F_V ( :, F % ENERGY_DENSITY_B ), &
!              S_F_1    =>  F_V ( :, F % MOMENTUM_DENSITY_D_1 ), &
!              S_F_2    =>  F_V ( :, F % MOMENTUM_DENSITY_D_2 ), &
!              S_F_3    =>  F_V ( :, F % MOMENTUM_DENSITY_D_3 ), &
!              N_F      =>  F_V ( :, F % ELECTRON_DENSITY_B ), &
!              E_R_0    =>  Y_I_R_V ( :, iEnergy_R ), &
!              S_R_1_0  =>  Y_I_R_V ( :, iMomentum_R ( 1 ) ), &
!              S_R_2_0  =>  Y_I_R_V ( :, iMomentum_R ( 2 ) ), &
!              S_R_3_0  =>  Y_I_R_V ( :, iMomentum_R ( 3 ) ), &
!              N_R_0    =>  Y_I_R_V ( :, iNumber_R ), &
!              E_F_0    =>  Y_I_F_V ( :, iEnergy_F ), &
!              S_F_1_0  =>  Y_I_F_V ( :, iMomentum_F ( 1 ) ), &
!              S_F_2_0  =>  Y_I_F_V ( :, iMomentum_F ( 2 ) ), &
!              S_F_3_0  =>  Y_I_F_V ( :, iMomentum_F ( 3 ) ), &
!              N_F_0    =>  Y_I_F_V ( :, iNumber_F ), &
!           KK_R_E      =>  KK_R_V ( :, iEnergy_R ), &
!           KK_R_S_1    =>  KK_R_V ( :, iMomentum_R ( 1 ) ), &
!           KK_R_S_2    =>  KK_R_V ( :, iMomentum_R ( 2 ) ), &
!           KK_R_S_3    =>  KK_R_V ( :, iMomentum_R ( 3 ) ), &
!           KK_R_N      =>  KK_R_V ( :, iNumber_R ), &
!           KK_F_E      =>  KK_F_V ( :, iEnergy_F ), &
!           KK_F_S_1    =>  KK_F_V ( :, iMomentum_F ( 1 ) ), &
!           KK_F_S_2    =>  KK_F_V ( :, iMomentum_F ( 2 ) ), &
!           KK_F_S_3    =>  KK_F_V ( :, iMomentum_F ( 3 ) ), &
!           KK_F_N      =>  KK_F_V ( :, iNumber_F ), &
!              Err      =>  ID_V ( :, ID % ERROR ), &
!              N_I      =>  ID_V ( :, ID % N_ITERATIONS ), &
!              R_Max    =>  ID_V ( :, ID % RESIDUAL_MAX ), &
!              R_R_E    =>  ID_V ( :, ID % RESIDUAL_RADIATION_ENERGY ), &
!              R_R_N    =>  ID_V ( :, ID % RESIDUAL_RADIATION_NUMBER ), &
!              R_F_E    =>  ID_V ( :, ID % RESIDUAL_FLUID_ENERGY ), &
!              R_F_N    =>  ID_V ( :, ID % RESIDUAL_FLUID_NUMBER ), &
!           ProperCell  =>  C % ProperCell )

!       select type ( G  =>  R % Geometry )
!       class is ( Gravitation_G_Form )

!         associate &
!           ( GSV  =>  G % Storage ( iC ) % Value )
!         associate &
!           ( M_DD_11  =>  GSV ( :, G % METRIC_F_DD_11 ), &
!             M_DD_22  =>  GSV ( :, G % METRIC_F_DD_22 ), &
!             M_DD_33  =>  GSV ( :, G % METRIC_F_DD_33 ) )

!         nV  =  size ( ProperCell )

! call Show ( '>>> Stage' )
! call Show ( iS, '>>> iS' )

! iShow = 3
! call Show ( iShow, '>>> iShow' )

!         do iV = 1, nV
!           if ( ProperCell ( iV ) ) then      

! if ( iV == iShow ) then
!   call Show ( E_R_0 ( iV ), '>>> E_R_0' )
!   call Show ( E_F_0 ( iV ), '>>> E_F_0' )
!   call Show ( N_R_0 ( iV ), '>>> N_R_0' )
!   call Show ( N_F_0 ( iV ), '>>> N_F_0' )
! end if

!             !-- Iterate radiation and fluid energy to convergence

!             iI  =  0
!             Err ( iV )  =  - huge ( 1.0_KDR )
!             E_R_P  =  E_R_0 ( iV )
!             N_R_P  =  N_R_0 ( iV )
!             E_F_P  =  E_F_0 ( iV )
!             N_F_P  =  N_F_0 ( iV )
!             Implicit: do 

!               iI  =  iI + 1

! if ( iV == iShow ) then
!   call Show ( '    >>> Iteration' )
!   call Show ( iI, '>>> iI' )
! end if

!               !-- Compute interactions

!               call I % Compute ( iC, iV )

!               !-- Compute energy and number updates

!               KK_R_E ( iV )  &
!                 =  ( Xi_J ( iV )  -  Chi_J ( iV )  *  J ( iV ) ) &
!                    /  ( 1.0_KDR  +  Chi_J ( iV ) * dT )
!               KK_R_N ( iV )  &
!                 =  ( Xi_N ( iV )  -  Chi_N ( iV )  *  N ( iV ) ) &
!                    /  ( 1.0_KDR  +  Chi_N ( iV ) * dT )
! if ( iV == iShow ) then
!   call Show ( dT * AA * KK_R_E ( iV ), '>>> dT * AA * KK_R_E' )
!   call Show ( dT * AA * KK_R_N ( iV ), '>>> dT * AA * KK_R_N' )
!   call Show ( 0.1 * E_R_P, '>>> 0.1 * E_R_P' )
!   call Show ( 0.1 * N_R_P, '>>> 0.1 * N_R_P' )
! end if

!               KK_R_E ( iV )  &
!                 =  sign ( min ( abs ( KK_R_E ( iV ) ),  &
!                                 0.1 * E_R_P  /  ( dT * AA ) ), &
!                           KK_R_E ( iV ) )
!               KK_R_N ( iV )  &
!                 =  sign ( min ( abs ( KK_R_N ( iV ) ),  &
!                                 0.1 * N_R_P  /  ( dT * AA ) ), &
!                           KK_R_N ( iV ) )
! if ( iV == iShow ) then
!   call Show ( dT * AA * KK_R_E ( iV ), '>>> dT * AA * KK_R_E (limited)' )
!   call Show ( dT * AA * KK_R_N ( iV ), '>>> dT * AA * KK_R_N (limited)' )
! end if

!               KK_F_E ( iV )  =  - KK_R_E ( iV )
!               KK_F_N ( iV )  =  - KK_R_N ( iV )

!               !-- Apply energy and number updates

!               E_R ( iV )  =  E_R_0 ( iV )  +  dT * AA * KK_R_E ( iV )
!               N_R ( iV )  =  N_R_0 ( iV )  +  dT * AA * KK_R_N ( iV )

!               E_F ( iV )  =  E_F_0 ( iV )  +  dT * AA * KK_F_E ( iV )
!               N_F ( iV )  =  N_F_0 ( iV )  +  dT * AA * KK_F_N ( iV )

! if ( iV == iShow ) then
!   call Show ( E_R ( iV ), '>>> E_R' )
!   call Show ( E_F ( iV ), '>>> E_F' )
!   call Show ( N_R ( iV ), '>>> N_R' )
!   call Show ( N_F ( iV ), '>>> N_F' )
! end if

! !               Edit_E  =  .false.
! !               if ( E_R ( iV )  <  0.0_KDR ) then
! !                 E_R ( iV )  =  E_R_0 ( iV )  &
! !                                +  ( 0.1 * E_R_P  -  E_R_0 ( iV ) ) 
! !                 E_F ( iV )  =  E_F_0 ( iV )  &
! !                                -  ( 0.1 * E_R_P  -  E_R_0 ( iV ) )
! !                 Edit_E  =  .true.
! ! if ( iV == iShow ) then
! !   call Show ( E_R ( iV ), '>>> E_R (edit)' )
! !   call Show ( E_F ( iV ), '>>> E_F (edit)' )
! ! end if
! !               end if

! !               Edit_N  =  .false.
! !               if ( N_R ( iV )  <  0.0_KDR ) then
! !                 N_R ( iV )  =  N_R_0 ( iV )  &
! !                                +  ( 0.1 * N_R_P  -  N_R_0 ( iV ) ) 
! !                 N_F ( iV )  =  N_F_0 ( iV )  &
! !                                -  ( 0.1 * N_R_P  -  N_R_0 ( iV ) ) 
! !                 Edit_N  =  .true.
! ! if ( iV == iShow ) then
! !   call Show ( N_R ( iV ), '>>> N_R (edit)' )
! !   call Show ( N_F ( iV ), '>>> N_F (edit)' )
! ! end if
! !               end if

! !               if ( Edit_E .or. Edit_N ) then
! !                 iI  =  0
! ! if ( iV == iShow ) then
! !   call Show ( N_R ( iV ), '>>> Restarting iteration' )
! ! end if
! !               end if

!               !-- Exit test

!               if ( iI  >  1 ) then

!                 associate &
!                   ( dE_R    =>  Res_R_E ( iI ), &
!                     dN_R    =>  Res_R_N ( iI ), &
!                     dE_F    =>  Res_F_E ( iI ), &
!                     dN_F    =>  Res_F_N ( iI ), &
!                     dE_R_P  =>  Res_R_E ( iI - 1 ), &
!                     dN_R_P  =>  Res_R_N ( iI - 1 ), &
!                     dE_F_P  =>  Res_F_E ( iI - 1 ), &
!                     dN_F_P  =>  Res_F_N ( iI - 1 ) )

!                 dE_R  =  abs ( E_R ( iV )  -  E_R_P )  &
!                          /  max ( abs ( E_R_0 ( iV ) ), SqrtTiny )
!                 dN_R  =  abs ( N_R ( iV )  -  N_R_P )  &
!                          /  max ( abs ( N_R_0 ( iV ) ), SqrtTiny )
!                 dE_F  =  abs ( E_F ( iV )  -  E_F_P )  &
!                          /  max ( abs ( E_F_0 ( iV ) ), SqrtTiny )
!                 dN_F  =  abs ( N_F ( iV )  -  N_F_P )  &
!                          /  max ( abs ( N_F_0 ( iV ) ), SqrtTiny )
! if ( iV == iShow ) then
!   call Show ( dE_R, '>>> dE_R' )
!   call Show ( dN_R, '>>> dN_R' )
! end if

!                 N_I   ( iV )  =  iI
!                 R_Max ( iV )  =  max ( dE_R, dN_R, dE_F, dN_F )
!                 R_R_E ( iV )  =  dE_R
!                 R_R_N ( iV )  =  dN_R
!                 R_F_E ( iV )  =  dE_F
!                 R_F_N ( iV )  =  dN_F

!                 if ( dE_R  <  Tol .and. dE_F  <  Tol  &
!                      .and. dN_R  <  Tol .and. dN_F  <  Tol ) &
!                 then
!                   Err ( iV )  =  0.0_KDR  !-- Converged
!                 else if ( iI  ==  Max_I ) then
!                   Err ( iV )  =  1.0_KDR  !-- Maximum iterations
!                   ErrorRank  =  C % Communicator % Rank
!                   call Show ( 'Implicit solve maximum iterations', &
!                               CONSOLE % ERROR )
!                   call Show ( ErrorRank, 'ErrorRank', CONSOLE % ERROR )
!                   call Show ( iS, 'iS', CONSOLE % ERROR )
!                   call Show ( iV, 'iV', CONSOLE % ERROR )
!                   call Show ( Res_R_E ( : iI ), 'Res_R_E', CONSOLE % ERROR )
!                   call Show ( Res_R_N ( : iI ), 'Res_R_N', CONSOLE % ERROR )
!                   call Show ( Res_F_E ( : iI ), 'Res_F_E', CONSOLE % ERROR )
!                   call Show ( Res_F_N ( : iI ), 'Res_F_N', CONSOLE % ERROR )
!                   call PROGRAM_HEADER % Abort ( )
!                 else if &
!                   ( iI  >  5  &
!                     ! .and. ( ( dE_R > Tol .and. dE_R_P > Tol &
!                     !           .and. dE_R > dE_R_P )  &
!                     !         .or. ( dE_F > Tol .and. dE_F_P > Tol &
!                     !                .and. dE_F > dE_F_P ) &
!                     !         .or. ( dN_R > Tol .and. dN_R_P > Tol &
!                     !                .and. dN_R > dN_R_P ) &
!                     !         .or. ( dN_F > Tol .and. dN_F_P > Tol &
!                     !                .and. dN_F > dN_F_P ) ) ) &
!                     .and. ( ( dE_R > Tol .and. dE_R_P > Tol &
!                               .and. dE_R > dE_R_P &
!                               .and. dN_R > Tol .and. dN_R_P > Tol &
!                               .and. dN_R > dN_R_P )  &
!                             .and. &
!                             ( dE_F > Tol .and. dE_F_P > Tol &
!                               .and. dE_F > dE_F_P &
!                               .and. dN_F > Tol .and. dN_F_P > Tol &
!                               .and. dN_F > dN_F_P ) ) ) &
!                 then
!                   Err ( iV )  =  2.0_KDR  !-- Diverging
!                   ErrorRank  =  C % Communicator % Rank
!                   call Show ( 'Implicit solve diverging', CONSOLE % ERROR )
!                   call Show ( ErrorRank, 'ErrorRank', CONSOLE % ERROR )
!                   call Show ( iS, 'iS', CONSOLE % ERROR )
!                   call Show ( iV, 'iV', CONSOLE % ERROR )
!                   call Show ( Res_R_E ( : iI ), 'Res_R_E', CONSOLE % ERROR )
!                   call Show ( Res_R_N ( : iI ), 'Res_R_N', CONSOLE % ERROR )
!                   call Show ( Res_F_E ( : iI ), 'Res_F_E', CONSOLE % ERROR )
!                   call Show ( Res_F_N ( : iI ), 'Res_F_N', CONSOLE % ERROR )
!                   call PROGRAM_HEADER % Abort ( )
!                 end if

!                 if ( Err ( iV )  >=  0.0_KDR ) &
!                   exit Implicit

!                 end associate !-- dE_R, etc.

!               end if !-- iI > 1 (exit test)

!               !-- Prepare for next iteration

!               E_R_P  =  E_R ( iV )
!               N_R_P  =  N_R ( iV )

!               E_F_P  =  E_F ( iV )
!               N_F_P  =  N_F ( iV )

!               call R % ComputeFromBalanced ( iC, iV )
!               call F % ComputeFromBalanced ( iC, iV )

!             end do Implicit

!             !--- Momentum update

!             if ( Err ( iV )  <  1.1_KDR ) then
!               KK_R_S_1 ( iV )  &
!                 =  ( Xi_H ( iV )  &
!                      -  Chi_H ( iV )  *  M_DD_11 ( iV ) * H_1 ( iV ) ) &
!                    /  ( 1.0_KDR  +  Chi_H ( iV ) * dT )
!               KK_R_S_2 ( iV )  &
!                 =  ( Xi_H ( iV )  &
!                      -  Chi_H ( iV )  *  M_DD_22 ( iV ) * H_2 ( iV ) ) &
!                    /  ( 1.0_KDR  +  Chi_H ( iV ) * dT )
!               KK_R_S_3 ( iV )  &
!                 =  ( Xi_H ( iV )  &
!                      -  Chi_H ( iV )  *  M_DD_33 ( iV ) * H_3 ( iV ) ) &
!                    /  ( 1.0_KDR  +  Chi_H ( iV ) * dT )
!             else
!               KK_R_S_1 ( iV )  =  0.0_KDR
!               KK_R_S_2 ( iV )  =  0.0_KDR
!               KK_R_S_3 ( iV )  =  0.0_KDR
!             end if

!             KK_F_S_1 ( iV )  =  - KK_R_S_1 ( iV )
!             KK_F_S_2 ( iV )  =  - KK_R_S_2 ( iV )
!             KK_F_S_3 ( iV )  =  - KK_R_S_3 ( iV )

!           else

!             KK_R_E   ( iV )  =  0.0_KDR
!             KK_R_S_1 ( iV )  =  0.0_KDR
!             KK_R_S_2 ( iV )  =  0.0_KDR
!             KK_R_S_3 ( iV )  =  0.0_KDR
!             KK_R_N   ( iV )  =  0.0_KDR

!             KK_F_E   ( iV )  =  0.0_KDR
!             KK_F_S_1 ( iV )  =  0.0_KDR
!             KK_F_S_2 ( iV )  =  0.0_KDR
!             KK_F_S_3 ( iV )  =  0.0_KDR
!             KK_F_N   ( iV )  =  0.0_KDR

!           end if !-- ProperCell

!         end do !-- iV

!         end associate !-- M_DD_11, etc.
!         end associate !-- GSV

!       class default
!         call Show ( 'Gravitation type not recognized', CONSOLE % ERROR )
!         call Show ( 'Step_RK_RM__Form', 'module', CONSOLE % ERROR )
!         call Show ( 'SolveUpdateImplicit', 'subroutine', CONSOLE % ERROR )
!         call PROGRAM_HEADER % Abort ( )
!       end select !-- G

!       end associate !-- Xi_J, etc.
!       end associate !-- I_V, etc.

!       class default
!         call Show ( 'Chart type not recognized', CONSOLE % ERROR )
!         call Show ( 'Step_RK_RM__Form', 'module', CONSOLE % ERROR )
!         call Show ( 'SolveUpdateImplicit', 'subroutine', CONSOLE % ERROR )
!         call PROGRAM_HEADER % Abort ( )
!       end select !-- C
!     end do !-- iC

!     end select !-- ID
!     end associate !-- KK_R, etc.
!     end select !-- I
!     end select !-- F
!     end select !-- R
!     end associate !-- S_1, etc.

!   end subroutine SolveUpdateImplicit


  subroutine SolveUpdateImplicit  ( S, T, dT, iS )

    class ( Step_RK_NM_G_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
       T, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iS

    integer ( KDI ) :: &
      iC, &  !-- iChart
      iV, &  !-- iValue
      iR, &  !-- iRelaxation
      iI, &  !-- iIteration
      nV
    integer ( KDI ) :: &
      iEnergy_R, &
      iEnergy_F, &
      iNumber_R, &
      iNumber_F, &
      ErrorRank
integer ( KDI ) :: &
  iShow
    integer ( KDI ), dimension ( 3 ) :: &
      iMomentum_R, &
      iMomentum_F
    real ( KDR ) :: &
      J_Eq_0, &  !-- upon entry
      J_Eq_P, &  !-- previous iteration
      N_Eq_0, &
      N_Eq_P, &
      E_R_P, E_R_N, &  !-- previous, new
      N_R_P, N_R_N, &
      dOmega, &
      SqrtTiny

    call Show ( 'SolveUpdateImplicit', CONSOLE % INFO_5 )
    call Show ( S % Name, 'Step', CONSOLE % INFO_5 )

    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

    associate &
      ( S_R  =>  S % Step_CS_1, &
        S_F  =>  S % Step_CS_2 )
    select type ( R  =>  S_R % CurrentSet )
      class is ( NeutrinoMoments_G_Form )
    select type ( F  =>  S_F % CurrentSet )
      class is ( Fluid_P_HN_Form )
    select type ( I  =>  R % Interactions )
      class is ( Interactions_NM_G_Form )
    associate &
      ( Y_I_R     =>  S_R % Intermediate, &
        Y_I_F     =>  S_F % Intermediate, &
         KK_R     =>  S_R % SlopeStageImplicit ( iS ) % Element, &
         KK_F     =>  S_F % SlopeStageImplicit ( iS ) % Element, &
         AA       =>  S_R % AA ( iS ) % Value ( iS ), &
         Tol      =>  S   % ImplicitTolerance, &
         Max_I    =>  S   % MaxImplicitIterations, &
         Max_R    =>  S   % MaxRelaxationIterations, &
        Res_J_Eq  =>  S   % Residual_J_Eq, & 
        Res_N_Eq  =>  S   % Residual_N_Eq ) 
    select type ( ID  =>  S % ImplicitDiagnostics ( iS ) % Element )
      class is ( ImplicitDiagnostics_NM_G_Form )

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

    do iC  =  1,  S % Atlas % nCharts
      select type ( C  =>  S % Atlas % Chart ( iC ) % Element )
        class is ( Chart_GS_Form )
      associate &
        (     I_V  =>      I % Storage ( iC ) % Value, &
              R_V  =>      R % Storage ( iC ) % Value, &
              F_V  =>      F % Storage ( iC ) % Value, &
          Y_I_R_V  =>  Y_I_R % Storage ( iC ) % Value, &
          Y_I_F_V  =>  Y_I_F % Storage ( iC ) % Value, &
           KK_R_V  =>   KK_R % Storage ( iC ) % Value, &
           KK_F_V  =>   KK_F % Storage ( iC ) % Value, &
             ID_V  =>     ID % Storage ( iC ) % Value )
      associate &
        (   Xi_J      =>  I_V ( :, I % EMISSIVITY_J ), &
            Xi_H      =>  I_V ( :, I % EMISSIVITY_H ), &
            Xi_N      =>  I_V ( :, I % EMISSIVITY_N ), &
           Chi_J      =>  I_V ( :, I % OPACITY_J ), &
           Chi_H      =>  I_V ( :, I % OPACITY_H ), &
           Chi_N      =>  I_V ( :, I % OPACITY_N ), &
             J        =>  R_V ( :, R % ENERGY_DENSITY_C ), &
             J_Eq     =>  R_V ( :, R % ENERGY_DENSITY_C_EQ ), &
             H_1      =>  R_V ( :, R % MOMENTUM_DENSITY_C_U_1 ), &
             H_2      =>  R_V ( :, R % MOMENTUM_DENSITY_C_U_2 ), &
             H_3      =>  R_V ( :, R % MOMENTUM_DENSITY_C_U_3 ), &
             N        =>  R_V ( :, R % NUMBER_DENSITY_C ), &
             N_Eq     =>  R_V ( :, R % NUMBER_DENSITY_C_EQ ), &
             E_R      =>  R_V ( :, R % ENERGY_DENSITY_B ), &
             S_R_1    =>  R_V ( :, R % MOMENTUM_DENSITY_B_D_1 ), &
             S_R_2    =>  R_V ( :, R % MOMENTUM_DENSITY_B_D_2 ), &
             S_R_3    =>  R_V ( :, R % MOMENTUM_DENSITY_B_D_3 ), &
             N_R      =>  R_V ( :, R % NUMBER_DENSITY_B ), &
             E_F      =>  F_V ( :, F % ENERGY_DENSITY_B ), &
             S_F_1    =>  F_V ( :, F % MOMENTUM_DENSITY_D_1 ), &
             S_F_2    =>  F_V ( :, F % MOMENTUM_DENSITY_D_2 ), &
             S_F_3    =>  F_V ( :, F % MOMENTUM_DENSITY_D_3 ), &
             N_F      =>  F_V ( :, F % ELECTRON_DENSITY_B ), &
             E_R_0    =>  Y_I_R_V ( :, iEnergy_R ), &
             S_R_1_0  =>  Y_I_R_V ( :, iMomentum_R ( 1 ) ), &
             S_R_2_0  =>  Y_I_R_V ( :, iMomentum_R ( 2 ) ), &
             S_R_3_0  =>  Y_I_R_V ( :, iMomentum_R ( 3 ) ), &
             N_R_0    =>  Y_I_R_V ( :, iNumber_R ), &
             E_F_0    =>  Y_I_F_V ( :, iEnergy_F ), &
             S_F_1_0  =>  Y_I_F_V ( :, iMomentum_F ( 1 ) ), &
             S_F_2_0  =>  Y_I_F_V ( :, iMomentum_F ( 2 ) ), &
             S_F_3_0  =>  Y_I_F_V ( :, iMomentum_F ( 3 ) ), &
             N_F_0    =>  Y_I_F_V ( :, iNumber_F ), &
          KK_R_E      =>  KK_R_V ( :, iEnergy_R ), &
          KK_R_S_1    =>  KK_R_V ( :, iMomentum_R ( 1 ) ), &
          KK_R_S_2    =>  KK_R_V ( :, iMomentum_R ( 2 ) ), &
          KK_R_S_3    =>  KK_R_V ( :, iMomentum_R ( 3 ) ), &
          KK_R_N      =>  KK_R_V ( :, iNumber_R ), &
          KK_F_E      =>  KK_F_V ( :, iEnergy_F ), &
          KK_F_S_1    =>  KK_F_V ( :, iMomentum_F ( 1 ) ), &
          KK_F_S_2    =>  KK_F_V ( :, iMomentum_F ( 2 ) ), &
          KK_F_S_3    =>  KK_F_V ( :, iMomentum_F ( 3 ) ), &
          KK_F_N      =>  KK_F_V ( :, iNumber_F ), &
             N_I      =>  ID_V ( :, ID % N_ITERATIONS ), &
          Omega       =>  ID_V ( :, ID % RELAXATION ), &   
             R_Max    =>  ID_V ( :, ID % RESIDUAL ), &
             R_J_Eq   =>  ID_V ( :, ID % RESIDUAL_ENERGY_EQ ), &
             R_N_Eq   =>  ID_V ( :, ID % RESIDUAL_NUMBER_EQ ), &
          ProperCell  =>  C % ProperCell )

      select type ( G  =>  R % Geometry )
      class is ( Gravitation_G_Form )

        associate &
          ( GSV  =>  G % Storage ( iC ) % Value )
        associate &
          ( M_DD_11  =>  GSV ( :, G % METRIC_F_DD_11 ), &
            M_DD_22  =>  GSV ( :, G % METRIC_F_DD_22 ), &
            M_DD_33  =>  GSV ( :, G % METRIC_F_DD_33 ) )

        nV  =  size ( ProperCell )

! call Show ( '>>> Stage' )
! call Show ( iS, '>>> iS' )

! iShow = 31
! call Show ( iShow, '>>> iShow' )

        do iV = 1, nV
          if ( ProperCell ( iV ) ) then      

            !-- Iterate radiation and fluid energy to convergence

            J_Eq_0  =  J_Eq ( iV )
            N_Eq_0  =  N_Eq ( iV )

! if ( iV == iShow ) then
!   call Show ( E_R_0 ( iV ), '>>> E_R_0' )
!   call Show ( E_F_0 ( iV ), '>>> E_F_0' )
!   call Show ( N_R_0 ( iV ), '>>> N_R_0' )
!   call Show ( N_F_0 ( iV ), '>>> N_F_0' )
!   call Show ( J_Eq_0, '>>> J_Eq_0' )
!   call Show ( N_Eq_0, '>>> N_Eq_0' )
! end if
     
            dOmega  =  1.0_KDR  /  Max_R
            if ( Omega ( iV )  ==  0.0_KDR ) then
              Omega ( iV )  =  1.0_KDR
            else
              Omega ( iV )  =  min ( Omega ( iV )  +  dOmega,  1.0_KDR )
            end if

            iR  =  0
            Relaxation: do 

              iR  =  iR + 1

              if ( iR  >  1 ) &
                Omega ( iV )  =  Omega ( iV )  -  dOmega

              if ( Omega ( iV )  <  0.99 * dOmega ) then
                ! ErrorRank  =  C % Communicator % Rank
                ! call Show ( 'Relaxation parameter too small', CONSOLE % ERROR )
                ! call Show ( iR, 'iRelaxation', CONSOLE % ERROR )
                ! call Show ( Omega ( iV ), 'Relaxation', CONSOLE % ERROR )
                ! call Show ( ErrorRank, 'ErrorRank', CONSOLE % ERROR )
                ! call Show ( iS, 'iStage', CONSOLE % ERROR )
                ! call Show ( iV, 'iValue', CONSOLE % ERROR )
!                call PROGRAM_HEADER % Abort ( )
                E_R ( iV )  =  E_R_0 ( iV )
                N_R ( iV )  =  N_R_0 ( iV )
                E_F ( iV )  =  E_F_0 ( iV )
                N_F ( iV )  =  N_F_0 ( iV )
                call R % ComputeFromBalanced ( iC, iV )
                call F % ComputeFromBalanced ( iC, iV )
                exit Relaxation
              end if

! if ( iV == iShow ) then
!   call Show ( '  >>> Relaxation loop' )
!   call Show ( iR, '>>> iR' )
!   call Show ( Omega ( iV ), '>>> Relaxation' )
! end if

              iI  =  0
              Implicit: do 

                iI  =  iI + 1

! if ( iV == iShow ) then
!   call Show ( '    >>> Iteration' )
!   call Show ( iI, '>>> iI' )
! end if

                !-- Compute interactions

                call I % Compute ( iC, iV )

                !-- For vanishing radiation initial conditions
                if ( J_Eq_0  ==  0.0_KDR )  &
                  J_Eq_0  =  J_Eq ( iV )
                if ( N_Eq_0  ==  0.0_KDR )  &
                  N_Eq_0  =  N_Eq ( iV )

! if ( iV == iShow ) then
!   call Show ( J_Eq ( iV ), '>>> J_Eq' )
!   call Show ( N_Eq ( iV ), '>>> N_Eq' )
! end if

                !-- Compute energy and number updates

                KK_R_E ( iV )  &
                  =  ( Xi_J ( iV )  -  Chi_J ( iV )  *  J ( iV ) ) &
                     /  ( 1.0_KDR  +  Chi_J ( iV ) * dT )
                KK_R_N ( iV )  &
                  =  ( Xi_N ( iV )  -  Chi_N ( iV )  *  N ( iV ) ) &
                     /  ( 1.0_KDR  +  Chi_N ( iV ) * dT )

! if ( iV == iShow ) then
!   call Show ( dT * AA * KK_R_E ( iV ), '>>> dT * AA * KK_R_E (raw)' )
!   call Show ( dT * AA * KK_R_N ( iV ), '>>> dT * AA * KK_R_N (raw)' )
! end if

                !-- Adjust and apply energy and number updates

                E_R_P  =  E_R ( iV )
                N_R_P  =  N_R ( iV )

                E_R_N  =  E_R_0 ( iV )  +  dT * AA * KK_R_E ( iV )  
                N_R_N  =  N_R_0 ( iV )  +  dT * AA * KK_R_N ( iV )

                E_R ( iV )  =  ( 1.0_KDR  -  Omega ( iV ) )  *  E_R_P  &
                                          +  Omega ( iV )    *  E_R_N
                N_R ( iV )  =  ( 1.0_KDR  -  Omega ( iV ) )  *  N_R_P  &
                                          +  Omega ( iV )    *  N_R_N

                KK_R_E ( iV )  =  ( E_R ( iV )  -  E_R_0 ( iV ) )  &
                                  /  ( dT * AA )
                KK_R_N ( iV )  =  ( N_R ( iV )  -  N_R_0 ( iV ) )  &
                                  /  ( dT * AA )

                KK_F_E ( iV )  =  - KK_R_E ( iV )
                KK_F_N ( iV )  =  - KK_R_N ( iV )

                E_F ( iV )  =  E_F_0 ( iV )  +  dT * AA * KK_F_E ( iV )
                N_F ( iV )  =  N_F_0 ( iV )  +  dT * AA * KK_F_N ( iV )

! if ( iV == iShow ) then
!   call Show ( dT * AA * KK_R_E ( iV ), '>>> dT * AA * KK_R_E (adjusted)' )
!   call Show ( dT * AA * KK_R_N ( iV ), '>>> dT * AA * KK_R_N (adjusted)' )
! end if

! if ( iV == iShow ) then
!   call Show ( E_R ( iV ), '>>> E_R' )
!   call Show ( E_F ( iV ), '>>> E_F' )
!   call Show ( N_R ( iV ), '>>> N_R' )
!   call Show ( N_F ( iV ), '>>> N_F' )
! end if

                if ( E_R ( iV )  <  0.0_KDR .or. N_R ( iV )  <  0.0_KDR ) then
                  ! if ( iV  ==  iShow ) &
                  !   call Show ( 'Implicit solve negative radiation density', &
                  !               CONSOLE % ERROR )
                  exit Implicit
                end if

                !-- Exit test

                if ( iI  >  1 ) then

                  associate &
                    ( dJ_Eq    =>  Res_J_Eq ( iI ), &
                      dN_Eq    =>  Res_N_Eq ( iI ), &
                      dJ_Eq_P  =>  Res_J_Eq ( iI - 1 ), &
                      dN_Eq_P  =>  Res_N_Eq ( iI - 1 ) )

                  dJ_Eq  =  abs ( J_Eq ( iV )  -  J_Eq_P )  &
                            /  max ( abs ( J_Eq_0 ), SqrtTiny )
                  dN_Eq  =  abs ( N_Eq ( iV )  -  N_Eq_P )  &
                            /  max ( abs ( N_Eq_0 ), SqrtTiny )

! if ( iV == iShow ) then
!   call Show ( dJ_Eq, '>>> dJ_Eq' )
!   call Show ( dN_Eq, '>>> dN_Eq' )
! end if

                  N_I    ( iV )  =  iI
                  R_Max  ( iV )  =  max ( dJ_Eq, dN_Eq )
                  R_J_Eq ( iV )  =  dJ_Eq
                  R_N_Eq ( iV )  =  dN_Eq

                  if ( dJ_Eq  <  Tol .and. dN_Eq  <  Tol ) then
                    exit Relaxation
                  ! else if &
                  !   ( iI  >  10  &
                  !     .and.  &
                  !     dJ_Eq > Tol .and. dJ_Eq_P > Tol .and. dJ_Eq > dJ_Eq_P &
                  !     .and. &
                  !     dN_Eq > Tol .and. dN_Eq_P > Tol .and. dN_Eq > dN_Eq_P ) &
                  ! then
                  !   if ( iV  ==  iShow ) &
                  !     call Show ( '>>> Implicit solve diverging', &
                  !                 CONSOLE % ERROR )
                  !   exit Implicit
                  else if ( iI  ==  Max_I ) then
!                    if ( iV  ==  iShow ) &
!                      call Show ( 'Implicit solve maximum iterations', &
!                                  CONSOLE % ERROR )
                    exit Implicit
                  end if

                  end associate !-- dJ_Eq, etc.

                end if !-- iI > 1 (exit test)

                !-- Prepare for next iteration

                J_Eq_P  =  J_Eq ( iV )
                N_Eq_P  =  N_Eq ( iV )

                call R % ComputeFromBalanced ( iC, iV )
                call F % ComputeFromBalanced ( iC, iV )

              end do Implicit

              E_R ( iV )  =  E_R_0 ( iV )
              N_R ( iV )  =  N_R_0 ( iV )
              E_F ( iV )  =  E_F_0 ( iV )
              N_F ( iV )  =  N_F_0 ( iV )
              call R % ComputeFromBalanced ( iC, iV )
              call F % ComputeFromBalanced ( iC, iV )

              ! if ( iV  ==  iShow ) then
              !   ErrorRank  =  C % Communicator % Rank
              !   call Show ( Res_J_Eq ( : iI ), '>>> Res_J_Eq', &
              !               CONSOLE % ERROR )
              !   call Show ( Res_N_Eq ( : iI ), '>>> Res_N_Eq', &
              !               CONSOLE % ERROR )
              !   call Show ( ErrorRank, '>>> ErrorRank', CONSOLE % ERROR )
              !   call Show ( iS, '>>> iS', CONSOLE % ERROR )
              !   call Show ( iV, '>>> iV', CONSOLE % ERROR )
              !   call Show ( '>>> Adjusting relaxation', CONSOLE % ERROR )
              ! end if
              ! call PROGRAM_HEADER % Abort ( )

            end do Relaxation

            !--- Momentum update

            KK_R_S_1 ( iV )  &
              =  ( Xi_H ( iV )  &
                      -  Chi_H ( iV )  *  M_DD_11 ( iV ) * H_1 ( iV ) ) &
                 /  ( 1.0_KDR  +  Chi_H ( iV ) * dT )
            KK_R_S_2 ( iV )  &
              =  ( Xi_H ( iV )  &
                      -  Chi_H ( iV )  *  M_DD_22 ( iV ) * H_2 ( iV ) ) &
                 /  ( 1.0_KDR  +  Chi_H ( iV ) * dT )
            KK_R_S_3 ( iV )  &
              =  ( Xi_H ( iV )  &
                      -  Chi_H ( iV )  *  M_DD_33 ( iV ) * H_3 ( iV ) ) &
                 /  ( 1.0_KDR  +  Chi_H ( iV ) * dT )

            KK_F_S_1 ( iV )  =  - KK_R_S_1 ( iV )
            KK_F_S_2 ( iV )  =  - KK_R_S_2 ( iV )
            KK_F_S_3 ( iV )  =  - KK_R_S_3 ( iV )

          else

            KK_R_E   ( iV )  =  0.0_KDR
            KK_R_S_1 ( iV )  =  0.0_KDR
            KK_R_S_2 ( iV )  =  0.0_KDR
            KK_R_S_3 ( iV )  =  0.0_KDR
            KK_R_N   ( iV )  =  0.0_KDR

            KK_F_E   ( iV )  =  0.0_KDR
            KK_F_S_1 ( iV )  =  0.0_KDR
            KK_F_S_2 ( iV )  =  0.0_KDR
            KK_F_S_3 ( iV )  =  0.0_KDR
            KK_F_N   ( iV )  =  0.0_KDR

          end if !-- ProperCell
        end do !-- iV

        end associate !-- M_DD_11, etc.
        end associate !-- GSV

      class default
        call Show ( 'Gravitation type not recognized', CONSOLE % ERROR )
        call Show ( 'Step_RK_RM__Form', 'module', CONSOLE % ERROR )
        call Show ( 'SolveUpdateImplicit', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- G

      end associate !-- Xi_J, etc.
      end associate !-- I_V, etc.

      class default
        call Show ( 'Chart type not recognized', CONSOLE % ERROR )
        call Show ( 'Step_RK_NM_G__Form', 'module', CONSOLE % ERROR )
        call Show ( 'SolveUpdateImplicit', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- C
    end do !-- iC

    end select !-- ID
    end associate !-- KK_R, etc.
    end select !-- I
    end select !-- F
    end select !-- R
    end associate !-- S_1, etc.

  end subroutine SolveUpdateImplicit


end module Step_RK_NM_G__Form
