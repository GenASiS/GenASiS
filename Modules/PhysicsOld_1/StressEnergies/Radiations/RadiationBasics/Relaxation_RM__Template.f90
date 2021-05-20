module Relaxation_RM__Template

  use Basics
  use Mathematics
  use RadiationMoments_Form
  use Sources_RM__Form
  use RadiationMoments_ASC__Form

  implicit none
  private
  
  type, public, abstract :: Relaxation_RM_Template
    integer ( KDI ) :: &
      IGNORABILITY = 0
    integer ( KDI ) :: &
      iEnergy, &
      iMomentum_1, iMomentum_2, iMomentum_3
    character ( LDL ) :: &
      Name = '', &
      Type = ''
    type ( LinearEquations_LAPACK_Form ), dimension ( : ), allocatable :: &
      LinearEquations
    procedure ( AS ), public, nopass, pointer :: &
      Apply => null ( )
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure ( AS ), private, nopass, deferred :: &
      ApplySubroutine
    procedure, public, pass :: &
      FinalizeTemplate
    procedure, public, pass :: &
      SetIndices
    procedure, public, pass ( R ) :: &
      ComputeIncrements
  end type Relaxation_RM_Template

  abstract interface
    subroutine AS ( S, RadiationMoments, Sources_RM, Increment, Chart, &
                   TimeStep, iStage, GeometryOption, iStrgeometryValueOption )
      use Basics
      use Mathematics
      class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
        S
      class ( CurrentTemplate ), intent ( inout ) :: &
        RadiationMoments
      class ( Sources_C_Form ), intent ( inout ) :: &
        Sources_RM
      type ( StorageForm ), intent ( inout ) :: &
        Increment
      class ( ChartTemplate ), intent ( in ) :: &
        Chart
      real ( KDR ), intent ( in ) :: &
        TimeStep
      integer ( KDI ), intent ( in ) :: &
        iStage
      class ( GeometryFlatForm ), intent ( in ), optional :: &
        GeometryOption
      integer ( KDI ), intent ( in ), optional :: &
        iStrgeometryValueOption
    end subroutine AS
  end interface

    private :: &
      ComputeCoefficients, &
      ComputeComovingIncrements, &
      ComputeConservedIncrements, &
      ComputeSources

contains

  subroutine InitializeTemplate ( R, Name )

    class ( Relaxation_RM_Template ), intent ( inout ) :: &
      R
    character ( * ), intent ( in ) :: &
      Name

    if ( R % Type == '' ) &
      R % Type  =  'a Relaxation_RM'

    R % Name  =  Name

    R % IGNORABILITY = CONSOLE % INFO_1 
    call Show ( 'Initializing ' // trim ( R % Type ), R % IGNORABILITY )
    call Show ( R % Name, 'Name', R % IGNORABILITY )

    if ( .not. associated ( R % Apply ) ) then
      call Show ( 'Apply procedure pointer member not set', CONSOLE % ERROR )
      call Show ( 'Set before calling InitializeTemplate', CONSOLE % ERROR )
      call Show ( 'InitializeTemplate', 'subroutine', CONSOLE % ERROR )
      call Show ( 'Relaxation_RM__Template', 'module', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

  end subroutine InitializeTemplate


  subroutine FinalizeTemplate ( R )

    class ( Relaxation_RM_Template ), intent ( inout ) :: &
      R

    if ( allocated ( R % LinearEquations ) ) &
      deallocate ( R % LinearEquations )

    call Show ( 'Finalizing ' // trim ( R % Type ), R % IGNORABILITY )
    call Show ( R % Name, 'Name', R % IGNORABILITY )

  end subroutine FinalizeTemplate


  subroutine SetIndices ( R, RMA )

    class ( Relaxation_RM_Template ), intent ( inout ) :: &
      R
    class ( RadiationMoments_ASC_Form ), intent ( in ), target :: &
      RMA

    class ( RadiationMomentsForm ), pointer :: &
      RM

    RM  =>  RMA % RadiationMoments ( )

    call Search ( RM % iaConserved, RM % CONSERVED_ENERGY, &
                  R % iEnergy )
    call Search ( RM % iaConserved, RM % CONSERVED_MOMENTUM_D ( 1 ), &
                  R % iMomentum_1 )
    call Search ( RM % iaConserved, RM % CONSERVED_MOMENTUM_D ( 2 ), &
                  R % iMomentum_2 )
    call Search ( RM % iaConserved, RM % CONSERVED_MOMENTUM_D ( 3 ), &
                  R % iMomentum_3 )
!    if ( associated ( NM ) ) &
!      call Search ( NM % iaConserved, NM % CONSERVED_NUMBER, iNumber )

    nullify ( RM )

  end subroutine SetIndices


  subroutine ComputeIncrements &
               ( S, Sources_RM, LE, Increment, R, RadiationMoments, &
                 TimeStep, iStage, iV, G, iStrgeometryValueOption )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    class ( Sources_C_Form ), intent ( inout ) :: &
      Sources_RM
    type ( LinearEquations_LAPACK_Form ), intent ( inout ) :: &
      LE
    type ( StorageForm ), intent ( inout ) :: &
      Increment
    class ( Relaxation_RM_Template ), intent ( in ) :: &
      R
    class ( CurrentTemplate ), intent ( in ) :: &
      RadiationMoments
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage, &
      iV  !-- iValue
    integer ( KDI ), intent ( in ), optional :: &
      iStrgeometryValueOption

    integer ( KDI ) :: &
      iS
    real ( KDR ) :: &
      dJ, dH_D_1, dH_D_2, dH_D_3
    real ( KDR ), dimension ( 3, 3 ) :: &
      k_DD  !-- Eddington tensor components ( K / J, Stress / EnergyDensity )

    select type ( RM => RadiationMoments )
    class is ( RadiationMomentsForm )

    associate ( I => RM % Interactions )

    select type ( SRM => Sources_RM )
    class is ( Sources_RM_Form )

    if ( present ( iStrgeometryValueOption ) ) then
      iS = iStrgeometryValueOption
    else
      iS = iV
    end if

    call ComputeCoefficients &
           ( RM, &
             Increment % Value ( iV, R % iEnergy ), &
             Increment % Value ( iV, R % iMomentum_1 ), &
             Increment % Value ( iV, R % iMomentum_2 ), &
             Increment % Value ( iV, R % iMomentum_3 ), &
             I  % Value ( iV, I % EMISSIVITY_J ), &
             I  % Value ( iV, I % OPACITY_J ), &
             I  % Value ( iV, I % OPACITY_H ), &
             RM % Value ( iV, RM % COMOVING_ENERGY ), &
             RM % Value ( iV, RM % COMOVING_MOMENTUM_U ( 1 ) ), &
             RM % Value ( iV, RM % COMOVING_MOMENTUM_U ( 2 ) ), &
             RM % Value ( iV, RM % COMOVING_MOMENTUM_U ( 3 ) ), &
             RM % Value ( iV, RM % STRESS_FACTOR ), &
             RM % Value ( iV, RM % FLUID_VELOCITY_U ( 1 ) ), &
             RM % Value ( iV, RM % FLUID_VELOCITY_U ( 2 ) ), &
             RM % Value ( iV, RM % FLUID_VELOCITY_U ( 3 ) ), &
             G  % Value ( iS, G % METRIC_DD_22 ), &
             G  % Value ( iS, G % METRIC_DD_33 ), &
             TimeStep, LE % Matrix, &
             LE % RightHandSide, k_DD )

    call ComputeComovingIncrements &
           ( LE, dJ, dH_D_1, dH_D_2, dH_D_3 )

    call ComputeConservedIncrements &
           ( k_DD, dJ, dH_D_1, dH_D_2, dH_D_3, &
             RM % Value ( iV, RM % FLUID_VELOCITY_U ( 1 ) ), &
             RM % Value ( iV, RM % FLUID_VELOCITY_U ( 2 ) ), &
             RM % Value ( iV, RM % FLUID_VELOCITY_U ( 3 ) ), &
             G  % Value ( iS, G % METRIC_DD_22 ), &
             G  % Value ( iS, G % METRIC_DD_33 ), &
             Increment % Value ( iV, R % iEnergy ), &
             Increment % Value ( iV, R % iMomentum_1 ), &
             Increment % Value ( iV, R % iMomentum_2 ), &
             Increment % Value ( iV, R % iMomentum_3 ) )

    call ComputeSources &
           ( SRM % Value ( iV, SRM % INTERACTIONS_J ), &
             SRM % Value ( iV, SRM % INTERACTIONS_H_D ( 1 ) ), &
             SRM % Value ( iV, SRM % INTERACTIONS_H_D ( 2 ) ), &
             SRM % Value ( iV, SRM % INTERACTIONS_H_D ( 3 ) ), &
               I % Value ( iV, I % EMISSIVITY_J ), &
               I % Value ( iV, I % OPACITY_J ), &
               I % Value ( iV, I % OPACITY_H ), &
              RM % Value ( iV, RM % COMOVING_ENERGY ), &
              RM % Value ( iV, RM % COMOVING_MOMENTUM_U ( 1 ) ), &
              RM % Value ( iV, RM % COMOVING_MOMENTUM_U ( 2 ) ), &
              RM % Value ( iV, RM % COMOVING_MOMENTUM_U ( 3 ) ), &
              dJ, dH_D_1, dH_D_2, dH_D_3, & 
               G % Value ( iS, G % METRIC_DD_22 ), &
               G % Value ( iS, G % METRIC_DD_33 ), &
               S % B ( iStage ) )

    end select !-- SRM
    end associate !-- I
    end select !-- RM

  end subroutine ComputeIncrements


  subroutine ComputeCoefficients &
               ( RM, dE, dS_D_1, dS_D_2, dS_D_3, Xi_J, Chi_J, Chi_H, &
                 J, H_1, H_2, H_3, SF, V_1, V_2, V_3, M_DD_22, M_DD_33, dt, &
                 A, B, k_DD )

    class ( RadiationMomentsForm ), intent ( in ) :: &
      RM
    real ( KDR ), intent ( in ) :: &
      dE, dS_D_1, dS_D_2, dS_D_3, &
      Xi_J, Chi_J, Chi_H, &
      J, H_1, H_2, H_3, SF, &
      V_1, V_2, V_3, &
      M_DD_22, M_DD_33, &
      dt
    real ( KDR ), dimension ( 4, 4 ), intent ( out ) :: &
      A
    real ( KDR ), dimension ( 4, 1 ), intent ( out ) :: &
      B
    real ( KDR ), dimension ( 3, 3 ), intent ( out ) :: &
      k_DD
    
    call RM % ComputeComovingStress_D &
           ( k_DD ( 1, : ), [ H_1, H_2, H_3 ], J, SF, M_DD_22, M_DD_33, iD = 1 )
    call RM % ComputeComovingStress_D &
           ( k_DD ( 2, : ), [ H_1, H_2, H_3 ], J, SF, M_DD_22, M_DD_33, iD = 2 )
    call RM % ComputeComovingStress_D &
           ( k_DD ( 3, : ), [ H_1, H_2, H_3 ], J, SF, M_DD_22, M_DD_33, iD = 3 )

    if ( J > 0.0_KDR ) then
      k_DD  =  k_DD / J
    else
      k_DD  =  0.0_KDR
    end if

    A ( 1, 1 )  =  1.0_KDR  +  Chi_J * dt
    A ( 2, 1 )  =            V_1 * ( 1.0_KDR  +  Chi_J * dt )  &
                   +  dot_product ( k_DD ( 1, : ), [ V_1, V_2, V_3 ] )
    A ( 3, 1 )  =  M_DD_22 * V_2 * ( 1.0_KDR  +  Chi_J * dt )  &
                   +  dot_product ( k_DD ( 2, : ), [ V_1, V_2, V_3 ] )
    A ( 4, 1 )  =  M_DD_33 * V_3 * ( 1.0_KDR  +  Chi_J * dt )  &
                   +  dot_product ( k_DD ( 3, : ), [ V_1, V_2, V_3 ] )

    A ( 1, 2 )  =  ( 2.0_KDR  +  Chi_H * dt )  *  V_1
    A ( 2, 2 )  =  1.0_KDR  +  Chi_H * dt
    A ( 3, 2 )  =  0.0_KDR
    A ( 4, 2 )  =  0.0_KDR

    A ( 1, 3 )  =  ( 2.0_KDR  +  Chi_H * dt )  *  V_2
    A ( 2, 3 )  =  0.0_KDR
    A ( 3, 3 )  =  1.0_KDR  +  Chi_H * dt
    A ( 4, 3 )  =  0.0_KDR

    A ( 1, 4 )  =  ( 2.0_KDR  +  Chi_H * dt )  *  V_3
    A ( 2, 4 )  =  0.0_KDR
    A ( 3, 4 )  =  0.0_KDR
    A ( 4, 4 )  =  1.0_KDR  +  Chi_H * dt

    B ( 1, 1 )  =  dE    +  (    Xi_J  -  Chi_J * J  &
                              -  Chi_H * (             V_1 * H_1  &
                                           * M_DD_22 * V_2 * H_2  &
                                           * M_DD_33 * V_3 * H_3 ) ) * dt
    B ( 2, 1 )  =  dS_D_1  +            ( -  Chi_H * H_1  &
                                        +  V_1 * ( Xi_J  -  Chi_J * J ) ) * dt
    B ( 3, 1 )  =  dS_D_2  +  M_DD_22 * ( -  Chi_H * H_2  &
                                        +  V_2 * ( Xi_J  -  Chi_J * J ) ) * dt
    B ( 4, 1 )  =  dS_D_3  +  M_DD_33 * ( -  Chi_H * H_3  &
                                        +  V_3 * ( Xi_J  -  Chi_J * J ) ) * dt

  end subroutine ComputeCoefficients


  subroutine ComputeComovingIncrements ( LE, dJ, dH_D_1, dH_D_2, dH_D_3 )

    type ( LinearEquations_LAPACK_Form ), intent ( inout ) :: &
      LE
    real ( KDR ), intent ( out ) :: &
      dJ, dH_D_1, dH_D_2, dH_D_3

    call LE % Solve ( )

    dJ      =  LE % Solution ( 1, 1 )
    dH_D_1  =  LE % Solution ( 2, 1 )
    dH_D_2  =  LE % Solution ( 3, 1 )
    dH_D_3  =  LE % Solution ( 4, 1 )

  end subroutine ComputeComovingIncrements


  subroutine ComputeConservedIncrements &
               ( k_DD, dJ, dH_D_1, dH_D_2, dH_D_3, V_1, V_2, V_3, &
                 M_DD_22, M_DD_33, dE, dS_D_1, dS_D_2, dS_D_3 )

    real ( KDR ), dimension ( 3, 3 ), intent ( in ) :: &
      k_DD
    real ( KDR ), intent ( in ) :: &
      dJ, dH_D_1, dH_D_2, dH_D_3, &
      V_1, V_2, V_3, &
      M_DD_22, M_DD_33
    real ( KDR ), intent ( out ) :: &
      dE, dS_D_1, dS_D_2, dS_D_3

    dE    =  dJ    +  2.0_KDR  *  dot_product ( [ V_1, V_2, V_3 ], &
                                                [ dH_D_1, dH_D_2, dH_D_3 ] )
    dS_D_1  =  dH_D_1  &
             +  ( V_1  &
                  +  dot_product ( k_DD ( 1, : ), [ V_1, V_2, V_3 ] ) )  *  dJ  

    dS_D_2  =  dH_D_2  &
             +  ( M_DD_22 * V_2  &
                  +  dot_product ( k_DD ( 2, : ), [ V_1, V_2, V_3 ] ) )  *  dJ  

    dS_D_3  =  dH_D_3  &
             +  ( M_DD_33 * V_3  &
                  +  dot_product ( k_DD ( 3, : ), [ V_1, V_2, V_3 ] ) )  *  dJ  

  end subroutine ComputeConservedIncrements


  subroutine ComputeSources &
               ( Q, A_D_1, A_D_2, A_D_3, Xi_J, Chi_J, Chi_H, &
                 J, H_U_1, H_U_2, H_U_3, dJ, dH_D_1, dH_D_2, dH_D_3, &
                 M_DD_22, M_DD_33, Weight_RK )

    real ( KDR ), intent ( inout ) :: &
      Q, A_D_1, A_D_2, A_D_3
    real ( KDR ), intent ( in ) :: &
      Xi_J, Chi_J, Chi_H, &
      J, H_U_1, H_U_2, H_U_3, &
      dJ, dH_D_1, dH_D_2, dH_D_3, &
      M_DD_22, M_DD_33, &
      Weight_RK

    Q      =  Q      &
              +  Weight_RK * ( Xi_J  &
                               -  Chi_J  *  (               J  +      dJ ) )
    A_D_1  =  A_D_1  &
              +  Weight_RK * (  &
                               -  Chi_H  *  (           H_U_1  +  dH_D_1 ) )
    A_D_2  =  A_D_2  &
              +  Weight_RK * (  &
                               -  Chi_H  *  ( M_DD_22 * H_U_2  +  dH_D_2 ) )
    A_D_3  =  A_D_3  &
              +  Weight_RK * (  &
                               -  Chi_H  *  ( M_DD_33 * H_U_3  +  dH_D_3 ) )

  end subroutine ComputeSources


end module Relaxation_RM__Template
