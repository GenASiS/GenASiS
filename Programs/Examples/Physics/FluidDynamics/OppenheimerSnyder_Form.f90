#include "Preprocessor"

module OppenheimerSnyder_Form

  !-- For example, Misner, Thorne, Wheeler p. 663, Eqs. (25.28)-(25.29)

  use GenASiS

  implicit none
  private

  type, public, extends ( Universe_F_CC_Form ) :: OppenheimerSnyderForm
    real ( KDR ) :: &
      Mass, &
      DensityInitial, &
      RadiusInitial, &
      TimeScale, &
      AtmosphereParameter
  !   type ( RootFinderForm ), allocatable :: &
  !     RootFinder
    type ( Fluid_D_Form ), allocatable :: &
      Reference, &
      Difference
  contains
    procedure, private, pass :: &
      Initialize_H
    procedure, public, pass :: &
      ComputeError
    final :: &
      Finalize
  end type OppenheimerSnyderForm

    private :: &
      InitializeUniverse, &
      ! InitializeDiagnostics, &
      SetInitial!, &
      ! SetReference

      private :: &
        SetFluid

        private :: &
          SetFluidKernel

contains

 
  subroutine Initialize_H ( U, NameOption )

    class ( OppenheimerSnyderForm ), intent ( inout ), target :: &
      U
    character ( * ), intent ( in ), optional  :: &
      NameOption

    character ( LDL ) :: &
      Name

    if ( U % Type  ==  '' ) &
      U % Type  =  'an OppenheimerSnyder'

    Name  =  'OppenheimerSnyder'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    call InitializeUniverse ( U, Name )
!    call InitializeDiagnostics ( U )

   if ( .not. associated ( U % Integrator % SetInitial ) ) &
     U % Integrator % SetInitial  =>  SetInitial

    U % Integrator % System  =>  U

  end subroutine Initialize_H


  subroutine ComputeError ( OS )

    class ( OppenheimerSnyderForm ), intent ( inout ) :: &
      OS

    real ( KDR ) :: &
      L1
    type ( CollectiveOperation_R_Form ) :: &
      CO
    
    select type ( A  =>  OS % Reference % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C    =>  A % Chart_GS, &
        F_R  =>  OS % Reference, &
        F_D  =>  OS % Difference )
    associate &
      ( FV_R  =>  F_R % Storage_GS % Value, &
        FV_D  =>  F_D % Storage_GS % Value )

    call CO % Initialize ( C % Communicator, [ 2 ], [ 2 ] )

    associate &
      ( D  =>  FV_D ( :, F_D % BARYON_DENSITY_C ), &
        R  =>  FV_R ( :, F_R % BARYON_DENSITY_C ), &
        Norm_D  =>  CO % Incoming % Value ( 1 ), &
        Norm_R  =>  CO % Incoming % Value ( 2 ) )

    CO % Outgoing % Value ( 1 ) &
      =  sum ( abs ( D ), mask = C % ProperCell )
    CO % Outgoing % Value ( 2 ) &
      =  sum ( abs ( R ), mask = C % ProperCell )

    call CO % Reduce ( REDUCTION % SUM )

    L1  =  Norm_D / Norm_R
    call Show ( L1, '*** L1 error', &
                nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )

    end associate !-- D, etc.
    end associate !-- FV_R, etc.
    end associate !-- C, etc.
    end select !-- A

  end subroutine ComputeError


  impure elemental subroutine Finalize ( OS )
    
    type ( OppenheimerSnyderForm ), intent ( inout ) :: &
      OS

    if ( allocated ( OS % Difference ) ) &
      deallocate ( OS % Difference )
    if ( allocated ( OS % Reference ) ) &
      deallocate ( OS % Reference )
    ! if ( allocated ( OS % RootFinder ) ) &
    !   deallocate ( OS % RootFinder )

  end subroutine Finalize


  subroutine InitializeUniverse ( OS, Name )

    class ( OppenheimerSnyderForm ), intent ( inout ) :: &
      OS
    character ( * ), intent ( in )  :: &
      Name

    ! integer ( KDI ) :: &
    !   iD

    call OS % Initialize &
           ( FluidType = 'DUST', &
             GravitationType = 'NEWTON_SG', &
             NameOption = Name, &
             nCellsPolarOption = 128 )

    !-- Modify from default OUTFLOW to INFLOW outer radial boundary condition
    select type ( I  =>  OS % Integrator )
      class is ( Integrator_CS_Form )
    associate &
      ( F  =>  I % CurrentSet_X )
    call F % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'INFLOW    ' ], iC = 1, iD = 1 )
    end associate !-- F
    end select !-- I

    ! OS % Integrator % SetReference  =>  SetReference

  end subroutine InitializeUniverse


  subroutine SetInitial ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    real ( KDR ) :: &
      Pi, &
      DensityFactor, &
      RadiusFactor, &
      Eta

    select type ( OS  =>  I % System )
      class is ( OppenheimerSnyderForm )
    select type ( I )
      class is ( Integrator_CS_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_D_Form )
    select type ( A  =>  F % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS )

    call Show ( 'Setting OppenheimerSnyder' )

    associate &
      ( R_Max => C % MaxCoordinate ( 1 ), &
          R_0 => OS % RadiusInitial, &
          D_0 => OS % DensityInitial, &
          Tau => OS % TimeScale, &
           AP => OS % AtmosphereParameter, &
            M => OS % Mass, &
           DF => DensityFactor, &
           RF => RadiusFactor )

     Pi  =  CONSTANT % PI
      M  =  1.0_KDR
    D_0  =  1.0e-3_KDR
     DF  =  1.0e2_KDR
     AP  =  1.0e-6_KDR
    call PROGRAM_HEADER % GetParameter (   M, 'Mass' )
    call PROGRAM_HEADER % GetParameter ( D_0, 'DensityInitial' )
    call PROGRAM_HEADER % GetParameter (  DF, 'DensityFactor' )
    call PROGRAM_HEADER % GetParameter (  AP, 'AtmosphereParameter' )

    Tau    =  sqrt ( 3.0 / ( 8.0 * Pi * D_0 ) )
      R_0  =  ( 3.0 * M / ( 4.0 * Pi * D_0 ) ) ** ( 1.0_KDR / 3.0_KDR )
     RF    =  DF ** ( - 1.0_KDR / 3.0_KDR ) 
    Eta    =  acos ( 2.0 * RF  -  1.0 )

    I % T_Finish  =  0.5 * Tau * ( Eta  +  sin ( Eta ) )

    call Show ( M, 'Mass' )
    call Show ( D_0, 'DensityInitial' )
    call Show ( R_0, 'RadiusInitial' )
    call Show ( DF, 'DensityFactor' )
    call Show ( RF, 'RadiusFactor' )
    call Show ( Pi / 2  *  Tau, 'CollapseTime' )
    call Show ( AP, 'AtmosphereParameter' )

    if ( R_0 > R_Max  ) then
      call Show ( 'RadiusInitial too large', CONSOLE % ERROR )
      call Show ( R_0, 'RadiusInitial', CONSOLE % ERROR )
      call Show ( C % MaxCoordinate ( 1 ), 'RadiusMax', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

!    allocate ( OS % RootFinder )
!    associate ( RF => OS % RootFinder )
!    call RF % Initialize ( OS )
!    RF % EvaluateZero  =>  EvaluateZeroEta

    call SetFluid ( OS, F )

    end associate !-- R_Max, etc.

    end associate !-- C
    end select !-- A
    end select !-- F
    end select !-- I
    end select !-- OS

  end subroutine SetInitial


  subroutine SetFluid ( OS, F )

    class ( OppenheimerSnyderForm ), intent ( inout ) :: &
      OS
    class ( Fluid_D_Form ), intent ( inout ) :: &
      F

    real ( KDR ) :: &
      Pi, &
      Eta, &
      Radius, &
      Density

    Pi  =  CONSTANT % PI

    associate &
!      ( RF   => OS % RootFinder, &
      (  D_0 => OS % DensityInitial, &
         R_0 => OS % RadiusInitial )

!    call RF % Solve ( [ 0.0_KDR, Pi ], Eta )

!    Radius   =  0.5 * R_0 * ( 1 + cos ( Eta ) )
!    Density  =  D_0 * ( R_0 / Radius ) ** 3
    Radius   =  R_0
    Density  =  D_0

    end associate !-- RF, etc.

    select type ( A  =>  F % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS, &
        G  =>  F % Geometry )
    associate &
      ( FV  =>  F % Storage_GS % Value, &
        GV  =>  G % Storage_GS % Value )

!  subroutine SetFluidKernel &
!               ( ProperCell, R_E, R_W, R_D, Density, N, V_1, V_2, V_3 )
   call SetFluidKernel &
          ( ProperCell = C % ProperCell, &
            R_E = GV ( :, G % EDGE_I_U_1 ), &
            R_W = GV ( :, G % WIDTH_U_1 ), &
            R_C = GV ( :, G % CENTER_U_1 ), &
            R_D = Radius, &
            D = Density, &
            M = OS % Mass, &
            AP = OS % AtmosphereParameter, &
            N = FV ( :, F % BARYON_DENSITY_C ), &
            V_1 = FV ( :, F % VELOCITY_U_1 ), &
            V_2 = FV ( :, F % VELOCITY_U_2 ), &
            V_3 = FV ( :, F % VELOCITY_U_3 ) )

    end associate !-- FV, etc.
    end associate !-- C, etc.
    end select !-- A

  end subroutine SetFluid


  subroutine SetFluidKernel &
               ( ProperCell, R_E, R_W, R_C, R_D, D, M, AP, N, V_1, V_2, V_3 )

    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      ProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      R_E, &
      R_W, &
      R_C
    real ( KDR ), intent ( in ) :: &
      R_D, &
      D, &
      M, &
      AP
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      N, &
      V_1, V_2, V_3

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nV
    real ( KDR ) :: &
      R_I, R_O

    nV  =  size ( N )

    !$OMP parallel do &
    !$OMP schedule ( OMP_SCHEDULE_HOST )
    do iV  =  1,  nV

!-- Establish atmosphere at outer radial boundary   
!      if ( .not. ProperCell ( iV ) ) &
!        cycle

      R_I  =  R_E ( iV )
      R_O  =  R_E ( iV )  +  R_W ( iV )
      if ( R_O  <=  R_D ) then
        N   ( iV )  =  D
        V_1 ( iV )  =  0.0_KDR
      else if ( R_I  <  R_D .and. R_O  >  R_D ) then
        N   ( iV )  =  D * ( R_D ** 3  -  R_I ** 3 ) &
                       / ( R_O ** 3  -  R_I ** 3 )
        V_1 ( iV )  =  0.0_KDR
      else
        N   ( iV )  =  AP  *  D  *  ( R_C ( iV ) / R_D ) ** ( -1.5_KDR )
        V_1 ( iV )  =  - sqrt ( 2.0_KDR * M / R_C ( iV ) )
      end if

      V_2 ( iV )  =  0.0_KDR
      V_3 ( iV )  =  0.0_KDR

    end do !-- iV
    !$OMP end parallel do

  end subroutine SetFluidKernel


end module OppenheimerSnyder_Form
