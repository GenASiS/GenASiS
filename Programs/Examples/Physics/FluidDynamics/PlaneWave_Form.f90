#include "Preprocessor"

module PlaneWave_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( Universe_F_B_Form ) :: PlaneWaveForm
    integer ( KDI ) :: &
      nPeriods
    integer ( KDI ), dimension ( 3 ) :: &
      nWavelengths
    real ( KDR ) :: &
      Speed, &
      Period
    real ( KDR ), dimension ( 3 ) :: &
      Wavenumber
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
    procedure, public, nopass :: &
      ShowSystem
    procedure, private, pass :: &
      Waveform
  end type PlaneWaveForm

    private :: &
      InitializeUniverse, &
      InitializeDiagnostics, &
      SetInitial, &
      SetReference

      private :: &
        SetFluid

        private :: &
          SetFluidKernel

contains


  subroutine Initialize_H ( U, NameOption )

    class ( PlaneWaveForm ), intent ( inout ), target :: &
      U
    character ( * ), intent ( in ), optional  :: &
      NameOption

    character ( LDL ) :: &
      Name

    if ( U % Type  ==  '' ) &
      U % Type  =  'a PlaneWave'

    Name  =  'PlaneWave'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    call InitializeUniverse ( U, Name )
    call InitializeDiagnostics ( U )

  end subroutine Initialize_H


  subroutine ComputeError ( PW )

    class ( PlaneWaveForm ), intent ( in ) :: &
      PW

    real ( KDR ) :: &
      L1
    type ( CollectiveOperation_R_Form ) :: &
      CO

    select type ( A  =>  PW % Reference % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C    =>  A % Chart_GS, &
        F_R  =>  PW % Reference, &
        F_D  =>  PW % Difference )
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


  impure elemental subroutine Finalize ( PW )

    type ( PlaneWaveForm ), intent ( inout ) :: &
      PW

    if ( allocated ( PW % Difference ) ) &
      deallocate ( PW % Difference )
    if ( allocated ( PW % Reference ) ) &
      deallocate ( PW % Reference )

  end subroutine Finalize


  subroutine ShowSystem ( I )

    class ( Integrator_H_Form ), intent ( in ) :: &
      I

    select type ( PW  =>  I % System )
      class is ( PlaneWaveForm )

    call PW % Universe_H_Form % ShowSystem ( I )

    call Show ( PW % nWavelengths, 'nWavelengths' )
    call Show ( PW % nPeriods,     'nPeriods' )
    call Show ( PW % Period,       'Period' )

    end select !-- PW

  end subroutine ShowSystem


  function Waveform ( PW, X ) result ( W )

    !-- Waveform with a full period in the range 0 < X < 1

    class ( PlaneWaveForm ), intent ( in ) :: &
      PW
    real ( KDR ), intent ( in ) :: &
      X
    real ( KDR ) :: &
      W

    call Show ( 'Waveform should be overridden', CONSOLE % WARNING )
    call Show ( 'PlaneWave_Form', 'module', CONSOLE % WARNING )
    call Show ( 'Waveform', 'function', CONSOLE % WARNING )

  end function Waveform


  subroutine InitializeUniverse ( PW, Name )

    class ( PlaneWaveForm ), intent ( inout ) :: &
      PW
    character ( * ), intent ( in )  :: &
      Name

    integer ( KDI ) :: &
      iD

    call PW % Initialize &
           ( FluidType = 'DUST', &
             GravitationType = 'GALILEO', &
             NameOption = Name, &
             nCellsOption = [ 128, 128, 128 ] )

    select type ( I  =>  PW % Integrator )
      class is ( Integrator_CS_Form )
    associate &
      ( F  =>  I % CurrentSet_X )
    do iD  =  1, 3
      call F % SetBoundaryConditionsFace &
             ( [ 'PERIODIC', 'PERIODIC' ], iC = 1, iD = iD )
    end do !-- iD
    end associate !-- F
    end select !-- I
             
    if ( .not. associated ( PW % Integrator % SetInitial ) ) &
      PW % Integrator % SetInitial  =>  SetInitial
    if ( .not. associated ( PW % Integrator % ShowSystem ) ) &
      PW % Integrator % ShowSystem  =>  ShowSystem

    PW % Integrator % SetReference  =>  SetReference

  end subroutine InitializeUniverse


  subroutine InitializeDiagnostics ( PW )

    class ( PlaneWaveForm ), intent ( inout ) :: &
      PW

    allocate &
      ( PW % Reference, &
        PW % Difference )
    associate &
      ( F_R  =>  PW % Reference, &
        F_D  =>  PW % Difference, &
        G    =>  PW % Integrator % Geometry_X, &
        S    =>  PW % Integrator % Checkpoint_X )

    call F_R % Initialize ( G, PW % Units_F, NameOption = 'Reference' )
    call F_D % Initialize ( G, PW % Units_F, NameOption = 'Difference' )
    call F_R % SetStream ( S )
    call F_D % SetStream ( S )

    end associate !-- FA_R, etc.

  end subroutine InitializeDiagnostics


  subroutine SetInitial ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    select type ( PW  =>  I % System )
      class is ( PlaneWaveForm )
    select type ( I )
      class is ( Integrator_CS_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_D_Form )
    select type ( A  =>  F % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS )

    PW % nWavelengths = 0
    PW % nWavelengths ( 1  :  C % nDimensions )  =  1
    call PROGRAM_HEADER % GetParameter ( PW % nWavelengths, 'nWavelengths' )

    associate ( BoxSize  =>  C % MaxCoordinate  -  C % MinCoordinate )
    where ( BoxSize  >  0.0_KDR )
      PW % Wavenumber  =  PW % nWavelengths / BoxSize
    elsewhere
      PW % Wavenumber  =  0.0_KDR
    end where

    PW % Speed  =  1.0_KDR
    call PROGRAM_HEADER % GetParameter ( PW % Speed, 'Speed' )

    associate &
      ( K      =>  PW % Wavenumber, &
        Abs_K  =>  sqrt ( dot_product ( PW % Wavenumber, PW % Wavenumber ) ), &
        V      =>  PW % Speed )
    PW % Period  =  1.0_KDR / ( Abs_K * V )

    PW % nPeriods  =  1
    call PROGRAM_HEADER % GetParameter ( PW % nPeriods, 'nPeriods' )

    I % T_Finish  =  PW % nPeriods  *  PW % Period

    call SetFluid ( PW, F )

    end associate !-- K, etc.
    end associate !-- BoxSize
    end associate !-- C
    end select !-- A
    end select !-- F
    end select !-- I
    end select !-- PW

  end subroutine SetInitial


  subroutine SetReference ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    select type ( PW  =>  I % System )
      class is ( PlaneWaveForm )
    select type ( I  =>  PW % Integrator )
      class is ( Integrator_CS_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_D_Form )
    associate &
      ( F_R  =>  PW % Reference, &
        F_D  =>  PW % Difference )

    call SetFluid ( PW, F_R )

    call F_D % MultiplyAdd ( F, F_R, -1.0_KDR )

    end associate !-- F_R, etc.
    end select !-- F
    end select !-- I
    end select !-- PW

  end subroutine SetReference


  subroutine SetFluid ( PW, F )

    class ( PlaneWaveForm ), intent ( inout ) :: &
      PW
    class ( Fluid_D_Form ), intent ( inout ) :: &
      F

    select type ( A  =>  F % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS, &
        G  =>  F % Geometry )
    associate &
      ( FV  =>  F % Storage_GS % Value, &
        GV  =>  G % Storage_GS % Value )

    call SetFluidKernel &
           (  N = FV ( :, F % BARYON_DENSITY_C ), &
             VX = FV ( :, F % VELOCITY_U_1 ), &
             VY = FV ( :, F % VELOCITY_U_2 ), &
             VZ = FV ( :, F % VELOCITY_U_3 ), &
             PW = PW, &
             ProperCell = C % ProperCell, &
             X = GV ( :, G % CENTER_U_1 ), &
             Y = GV ( :, G % CENTER_U_2 ), &
             Z = GV ( :, G % CENTER_U_3 ), &
             K = PW % Wavenumber, &
             V = PW % Speed, &
             T = PW % Integrator % T )

    end associate !-- FV, etc.
    end associate !-- C, etc.
    end select !-- A

  end subroutine SetFluid


  subroutine SetFluidKernel &
               ( N, VX, VY, VZ, PW, ProperCell, X, Y, Z, K, V, T )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      N, &
      VX, VY, VZ
    class ( PlaneWaveForm ), intent ( in ) :: &
      PW
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      ProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      X, Y, Z
    real ( KDR ), dimension ( 3 ), intent ( in ) :: &
      K
    real ( KDR ), intent ( in ) :: &
      V, &
      T

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nV
    real ( KDR ) :: &
      Abs_K

    nV  =  size ( X )
    
    Abs_K  =  sqrt ( dot_product ( K, K ) )

    !$OMP parallel do &
    !$OMP schedule ( OMP_SCHEDULE_HOST ) firstprivate ( Abs_K )
    do iV  =  1,  nV

      if ( .not. ProperCell ( iV ) ) &
        cycle

      VX ( iV )  =  V  *  K ( 1 )  /  Abs_K
      VY ( iV )  =  V  *  K ( 2 )  /  Abs_K
      VZ ( iV )  =  V  *  K ( 3 )  /  Abs_K

      N ( iV )  =  PW % Waveform &
                     (    K ( 1 )  *  ( X ( iV )  -  VX ( iV )  *  T ) &
                       +  K ( 2 )  *  ( Y ( iV )  -  VY ( iV )  *  T ) &
                       +  K ( 3 )  *  ( Z ( iV )  -  VZ ( iV )  *  T ) )

    end do !-- iV
    !$OMP end parallel do

  end subroutine SetFluidKernel


end module PlaneWave_Form
