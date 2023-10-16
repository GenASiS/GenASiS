#include "Preprocessor"

module PlaneWaveAdvection_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( Universe_F_B_Form ) :: PlaneWaveAdvectionForm
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
    procedure, public, pass :: &
      ShowParameters
    procedure, private, pass :: &
      Waveform
  end type PlaneWaveAdvectionForm

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


  subroutine Initialize_H ( U, Name, CommunicatorOption )

    class ( PlaneWaveAdvectionForm ), intent ( inout ), target :: &
      U
    character ( * ), intent ( in ) :: &
      Name
    type ( CommunicatorForm ), intent ( in ), target, optional :: &
      CommunicatorOption

    if ( U % Type  ==  '' ) &
      U % Type  =  'a PlaneWaveAdvection'

    call InitializeUniverse ( U, Name )
    call InitializeDiagnostics ( U )

  end subroutine Initialize_H


  subroutine ComputeError ( PWA )

    class ( PlaneWaveAdvectionForm ), intent ( in ) :: &
      PWA

    real ( KDR ) :: &
      L1
    type ( CollectiveOperation_R_Form ) :: &
      CO

    select type ( A  =>  PWA % Reference % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C    =>  A % Chart_GS, &
        F_R  =>  PWA % Reference, &
        F_D  =>  PWA % Difference )
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


  impure elemental subroutine Finalize ( PWA )

    type ( PlaneWaveAdvectionForm ), intent ( inout ) :: &
      PWA

    if ( allocated ( PWA % Difference ) ) &
      deallocate ( PWA % Difference )
    if ( allocated ( PWA % Reference ) ) &
      deallocate ( PWA % Reference )

  end subroutine Finalize


  subroutine ShowParameters ( U )

    class ( PlaneWaveAdvectionForm ), intent ( in ) :: &
      U

    call U % Universe_F_B_Form % ShowParameters ( )

    call Show ( U % nPeriods,     'nPeriods' )
    call Show ( U % nWavelengths, 'nWavelengths' )
    call Show ( U % Speed,        'Speed' )
    call Show ( U % Period,       'Period' )
    call Show ( U % Wavenumber,   'Wavenumber' )

  end subroutine ShowParameters


  function Waveform ( PWA, X ) result ( W )

    !-- Waveform with a full period in the range 0 < X < 1

    class ( PlaneWaveAdvectionForm ), intent ( in ) :: &
      PWA
    real ( KDR ), intent ( in ) :: &
      X
    real ( KDR ) :: &
      W
    
    W = huge ( 1.0_KDR ) 
    call Show ( 'Waveform should be overridden', CONSOLE % WARNING )
    call Show ( 'PlaneWaveAdvection_Form', 'module', CONSOLE % WARNING )
    call Show ( 'Waveform', 'function', CONSOLE % WARNING )

  end function Waveform


  subroutine InitializeUniverse ( PWA, Name )

    class ( PlaneWaveAdvectionForm ), intent ( inout ) :: &
      PWA
    character ( * ), intent ( in )  :: &
      Name

    integer ( KDI ) :: &
      iD

    call PWA % Initialize &
           ( FluidType = 'DUST', &
             GravitationType = 'GALILEO', &
             Name = Name, &
             nCellsOption = [ 128, 128, 128 ] )

    select type ( I  =>  PWA % Integrator )
      class is ( Integrator_CS_Form )
    associate &
      ( F  =>  I % CurrentSet_X )
    do iD  =  1, 3
      call F % SetBoundaryConditionsFace &
             ( [ 'PERIODIC', 'PERIODIC' ], iC = 1, iD = iD )
    end do !-- iD
    end associate !-- F
    end select !-- I
             
    PWA % Integrator % SetInitial    =>  SetInitial
    PWA % Integrator % SetReference  =>  SetReference

  end subroutine InitializeUniverse


  subroutine InitializeDiagnostics ( PWA )

    class ( PlaneWaveAdvectionForm ), intent ( inout ) :: &
      PWA

    allocate &
      ( PWA % Reference, &
        PWA % Difference )
    associate &
      ( F_R  =>  PWA % Reference, &
        F_D  =>  PWA % Difference, &
        G    =>  PWA % Integrator % Geometry_X, &
        S    =>  PWA % Integrator % Checkpoint_X )

    call F_R % Initialize ( G, PWA % Units_F, NameOption = 'Reference' )
    call F_D % Initialize ( G, PWA % Units_F, NameOption = 'Difference' )
    call F_R % SetStream ( S )
    call F_D % SetStream ( S )

    end associate !-- FA_R, etc.

  end subroutine InitializeDiagnostics


  subroutine SetInitial ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    select type ( PWA  =>  I % System )
      class is ( PlaneWaveAdvectionForm )
    select type ( I )
      class is ( Integrator_CS_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_D_Form )
    select type ( A  =>  F % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS )

    PWA % nWavelengths = 0
    PWA % nWavelengths ( 1  :  C % nDimensions )  =  1
    call PROGRAM_HEADER % GetParameter ( PWA % nWavelengths, 'nWavelengths' )

    associate ( BoxSize  =>  C % MaxCoordinate  -  C % MinCoordinate )
    where ( BoxSize  >  0.0_KDR )
      PWA % Wavenumber  =  PWA % nWavelengths  /  BoxSize
    elsewhere
      PWA % Wavenumber  =  0.0_KDR
    end where
    end associate !-- BoxSize

    PWA % Speed  =  1.0_KDR
    call PROGRAM_HEADER % GetParameter ( PWA % Speed, 'Speed' )

    associate &
      ( K      =>  PWA % Wavenumber, &
        Abs_K  =>  sqrt ( dot_product &
                            ( PWA % Wavenumber, PWA % Wavenumber ) ), &
        V      =>  PWA % Speed )
    PWA % Period  =  1.0_KDR / ( Abs_K * V )
    end associate !-- K, etc.

    PWA % nPeriods  =  1
    call PROGRAM_HEADER % GetParameter ( PWA % nPeriods, 'nPeriods' )

    I % T_Finish  =  PWA % nPeriods  *  PWA % Period

    call SetFluid ( PWA, F )

    end associate !-- C
    end select !-- A
    end select !-- F
    end select !-- I
    end select !-- PWA

  end subroutine SetInitial


  subroutine SetReference ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    select type ( PWA  =>  I % System )
      class is ( PlaneWaveAdvectionForm )
    select type ( I  =>  PWA % Integrator )
      class is ( Integrator_CS_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_D_Form )
    associate &
      ( F_R  =>  PWA % Reference, &
        F_D  =>  PWA % Difference )

    call SetFluid ( PWA, F_R )

    call F_D % MultiplyAdd ( F, F_R, -1.0_KDR )

    end associate !-- F_R, etc.
    end select !-- F
    end select !-- I
    end select !-- PWA

  end subroutine SetReference


  subroutine SetFluid ( PWA, F )

    class ( PlaneWaveAdvectionForm ), intent ( inout ) :: &
      PWA
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
             PWA = PWA, &
             ProperCell = C % ProperCell, &
             X = GV ( :, G % CENTER_U_1 ), &
             Y = GV ( :, G % CENTER_U_2 ), &
             Z = GV ( :, G % CENTER_U_3 ), &
             K = PWA % Wavenumber, &
             V = PWA % Speed, &
             T = PWA % Integrator % T )

    end associate !-- FV, etc.
    end associate !-- C, etc.
    end select !-- A

  end subroutine SetFluid


  subroutine SetFluidKernel &
               ( N, VX, VY, VZ, PWA, ProperCell, X, Y, Z, K, V, T )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      N, &
      VX, VY, VZ
    class ( PlaneWaveAdvectionForm ), intent ( in ) :: &
      PWA
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

      N ( iV )  =  PWA % Waveform &
                     (    K ( 1 )  *  ( X ( iV )  -  VX ( iV )  *  T ) &
                       +  K ( 2 )  *  ( Y ( iV )  -  VY ( iV )  *  T ) &
                       +  K ( 3 )  *  ( Z ( iV )  -  VZ ( iV )  *  T ) )

    end do !-- iV
    !$OMP end parallel do

  end subroutine SetFluidKernel


end module PlaneWaveAdvection_Form
