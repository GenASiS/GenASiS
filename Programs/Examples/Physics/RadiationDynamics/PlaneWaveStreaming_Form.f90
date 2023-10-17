#include "Preprocessor"

module PlaneWaveStreaming_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( Universe_R_B_Form ) :: PlaneWaveStreamingForm
    integer ( KDI ) :: &
      nPeriods
    integer ( KDI ), dimension ( 3 ) :: &
      nWavelengths
    real ( KDR ) :: &
      Speed, &
      Period
    real ( KDR ), dimension ( 3 ) :: &
      Wavenumber
!    type ( Fluid_D_Form ), allocatable :: &
!      Reference, &
!      Difference
  contains
    procedure, private, pass :: &
      Initialize_PWS
    generic, public :: &
      Initialize => Initialize_PWS
    final :: &
      Finalize
    procedure, public, pass :: &
      ShowParameters
    procedure, private, pass :: &
      Waveform
  end type PlaneWaveStreamingForm

    private :: &
      InitializeUniverse, &
    !   InitializeDiagnostics, &
      SetInitial!, &
    !   SetReference

      private :: &
        SetRadiation

        private :: &
          SetRadiationKernel


contains


  subroutine Initialize_PWS ( PWS, FormalismType, Name )

    class ( PlaneWaveStreamingForm ), intent ( inout ), target :: &
      PWS
    character ( * ), intent ( in )  :: &
      FormalismType, &
      Name

    if ( PWS % Type  ==  '' ) &
      PWS % Type  =  'a PlaneWaveStreaming'

    call InitializeUniverse ( PWS, FormalismType, Name )
    ! call InitializeDiagnostics ( PWS, FormalismType )
 
  end subroutine Initialize_PWS


  impure elemental subroutine Finalize ( PWS )

    type ( PlaneWaveStreamingForm ), intent ( inout ) :: &
      PWS

    ! if ( allocated ( PW % Difference ) ) &
    !   deallocate ( PW % Difference )
    ! if ( allocated ( PW % Reference ) ) &
    !   deallocate ( PW % Reference )

  end subroutine Finalize


  subroutine ShowParameters ( U )

    class ( PlaneWaveStreamingForm ), intent ( in ) :: &
      U

    call U % Universe_R_B_Form % ShowParameters ( )

    call Show ( U % nPeriods,     'nPeriods' )
    call Show ( U % nWavelengths, 'nWavelengths' )
    call Show ( U % Speed,        'Speed' )
    call Show ( U % Period,       'Period' )
    call Show ( U % Wavenumber,   'Wavenumber' )

  end subroutine ShowParameters


  function Waveform ( PWA, X ) result ( W )

    !-- Waveform with a full period in the range 0 < X < 1

    class ( PlaneWaveStreamingForm ), intent ( in ) :: &
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


  subroutine InitializeUniverse ( PWS, FormalismType, Name )

    class ( PlaneWaveStreamingForm ), intent ( inout ) :: &
      PWS
    character ( * ), intent ( in )  :: &
      FormalismType, &
      Name

    integer ( KDI ) :: &
      iD

    call PWS % Initialize &
           ( RadiationName = [ 'Radiation_1', 'Radiation_2' ], &
             RadiationType = [ 'GENERIC', 'GENERIC' ], &
             FormalismType = FormalismType, &
             Name = Name, &
             ApplyInteractionsOption = .false., &
             EvolveFluidOption = .false., &
             nCellsPositionOption = [ 128, 128, 128 ] )
             ! EnergySpacingOption = 'COMPACTIFIED', &
             ! nCellsEnergyOption = 4 )

    select type ( I  =>  PWS % Integrator )
      class is ( Integrator_CS_Form )
    associate &
      ( F  =>  I % CurrentSet_X )
    do iD  =  1, 3
      call F % SetBoundaryConditionsFace &
             ( [ 'PERIODIC', 'PERIODIC' ], iC = 1, iD = iD )
    end do !-- iD
    end associate !-- F
    end select !-- I
             
    select type ( I  =>  PWS % Integrator )
      class is ( Integrator_CS_1D_BM_CS_Form )
    associate &
      ( R  =>  I % CurrentSet_X_1D )
    do iD  =  1, 3
      call R % SetBoundaryConditionsFace &
             ( [ 'PERIODIC', 'PERIODIC' ], iC = 1, iD = iD )
    end do !-- iD
    end associate !-- R
    end select !-- I
             
    PWS % Integrator % SetInitial    =>  SetInitial
    ! PW % Integrator % SetReference  =>  SetReference

  end subroutine InitializeUniverse


  subroutine SetInitial ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    integer ( KDI ) :: &
      Direction

    select type ( PWS  =>  I % System )
      class is ( PlaneWaveStreamingForm )
    select type ( I )
      class is ( Integrator_CS_1D_CS_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_D_Form )
    select type ( A  =>  F % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS )

    Direction  =  ( -1 ) ** ( PWS % iRadiation  -  1 )

    PWS % nWavelengths = 0
    PWS % nWavelengths ( 1  :  C % nDimensions )  =  1
    call PROGRAM_HEADER % GetParameter ( PWS % nWavelengths, 'nWavelengths' )

    associate ( BoxSize  =>  C % MaxCoordinate  -  C % MinCoordinate )
    where ( BoxSize  >  0.0_KDR )
      PWS % Wavenumber  =  PWS % nWavelengths  /  BoxSize
    elsewhere
      PWS % Wavenumber  =  0.0_KDR
    end where
    end associate !-- BoxSize

    PWS % Speed  =  1.0_KDR
    call PROGRAM_HEADER % GetParameter ( PWS % Speed, 'Speed' )
    PWS % Speed  =  Direction  *  PWS % Speed

    associate &
      ( K      =>  PWS % Wavenumber, &
        Abs_K  =>  sqrt ( dot_product &
                            ( PWS % Wavenumber, PWS % Wavenumber ) ), &
        V      =>  PWS % Speed )
    PWS % Period  =  1.0_KDR / ( Abs_K * abs ( V ) )
    end associate !-- K, etc.

    PWS % nPeriods  =  1
    call PROGRAM_HEADER % GetParameter ( PWS % nPeriods, 'nPeriods' )

    I % T_Finish  =  PWS % nPeriods  *  PWS % Period

    select type ( I )
    class is ( Integrator_CS_1D_BM_CS_Form )

      select type ( R  =>  I % CurrentSet_X_1D )
        class is ( RadiationMoments_BM_Form )
    
      call SetRadiation ( PWS, R )

      end select !-- R

    class default
      call Show ( 'Integrator type not recognized', CONSOLE % ERROR )
      call Show ( 'PlaneWaveStreamingForm', 'module', CONSOLE % ERROR )
      call Show ( 'SetInitial', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- I

    end associate !-- C
    end select !-- A
    end select !-- F
    end select !-- I
    end select !-- PWS

  end subroutine SetInitial


  subroutine SetRadiation ( PWS, R )

    class ( PlaneWaveStreamingForm ), intent ( in ) :: &
      PWS
    type ( RadiationMoments_BM_Form ), intent ( inout ) :: &
      R  !-- RadiationSection

    select type ( A  =>  R % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS, &
        G  =>  R % Geometry )
    associate &
      ( RV  =>  R % Storage_GS % Value, &
        GV  =>  G % Storage_GS % Value )

    call SetRadiationKernel &
           ( J  = RV ( :, R % ENERGY_DENSITY_C ), &
             HX = RV ( :, R % MOMENTUM_DENSITY_C_U_1 ), &
             HY = RV ( :, R % MOMENTUM_DENSITY_C_U_2 ), &
             HZ = RV ( :, R % MOMENTUM_DENSITY_C_U_3 ), &
             PWS = PWS, &
             ProperCell = C % ProperCell, &
             X  = GV ( :, G % CENTER_U_1 ), &
             Y  = GV ( :, G % CENTER_U_2 ), &
             Z  = GV ( :, G % CENTER_U_3 ), &
             K  = PWS % Wavenumber, &
             V  = PWS % Speed, &
             T  = PWS % Integrator % T )

    end associate !-- RV, etc.
    end associate !-- C, etc.
    end select !-- A

  end subroutine SetRadiation


  subroutine SetRadiationKernel &
               ( J, HX, HY, HZ, PWS, ProperCell, X, Y, Z, K, V, T )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      J, &
      HX, HY, HZ
    class ( PlaneWaveStreamingForm ), intent ( in ) :: &
      PWS
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
      Abs_K, &
      VX, VY, VZ

    nV = size ( X )
    
    Abs_K  =  sqrt ( dot_product ( K, K ) )
       VX  =  V * K ( 1 ) / Abs_K
       VY  =  V * K ( 2 ) / Abs_K
       VZ  =  V * K ( 3 ) / Abs_K

    !$OMP parallel do &
    !$OMP schedule ( OMP_SCHEDULE_HOST ) firstprivate ( Abs_K, VX, VY, VZ )
    do iV = 1, nV

      if ( .not. ProperCell ( iV ) ) &
        cycle

      J ( iV )  =  PWS % Waveform &
                     (    K ( 1 ) * ( X ( iV )  -  VX * T ) &
                       +  K ( 2 ) * ( Y ( iV )  -  VY * T ) &
                       +  K ( 3 ) * ( Z ( iV )  -  VZ * T ) )

      HX ( iV )  =  VX  *  J ( iV )
      HY ( iV )  =  VY  *  J ( iV )
      HZ ( iV )  =  VZ  *  J ( iV )

    end do !-- iV
    !$OMP end parallel do

  end subroutine SetRadiationKernel


end module PlaneWaveStreaming_Form
