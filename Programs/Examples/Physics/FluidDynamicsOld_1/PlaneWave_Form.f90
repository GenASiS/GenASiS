#include "Preprocessor"

module PlaneWave_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( FluidBoxForm ) :: PlaneWaveForm
    real ( KDR ) :: &
      Speed
    real ( KDR ), dimension ( 3 ) :: &
      Wavenumber
    type ( Fluid_D_A_Form ), allocatable :: &
      Reference, &
      Difference
  contains
    procedure, private, pass :: &
      Initialize_H
    procedure, public, pass :: &
      ComputeError
    final :: &
      Finalize
    procedure, private, pass :: &
      Waveform
  end type PlaneWaveForm

  class ( PlaneWaveForm ), public, pointer :: &
    PLANE_WAVE => null ( )  !-- Makes instance of PlaneWave accessible to 
                            !   SetInitial and SetReference

    private :: &
      InitializeFluidBox, &
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

    if ( .not. associated ( PLANE_WAVE ) ) &
      PLANE_WAVE  =>  U

    if ( U % Type  ==  '' ) &
      U % Type  =  'a PlaneWave'

    Name  =  'PlaneWave'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    call InitializeFluidBox ( U, Name )
    call InitializeDiagnostics ( U )

    U % Integrator % SetInitial  =>  SetInitial

  end subroutine Initialize_H


  subroutine ComputeError ( PW )

    class ( PlaneWaveForm ), intent ( in ) :: &
      PW

    real ( KDR ) :: &
      L1
    type ( CollectiveOperation_R_Form ) :: &
      CO

    select type ( FC_R  =>  PW % Reference % FieldSet_C ( 1 ) % Element )
      class is ( Fluid_D_C_Form )
    select type ( FC_D  =>  PW % Difference % FieldSet_C ( 1 ) % Element )
      class is ( Fluid_D_C_Form )
    select type ( C  =>  FC_R % Chart )
      class is ( Chart_GS_Form )
    associate &
      ( FV_R  =>  FC_R % Storage_FSC % Storage % Value, &
        FV_D  =>  FC_D % Storage_FSC % Storage % Value )

    call CO % Initialize ( C % Communicator, [ 2 ], [ 2 ] )

    associate &
      ( D  =>  FV_D ( :, FC_D % BARYON_DENSITY_C ), &
        R  =>  FV_R ( :, FC_R % BARYON_DENSITY_C ), &
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

    end associate !-- FV_R, etc.
    end associate !-- D, etc.
    end select !-- C
    end select !-- FC_D
    end select !-- FC_R

  end subroutine ComputeError


  impure elemental subroutine Finalize ( PW )

    type ( PlaneWaveForm ), intent ( inout ) :: &
      PW

    if ( allocated ( PW % Difference ) ) &
      deallocate ( PW % Difference )
    if ( allocated ( PW % Reference ) ) &
      deallocate ( PW % Reference )

  end subroutine Finalize


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


  subroutine InitializeFluidBox ( PW, Name )

    class ( PlaneWaveForm ), intent ( inout ) :: &
      PW
    character ( * ), intent ( in )  :: &
      Name
    
    call PW % Initialize &
           ( FluidType = 'DUST', &
             GravitationType = 'GALILEO', &
             NameOption = Name, &
             nCellsOption = [ 128, 128, 128 ] )
             
    PW % Integrator % SetReference  =>  SetReference

  end subroutine InitializeFluidBox


  subroutine InitializeDiagnostics ( PW )

    class ( PlaneWaveForm ), intent ( inout ) :: &
      PW

    allocate &
      ( PW % Reference, &
        PW % Difference )
    associate &
      ( FA_R  =>  PW % Reference, &
        FA_D  =>  PW % Difference, &
        GA    =>  PW % Integrator % Geometry_X_A, &
        SA    =>  PW % Integrator % Checkpoint_X_A )

    call FA_R % Initialize ( GA, PW % Units_F, NameOption = 'Reference' )
    call FA_D % Initialize ( GA, PW % Units_F, NameOption = 'Difference' )
    call FA_R % SetStream ( SA )
    call FA_D % SetStream ( SA )

    end associate !-- FA_R, etc.

  end subroutine InitializeDiagnostics


  subroutine SetInitial ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    integer ( KDI ) :: &
      nPeriods
    integer ( KDI ), dimension ( 3 ) :: &
      nWavelengths
    real ( KDR ) :: &
      Period

    associate &
      ( PW  =>  PLANE_WAVE )
    select type ( I )
      class is ( Integrator_CSA_Form )
    select type ( FA  =>  I % CurrentSet_X_A )
      class is ( Fluid_D_A_Form )
    select type ( FC  =>  FA % FieldSet_C ( 1 ) % Element )
      class is ( Fluid_D_C_Form )
    select type ( C  =>  FC % Chart )
      class is ( Chart_GS_Form )

    call Show ( 'Setting PlaneWave' )

    nWavelengths = 0
    nWavelengths ( 1  :  C % nDimensions )  =  1
    call PROGRAM_HEADER % GetParameter ( nWavelengths, 'nWavelengths' )

    associate ( BoxSize  =>  C % MaxCoordinate  -  C % MinCoordinate )
    where ( BoxSize  >  0.0_KDR )
      PW % Wavenumber  =  nWavelengths / BoxSize
    elsewhere
      PW % Wavenumber  =  0.0_KDR
    end where

    PW % Speed = 1.0_KDR
    call PROGRAM_HEADER % GetParameter ( PW % Speed, 'Speed' )

    associate &
      ( K      =>  PW % Wavenumber, &
        Abs_K  =>  sqrt ( dot_product ( PW % Wavenumber, PW % Wavenumber ) ), &
        V      =>  PW % Speed )
    Period  =  1.0_KDR / ( Abs_K * V )
    call Show ( Period, 'Period' )

    nPeriods  =  1
    call PROGRAM_HEADER % GetParameter ( nPeriods, 'nPeriods' )
    call Show ( nPeriods, 'nPeriods' )

    I % T_Finish  =  nPeriods * Period

    call SetFluid ( PW, FC )

    end associate !-- K, etc.
    end associate !-- BoxSize
    end select !-- C
    end select !-- FC
    end select !-- FA
    end select !-- I
    end associate !-- PW

  end subroutine SetInitial


  subroutine SetReference ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    associate &
      ( PW  =>  PLANE_WAVE )
    select type ( I  =>  PW % Integrator )
      class is ( Integrator_CSA_Form )
    select type ( FC    =>  I % CurrentSet_X_A % FieldSet_C ( 1 ) % Element )
      class is ( Fluid_D_C_Form )
    select type ( FC_R  =>  PW % Reference % FieldSet_C ( 1 ) % Element )
      class is ( Fluid_D_C_Form )
    select type ( FC_D  =>  PW % Difference % FieldSet_C ( 1 ) % Element )
      class is ( Fluid_D_C_Form )

    call SetFluid ( PW, FC_R )

    associate &
      ( FV    =>  FC   % Storage_FSC % Storage % Value, &
        FV_R  =>  FC_R % Storage_FSC % Storage % Value, &
        FV_D  =>  FC_D % Storage_FSC % Storage % Value )

    call MultiplyAdd ( FV, FV_R, -1.0_KDR, FV_D )

    end associate !-- FV, etc.
    end select !-- FC_D
    end select !-- FC_R
    end select !-- FC
    end select !-- I
    end associate !-- PW

  end subroutine SetReference


  subroutine SetFluid ( PW, FC )

    class ( PlaneWaveForm ), intent ( inout ) :: &
      PW
    class ( Fluid_D_C_Form ), intent ( inout ) :: &
      FC

    associate &
      ( GC  =>  FC % Geometry_C )
    select type ( C  =>  FC % Chart )
      class is ( Chart_GS_Form )
    associate &
      ( FS  =>  FC % Storage_FSC % Storage, &
        GS  =>  GC % Storage_FSC % Storage )

    call SetFluidKernel &
           (  N = FS % Value ( :, FC % BARYON_DENSITY_C ), &
             VX = FS % Value ( :, FC % VELOCITY_U_1 ), &
             VY = FS % Value ( :, FC % VELOCITY_U_2 ), &
             VZ = FS % Value ( :, FC % VELOCITY_U_3 ), &
             PW = PW, &
             ProperCell = C % ProperCell, &
             X = GS % Value ( :, GC % CENTER_U_1 ), &
             Y = GS % Value ( :, GC % CENTER_U_2 ), &
             Z = GS % Value ( :, GC % CENTER_U_3 ), &
             K = PW % Wavenumber, &
             V = PW % Speed, &
             T = PW % Integrator % T )

    end associate !-- FS, etc.
    end select !-- C
    end associate !-- GC

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

    nV = size ( X )
    
    Abs_K = sqrt ( dot_product ( K, K ) )

    !$OMP parallel do &
    !$OMP schedule ( OMP_SCHEDULE_HOST )
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
