module LaplacianMultipole_Template

  use Basics
  use Manifolds

  implicit none
  private

  type, public, abstract :: LaplacianMultipoleTemplate
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      nRadialCells = 0, &
      nAngularMoments = 0, &
      nEquations = 0, &
      MaxDegree = 0, &  !-- Max L
      MaxOrder  = 0     !-- Max M
    integer ( KDI ) :: &
      iTimerMoments = 0, &
      iTimerClearMoments = 0, &
      iTimerLocalMoments = 0, &
      iTimerReduceMoments = 0, &
      iTimerAddMoments = 0
    integer ( KDI ) :: &
        REGULAR_COS = 1, &
      IRREGULAR_COS = 2, &
        REGULAR_SIN = 3, &
      IRREGULAR_SIN = 4
    integer ( KDI ) :: &
      iSolidHarmonics, &
      iSolidHarmonics_P1, &
      iSolidHarmonics_P2
    real ( KDR ), dimension ( 3 ) :: &
      Origin = 0.0_KDR
    real ( KDR ), dimension ( :, :, : ), pointer :: &
        M_RC => null ( ),   M_IC => null ( ), &
        M_RS => null ( ),   M_IS => null ( ), &
      MyM_RC => null ( ), MyM_IC => null ( ), &
      MyM_RS => null ( ), MyM_IS => null ( )
    logical ( KDL ) :: &
      UseDevice = .false.
    character ( LDF ) :: &
      Type = '', &
      Name = ''
    type ( StorageForm ), allocatable :: &
        Moments, &
      MyMoments
    type ( StorageForm ), dimension ( : ), allocatable :: &
      SolidHarmonics_1D
    type ( CollectiveOperation_R_Form ), allocatable :: &
      ReductionMoments
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure, public, pass :: &
      InitializeTimers
    procedure, public, pass :: &
      ComputeMoments
    procedure, public, pass :: &
      FinalizeTemplate
    procedure, private, pass :: &
      SetParameters
    procedure ( SPA ), private, pass, deferred :: &
      SetParametersAtlas
    procedure ( ASH ), private, pass, deferred :: &
      AllocateSolidHarmonics
    procedure, private, pass :: &
      AllocateMoments
    procedure ( CML ), private, pass, deferred :: &
      ComputeMomentContributions
  end type LaplacianMultipoleTemplate


  abstract interface

    subroutine SPA ( L, A )
      use Manifolds
      import LaplacianMultipoleTemplate
      implicit none
      class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
        L
      class ( AtlasHeaderForm ), intent ( in ), target :: &
        A
    end subroutine SPA

    subroutine ASH ( L )
      import LaplacianMultipoleTemplate
      implicit none
      class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
        L
    end subroutine ASH

    subroutine CML ( L, Source )
      use Basics
      import LaplacianMultipoleTemplate
      implicit none
      class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
        L
      type ( StorageForm ), intent ( in ) :: &
        Source  
    end subroutine CML

  end interface


    private :: &
      AllocateReduction, &
      AssignMomentPointers


contains


  subroutine InitializeTemplate ( L, A, MaxDegree, nEquations )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      L
    class ( AtlasHeaderForm ), intent ( in ) :: &
      A
    integer ( KDI ), intent ( in ) :: &
      MaxDegree, &
      nEquations

    L % IGNORABILITY = A % IGNORABILITY

    if ( L % Type == '' ) &
      L % Type = 'a LaplacianMultipole' 

    L % Name = 'Laplacian_' // trim ( A % Name )

    call Show ( 'Initializing ' // trim ( L % Type ), L % IGNORABILITY )
    call Show ( L % Name, 'Name', L % IGNORABILITY )

    call L % SetParameters ( A, MaxDegree, nEquations )
    call L % AllocateSolidHarmonics ( )
    call L % AllocateMoments ( )

  end subroutine InitializeTemplate


  subroutine InitializeTimers ( L, BaseLevel )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      L
    integer ( KDI ), intent ( in ) :: &
      BaseLevel

    call PROGRAM_HEADER % AddTimer &
           ( 'Moments', L % iTimerMoments, Level = BaseLevel )
      call PROGRAM_HEADER % AddTimer &
             ( 'ClearMoments', L % iTimerClearMoments, Level = BaseLevel + 1 )
      call PROGRAM_HEADER % AddTimer &
             ( 'LocalMoments', L % iTimerLocalMoments, Level = BaseLevel + 1 )
      call PROGRAM_HEADER % AddTimer &
             ( 'ReduceMoments', L % iTimerReduceMoments, &
               Level = BaseLevel + 1 )
      call PROGRAM_HEADER % AddTimer &
             ( 'AddMoments', L % iTimerAddMoments, Level = BaseLevel + 1 )

  end subroutine InitializeTimers


  subroutine ComputeMoments ( L, Source )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      L
    type ( StorageForm ), intent ( in ) :: &
      Source !-- array over levels    

    type ( TimerForm ), pointer :: &
      Timer, &
      Timer_CM, &
      Timer_LM, &
      Timer_RM, &
      Timer_AM

    Timer     =>  PROGRAM_HEADER % TimerPointer ( L % iTimerMoments )
    Timer_CM  =>  PROGRAM_HEADER % TimerPointer ( L % iTimerClearMoments )
    Timer_LM  =>  PROGRAM_HEADER % TimerPointer ( L % iTimerLocalMoments )
    Timer_RM  =>  PROGRAM_HEADER % TimerPointer ( L % iTimerReduceMoments )
    Timer_AM  =>  PROGRAM_HEADER % TimerPointer ( L % iTimerAddMoments )

    if ( associated ( Timer ) ) call Timer % Start ( )

    call Show ( 'Computing Moments', L % IGNORABILITY + 2 )

    associate ( MyM  =>  L % MyMoments )

    if ( associated ( Timer_CM ) ) call Timer_CM % Start ( )
      call Clear ( MyM % Value, UseDeviceOption = MyM % AllocatedDevice )
    if ( associated ( Timer_CM ) ) call Timer_CM % Stop ( )

    if ( associated ( Timer_LM ) ) call Timer_LM % Start ( )
    call L % ComputeMomentContributions ( Source )
    if ( associated ( Timer_LM ) ) call Timer_LM % Stop ( )

    if ( associated ( Timer_RM ) ) call Timer_RM % Start ( )
    call MyM % UpdateHost ( )
    call L % ReductionMoments % Reduce ( REDUCTION % SUM )
    if ( associated ( Timer_RM ) ) call Timer_RM % Stop ( )

call Show ( L % M_RC, '>>> L % M_RC' )

call PROGRAM_HEADER % ShowStatistics &
       ( CONSOLE % INFO_1, &
         CommunicatorOption = PROGRAM_HEADER % Communicator )
call Show ( '>>> Aborting during development', CONSOLE % ERROR )
call Show ( 'LaplacianMultipole_Template', 'module', CONSOLE % ERROR )
call Show ( 'ComputeMoments', 'subroutine', CONSOLE % ERROR )
call PROGRAM_HEADER % Abort ( )

    end associate !-- MyM

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ComputeMoments


  impure elemental subroutine FinalizeTemplate ( L )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      L

    if ( allocated ( L % ReductionMoments ) ) &
      deallocate ( L % ReductionMoments )

    if ( allocated ( L % SolidHarmonics_1D ) ) &
      deallocate ( L % SolidHarmonics_1D )

    if ( allocated ( L % MyMoments ) ) &
      deallocate ( L % MyMoments )
    if ( allocated ( L % Moments ) ) &
      deallocate ( L % Moments )

    nullify ( L % MyM_IS )
    nullify ( L % MyM_RS )
    nullify ( L % MyM_IC )
    nullify ( L % MyM_RC )
    nullify ( L % M_IS )
    nullify ( L % M_RS )
    nullify ( L % M_IC )
    nullify ( L % M_RC )

    if ( L % Name == '' ) return

    call Show ( 'Finalizing ' // trim ( L % Type ), L % IGNORABILITY )
    call Show ( L % Name, 'Name', L % IGNORABILITY )
    
  end subroutine FinalizeTemplate


  subroutine SetParameters ( L, A, MaxDegree, nEquations )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      L
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      A
    integer ( KDI ), intent ( in ) :: &
      MaxDegree, &
      nEquations

    integer ( KDI ) :: &
      iL, &
      iM

    L % MaxDegree  =  MaxDegree
    L % MaxOrder   =  MaxDegree

    call L % SetParametersAtlas ( A )

    associate &
      (  L  =>  L % MaxDegree, &
         M  =>  L % MaxOrder, &
        nA  =>  L % nAngularMoments )
    nA = 0
    do iM  =  0, M
      do iL  =  iM, L
        nA  =  nA + 1
      end do
    end do
    end associate !-- L, etc.

    L % nEquations  =  nEquations

    call Show ( L % MaxDegree, 'MaxDegree (l)', L % IGNORABILITY )
    call Show ( L % MaxOrder, 'MaxOrder (m)', L % IGNORABILITY )
    call Show ( L % nRadialCells, 'nRadialCells', L % IGNORABILITY )
    call Show ( L % nAngularMoments, 'nAngularMoments', L % IGNORABILITY )
    call Show ( L % nEquations, 'nEquations', L % IGNORABILITY )
    call Show ( L % UseDevice, 'UseDevice', L % IGNORABILITY )

  end subroutine SetParameters


  subroutine AllocateMoments ( L )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      L

    if ( allocated ( L % Moments ) ) &
      deallocate ( L % Moments )
    if ( allocated ( L % MyMoments ) ) &
      deallocate ( L % MyMoments )

    allocate ( L % Moments )
    allocate ( L % MyMoments )
    associate &
      (         M  =>  L % Moments, &
              MyM  =>  L % MyMoments, &
               nA  =>  L % nAngularMoments, &
               nR  =>  L % nRadialCells, &
               nE  =>  L % nEquations )

    call   M % Initialize ( [ nA * nR * nE, 4 ] )
    call MyM % Initialize ( [ nA * nR * nE, 4 ], PinnedOption = .true. )
      !-- 4: RegularCos, IrregularCos, RegularSin, IrregularSin
    if ( L % UseDevice ) then
      call   M % AllocateDevice ( )
      call MyM % AllocateDevice ( )
    end if

    call AllocateReduction ( L, M % Value, MyM % Value )

    call AssignMomentPointers &
           ( L, L %   Moments % Value ( :, L %   REGULAR_COS ), &
                L %   Moments % Value ( :, L % IRREGULAR_COS ), &
                L %   Moments % Value ( :, L %   REGULAR_SIN ), &
                L %   Moments % Value ( :, L % IRREGULAR_SIN ), &
                L % MyMoments % Value ( :, L %   REGULAR_COS ), &
                L % MyMoments % Value ( :, L % IRREGULAR_COS ), &
                L % MyMoments % Value ( :, L %   REGULAR_SIN ), &
                L % MyMoments % Value ( :, L % IRREGULAR_SIN ), &
                L %   M_RC, L %   M_IC, L %   M_RS, L %   M_IS, &
                L % MyM_RC, L % MyM_IC, L % MyM_RS, L % MyM_IS )

    end associate !-- M, etc.

  end subroutine AllocateMoments


  subroutine AllocateReduction ( L, M_Value, MyM_Value )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      L
    real ( KDR ), dimension ( :, : ), intent ( in ), target, contiguous :: &
        M_Value, &
      MyM_Value

    real ( KDR ), dimension ( : ), pointer :: &
        Moment_1D, &
      MyMoment_1D

      Moment_1D ( 1 : size (   M_Value ) )  =>    M_Value
    MyMoment_1D ( 1 : size ( MyM_Value ) )  =>  MyM_Value

    if ( allocated ( L % ReductionMoments ) ) &
      deallocate ( L % ReductionMoments )
    allocate ( L % ReductionMoments )
    associate &
      (  RM  =>  L % ReductionMoments, &
        PHC  =>  PROGRAM_HEADER % Communicator )
      call RM % Initialize &
             ( PHC, OutgoingValue = MyMoment_1D, IncomingValue = Moment_1D )
    end associate !-- RM, etc.

    nullify ( Moment_1D, MyMoment_1D )

  end subroutine AllocateReduction


  subroutine AssignMomentPointers &
               ( L,   M_RC_1D,   M_IC_1D,   M_RS_1D,   M_IS_1D, &
                    MyM_RC_1D, MyM_IC_1D, MyM_RS_1D, MyM_IS_1D, &
                      M_RC_3D,   M_IC_3D,   M_RS_3D,   M_IS_3D, &
                    MyM_RC_3D, MyM_IC_3D, MyM_RS_3D, MyM_IS_3D )

    class ( LaplacianMultipoleTemplate ), intent ( in ) :: &
      L
    real ( KDR ), dimension ( : ), intent ( in ), target :: &
        M_RC_1D,   M_IC_1D, &
        M_RS_1D,   M_IS_1D, &
      MyM_RC_1D, MyM_IC_1D, &
      MyM_RS_1D, MyM_IS_1D
    real ( KDR ), dimension ( :, :, : ), intent ( out ), pointer :: &
        M_RC_3D,   M_IC_3D, &
        M_RS_3D,   M_IS_3D, &
      MyM_RC_3D, MyM_IC_3D, &
      MyM_RS_3D, MyM_IS_3D

    associate &
      ( nE => L % nEquations, &
        nR => L % nRadialCells, &
        nA => L % nAngularMoments )

      M_RC_3D ( 1 : nA, 1 : nR, 1 : nE )  =>    M_RC_1D
      M_IC_3D ( 1 : nA, 1 : nR, 1 : nE )  =>    M_IC_1D
      M_RS_3D ( 1 : nA, 1 : nR, 1 : nE )  =>    M_RS_1D
      M_IS_3D ( 1 : nA, 1 : nR, 1 : nE )  =>    M_IS_1D
    MyM_RC_3D ( 1 : nA, 1 : nR, 1 : nE )  =>  MyM_RC_1D
    MyM_IC_3D ( 1 : nA, 1 : nR, 1 : nE )  =>  MyM_IC_1D
    MyM_RS_3D ( 1 : nA, 1 : nR, 1 : nE )  =>  MyM_RS_1D
    MyM_IS_3D ( 1 : nA, 1 : nR, 1 : nE )  =>  MyM_IS_1D

    end associate !-- nR, nA

  end subroutine AssignMomentPointers


end module LaplacianMultipole_Template
