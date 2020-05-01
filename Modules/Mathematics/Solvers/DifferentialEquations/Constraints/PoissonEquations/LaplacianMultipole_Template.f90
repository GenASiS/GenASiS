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
        Moment_RC => null ( ),   Moment_IC => null ( ), &
        Moment_RS => null ( ),   Moment_IS => null ( ), &
      MyMoment_RC => null ( ), MyMoment_IC => null ( ), &
      MyMoment_RS => null ( ), MyMoment_IS => null ( )
    logical ( KDL ) :: &
      UseDevice = .false.
    character ( LDF ) :: &
      Type = '', &
      Name = ''
    type ( StorageForm ), allocatable :: &
        Moments, &
      MyMoments
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
    procedure ( ARC ), private, pass, deferred :: &
      AllocateRectangularCoordinates
    procedure ( ASH ), private, pass, deferred :: &
      AllocateSolidHarmonics
    procedure, private, pass :: &
      AllocateMoments
    procedure, private, pass :: &
      ComputeMomentContributions
    procedure ( CSH_0_0 ), private, pass, deferred :: &
      ComputeSolidHarmonics_0_0
    procedure ( CSH_iM_iM ), private, pass, deferred :: &
      ComputeSolidHarmonics_iM_iM
    procedure ( CSH_iL_iM_1 ), private, pass, deferred :: &
      ComputeSolidHarmonics_iL_iM_1
    procedure ( CSH_iL_iM_2 ), private, pass, deferred :: &
      ComputeSolidHarmonics_iL_iM_2
    procedure ( SMC ), private, pass, deferred :: &
      SumMomentContributions
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

    subroutine ARC ( L )
      import LaplacianMultipoleTemplate
      implicit none
      class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
        L
    end subroutine ARC

    subroutine ASH ( L )
      import LaplacianMultipoleTemplate
      implicit none
      class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
        L
    end subroutine ASH

    subroutine CSH_0_0 ( L, iSH_0, iSH_PD )
      use Basics
      import LaplacianMultipoleTemplate
      implicit none
      class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
        L
      integer ( KDI ), intent ( in ) :: &
        iSH_0, iSH_PD
    end subroutine CSH_0_0

    subroutine CSH_iM_iM ( L, iM, iSH_0, iSH_PD )
      use Basics
      import LaplacianMultipoleTemplate
      implicit none
      class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
        L
      integer ( KDI ), intent ( in ) :: &
        iM, &
        iSH_0, iSH_PD
    end subroutine CSH_iM_iM

    subroutine CSH_iL_iM_1 ( L, iM, iSH_0, iSH_1 )
      use Basics
      import LaplacianMultipoleTemplate
      implicit none
      class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
        L
      integer ( KDI ), intent ( in ) :: &
        iM, &
        iSH_0, iSH_1
    end subroutine CSH_iL_iM_1

    subroutine CSH_iL_iM_2 ( L, iL, iM, iSH_0, iSH_1, iSH_2 )
      use Basics
      import LaplacianMultipoleTemplate
      implicit none
      class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
        L
      integer ( KDI ), intent ( in ) :: &
        iL, iM, &
        iSH_0, iSH_1, iSH_2
    end subroutine CSH_iL_iM_2

    subroutine SMC ( L, Source, iA, iSH_0 )
      use Basics
      import LaplacianMultipoleTemplate
      implicit none
      class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
        L
      class ( * ), intent ( in ) :: &
        Source
      integer ( KDI ), intent ( in ) :: &
        iA, &  
        iSH_0  
    end subroutine SMC

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
    call L % AllocateRectangularCoordinates ( )
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

    associate &
      (   M  =>  L %   Moments, &
        MyM  =>  L % MyMoments )

    if ( associated ( Timer_CM ) ) call Timer_CM % Start ( )
      call Clear ( MyM % Value, UseDeviceOption = MyM % AllocatedDevice )
    if ( associated ( Timer_CM ) ) call Timer_CM % Stop ( )

    if ( associated ( Timer_LM ) ) call Timer_LM % Start ( )
    call L % ComputeMomentContributions ( Source )
    if ( associated ( Timer_LM ) ) call Timer_LM % Stop ( )

    if ( associated ( Timer_RM ) ) call Timer_RM % Start ( )
    call MyM % UpdateHost ( )
    call L % ReductionMoments % Reduce ( REDUCTION % SUM )
    call M % UpdateDevice ( ) 
    if ( associated ( Timer_RM ) ) call Timer_RM % Stop ( )

call Show ( M % Value ( 1 : 10, 1 ), '>>> M % Value ( 1 : 10, 1 )' )

call PROGRAM_HEADER % ShowStatistics &
       ( CONSOLE % INFO_1, &
         CommunicatorOption = PROGRAM_HEADER % Communicator )
call Show ( '>>> Aborting during development', CONSOLE % ERROR )
call Show ( 'LaplacianMultipole_Template', 'module', CONSOLE % ERROR )
call Show ( 'ComputeMoments', 'subroutine', CONSOLE % ERROR )
call PROGRAM_HEADER % Abort ( )

    end associate !-- M, etc.

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ComputeMoments


  impure elemental subroutine FinalizeTemplate ( L )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      L

    if ( allocated ( L % ReductionMoments ) ) &
      deallocate ( L % ReductionMoments )

    if ( allocated ( L % MyMoments ) ) &
      deallocate ( L % MyMoments )
    if ( allocated ( L % Moments ) ) &
      deallocate ( L % Moments )

    nullify ( L % MyMoment_IS )
    nullify ( L % MyMoment_RS )
    nullify ( L % MyMoment_IC )
    nullify ( L % MyMoment_RC )
    nullify ( L % Moment_IS )
    nullify ( L % Moment_RS )
    nullify ( L % Moment_IC )
    nullify ( L % Moment_RC )

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
                L %   Moment_RC, L %   Moment_IC, &
                L %   Moment_RS, L %   Moment_IS, &
                L % MyMoment_RC, L % MyMoment_IC, &
                L % MyMoment_RS, L % MyMoment_IS )

    end associate !-- M, etc.

  end subroutine AllocateMoments


  subroutine ComputeMomentContributions ( L, Source )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      L
    class ( * ), intent ( in ) :: &
      Source  

    integer ( KDI ) :: &
      iA, &   !-- iAngularMoment
      iM, &   !-- iOrder
      iL, &   !-- iDegree
      iSH_PD  !-- iSolidHarmonic_PreviousDiagonal
    integer ( KDI ), pointer :: &
      iSH_0, &  !-- iSolidHarmonic_Current
      iSH_1, &  !-- iSolidHarmonic_Previous_1
      iSH_2     !-- iSolidHarmonic_Previous_2
    integer ( KDI ), dimension ( 3 ), target :: &
      iSH

    iSH_0  =>  iSH ( 1 )
    iSH_1  =>  iSH ( 2 )
    iSH_2  =>  iSH ( 3 )

        iA  =  1
       iSH  =  [ 1, 2, 3 ]
    iSH_PD  =  4

    do iM  =  0, L % MaxOrder

      !-- ( L, M ) = ( iM, iM )
      !-- Note iL = iM
      if ( iM  ==  0 ) then
        call L % ComputeSolidHarmonics_0_0 ( iSH_0, iSH_PD )
      else
        call L % ComputeSolidHarmonics_iM_iM ( iM, iSH_0, iSH_PD )
      end if

      do iL  =  iM, L % MaxDegree

        if ( iL  ==  iM + 1 ) then
          !-- ( L, M ) = ( iM + 1, iM )
          !-- Note iL = iM + 1
          call L % ComputeSolidHarmonics_iL_iM_1 &
                 ( iM, iSH_0, iSH_1 )
        else if ( iL  >=  iM  +  2 ) then
          !-- ( L, M ) = ( iL, iM )
          call L % ComputeSolidHarmonics_iL_iM_2 &
                 ( iL, iM, iSH_0, iSH_1, iSH_2 )
          end if

        call L % SumMomentContributions ( Source, iA, iSH_0 )

         iA  =  iA + 1
        iSH  =  cshift ( iSH, -1 )

      end do !-- iL
    end do !-- iM

  end subroutine ComputeMomentContributions


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
