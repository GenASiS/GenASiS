module LaplacianMultipoleOld_2__Template

  use Basics
  use Manifolds

  implicit none
  private

  type, public, abstract :: LaplacianMultipoleOld_2_Template
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      nRadialCells = 0, &
      nAngularMoments = 0, &
      nEquations = 0, &
      MaxDegree = 0, &  !-- Max L
      MaxOrder  = 0     !-- Max M
    integer ( KDI ) :: &
      iTimerComputeMoments = 0, &
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
      UseDevice = .false., &
      ReductionUseDevice = .false.
    character ( LDF ) :: &
      Type = '', &
      Name = ''
    type ( StorageForm ), allocatable :: &
        Moments, &
      MyMoments, &
      RadialEdges
    type ( CollectiveOperation_R_Form ), allocatable :: &
      ReductionMoments
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure, public, pass :: &
      InitializeTimers
    procedure, public, pass :: &
      ComputeMoments
    procedure ( CSH_0_0 ), public, pass, deferred :: &
      ComputeSolidHarmonics_0_0
    procedure ( CSH_iM_iM ), public, pass, deferred :: &
      ComputeSolidHarmonics_iM_iM
    procedure ( CSH_iL_iM_1 ), public, pass, deferred :: &
      ComputeSolidHarmonics_iL_iM_1
    procedure ( CSH_iL_iM_2 ), public, pass, deferred :: &
      ComputeSolidHarmonics_iL_iM_2
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
      ComputeMomentsLocal
    procedure ( CMLA ), private, pass, deferred :: &
      ComputeMomentLocalAtlas
  end type LaplacianMultipoleOld_2_Template


  abstract interface

    subroutine CSH_0_0 ( L, iSH_0, iSH_PD )
      use Basics
      import LaplacianMultipoleOld_2_Template
      implicit none
      class ( LaplacianMultipoleOld_2_Template ), intent ( inout ) :: &
        L
      integer ( KDI ), intent ( in ) :: &
        iSH_0, iSH_PD
    end subroutine CSH_0_0

    subroutine CSH_iM_iM ( L, iM, iSH_0, iSH_PD )
      use Basics
      import LaplacianMultipoleOld_2_Template
      implicit none
      class ( LaplacianMultipoleOld_2_Template ), intent ( inout ) :: &
        L
      integer ( KDI ), intent ( in ) :: &
        iM, &
        iSH_0, iSH_PD
    end subroutine CSH_iM_iM

    subroutine CSH_iL_iM_1 ( L, iM, iSH_0, iSH_1 )
      use Basics
      import LaplacianMultipoleOld_2_Template
      implicit none
      class ( LaplacianMultipoleOld_2_Template ), intent ( inout ) :: &
        L
      integer ( KDI ), intent ( in ) :: &
        iM, &
        iSH_0, iSH_1
    end subroutine CSH_iL_iM_1

    subroutine CSH_iL_iM_2 ( L, iL, iM, iSH_0, iSH_1, iSH_2 )
      use Basics
      import LaplacianMultipoleOld_2_Template
      implicit none
      class ( LaplacianMultipoleOld_2_Template ), intent ( inout ) :: &
        L
      integer ( KDI ), intent ( in ) :: &
        iL, iM, &
        iSH_0, iSH_1, iSH_2
    end subroutine CSH_iL_iM_2

    subroutine SPA ( L, A )
      use Manifolds
      import LaplacianMultipoleOld_2_Template
      implicit none
      class ( LaplacianMultipoleOld_2_Template ), intent ( inout ) :: &
        L
      class ( AtlasHeaderForm ), intent ( in ), target :: &
        A
    end subroutine SPA

    subroutine ARC ( L )
      import LaplacianMultipoleOld_2_Template
      implicit none
      class ( LaplacianMultipoleOld_2_Template ), intent ( inout ) :: &
        L
    end subroutine ARC

    subroutine ASH ( L )
      import LaplacianMultipoleOld_2_Template
      implicit none
      class ( LaplacianMultipoleOld_2_Template ), intent ( inout ) :: &
        L
    end subroutine ASH

    subroutine CMLA ( L, Source, iA, iSH_0 )
      use Basics
      use Manifolds
      import LaplacianMultipoleOld_2_Template
      implicit none
      class ( LaplacianMultipoleOld_2_Template ), intent ( inout ) :: &
        L
      class ( FieldAtlasTemplate ), intent ( in ) :: &
        Source
      integer ( KDI ), intent ( in ) :: &
        iA, &  
        iSH_0  
    end subroutine CMLA

  end interface

!-- FIXME: With GCC 6.1.0, must be public to trigger .smod generation
!    private :: &
    public :: &
      AddMomentShellsKernel

    interface

      module subroutine AddMomentShellsKernel &
                          ( M_RC, M_IC, M_RS, M_IS, nA, nE, nR, &
                            UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, : ) :: &
          M_RC, M_IC, M_RS, M_IS
        integer ( KDI ), intent ( in ) :: &
          nA, nE, nR
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine AddMomentShellsKernel

    end interface


    private :: &
      AllocateReduction, &
      AssignMomentPointers


contains


  subroutine InitializeTemplate ( L, A, MaxDegree, nEquations )

    class ( LaplacianMultipoleOld_2_Template ), intent ( inout ) :: &
      L
    class ( AtlasHeaderForm ), intent ( in ) :: &
      A
    integer ( KDI ), intent ( in ) :: &
      MaxDegree, &
      nEquations

    L % IGNORABILITY = A % IGNORABILITY

    if ( L % Type == '' ) &
      L % Type = 'a LaplacianMultipoleOld_2' 

    L % Name = 'Laplacian_' // trim ( A % Name )

    call Show ( 'Initializing ' // trim ( L % Type ), L % IGNORABILITY )
    call Show ( L % Name, 'Name', L % IGNORABILITY )

    call L % SetParameters ( A, MaxDegree, nEquations )
    call L % AllocateRectangularCoordinates ( )
    call L % AllocateSolidHarmonics ( )
    call L % AllocateMoments ( )

  end subroutine InitializeTemplate


  subroutine InitializeTimers ( L, BaseLevel )

    class ( LaplacianMultipoleOld_2_Template ), intent ( inout ) :: &
      L
    integer ( KDI ), intent ( in ) :: &
      BaseLevel

    call PROGRAM_HEADER % AddTimer &
           ( 'ComputeMoments', L % iTimerComputeMoments, Level = BaseLevel )
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

    class ( LaplacianMultipoleOld_2_Template ), intent ( inout ) :: &
      L
    class ( FieldAtlasTemplate ), intent ( in ) :: &
      Source

    type ( TimerForm ), pointer :: &
      Timer, &
      Timer_CM, &
      Timer_LM, &
      Timer_RM, &
      Timer_AM

    Timer     =>  PROGRAM_HEADER % TimerPointer ( L % iTimerComputeMoments )
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
    call MyM % Clear ( )
    if ( associated ( Timer_CM ) ) call Timer_CM % Stop ( )

    if ( associated ( Timer_LM ) ) call Timer_LM % Start ( )
    call L % ComputeMomentsLocal ( Source )
    if ( associated ( Timer_LM ) ) call Timer_LM % Stop ( )

    if ( associated ( Timer_RM ) ) call Timer_RM % Start ( )
    if ( .not. L % ReductionUseDevice ) call MyM % UpdateHost ( ) 
    call L % ReductionMoments % Reduce ( REDUCTION % SUM )
    if ( .not. L % ReductionUseDevice ) call M % UpdateDevice ( ) 
    if ( associated ( Timer_RM ) ) call Timer_RM % Stop ( )

    if ( associated ( Timer_AM ) ) call Timer_AM % Start ( )
      call AddMomentShellsKernel &
             ( L % Moment_RC, L % Moment_IC, L % Moment_RS, L % Moment_IS, &
               L % nAngularMoments, L % nEquations, L % nRadialCells, &
               UseDeviceOption = L % UseDevice )
    if ( associated ( Timer_AM ) ) call Timer_AM % Stop ( )
    
    end associate !-- M, etc.

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ComputeMoments


  impure elemental subroutine FinalizeTemplate ( L )

    class ( LaplacianMultipoleOld_2_Template ), intent ( inout ) :: &
      L

    if ( allocated ( L % ReductionMoments ) ) &
      deallocate ( L % ReductionMoments )

    if ( allocated ( L % RadialEdges ) ) &
      deallocate ( L % RadialEdges )
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

    class ( LaplacianMultipoleOld_2_Template ), intent ( inout ) :: &
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
    call Show ( L % ReductionUseDevice, 'ReductionUseDevice', L % IGNORABILITY )

  end subroutine SetParameters


  subroutine AllocateMoments ( L )

    class ( LaplacianMultipoleOld_2_Template ), intent ( inout ) :: &
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

    call   M % Initialize &
             ( [ nA * nR * nE, 4 ], PinnedOption = L % UseDevice )
    call MyM % Initialize &
             ( [ nA * nR * nE, 4 ], PinnedOption = L % UseDevice )
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


  subroutine ComputeMomentsLocal ( L, Source )

    class ( LaplacianMultipoleOld_2_Template ), intent ( inout ) :: &
      L
    class ( FieldAtlasTemplate ), intent ( in ) :: &
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

        call L % ComputeMomentLocalAtlas ( Source, iA, iSH_0 )

         iA  =  iA + 1
        iSH  =  cshift ( iSH, -1 )

      end do !-- iL
    end do !-- iM

  end subroutine ComputeMomentsLocal


  subroutine AllocateReduction ( L, M_Value, MyM_Value )

    class ( LaplacianMultipoleOld_2_Template ), intent ( inout ) :: &
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
      if ( L % ReductionUseDevice ) &
        call RM % AllocateDevice ( )
    end associate !-- RM, etc.

    nullify ( Moment_1D, MyMoment_1D )

  end subroutine AllocateReduction


  subroutine AssignMomentPointers &
               ( L,   M_RC_1D,   M_IC_1D,   M_RS_1D,   M_IS_1D, &
                    MyM_RC_1D, MyM_IC_1D, MyM_RS_1D, MyM_IS_1D, &
                      M_RC_3D,   M_IC_3D,   M_RS_3D,   M_IS_3D, &
                    MyM_RC_3D, MyM_IC_3D, MyM_RS_3D, MyM_IS_3D )

    class ( LaplacianMultipoleOld_2_Template ), intent ( in ) :: &
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
      ( nA => L % nAngularMoments, &
        nE => L % nEquations, &
        nR => L % nRadialCells )

      M_RC_3D ( 1 : nR, 1 : nE, 1 : nA )  =>    M_RC_1D
      M_IC_3D ( 1 : nR, 1 : nE, 1 : nA )  =>    M_IC_1D
      M_RS_3D ( 1 : nR, 1 : nE, 1 : nA )  =>    M_RS_1D
      M_IS_3D ( 1 : nR, 1 : nE, 1 : nA )  =>    M_IS_1D
    MyM_RC_3D ( 1 : nR, 1 : nE, 1 : nA )  =>  MyM_RC_1D
    MyM_IC_3D ( 1 : nR, 1 : nE, 1 : nA )  =>  MyM_IC_1D
    MyM_RS_3D ( 1 : nR, 1 : nE, 1 : nA )  =>  MyM_RS_1D
    MyM_IS_3D ( 1 : nR, 1 : nE, 1 : nA )  =>  MyM_IS_1D

    end associate !-- nR, nA

  end subroutine AssignMomentPointers


end module LaplacianMultipoleOld_2__Template
