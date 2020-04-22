module LaplacianMultipole_Template

  use Basics
  use Manifolds

  implicit none
  private

  type, public, abstract :: LaplacianMultipoleTemplate
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      iTimerMoments = 0, &
      iTimerClearMoments = 0, &
      iTimerLocalMoments = 0, &
      iTimerReduceMoments = 0, &
      iTimerAddMoments = 0, &
      nRadialCells = 0, &
      nAngularMomentCells = 0, &
      nEquations = 0, &
      MaxDegree = 0, &  !-- Max L
      MaxOrder  = 0     !-- Max M
    integer ( KDI ) :: &
      REGULAR_COSINE    = 1, &
      IRREGULAR_COSINE  = 2, &
      REGULAR_SINE      = 3, &
      IRREGULAR_SINE    = 4, &
      N_SOLID_HARMONICS = 0, & 
      DELTA_FACTOR      = 0
    real ( KDR ), dimension ( 3 ) :: &
      Origin = 0.0_KDR
    real ( KDR ), dimension ( : ), pointer :: &
      RadialEdge       => null ( ), &
      SolidHarmonic_RC => null ( ), &  !-- Regular Cos
      SolidHarmonic_IC => null ( ), &  !-- Irregular Cos
      SolidHarmonic_RS => null ( ), &  !-- Regular Sin
      SolidHarmonic_IS => null ( ), &  !-- Irregular Sin
      Delta            => null ( ), &
      MyMoment_RC_1D   => null ( ), &
      MyMoment_IC_1D   => null ( ), &
      MyMoment_RS_1D   => null ( ), &
      MyMoment_IS_1D   => null ( ), &
      Moment_RC_1D     => null ( ), &    
      Moment_IC_1D     => null ( ), &
      Moment_RS_1D     => null ( ), &
      Moment_IS_1D     => null ( )
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      MyM_RC => null ( ), MyM_IC => null ( ), &
      MyM_RS => null ( ), MyM_IS => null ( ), &
      M_RC   => null ( ), M_IC   => null ( ), &
      M_RS   => null ( ), M_IS   => null ( )
    character ( LDF ) :: &
      Type = '', &
      Name = ''
    type ( StorageForm ), allocatable :: &
      RadialEdges, &
      SolidHarmonics, &
      MyMoments, &
      Moments
    type ( CollectiveOperation_R_Form ), allocatable :: &
      ReductionMoments
    !-- FIXME: This does not work with GCC 6.1.0
!    procedure ( ), pointer, nopass :: &
!      ComputeSolidHarmonicsKernel => ComputeSolidHarmonicsKernel
!    procedure ( ), pointer, nopass :: &
!      ComputeMomentContributionsKernel => ComputeMomentContributionsKernel
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure, public, pass :: &
      InitializeTimers
    procedure ( SP ), private, pass, deferred :: &
      SetParameters   
    procedure, public, pass :: &
      ComputeMoments
    procedure, public, pass :: &
      FinalizeTemplate
    procedure, public, pass :: &
      SetMomentStorage
    procedure ( CML ), private, pass, deferred :: &
      ComputeMomentsLocal
  end type LaplacianMultipoleTemplate

  abstract interface

    subroutine SP ( LM, A, MaxDegree, nEquations )
      use Basics
      use Manifolds
      import LaplacianMultipoleTemplate
      implicit none
      class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
        LM
      class ( AtlasHeaderForm ), intent ( in ), target :: &
        A
      integer ( KDI ), intent ( in ) :: &
        MaxDegree, &
        nEquations
    end subroutine SP

    subroutine CML ( LM, Source )
      use Basics
      import LaplacianMultipoleTemplate
      implicit none
      class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
        LM
      type ( StorageForm ), intent ( in ) :: &
        Source  
    end subroutine CML

  end interface

!-- FIXME: With GCC 6.1.0, must be public to trigger .smod generation
!    private :: &
    public :: &
      ComputeMomentContributionsKernel, &
      AddMomentShellsKernel, &
      ComputeSolidHarmonicsKernel

!-- FIXME: With GCC 6.1.0, must be public to trigger .smod generation
!      private :: &
      public :: &
        ComputeMomentContributions_MR_MI_Kernel, &
        ComputeSolidHarmonics_C_M_0_Kernel, &
        ComputeSolidHarmonics_C_S_Kernel

    interface 

      module subroutine ComputeMomentContributionsKernel &
                          ( MyM_RC, MyM_IC, MyM_RS, MyM_IS, &
                            SH_RC, SH_IC, SH_RS, SH_IS, &
                            Source, Volume, iaSource, L, nE, nA, iR )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
          MyM_RC, MyM_IC, &  !-- MyMoment_RealCosine, etc.
          MyM_RS, MyM_IS
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          SH_RC, SH_IC, &  !-- SolidHarmonic_RealCosine, etc.
          SH_RS, SH_IS, &
          Source
        real ( KDR ), intent ( in ) :: &
          Volume
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          iaSource
        integer ( KDI ), intent ( in ) :: &
          L,  &  !-- MaxOrder
          nE, &  !-- nEquations
          nA, &  !-- nAngular
          iR     !-- iRadius
      end subroutine ComputeMomentContributionsKernel

      module subroutine AddMomentShellsKernel ( MR, MI, nE, nR, nA )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
          MR, MI
        integer ( KDI ), intent ( in ) :: &
          nE, nR, nA
      end subroutine AddMomentShellsKernel
        
      module subroutine ComputeSolidHarmonicsKernel &
                          ( CoordinateSystem, Position, Origin, RadialEdge, &
                            L, nDimensions, GridError, R_C, I_C, R_S, I_S, &
                            R, iR )
        use Basics
        implicit none
        character ( * ), intent ( in ) :: &
          CoordinateSystem
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          Position, &
          Origin, &
          RadialEdge
        integer ( KDI ), intent ( in ) :: &
          nDimensions, &
          L
        logical ( KDL ), intent ( out ) :: &
          GridError
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          R_C, I_C, R_S, I_S
        real ( KDR ), intent ( out ) :: &
          R
        integer ( KDI ), intent ( out ) :: &
          iR
      end subroutine ComputeSolidHarmonicsKernel

      module subroutine ComputeMomentContributions_MR_MI_Kernel &
                          ( MyMR, MyMI, SH_R, SH_I, Source_dV, nE, nA, iRS )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
          MyMR, MyMI
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          SH_R, SH_I
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          Source_dV
        integer ( KDI ), intent ( in ) :: &
          nE, &
          nA, &
          iRS
      end subroutine ComputeMomentContributions_MR_MI_Kernel

      module subroutine ComputeSolidHarmonics_C_M_0_Kernel &
                          ( X, Z, L, R_C, I_C )
        use Basics
        implicit none
        real ( KDR ), intent ( in ) :: &
          X, Z
        integer ( KDI ), intent ( in ) :: &
          L
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          R_C, I_C
      end subroutine ComputeSolidHarmonics_C_M_0_Kernel

      module subroutine ComputeSolidHarmonics_C_S_Kernel &
               ( X, Y, Z, L, R_C, I_C, R_S, I_S )
        use Basics
        implicit none
        real ( KDR ), intent ( in ) :: &
          X, Y, Z
        integer ( KDI ), intent ( in ) :: &
          L
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          R_C, I_C, &
          R_S, I_S
      end subroutine ComputeSolidHarmonics_C_S_Kernel

    end interface

    private :: &
      SetReduction, &
      AssignPointers

contains


  subroutine InitializeTemplate ( LM, A, MaxDegree, nEquations )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      LM
    class ( AtlasHeaderForm ), intent ( in ) :: &
      A
    integer ( KDI ), intent ( in ) :: &
      MaxDegree, &
      nEquations

    integer ( KDI ) :: &
      iV, &  !-- iValue
      iL, &
      iM

    LM % IGNORABILITY = A % IGNORABILITY

    if ( LM % Type == '' ) &
      LM % Type = 'a LaplacianMultipole' 

    LM % Name = 'Laplacian_' // trim ( A % Name )

    call Show ( 'Initializing ' // trim ( LM % Type ), LM % IGNORABILITY )
    call Show ( LM % Name, 'Name', LM % IGNORABILITY )

    call LM % SetParameters ( A, MaxDegree, nEquations )

    if ( MaxDegree > 0 ) then
      LM % N_SOLID_HARMONICS  =  4  !-- Regular and Irregular Cosine and Sine
    else
      LM % N_SOLID_HARMONICS  =  2  !-- Regular and Irregular Cosine
    end if
    LM % DELTA_FACTOR  =  LM % N_SOLID_HARMONICS + 1

    allocate ( LM % SolidHarmonics )
    associate ( SH  =>  LM % SolidHarmonics )

      call SH % Initialize &
             ( [ LM % nAngularMomentCells, LM % N_SOLID_HARMONICS + 1 ] )

      LM % SolidHarmonic_RC  =>  SH % Value ( :, LM % REGULAR_COSINE ) 
      LM % SolidHarmonic_IC  =>  SH % Value ( :, LM % IRREGULAR_COSINE ) 

      if ( LM % MaxOrder > 0 ) then
        LM % SolidHarmonic_RS  =>  SH % Value ( :, LM % REGULAR_SINE ) 
        LM % SolidHarmonic_IS  =>  SH % Value ( :, LM % IRREGULAR_SINE ) 
      end if

      LM % Delta  =>  SH % Value ( :, LM % DELTA_FACTOR )

    end associate !-- SH
    
    associate &
      ( L => LM % MaxDegree, &
        M => LM % MaxOrder )
    iV = 0
    do iM = 0, M
      do iL = iM, L
        iV = iV + 1
        if ( iM == 0 ) then
          LM % Delta ( iV ) = 1.0_KDR
        else
          LM % Delta ( iV ) = 2.0_KDR
        end if
      end do
    end do
    end associate !-- L, etc.

  end subroutine InitializeTemplate


  subroutine InitializeTimers ( LM, BaseLevel )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      LM
    integer ( KDI ), intent ( in ) :: &
      BaseLevel

    call PROGRAM_HEADER % AddTimer &
           ( 'Moments', LM % iTimerMoments, Level = BaseLevel )
      call PROGRAM_HEADER % AddTimer &
             ( 'ClearMoments', LM % iTimerClearMoments, Level = BaseLevel + 1 )
      call PROGRAM_HEADER % AddTimer &
             ( 'LocalMoments', LM % iTimerLocalMoments, Level = BaseLevel + 1 )
      call PROGRAM_HEADER % AddTimer &
             ( 'ReduceMoments', LM % iTimerReduceMoments, &
               Level = BaseLevel + 1 )
      call PROGRAM_HEADER % AddTimer &
             ( 'AddMoments', LM % iTimerAddMoments, Level = BaseLevel + 1 )

  end subroutine InitializeTimers


  subroutine ComputeMoments ( LM, Source )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      LM
    type ( StorageForm ), intent ( in ) :: &
      Source !-- array over levels    
    
    ! integer ( KDI ) :: &
    !   iL, &  !-- iLevel
    !   iC     !-- iCell
    type ( TimerForm ), pointer :: &
      Timer, &
      Timer_CM, &
      Timer_LM, &
      Timer_RM, &
      Timer_AM

    Timer     =>  PROGRAM_HEADER % TimerPointer ( LM % iTimerMoments )
    Timer_CM  =>  PROGRAM_HEADER % TimerPointer ( LM % iTimerClearMoments )
    Timer_LM  =>  PROGRAM_HEADER % TimerPointer ( LM % iTimerLocalMoments )
    Timer_RM  =>  PROGRAM_HEADER % TimerPointer ( LM % iTimerReduceMoments )
    Timer_AM  =>  PROGRAM_HEADER % TimerPointer ( LM % iTimerAddMoments )

    if ( associated ( Timer ) ) call Timer % Start ( )

    call Show ( 'Computing Moments', CONSOLE % INFO_3 )

    if ( associated ( Timer_CM ) ) call Timer_CM % Start ( )
    call Clear ( LM % MyM_RC )
    call Clear ( LM % MyM_IC )
    if ( LM % MaxOrder > 0 ) then
      call Clear ( LM % MyM_RS )
      call Clear ( LM % MyM_IS )
    end if
    if ( associated ( Timer_CM ) ) call Timer_CM % Stop ( )

    if ( associated ( Timer_LM ) ) call Timer_LM % Start ( )
    call LM % ComputeMomentsLocal ( Source )
    if ( associated ( Timer_LM ) ) call Timer_LM % Stop ( )

    if ( associated ( Timer_RM ) ) call Timer_RM % Start ( )
    call LM % ReductionMoments % Reduce ( REDUCTION % SUM )
    if ( associated ( Timer_RM ) ) call Timer_RM % Stop ( )

    if ( associated ( Timer_AM ) ) call Timer_AM % Start ( )
    call AddMomentShellsKernel &
           ( LM % M_RC, LM % M_IC, &
             LM % nEquations, LM % nRadialCells, LM % nAngularMomentCells )
    if ( LM % MaxOrder > 0 ) &
      call AddMomentShellsKernel &
             ( LM % M_RS, LM % M_IS, &
               LM % nEquations, LM % nRadialCells, LM % nAngularMomentCells )
    if ( associated ( Timer_AM ) ) call Timer_AM % Stop ( )

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ComputeMoments


  impure elemental subroutine FinalizeTemplate ( LM )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      LM

    if ( allocated ( LM % ReductionMoments ) ) &
      deallocate ( LM % ReductionMoments )

    if ( allocated ( LM % MyMoments ) ) deallocate ( LM % MyMoments )
    if ( allocated ( LM % Moments ) ) deallocate ( LM % Moments )
    if ( allocated ( LM % SolidHarmonics ) ) deallocate ( LM % SolidHarmonics )
    if ( allocated ( LM % RadialEdges ) ) deallocate ( LM % RadialEdges )

    nullify ( LM % M_IS )
    nullify ( LM % M_RS )
    nullify ( LM % M_IC )
    nullify ( LM % M_RC )
    nullify ( LM % MyM_IS )
    nullify ( LM % MyM_RS )
    nullify ( LM % MyM_IC )
    nullify ( LM % MyM_RC )

    nullify ( LM % Moment_IS_1D )
    nullify ( LM % Moment_IC_1D )
    nullify ( LM % Moment_RS_1D )
    nullify ( LM % Moment_RC_1D )
    nullify ( LM % MyMoment_IS_1D )
    nullify ( LM % MyMoment_IC_1D )
    nullify ( LM % MyMoment_RS_1D )
    nullify ( LM % MyMoment_RC_1D )

    nullify ( LM % Delta )
    nullify ( LM % SolidHarmonic_IS )
    nullify ( LM % SolidHarmonic_RS )
    nullify ( LM % SolidHarmonic_IC )
    nullify ( LM % SolidHarmonic_RC )

    nullify ( LM % RadialEdge )

    if ( LM % Name == '' ) return

    call Show ( 'Finalizing ' // trim ( LM % Type ), LM % IGNORABILITY )
    call Show ( LM % Name, 'Name', LM % IGNORABILITY )
    
  end subroutine FinalizeTemplate


  subroutine SetMomentStorage ( LM )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      LM

    if ( allocated ( LM % MyMoments ) ) deallocate ( LM % MyMoments )
    if ( allocated ( LM % Moments ) ) deallocate ( LM % Moments )

    allocate ( LM % Moments )
    allocate ( LM % MyMoments )
    associate &
      (        M  =>  LM % Moments, &
             MyM  =>  LM % MyMoments, &
        nA_nR_nE  =>  LM % nAngularMomentCells  *  LM % nRadialCells &
                      *  LM % nEquations )
    call   M % Initialize ( [ nA_nR_nE, LM % N_SOLID_HARMONICS ] )
    call MyM % Initialize ( [ nA_nR_nE, LM % N_SOLID_HARMONICS ] )

      LM % Moment_RC_1D  =>    M % Value ( :, LM % REGULAR_COSINE ) 
      LM % Moment_IC_1D  =>    M % Value ( :, LM % IRREGULAR_COSINE ) 
    LM % MyMoment_RC_1D  =>  MyM % Value ( :, LM % REGULAR_COSINE ) 
    LM % MyMoment_IC_1D  =>  MyM % Value ( :, LM % IRREGULAR_COSINE ) 
    if ( LM % MaxOrder > 0 ) then
        LM % Moment_RS_1D  =>    M % Value ( :, LM % REGULAR_SINE ) 
        LM % Moment_IS_1D  =>    M % Value ( :, LM % IRREGULAR_SINE ) 
      LM % MyMoment_RS_1D  =>  MyM % Value ( :, LM % REGULAR_SINE ) 
      LM % MyMoment_IS_1D  =>  MyM % Value ( :, LM % IRREGULAR_SINE ) 
    end if

    call SetReduction ( LM, MyM % Value, M % Value )

    end associate !-- M, etc.

    !-- FIXME: NAG has trouble with pointer rank reassignment when the 
    !          left-hand side is a member
    call AssignPointers &
           ( LM, LM % MyM_RC, LM % MyM_IC, LM % MyM_RS, LM % MyM_IS, &
             LM % M_RC, LM % M_IC, LM % M_RS, LM % M_IS )

  end subroutine SetMomentStorage


  subroutine SetReduction ( LM, MyM_Value, M_Value )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      LM
    real ( KDR ), dimension ( :, : ), intent ( in ), target, contiguous :: &
      MyM_Value, &
        M_Value

    real ( KDR ), dimension ( : ), pointer :: &
      MyMoment_1D, &
        Moment_1D

    MyMoment_1D ( 1 : size ( MyM_Value ) )  =>  MyM_Value
      Moment_1D ( 1 : size (   M_Value ) )  =>    M_Value

    if ( allocated ( LM % ReductionMoments ) ) &
      deallocate ( LM % ReductionMoments )
    allocate ( LM % ReductionMoments )
    associate &
      (  RM  =>  LM % ReductionMoments, &
        PHC  =>  PROGRAM_HEADER % Communicator )
      call RM % Initialize &
             ( PHC, OutgoingValue = MyMoment_1D, &
               IncomingValue = Moment_1D )
    end associate !-- RM, etc.

    nullify ( Moment_1D, MyMoment_1D )

  end subroutine SetReduction


  subroutine AssignPointers &
               ( LM, MyM_RC, MyM_IC, MyM_RS, MyM_IS, M_RC, M_IC, M_RS, M_IS )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      LM
    real ( KDR ), dimension ( :, :, : ), pointer, intent ( out ) :: &
      MyM_RC, MyM_IC, &
      MyM_RS, MyM_IS, &
      M_RC, M_IC, &
      M_RS, M_IS

    associate &
      ( nE => LM % nEquations, &
        nR => LM % nRadialCells, &
        nA => LM % nAngularMomentCells )

    MyM_RC ( 1 : nA, 1 : nR, 1 : nE ) => LM % MyMoment_RC_1D
    MyM_IC ( 1 : nA, 1 : nR, 1 : nE ) => LM % MyMoment_IC_1D
      M_RC ( 1 : nA, 1 : nR, 1 : nE ) => LM % Moment_RC_1D
      M_IC ( 1 : nA, 1 : nR, 1 : nE ) => LM % Moment_IC_1D

    if ( LM % MaxOrder > 0 ) then
      MyM_RS ( 1 : nA, 1 : nR, 1 : nE ) => LM % MyMoment_RS_1D
      MyM_IS ( 1 : nA, 1 : nR, 1 : nE ) => LM % MyMoment_IS_1D
        M_RS ( 1 : nA, 1 : nR, 1 : nE ) => LM % Moment_RS_1D
        M_IS ( 1 : nA, 1 : nR, 1 : nE ) => LM % Moment_IS_1D
    end if

    end associate !-- nR, nA

  end subroutine AssignPointers


end module LaplacianMultipole_Template
