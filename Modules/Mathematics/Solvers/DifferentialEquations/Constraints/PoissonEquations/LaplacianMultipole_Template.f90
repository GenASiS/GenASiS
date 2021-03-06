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
      iTimerComputeMoments = 0, &
      iTimerClearMoments = 0, &
      iTimerLocalMoments = 0, &
      iTimerReduceMoments = 0, &
      iTimerAddMoments = 0
    logical ( KDL ) :: &
      UseDevice = .false., &
      ReductionUseDevice = .false.
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      MyAngularMoment_3D => null ( ), &
        AngularMoment_3D => null ( ), &
       RadialMoment_R_3D => null ( ), &
       RadialMoment_I_3D => null ( )
    character ( LDF ) :: &
      Type = '', &
      Name = ''
    character ( LDL ), dimension ( : ), allocatable :: &
      AngularFunctionName, &
      MomentName
    type ( StorageForm ), allocatable :: &
      d_Radius_3_3, &
      CellFraction, &
      RadialFunctions_R, &
      RadialFunctions_I, &
      MyAngularMoments, &
        AngularMoments, &
      RadialMoments_R, &
      RadialMoments_I, &
      DeltaFactor
    type ( CollectiveOperation_R_Form ), allocatable :: &
      CO_AngularMoments
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
    procedure ( SKF ), private, pass, deferred :: &
      SetKernelFunctions
    procedure, private, pass :: &
      AllocateMoments
    procedure ( CAML ), private, pass, deferred :: &
      ComputeAngularMomentsLocal
    procedure, private, pass :: &
      ComputeRadialMoments
    procedure, public, pass :: &
      ShowMoments
    procedure, public, nopass :: &
      AssociatedLegendre
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

    subroutine SKF ( L )
      import LaplacianMultipoleTemplate
      implicit none
      class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
        L
    end subroutine SKF

    subroutine CAML ( L, Source )
      use Manifolds
      import LaplacianMultipoleTemplate
      implicit none
      class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
        L
      class ( FieldAtlasTemplate ), intent ( in ) :: &
        Source
    end subroutine CAML

  end interface


    private :: &
      AllocateReduction, &
      AssignMomentPointers

    private :: &
      ComputeRadialMomentsKernel


    interface

      module subroutine ComputeRadialMomentsKernel &
                          ( RM_R, RM_I, AM, RF_R, RF_I, dR33, nE, nAM, nR, &
                            UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
          RM_R, RM_I  !-- RadialMoment_Regular, _Irregular
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          AM  !-- AngularMoment
        real ( KDR ), dimension ( :, : ), intent ( in ) :: &
          RF_R, RF_I  !-- RadialFunction_Regular, _Irregular
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          dR33  !-- dRadius_3_3
        integer ( KDI ), intent ( in ) :: &
          nE, &   !-- nEquations
          nAM, &  !-- nAngularMoments
          nR      !-- nRadial
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine

    end interface


contains


  subroutine InitializeTemplate ( L, A, MaxDegree, nEquations )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      L
    class ( AtlasHeaderForm ), intent ( in ) :: &
      A
    integer ( KDI ), intent ( in ) :: &
      MaxDegree, &
      nEquations

    L % IGNORABILITY  =  A % IGNORABILITY

    if ( L % Type == '' ) &
      L % Type = 'a LaplacianMultipole' 

    L % Name = 'Laplacian_' // trim ( A % Name )

    call Show ( 'Initializing ' // trim ( L % Type ), L % IGNORABILITY )
    call Show ( L % Name, 'Name', L % IGNORABILITY )

    call L % SetParameters ( A, MaxDegree, nEquations )
    call L % SetKernelFunctions ( )
    call L % AllocateMoments ( )

  end subroutine InitializeTemplate


  subroutine InitializeTimers ( L, BaseLevel )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
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

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
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
      (   AM  =>  L %   AngularMoments, &
        MyAM  =>  L % MyAngularMoments )

    if ( associated ( Timer_CM ) ) call Timer_CM % Start ( )
    call MyAM % Clear ( )
    if ( associated ( Timer_CM ) ) call Timer_CM % Stop ( )

    if ( associated ( Timer_LM ) ) call Timer_LM % Start ( )
    call L % ComputeAngularMomentsLocal ( Source )
    if ( associated ( Timer_LM ) ) call Timer_LM % Stop ( )

    if ( associated ( Timer_RM ) ) call Timer_RM % Start ( )
    if ( .not. L % ReductionUseDevice ) call MyAM % UpdateHost ( ) 
    call L % CO_AngularMoments % Reduce ( REDUCTION % SUM )
    if ( .not. L % ReductionUseDevice ) call AM % UpdateDevice ( ) 
    if ( associated ( Timer_RM ) ) call Timer_RM % Stop ( )

    if ( associated ( Timer_AM ) ) call Timer_AM % Start ( )
      call L % ComputeRadialMoments ( )
    if ( associated ( Timer_AM ) ) call Timer_AM % Stop ( )
    
    end associate !-- M, etc.

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ComputeMoments


  impure elemental subroutine FinalizeTemplate ( L )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      L

    if ( allocated ( L % CO_AngularMoments ) ) &
      deallocate ( L % CO_AngularMoments )

    if ( allocated ( L % DeltaFactor ) ) &
      deallocate ( L % DeltaFactor )
    if ( allocated ( L % RadialMoments_I ) ) &
      deallocate ( L % RadialMoments_I )
    if ( allocated ( L % RadialMoments_R ) ) &
      deallocate ( L % RadialMoments_R )
    if ( allocated ( L % AngularMoments ) ) &
      deallocate ( L % AngularMoments )
    if ( allocated ( L % MyAngularMoments ) ) &
      deallocate ( L % MyAngularMoments )
    if ( allocated ( L % RadialFunctions_I ) ) &
      deallocate ( L % RadialFunctions_I )
    if ( allocated ( L % RadialFunctions_R ) ) &
      deallocate ( L % RadialFunctions_R )
    if ( allocated ( L % CellFraction ) ) &
      deallocate ( L % CellFraction )
    if ( allocated ( L % d_Radius_3_3 ) ) &
      deallocate ( L % d_Radius_3_3 )

    nullify ( L % RadialMoment_I_3D )
    nullify ( L % RadialMoment_R_3D )
    nullify ( L % AngularMoment_3D )
    nullify ( L % MyAngularMoment_3D )

    if ( allocated ( L % MomentName ) ) &
      deallocate ( L % MomentName )
    if ( allocated ( L % AngularFunctionName ) ) &
      deallocate ( L % AngularFunctionName )

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
      iL, &  !-- iDegree
      iM, &  !-- iOrder
      iE, &  !-- iEquation
      iA, &  !-- iAngular
      iAE    !-- iAngularEquations  
    character ( 1 ) :: &
      iE_Label
    character ( 2 ) :: &
      iL_Label, &
      iM_Label

    L % MaxDegree  =  MaxDegree
    L % MaxOrder   =  MaxDegree

    call L % SetParametersAtlas ( A )

    associate &
      (  L_Max  =>  L % MaxDegree, &
         M_Max  =>  L % MaxOrder, &
            nA  =>  L % nAngularMoments, &
            nE  =>  L % nEquations )

    nA = 0
    do iM  =  0, M_Max
      do iL  =  iM, L_Max
        nA  =  nA + 1
      end do
    end do
    nA  =  2 * nA  !-- Sine and Cosine parts 

    allocate ( L % AngularFunctionName ( nA ) )
    iA = 1
    do iM  =  0, M_Max
      do iL  =  iM, L_Max
        write ( iL_Label, fmt = '(i2.2)' ) iL
        write ( iM_Label, fmt = '(i2.2)' ) iM
        L % AngularFunctionName ( iA )  &
          =  'AngularFunction_' // iL_Label // '_' // iM_Label // '_Cos'
        L % AngularFunctionName ( iA + 1 )  &
          =  'AngularFunction_' // iL_Label // '_' // iM_Label // '_Sin'
        iA  =  iA + 2
      end do  !-- iL
    end do  !-- iM

    nE  =  nEquations

    allocate ( L % MomentName ( nA * nE ) )
    iAE = 1
    do iE  =  1, nE
      do iM  =  0, M_Max
        do iL  =  iM, L_Max
          write ( iL_Label, fmt = '(i2.2)' ) iL
          write ( iM_Label, fmt = '(i2.2)' ) iM
          write ( iE_Label, fmt = '(i1.1)' ) iE
          L % MomentName ( iAE )  &
            =  'Moment_' // iL_Label // '_' // iM_Label // '_Cos_Eq_' &
               // iE_Label
          L % MomentName ( iAE + 1 )  &
            =  'Moment_' // iL_Label // '_' // iM_Label // '_Sin_Eq_' &
               // iE_Label
          iAE  =  iAE + 2
        end do  !-- iL
      end do  !-- iM
    end do  !-- iE

    allocate ( L % DeltaFactor )
    associate ( DF => L % DeltaFactor )
    call DF % Initialize ( [ nA, 1 ] )
    iA = 1
    do iM  =  0, M_Max
      do iL  =  iM, L_Max
        if ( iM == 0 ) then
          DF % Value ( iA : iA + 1, 1 )  &
            =  1.0_KDR / ( 2.0_KDR * iL  +  1.0_KDR )
        else
          DF % Value ( iA : iA + 1, 1 )  &
            =  2.0_KDR / ( 2.0_KDR * iL  +  1.0_KDR )
        end if
        iA  =  iA + 2
      end do  !-- iL
    end do  !-- iM
    if ( L % UseDevice ) then
      call DF % AllocateDevice ( AssociateVariablesOption = .false. )
      call DF % UpdateDevice ( )
    end if
    end associate !-- DF

    end associate !-- L_Max, etc.

    call Show ( L % MaxDegree, 'MaxDegree (l)', L % IGNORABILITY )
    call Show ( L % MaxOrder, 'MaxOrder (m)', L % IGNORABILITY )
    call Show ( L % nRadialCells, 'nRadialCells', L % IGNORABILITY )
    call Show ( L % nAngularMoments, 'nAngularMoments', L % IGNORABILITY )
    call Show ( L % nEquations, 'nEquations', L % IGNORABILITY )
    call Show ( L % UseDevice, 'UseDevice', L % IGNORABILITY )
    call Show ( L % ReductionUseDevice, 'ReductionUseDevice', L % IGNORABILITY )

  end subroutine SetParameters


  subroutine AllocateMoments ( L )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      L

    if ( allocated ( L % AngularMoments ) ) &
      deallocate ( L % AngularMoments )
    if ( allocated ( L % MyAngularMoments ) ) &
      deallocate ( L % MyAngularMoments )

    allocate ( L % AngularMoments )
    allocate ( L % MyAngularMoments )
    allocate ( L % RadialMoments_R )
    allocate ( L % RadialMoments_I )
    associate &
      (         AM  =>  L % AngularMoments, &
              MyAM  =>  L % MyAngularMoments, &
              RM_R  =>  L % RadialMoments_R, &
              RM_I  =>  L % RadialMoments_I, &
               nAM  =>  L % nAngularMoments, &
                nR  =>  L % nRadialCells, &
                nE  =>  L % nEquations )

    call   AM % Initialize &
             ( [ nR, nAM * nE ], &
               NameOption = 'AngularMoments', &
               VariableOption = L % MomentName, &
               PinnedOption = L % UseDevice )
    call MyAM % Initialize &
             ( [ nR, nAM * nE ], &
               NameOption = 'MyAngularMoments', &
               VariableOption = L % MomentName, &
               PinnedOption = L % UseDevice )
    call RM_R % Initialize &
             ( [ nR + 1, nAM * nE ], &
               NameOption = 'RadialMoments_R', &
               VariableOption = L % MomentName, &
               PinnedOption = L % UseDevice )
    call RM_I % Initialize &
             ( [ nR + 1, nAM * nE ], &
               NameOption = 'RadialMoments_I', &
               VariableOption = L % MomentName, &
               PinnedOption = L % UseDevice )
    if ( L % UseDevice ) then
      call   AM % AllocateDevice ( AssociateVariablesOption = .false. )
      call MyAM % AllocateDevice ( AssociateVariablesOption = .false. )
      call RM_R % AllocateDevice ( AssociateVariablesOption = .false. )
      call RM_I % AllocateDevice ( AssociateVariablesOption = .false. )
    end if

    call Clear ( MyAM % Value, UseDeviceOption = L % UseDevice )

    call AllocateReduction ( L, AM % Value, MyAM % Value )

    call AssignMomentPointers &
           ( L, L %  AngularMoments % Value, L % MyAngularMoments % Value, &
                L % RadialMoments_R % Value, L %  RadialMoments_I % Value, &
                L % AngularMoment_3D,        L % MyAngularMoment_3D, &
                L % RadialMoment_R_3D,         L % RadialMoment_I_3D )

    end associate !-- M, etc.

  end subroutine AllocateMoments


  subroutine ComputeRadialMoments ( L )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      L

    call ComputeRadialMomentsKernel &
           ( L % RadialMoment_R_3D, L % RadialMoment_I_3D, &
             L % AngularMoment_3D, L % RadialFunctions_R % Value, &
             L % RadialFunctions_I % Value, L % d_Radius_3_3 % Value ( :, 1 ), &
             L % nEquations, L % nAngularMoments, L % nRadialCells, &
             UseDeviceOption = L % UseDevice )
    
  end subroutine ComputeRadialMoments


  subroutine ShowMoments ( L )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      L

    integer ( KDI ) :: &
      iL, &  !-- iDegree
      iM, &  !-- iOrder
      iE, &  !-- iEquation
      iA, &  !-- iAngular
      iAE    !-- iAngularEquations  

    associate &
      (  L_Max  =>  L % MaxDegree, &
         M_Max  =>  L % MaxOrder, &
            nE  =>  L % nEquations )

    call Show ( 'Displaying Angular Moments' )

    iAE = 1
    do iE  =  1, nE
      iA = 1
      do iM  =  0, M_Max
        do iL  =  iM, L_Max
          !-- Cos and Sin
          call Show ( L % AngularMoment_3D ( :, iA, iE ), &
                      L % MomentName ( iAE ) )
          call Show ( L % AngularMoment_3D ( :, iA + 1, iE ), &
                      L % MomentName ( iAE + 1) )
           iA  =   iA + 2
          iAE  =  iAE + 2
        end do  !-- iL
      end do  !-- iM
    end do  !-- iE

    call Show ( 'Displaying Regular Radial Moments' )

    iAE = 1
    do iE  =  1, nE
      iA = 1
      do iM  =  0, M_Max
        do iL  =  iM, L_Max
          !-- Cos and Sin
          call Show ( L % RadialMoment_R_3D ( :, iA, iE ), &
                      L % MomentName ( iAE ) )
          call Show ( L % RadialMoment_R_3D ( :, iA + 1, iE ), &
                      L % MomentName ( iAE + 1) )
           iA  =   iA + 2
          iAE  =  iAE + 2
        end do  !-- iL
      end do  !-- iM
    end do  !-- iE

    call Show ( 'Displaying Irregular Moments' )

    iAE = 1
    do iE  =  1, nE
      iA = 1
      do iM  =  0, M_Max
        do iL  =  iM, L_Max
          !-- Cos and Sin
          call Show ( L % RadialMoment_I_3D ( :, iA, iE ), &
                      L % MomentName ( iAE ) )
          call Show ( L % RadialMoment_I_3D ( :, iA + 1, iE ), &
                      L % MomentName ( iAE + 1) )
           iA  =   iA + 2
          iAE  =  iAE + 2
        end do  !-- iL
      end do  !-- iM
    end do  !-- iE

    end associate !-- L_Max, etc.

  end subroutine ShowMoments


  function AssociatedLegendre ( X, L, M ) result ( P_LM )
  
    !-- Normalized, see Numerical Recipes Third Edition, Section 6.7

    real ( KDR ), intent ( in ) :: &
      X
    integer ( KDI ), intent ( in ) :: &
      L, &
      M
    real ( KDR ) :: &
      P_LM

    integer ( KDI ) :: &
      iM, &
      iL
    real ( KDR ) :: &
      P_MM, &
      P_MP1_M, &  !-- P_M+1_M
      P_LL, &
      O_M_X2, &  !-- 1 - x ** 2
      Factor, &
      FactorOld, &
      FourPi

    if ( M < 0 .or. M > L .or. abs ( X ) > 1.0_KDR ) then
      call Show ( 'Arguments out of range', CONSOLE % ERROR )
      call Show ( M, 'M', CONSOLE % ERROR )
      call Show ( L, 'L', CONSOLE % ERROR )
      call Show ( X, 'X', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_Template', 'module', CONSOLE % ERROR )
      call Show ( 'AssociatedLegendre', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    FourPi  =  4.0_KDR  *  CONSTANT % PI

    !-- Compute P_mm

    P_MM  =  1.0_KDR

    if ( M > 0 ) then
      O_M_X2  =  ( 1.0_KDR - X ) * ( 1.0_KDR + X )
      Factor  =  1.0_KDR
      do iM  =  1, M
          P_MM  =  P_MM  *  Factor / ( Factor + 1.0_KDR )  *  O_M_X2  
        Factor  =  Factor + 2.0_KDR
      end do !-- iM
    end if !-- iM > 0

    P_MM  =  sqrt ( ( 2 * M + 1 ) * P_MM / FourPi )
    if ( mod ( M, 2 ) == 1 )  &
      P_MM  =  - P_MM

    if ( L == M ) then

      P_LM  =  P_MM
      return

    else

      !-- Compute P_lm

      P_MP1_M  =  sqrt ( 2.0_KDR * M  +  3.0_KDR )  *  X  *  P_MM

      if ( L == M + 1 ) then
        P_LM  =  P_MP1_M
        return

      else 

        !-- Compute P_lm, l > m + 1
        FactorOld  =  sqrt ( 2.0_KDR * M  +  3.0_KDR )
        do iL  =  M + 2, L
              Factor  =  sqrt ( ( 4.0_KDR * iL * iL  - 1.0_KDR )  &
                                /  ( iL * iL  -  M * M ) )
                P_LL  =  Factor * ( X * P_MP1_M  -  P_MM / FactorOld )
           FactorOld  =  Factor
                P_MM  =  P_MP1_M
             P_MP1_M  =  P_LL
        end do !-- iL

        P_LM  =  P_LL
        return

      end if  !-- L == M + 1
    end if !-- L == M

  end function AssociatedLegendre


  subroutine AllocateReduction ( L, AM_Value, MyAM_Value )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      L
    real ( KDR ), dimension ( :, : ), intent ( in ), target, contiguous :: &
        AM_Value, &
      MyAM_Value

    real ( KDR ), dimension ( : ), pointer :: &
        AngularMoment_1D, &
      MyAngularMoment_1D

      AngularMoment_1D ( 1 : size (   AM_Value ) )  =>    AM_Value
    MyAngularMoment_1D ( 1 : size ( MyAM_Value ) )  =>  MyAM_Value

    if ( allocated ( L % CO_AngularMoments ) ) &
      deallocate ( L % CO_AngularMoments )
    allocate ( L % CO_AngularMoments )
    associate &
      (  CO  =>  L % CO_AngularMoments, &
        PHC  =>  PROGRAM_HEADER % Communicator )
      call CO % Initialize &
             ( PHC, OutgoingValue = MyAngularMoment_1D, &
               IncomingValue = AngularMoment_1D )
      if ( L % ReductionUseDevice ) &
        call CO % AllocateDevice ( )
    end associate !-- RM, etc.

    nullify ( AngularMoment_1D, MyAngularMoment_1D )

  end subroutine AllocateReduction


  subroutine AssignMomentPointers &
               ( L, AM_2D, MyAM_2D, RM_R_2D, RM_I_2D, &
                    AM_3D, MyAM_3D, RM_R_3D, RM_I_3D )

    class ( LaplacianMultipoleTemplate ), intent ( in ) :: &
      L
    real ( KDR ), dimension ( :, : ), intent ( in ), target, contiguous :: &
        AM_2D, &
      MyAM_2D, &
      RM_R_2D, &
      RM_I_2D
    real ( KDR ), dimension ( :, :, : ), intent ( out ), pointer :: &
        AM_3D, &
      MyAM_3D, &
      RM_R_3D, &
      RM_I_3D

    associate &
      (  nE => L % nEquations, &
        nAM => L % nAngularMoments, &
         nR => L % nRadialCells )

      AM_3D ( 1 : nR,     1 : nAM, 1 : nE )  =>    AM_2D
    MyAM_3D ( 1 : nR,     1 : nAM, 1 : nE )  =>  MyAM_2D
    RM_R_3D ( 1 : nR + 1, 1 : nAM, 1 : nE )  =>  RM_R_2D
    RM_I_3D ( 1 : nR + 1, 1 : nAM, 1 : nE )  =>  RM_I_2D

    end associate  !-- nE, nA, nR

  end subroutine AssignMomentPointers


end module LaplacianMultipole_Template
