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
    logical ( KDL ) :: &
      UseDevice = .false., &
      ReductionUseDevice = .false.
    real ( KDR ), dimension ( :, :, : ), pointer :: &
        Moment_3D => null ( ), &
      MyMoment_3D => null ( )
    character ( LDF ) :: &
      Type = '', &
      Name = ''
    character ( LDL ), dimension ( : ), allocatable :: &
      AngularFunctionName, &
      MomentName
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
      FinalizeTemplate
    procedure, private, pass :: &
      SetParameters
    procedure ( SPA ), private, pass, deferred :: &
      SetParametersAtlas
    procedure ( SAF ), private, pass, deferred :: &
      SetAngularFunctions
    procedure, private, pass :: &
      AllocateMoments
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

    subroutine SAF ( L )
      import LaplacianMultipoleTemplate
      implicit none
      class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
        L
    end subroutine SAF

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

    L % IGNORABILITY  =  A % IGNORABILITY

    if ( L % Type == '' ) &
      L % Type = 'a LaplacianMultipole' 

    L % Name = 'Laplacian_' // trim ( A % Name )

    call Show ( 'Initializing ' // trim ( L % Type ), L % IGNORABILITY )
    call Show ( L % Name, 'Name', L % IGNORABILITY )

    call L % SetParameters ( A, MaxDegree, nEquations )
    call L % SetAngularFunctions ( )
    call L % AllocateMoments ( )

  end subroutine InitializeTemplate


  impure elemental subroutine FinalizeTemplate ( L )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      L

    if ( allocated ( L % ReductionMoments ) ) &
      deallocate ( L % ReductionMoments )

    if ( allocated ( L % RadialEdges ) ) &
      deallocate ( L % RadialEdges )
    if ( allocated ( L % MyMoments ) ) &
      deallocate ( L % MyMoments )
    if ( allocated ( L % Moments ) ) &
      deallocate ( L % Moments )

    nullify ( L % MyMoment_3D )
    nullify ( L % Moment_3D )

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
             ( [ nR, nA * nE ], &
               NameOption = 'Moments', &
               VariableOption = L % MomentName, &
               PinnedOption = L % UseDevice )
    call MyM % Initialize &
             ( [ nR, nA * nE ], &
               NameOption = 'MyMoments', &
               VariableOption = L % MomentName, &
               PinnedOption = L % UseDevice )
    if ( L % UseDevice ) then
      call   M % AllocateDevice ( )
      call MyM % AllocateDevice ( )
    end if

    call AllocateReduction ( L, M % Value, MyM % Value )

    call AssignMomentPointers &
           ( L, L % Moments % Value, L % MyMoments % Value, &
                L % Moment_3D,       L % MyMoment_3D )

    end associate !-- M, etc.

  end subroutine AllocateMoments


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
             P_MM  =  P_MP1_M
          P_MP1_M  =  P_LL
        end do !-- iL

        P_LM  =  P_LL
        return

      end if  !-- L == M + 1
    end if !-- L == M

  end function AssociatedLegendre


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
      if ( L % ReductionUseDevice ) &
        call RM % AllocateDevice ( )
    end associate !-- RM, etc.

    nullify ( Moment_1D, MyMoment_1D )

  end subroutine AllocateReduction


  subroutine AssignMomentPointers ( L, M_2D, MyM_2D, M_3D, MyM_3D )

    class ( LaplacianMultipoleTemplate ), intent ( in ) :: &
      L
    real ( KDR ), dimension ( :, : ), intent ( in ), target, contiguous :: &
        M_2D, &
      MyM_2D
    real ( KDR ), dimension ( :, :, : ), intent ( out ), pointer :: &
        M_3D, &
      MyM_3D

    associate &
      ( nE => L % nEquations, &
        nA => L % nAngularMoments, &
        nR => L % nRadialCells )

      M_3D ( 1 : nR, 1 : nA, 1 : nE )  =>    M_2D
    MyM_3D ( 1 : nR, 1 : nA, 1 : nE )  =>  MyM_2D

    end associate  !-- nE, nA, nR

  end subroutine AssignMomentPointers


end module LaplacianMultipole_Template
