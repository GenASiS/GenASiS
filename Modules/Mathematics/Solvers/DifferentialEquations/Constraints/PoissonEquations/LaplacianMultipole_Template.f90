module LaplacianMultipole_Template

  use Basics
  use Manifolds

  implicit none
  private

  type, public, abstract :: LaplacianMultipoleTemplate
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      nRadialCells = 0, &
      nAngularMomentCells = 0, &
      nEquations = 0, &
      MaxDegree = 0, &  !-- Max L
      MaxOrder  = 0     !-- Max M
    real ( KDR ), dimension ( 3 ) :: &
      Origin = 0.0_KDR
    real ( KDR ), dimension ( : ), allocatable :: &
      RadialEdge, &
      Delta
    real ( KDR ), dimension ( : ), pointer :: &
      SolidHarmonic_RC => null ( ), &  !-- Regular Cos
      SolidHarmonic_IC => null ( ), &  !-- Irregular Cos
      SolidHarmonic_RS => null ( ), &  !-- Regular Sin
      SolidHarmonic_IS => null ( ), &  !-- Irregular Sin
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
    type ( CollectiveOperation_R_Form ), allocatable :: &
      Reduction_RC, Reduction_IC, &
      Reduction_RS, Reduction_IS
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure ( SP ), private, pass, deferred :: &
      SetParameters   
    procedure, public, pass :: &
      ComputeSolidHarmonics
    procedure, public, pass :: &
      ComputeMoments
    procedure, public, pass :: &
      FinalizeTemplate
    procedure, public, pass :: &
      SetMomentStorage
    procedure ( CML ), private, pass, deferred :: &
      ComputeMomentsLocal
    procedure, public, pass :: &
      ComputeMomentContributions
  end type LaplacianMultipoleTemplate

  abstract interface

    subroutine SP ( LM, A, MaxDegree, nEquations )
      use Basics
      use Manifolds
      import LaplacianMultipoleTemplate
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
      class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
        LM
      type ( VariableGroupForm ), intent ( in ) :: &
        Source  
    end subroutine CML
  end interface

    private :: &
      AssignPointers, &
      ComputeSolidHarmonicsKernel_C_M_0, &
      ComputeSolidHarmonicsKernel_C_S, &
      AddMomentShells, &
      ComputeMomentContributionsKernel

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

    allocate ( LM % SolidHarmonic_RC ( LM % nAngularMomentCells ) )
    allocate ( LM % SolidHarmonic_IC ( LM % nAngularMomentCells ) )

    if ( LM % MaxOrder > 0 ) then
      allocate ( LM % SolidHarmonic_RS ( LM % nAngularMomentCells ) )
      allocate ( LM % SolidHarmonic_IS ( LM % nAngularMomentCells ) )
    end if

    allocate ( LM % Delta ( LM % nAngularMomentCells ) )
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


  subroutine ComputeSolidHarmonics &
               ( LM, CoordinateSystem, Position, nDimensions, R, iR )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      LM
    character ( * ), intent ( in ) :: &
      CoordinateSystem
    real ( KDR ), dimension ( 3 ), intent ( in ) :: &
      Position
    integer ( KDI ), intent ( in ) :: &
      nDimensions
    real ( KDR ), intent ( out ) :: &
      R
    integer ( KDI ), intent ( out ) :: &
      iR

    integer ( KDI ) :: &
      L
    real ( KDR ) :: &
      X, Y, Z
    real ( KDR ), dimension ( : ), pointer :: &
      R_C, I_C, &
      R_S, I_S

    L  =  LM % MaxDegree

    R_C  =>  LM % SolidHarmonic_RC ( : )
    I_C  =>  LM % SolidHarmonic_IC ( : )
    R_S  =>  LM % SolidHarmonic_RS ( : )
    I_S  =>  LM % SolidHarmonic_IS ( : )

    select case ( trim ( CoordinateSystem ) )
    case ( 'CARTESIAN' )
      X  =  Position ( 1 )  -  LM % Origin ( 1 )
      Y  =  Position ( 2 )  -  LM % Origin ( 2 )
      Z  =  Position ( 3 )  -  LM % Origin ( 3 )
    case ( 'CYLINDRICAL' )
      if ( nDimensions < 3 ) then
        X  =  Position ( 1 )
        Y  =  0.0_KDR
      else
        X  =  Position ( 1 )  *  cos ( Position ( 3 ) )
        Y  =  Position ( 1 )  *  sin ( Position ( 3 ) )
      end if
      Z  =  Position ( 2 )  -  LM % Origin ( 2 )
    case ( 'SPHERICAL' )
      if ( nDimensions < 3 ) then
        X  =  Position ( 1 )  *  sin ( Position ( 2 ) )
        Y  =  0.0_KDR
      else
        X  =  Position ( 1 )  *  sin ( Position ( 2 ) )  &
                              *  cos ( Position ( 3 ) )
        Y  =  Position ( 1 )  *  sin ( Position ( 2 ) )  &
                              *  sin ( Position ( 3 ) )
      end if
      Z  =  Position ( 1 )  *  cos ( Position ( 2 ) )
    end select !-- CoordinateSystem

    if ( nDimensions < 3 ) then
      call ComputeSolidHarmonicsKernel_C_M_0 ( X, Z, L, R_C, I_C )
    else
      call ComputeSolidHarmonicsKernel_C_S ( X, Y, Z, L, R_C, I_C, R_S, I_S )
    end if

    R  =  sqrt ( X * X  +  Y * Y  +  Z * Z )

    if ( R > LM % RadialEdge ( size ( LM % RadialEdge ) ) ) then
      call Show ( 'Radial grid not large enough', CONSOLE % ERROR )
      call Show ( R, 'R', CONSOLE % ERROR )
      call Show ( LM % RadialEdge ( size ( LM % RadialEdge ) ), 'R_Max', &
                  CONSOLE % ERROR )
      call Show ( 'ComputeMomentContributions', 'subroutine', &
                  CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_Template', 'module', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    call Search ( LM % RadialEdge, R, iR ) 
!call Show ( R, 'R', nLeadingLinesOption = 1 )
!call Show ( iR, 'iR' )
!call Show ( LM % RadialEdge ( iR ), 'R_in' )
!call Show ( LM % RadialEdge ( iR + 1 ), 'R_out' )

    nullify ( R_C, I_C, R_S, I_S )

  end subroutine ComputeSolidHarmonics


  subroutine ComputeMoments ( LM, Source )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      LM
    type ( VariableGroupForm ), intent ( in ) :: &
      Source !-- array over levels    
    
    ! integer ( KDI ) :: &
    !   iL, &  !-- iLevel
    !   iC     !-- iCell

    call Show ( 'Computing Moments', CONSOLE % INFO_3 )

    call Clear ( LM % MyM_RC )
    call Clear ( LM % MyM_IC )
    if ( LM % MaxOrder > 0 ) then
      call Clear ( LM % MyM_RS )
      call Clear ( LM % MyM_IS )
    end if

    call LM % ComputeMomentsLocal ( Source )

    call LM % Reduction_RC % Reduce ( REDUCTION % SUM )
    call LM % Reduction_IC % Reduce ( REDUCTION % SUM )
    if ( LM % MaxOrder > 0 ) then
      call LM % Reduction_RS % Reduce ( REDUCTION % SUM )
      call LM % Reduction_IS % Reduce ( REDUCTION % SUM )
    end if

    call AddMomentShells &
           ( LM % M_RC, LM % M_IC, &
             LM % nEquations, LM % nRadialCells, LM % nAngularMomentCells )
    if ( LM % MaxOrder > 0 ) &
      call AddMomentShells &
             ( LM % M_RS, LM % M_IS, &
               LM % nEquations, LM % nRadialCells, LM % nAngularMomentCells )

  end subroutine ComputeMoments


  impure elemental subroutine FinalizeTemplate ( LM )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      LM

    if ( allocated ( LM % Reduction_RC ) ) deallocate ( LM % Reduction_RC )
    if ( allocated ( LM % Reduction_RS ) ) deallocate ( LM % Reduction_RS )
    if ( allocated ( LM % Reduction_IC ) ) deallocate ( LM % Reduction_IC )
    if ( allocated ( LM % Reduction_IS ) ) deallocate ( LM % Reduction_IS )

    nullify ( LM % M_IS )
    nullify ( LM % M_RS )
    nullify ( LM % M_IC )
    nullify ( LM % M_RC )
    nullify ( LM % MyM_IS )
    nullify ( LM % MyM_RS )
    nullify ( LM % MyM_IC )
    nullify ( LM % MyM_RC )

    if ( associated ( LM % Moment_IS_1D ) ) &
       deallocate ( LM % Moment_IS_1D )
    if ( associated ( LM % Moment_IC_1D ) ) &
       deallocate ( LM % Moment_IC_1D )
    if ( associated ( LM % Moment_RS_1D ) ) &
       deallocate ( LM % Moment_RS_1D )
    if ( associated ( LM % Moment_RC_1D ) ) &
       deallocate ( LM % Moment_RC_1D )
    if ( associated ( LM % MyMoment_IS_1D ) ) &
       deallocate ( LM % MyMoment_IS_1D )
    if ( associated ( LM % MyMoment_IC_1D ) ) &
       deallocate ( LM % MyMoment_IC_1D )
    if ( associated ( LM % MyMoment_RS_1D ) ) &
       deallocate ( LM % MyMoment_RS_1D )
    if ( associated ( LM % MyMoment_RC_1D ) ) &
       deallocate ( LM % MyMoment_RC_1D )
    if ( associated ( LM % SolidHarmonic_IS ) ) &
      deallocate ( LM % SolidHarmonic_IS )
    if ( associated ( LM % SolidHarmonic_RS ) ) &
      deallocate ( LM % SolidHarmonic_RS )
    if ( associated ( LM % SolidHarmonic_IC ) ) &
      deallocate ( LM % SolidHarmonic_IC )
    if ( associated ( LM % SolidHarmonic_RC ) ) &
      deallocate ( LM % SolidHarmonic_RC )

    if ( allocated ( LM % Delta ) ) &
      deallocate ( LM % Delta )
    if ( allocated ( LM % RadialEdge ) ) &
      deallocate ( LM % RadialEdge )

    if ( LM % Name == '' ) return

    call Show ( 'Finalizing ' // trim ( LM % Type ), LM % IGNORABILITY )
    call Show ( LM % Name, 'Name', LM % IGNORABILITY )
    
  end subroutine FinalizeTemplate


  subroutine SetMomentStorage ( LM )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      LM

    if ( associated ( LM % Moment_IS_1D ) ) &
       deallocate ( LM % Moment_IS_1D )
    if ( associated ( LM % Moment_IC_1D ) ) &
       deallocate ( LM % Moment_IC_1D )
    if ( associated ( LM % Moment_RS_1D ) ) &
       deallocate ( LM % Moment_RS_1D )
    if ( associated ( LM % Moment_RC_1D ) ) &
       deallocate ( LM % Moment_RC_1D )
    if ( associated ( LM % MyMoment_IS_1D ) ) &
       deallocate ( LM % MyMoment_IS_1D )
    if ( associated ( LM % MyMoment_IC_1D ) ) &
       deallocate ( LM % MyMoment_IC_1D )
    if ( associated ( LM % MyMoment_RS_1D ) ) &
       deallocate ( LM % MyMoment_RS_1D )
    if ( associated ( LM % MyMoment_RC_1D ) ) &
       deallocate ( LM % MyMoment_RC_1D )

    associate &
      ( nR_nA_nE => LM % nRadialCells * LM % nAngularMomentCells &
                    * LM % nEquations )
    allocate ( LM % MyMoment_RC_1D ( nR_nA_nE ) )
    allocate ( LM % MyMoment_IC_1D ( nR_nA_nE ) )
    allocate ( LM % Moment_RC_1D ( nR_nA_nE ) )
    allocate ( LM % Moment_IC_1D ( nR_nA_nE ) )
    if ( LM % MaxOrder > 0 ) then
      allocate ( LM % MyMoment_RS_1D ( nR_nA_nE ) )
      allocate ( LM % MyMoment_IS_1D ( nR_nA_nE ) )
      allocate ( LM % Moment_RS_1D ( nR_nA_nE ) )
      allocate ( LM % Moment_IS_1D ( nR_nA_nE ) )
    end if
    end associate !-- nR_nA_nE

    !-- FIXME: NAG has trouble with pointer rank reassignment when the 
    !          left-hand side is a member
    call AssignPointers &
           ( LM, LM % MyM_RC, LM % MyM_IC, LM % MyM_RS, LM % MyM_IS, &
             LM % M_RC, LM % M_IC, LM % M_RS, LM % M_IS )

    if ( allocated ( LM % Reduction_RC ) ) deallocate ( LM % Reduction_RC )
    if ( allocated ( LM % Reduction_IC ) ) deallocate ( LM % Reduction_IC )
    if ( allocated ( LM % Reduction_RS ) ) deallocate ( LM % Reduction_RS )
    if ( allocated ( LM % Reduction_IS ) ) deallocate ( LM % Reduction_IS )

    associate ( PHC => PROGRAM_HEADER % Communicator )
    allocate ( LM % Reduction_RC )
    allocate ( LM % Reduction_IC )
    call LM % Reduction_RC % Initialize &
           ( PHC, OutgoingValue = LM % MyMoment_RC_1D, &
             IncomingValue = LM % Moment_RC_1D )
    call LM % Reduction_IC % Initialize &
           ( PHC, OutgoingValue = LM % MyMoment_IC_1D, &
             IncomingValue = LM % Moment_IC_1D )
    if ( LM % MaxOrder > 0 ) then
      allocate ( LM % Reduction_RS )
      allocate ( LM % Reduction_IS )
      call LM % Reduction_RS % Initialize &
             ( PHC, OutgoingValue = LM % MyMoment_RS_1D, &
               IncomingValue = LM % Moment_RS_1D )
      call LM % Reduction_IS % Initialize &
             ( PHC, OutgoingValue = LM % MyMoment_IS_1D, &
               IncomingValue = LM % Moment_IS_1D )
    end if
    end associate !-- PHC

  end subroutine SetMomentStorage


  subroutine ComputeMomentContributions &
               ( LM, Source, Width, VolumeJacobian, iaSource, nDimensions, &
                 iR )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      LM
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Source
    real ( KDR ), dimension ( 3 ), intent ( in ) :: &
      Width
    real ( KDR ), intent ( in ) :: &
      VolumeJacobian
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iaSource
    integer ( KDI ), intent ( in ) :: &
      nDimensions
    integer ( KDI ) :: &
      iR  !-- iRadius

    integer ( KDI ) :: &
      iE  !-- iEquation
    real ( KDR ) :: &
      dV
    real ( KDR ), dimension ( LM % nEquations ) :: &
      Source_dV

    dV  =  VolumeJacobian * product ( Width ( 1 : nDimensions ) )
    Source_dV  =  [ ( Source ( iaSource ( iE ) )  *  dV, &
                      iE = 1, LM % nEquations ) ] 
!call Show ( Source_dV, 'Source_dV' )

    call ComputeMomentContributionsKernel &
           ( LM % MyM_RC, LM % MyM_IC, &
             LM % SolidHarmonic_RC, LM % SolidHarmonic_IC, &
             Source_dV, LM % nEquations, LM % nAngularMomentCells, iR )
    if ( LM % MaxOrder > 0 ) &
      call ComputeMomentContributionsKernel &
             ( LM % MyM_RS, LM % MyM_IS, &
               LM % SolidHarmonic_RS, LM % SolidHarmonic_IS, &
               Source_dV, LM % nEquations, LM % nAngularMomentCells, iR )

  end subroutine ComputeMomentContributions


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


  subroutine ComputeSolidHarmonicsKernel_C_M_0 ( X, Z, L, R_C, I_C )

    real ( KDR ), intent ( in ) :: &
      X, Z
    integer ( KDI ), intent ( in ) :: &
      L
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      R_C, I_C

    integer ( KDI ) :: &
      iV, &  !-- iValue
      iL, &
      iM, &
      iPD   !-- iPreviousDiagonal
    real ( KDR ) :: &
      D_2

    D_2  =  X * X  +  Z * Z 

    iV = 0
    iM = 0

    !-- ( L, M ) = ( iM, iM )
    !-- Note iL = iM
    iV = iV + 1
    if ( iM == 0 ) then
      iV = 1
      R_C ( iV ) = 1.0_KDR
      I_C ( iV ) = 1.0_KDR / sqrt ( D_2 )
    else
      R_C ( iV ) = - ( X * R_C ( iPD ) ) / ( 2 * iM )
      I_C ( iV ) = - ( 2 * iM - 1 ) &
                     * ( X * I_C ( iPD ) ) / D_2
    end if
    iPD = iV

    if ( iM == L ) return

    !-- ( L, M ) = ( iM + 1, iM )
    !-- Note iL = iM + 1
    iV = iV + 1
    R_C ( iV ) = Z * R_C ( iV - 1 )  
    I_C ( iV ) = ( 2 * ( iM + 1 ) - 1 ) * Z * I_C ( iV - 1 ) / D_2  

    do iL = iM + 2, L
      !-- ( L, M ) = ( iL, iM )
      iV = iV + 1
      R_C ( iV ) &
        = ( ( 2 * iL - 1 ) * Z * R_C ( iV - 1 )  -  D_2 * R_C ( iV - 2 ) ) &
          / ( ( iL + iM ) * ( iL - iM ) )
      I_C ( iV ) &
        = ( ( 2 * iL - 1 ) * Z * I_C ( iV - 1 )  &
            -  ( ( iL - 1 )**2 - iM**2 ) * I_C ( iV - 2 ) ) &
          / D_2
    end do !-- iL

  end subroutine ComputeSolidHarmonicsKernel_C_M_0


  subroutine ComputeSolidHarmonicsKernel_C_S &
               ( X, Y, Z, L, R_C, I_C, R_S, I_S )

    real ( KDR ), intent ( in ) :: &
      X, Y, Z
    integer ( KDI ), intent ( in ) :: &
      L
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      R_C, I_C, &
      R_S, I_S

    integer ( KDI ) :: &
      iV, &  !-- iValue
      iL, &
      iM, &
      iPD   !-- iPreviousDiagonal
    real ( KDR ) :: &
      D_2

    D_2  =  X * X  +  Y * Y  +  Z * Z 

    iV = 0
    do iM = 0, L

      !-- ( L, M ) = ( iM, iM )
      !-- Note iL = iM
      iV = iV + 1
      if ( iM == 0 ) then
        iV = 1
        R_C ( iV ) = 1.0_KDR
        R_S ( iV ) = 0.0_KDR
        I_C ( iV ) = 1.0_KDR / sqrt ( D_2 )
        I_S ( iV ) = 0.0_KDR
      else
        R_C ( iV ) = - ( X * R_C ( iPD ) - Y * R_S ( iPD ) ) / ( 2 * iM )
        R_S ( iV ) = - ( Y * R_C ( iPD ) + X * R_S ( iPD ) ) / ( 2 * iM )
        I_C ( iV ) = - ( 2 * iM - 1 ) &
                       * ( X * I_C ( iPD ) - Y * I_S ( iPD ) ) / D_2
        I_S ( iV ) = - ( 2 * iM - 1 ) &
                       * ( Y * I_C ( iPD ) + X * I_S ( iPD ) ) / D_2
      end if
      iPD = iV

      if ( iM == L ) exit

      !-- ( L, M ) = ( iM + 1, iM )
      !-- Note iL = iM + 1
      iV = iV + 1
      R_C ( iV ) = Z * R_C ( iV - 1 )  
      R_S ( iV ) = Z * R_S ( iV - 1 )  
      I_C ( iV ) = ( 2 * ( iM + 1 ) - 1 ) * Z * I_C ( iV - 1 ) / D_2  
      I_S ( iV ) = ( 2 * ( iM + 1 ) - 1 ) * Z * I_S ( iV - 1 ) / D_2 

      do iL = iM + 2, L
        !-- (L,M) = (iL,iM)
        iV = iV + 1
        R_C ( iV ) &
          = ( ( 2 * iL - 1 ) * Z * R_C ( iV - 1 )  -  D_2 * R_C ( iV - 2 ) )&
            / ( ( iL + iM ) * ( iL - iM ) )
        R_S ( iV ) &
          = ( ( 2 * iL - 1 ) * Z * R_S ( iV - 1 )  -  D_2 * R_S ( iV - 2 ) )&
            / ( ( iL + iM ) * ( iL - iM ) )
        I_C ( iV ) &
          = ( ( 2 * iL - 1 ) * Z * I_C ( iV - 1 )  &
              -  ( ( iL - 1 )**2 - iM**2 ) * I_C ( iV - 2 ) ) &
            / D_2
        I_S ( iV ) &
          = ( ( 2 * iL - 1 ) * Z * I_S ( iV - 1 )  &
              -  ( ( iL - 1 )**2 - iM**2 ) * I_S ( iV - 2 ) ) &
            / D_2
      end do !-- iL

    end do !-- iM

  end subroutine ComputeSolidHarmonicsKernel_C_S


  subroutine AddMomentShells ( MR, MI, nE, nR, nA )

    real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
      MR, MI
    integer ( KDI ), intent ( in ) :: &
      nE, nR, nA

    integer ( KDI ) :: &
      iA, &  !-- angular index
      iR, &  !-- radial index
      iE     !-- equation index

    do iE = 1, nE
      !$OMP parallel do private ( iR, iA )
      do iR = 2, nR
        do iA = 1, nA
          MR ( iA, iR, iE )  =  MR ( iA, iR - 1, iE )  +  MR ( iA, iR, iE )
        end do
      end do
      !$OMP end parallel do
    end do

    do iE = 1, nE
      !$OMP parallel do private ( iR, iA )
      do iR = nR - 1, 1, -1
        do iA = 1, nA
          MI ( iA, iR, iE )  =  MI ( iA, iR + 1, iE )  +  MI ( iA, iR, iE )
        end do
      end do
      !$OMP end parallel do
    end do

  end subroutine AddMomentShells


  subroutine ComputeMomentContributionsKernel &
               ( MyMR, MyMI, SH_R, SH_I, Source_dV, nE, nA, iRS )

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

    integer ( KDI ) :: &
      iA, &  !-- iAngular
      iE     !-- iEquation   

    do iE = 1, nE
      do iA = 1, nA
        MyMR ( iA, iRS, iE )  &
          =  MyMR ( iA, iRS, iE )  +  SH_R ( iA )  *  Source_dV ( iE ) 
        MyMI ( iA, iRS, iE )  &
          =  MyMI ( iA, iRS, iE )  +  SH_I ( iA )  *  Source_dV ( iE )
      end do
    end do

  end subroutine ComputeMomentContributionsKernel


end module LaplacianMultipole_Template
