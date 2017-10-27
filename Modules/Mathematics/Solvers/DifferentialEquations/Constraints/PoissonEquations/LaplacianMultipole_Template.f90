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
    real ( KDR ), dimension ( :, : ), pointer :: &
      MyMRC => null ( ), MyMRS => null ( ), &
      MyMIC => null ( ), MyMIS => null ( ), &
      MRC   => null ( ), MRS   => null ( ), &
      MIC   => null ( ), MIS   => null ( )
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
      NameSolidHarmonics
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

    subroutine SP ( LM, A, MaxDegree )
      use Basics
      use Manifolds
      import LaplacianMultipoleTemplate
      class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
        LM
      class ( AtlasHeaderForm ), intent ( in ), target :: &
        A
      integer ( KDI ), intent ( in ) :: &
        MaxDegree
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


  subroutine InitializeTemplate ( LM, A, MaxDegree )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      LM
    class ( AtlasHeaderForm ), intent ( in ) :: &
      A
    integer ( KDI ), intent ( in ) :: &
      MaxDegree

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

    call LM % SetParameters ( A, MaxDegree )

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


  subroutine NameSolidHarmonics ( LM, R_C_Name, I_C_Name, R_S_Name, I_S_Name )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      LM
    character ( * ), dimension ( : ), intent ( out ), allocatable :: &
      R_C_Name, I_C_Name, &
      R_S_Name, I_S_Name

    integer ( KDI ) :: &
      iV, &  !-- iValue
      iL, &
      iM
    character ( 3 ) :: &
      Label_L, &
      Label_M

    associate &
      ( L => LM % MaxDegree, &
        M => LM % MaxOrder )

    allocate ( R_C_Name ( LM % nAngularMomentCells ) ) 
    allocate ( I_C_Name ( LM % nAngularMomentCells ) ) 
    if ( LM % MaxOrder > 0 ) then
      allocate ( R_S_Name ( LM % nAngularMomentCells ) ) 
      allocate ( I_S_Name ( LM % nAngularMomentCells ) ) 
    end if

    iV = 0
    do iM = 0, M
      do iL = iM, L
        iV = iV + 1
        write ( Label_L, fmt = '(i3.3)' ) iL
        write ( Label_M, fmt = '(i3.3)' ) iM
        R_C_Name ( iV ) = 'RegularCos_' // Label_L // '_' // Label_M
        I_C_Name ( iV ) = 'IrregularCos_' // Label_L // '_' // Label_M
        if ( LM % MaxOrder > 0 ) then
          R_S_Name ( iV ) = 'RegularSin_' // Label_L // '_' // Label_M
          I_S_Name ( iV ) = 'IrregularSin_' // Label_L // '_' // Label_M
        end if
      end do
    end do

    end associate !-- L, etc.

  end subroutine NameSolidHarmonics


  subroutine ComputeSolidHarmonics &
               ( LM, CoordinateSystem, Position, nDimensions )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      LM
    character ( * ), intent ( in ) :: &
      CoordinateSystem
    real ( KDR ), dimension ( 3 ), intent ( in ) :: &
      Position
    integer ( KDI ), intent ( in ) :: &
      nDimensions

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

    call Clear ( LM % MyMRC )
    call Clear ( LM % MyMIC )
    if ( LM % MaxOrder > 0 ) then
      call Clear ( LM % MyMRS )
      call Clear ( LM % MyMIS )
    end if

    call LM % ComputeMomentsLocal ( Source )

    call LM % Reduction_RC % Reduce ( REDUCTION % SUM )
    call LM % Reduction_IC % Reduce ( REDUCTION % SUM )
    if ( LM % MaxOrder > 0 ) then
      call LM % Reduction_RS % Reduce ( REDUCTION % SUM )
      call LM % Reduction_IS % Reduce ( REDUCTION % SUM )
    end if

    call AddMomentShells &
           ( LM % MRC, LM % MIC, &
             LM % nAngularMomentCells, LM % nRadialCells )
    if ( LM % MaxOrder > 0 ) &
      call AddMomentShells &
             ( LM % MRS, LM % MIS, &
               LM % nAngularMomentCells, LM % nRadialCells )

  end subroutine ComputeMoments


  impure elemental subroutine FinalizeTemplate ( LM )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      LM

    if ( allocated ( LM % Reduction_RC ) ) deallocate ( LM % Reduction_RC )
    if ( allocated ( LM % Reduction_RS ) ) deallocate ( LM % Reduction_RS )
    if ( allocated ( LM % Reduction_IC ) ) deallocate ( LM % Reduction_IC )
    if ( allocated ( LM % Reduction_IS ) ) deallocate ( LM % Reduction_IS )

    nullify ( LM % MIS )
    nullify ( LM % MIC )
    nullify ( LM % MRS )
    nullify ( LM % MRC )
    nullify ( LM % MyMIS )
    nullify ( LM % MyMIC )
    nullify ( LM % MyMRS )
    nullify ( LM % MyMRC )

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

    associate ( nR_nA => LM % nRadialCells * LM % nAngularMomentCells )
    allocate ( LM % MyMoment_RC_1D ( nR_nA ) )
    allocate ( LM % MyMoment_IC_1D ( nR_nA ) )
    allocate ( LM % Moment_RC_1D ( nR_nA ) )
    allocate ( LM % Moment_IC_1D ( nR_nA ) )
    if ( LM % MaxOrder > 0 ) then
      allocate ( LM % MyMoment_RS_1D ( nR_nA ) )
      allocate ( LM % MyMoment_IS_1D ( nR_nA ) )
      allocate ( LM % Moment_RS_1D ( nR_nA ) )
      allocate ( LM % Moment_IS_1D ( nR_nA ) )
    end if
    end associate !-- nR_nA

    !-- FIXME: NAG has trouble with pointer rank reassignment when the 
    !          left-hand side is a member
    call AssignPointers &
           ( LM, LM % MyMRC, LM % MyMRS, LM % MyMIC, LM % MyMIS, &
             LM % MRC, LM % MRS, LM % MIC, LM % MIS )

    if ( allocated ( LM % Reduction_RC ) ) deallocate ( LM % Reduction_RC )
    if ( allocated ( LM % Reduction_RS ) ) deallocate ( LM % Reduction_RS )
    if ( allocated ( LM % Reduction_IC ) ) deallocate ( LM % Reduction_IC )
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
               ( LM, CoordinateSystem, Position, Width, VolumeJacobian, &
                 Source, nDimensions )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      LM
    character ( * ), intent ( in ) :: &
      CoordinateSystem
    real ( KDR ), dimension ( 3 ), intent ( in ) :: &
      Position, &
      Width
    real ( KDR ), intent ( in ) :: &
      VolumeJacobian, &
      Source
    integer ( KDI ), intent ( in ) :: &
      nDimensions

    integer ( KDI ) :: &
      iRS  !-- iRadiusSplit
    real ( KDR ) :: &
      X, Y, Z, &
      R_P, &  !-- R_Perpendicular
      D, &
      Source_dV

    select case ( trim ( CoordinateSystem ) )
    case ( 'CARTESIAN' )
      X  =  Position ( 1 )  -  LM % Origin ( 1 )
      Y  =  Position ( 2 )  -  LM % Origin ( 2 )
      Z  =  Position ( 3 )  -  LM % Origin ( 3 )
      D  =  sqrt ( X * X  +  Y * Y  +  Z * Z )
    case ( 'CYLINDRICAL' )
      R_P  =  Position ( 1 )
        Z  =  Position ( 2 )  -  LM % Origin ( 2 )
        D  =  sqrt ( R_P * R_P  +  Z * Z )
    case ( 'SPHERICAL' )
      D  =  Position ( 1 )
    end select !-- CoordinateSystem

    if ( D > LM % RadialEdge ( size ( LM % RadialEdge ) ) ) then
      call Show ( 'Radial grid not large enough', CONSOLE % ERROR )
      call Show ( D, 'D', CONSOLE % ERROR )
      call Show ( LM % RadialEdge ( size ( LM % RadialEdge ) ), 'D_Max', &
                  CONSOLE % ERROR )
      call Show ( 'ComputeMomentContributions', 'subroutine', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_Template', 'module', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    call Search ( LM % RadialEdge, D, iRS ) 
!call Show ( D, 'D', nLeadingLinesOption = 1 )
!call Show ( iRS, 'iRS' )
!call Show ( LM % RadialEdge ( iRS ), 'R_in' )
!call Show ( LM % RadialEdge ( iRS + 1 ), 'R_out' )

    Source_dV  =  Source * VolumeJacobian &
                         * product ( Width ( 1 : nDimensions ) )
!call Show ( Source_dV, 'Source_dV' )

    call ComputeMomentContributionsKernel &
           ( LM % MyMRC, LM % MyMIC, &
             LM % SolidHarmonic_RC, LM % SolidHarmonic_IC, &
             Source_dV, LM % nAngularMomentCells, iRS )
    if ( LM % MaxOrder > 0 ) &
      call ComputeMomentContributionsKernel &
             ( LM % MyMRS, LM % MyMIS, &
               LM % SolidHarmonic_RS, LM % SolidHarmonic_IS, &
               Source_dV, LM % nAngularMomentCells, iRS )

  end subroutine ComputeMomentContributions


  subroutine AssignPointers &
               ( LM, MyMRC, MyMRS, MyMIC, MyMIS, MRC, MRS, MIC, MIS )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      LM
    real ( KDR ), dimension ( :, : ), pointer, intent ( out ) :: &
      MyMRC, MyMRS, &
      MyMIC, MyMIS, &
      MRC, MRS, &
      MIC, MIS

    associate &
      ( nR => LM % nRadialCells, &
        nA => LM % nAngularMomentCells )

    MyMRC ( 1 : nA, 1 : nR ) => LM % MyMoment_RC_1D
    MyMIC ( 1 : nA, 1 : nR ) => LM % MyMoment_IC_1D
      MRC ( 1 : nA, 1 : nR ) => LM % Moment_RC_1D
      MIC ( 1 : nA, 1 : nR ) => LM % Moment_IC_1D

    if ( LM % MaxOrder > 0 ) then
      MyMRS ( 1 : nA, 1 : nR ) => LM % MyMoment_RS_1D
      MyMIS ( 1 : nA, 1 : nR ) => LM % MyMoment_IS_1D
        MRS ( 1 : nA, 1 : nR ) => LM % Moment_RS_1D
        MIS ( 1 : nA, 1 : nR ) => LM % Moment_IS_1D
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


  subroutine ComputeSolidHarmonicsKernel_C_S ( X, Y, Z, L, R_C, I_C, R_S, I_S )

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
          = ( ( 2 * iL - 1 ) * Z * R_C ( iV - 1 )  -  D_2 * R_C ( iV - 2 ) ) &
            / ( ( iL + iM ) * ( iL - iM ) )
        R_S ( iV ) &
          = ( ( 2 * iL - 1 ) * Z * R_S ( iV - 1 )  -  D_2 * R_S ( iV - 2 ) ) &
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


  subroutine AddMomentShells ( MR, MI, nA, nR )

    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      MR, MI
    integer ( KDI ), intent ( in ) :: &
      nA, nR

    integer ( KDI ) :: &
      iA, &  !-- angular index
      iR     !-- radial index

    !$OMP parallel do private ( iR, iA )
    do iR = 2, nR
      do iA = 1, nA
        MR ( iA, iR )  =  MR ( iA, iR - 1 )  +  MR ( iA, iR )
      end do
    end do
    !$OMP end parallel do

    !$OMP parallel do private ( iR, iA )
    do iR = nR - 1, 1, -1
      do iA = 1, nA
        MI ( iA, iR )  =  MI ( iA, iR + 1 )  +  MI ( iA, iR )
      end do
    end do
    !$OMP end parallel do

  end subroutine AddMomentShells


  subroutine ComputeMomentContributionsKernel &
               ( MyMR, MyMI, R, I, Source_dV, nA, iRS )

    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      MyMR, MyMI
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      R, I
    real ( KDR ), intent ( in ) :: &
      Source_dV
    integer ( KDI ), intent ( in ) :: &
      nA, &
      iRS

    integer ( KDI ) :: &
      iA      !-- iAngular
    real ( KDR ), dimension ( nA ) :: &
      R_S, I_S

    do iA = 1, nA 
      R_S ( iA )  =  R ( iA )  *  Source_dV
      I_S ( iA )  =  I ( iA )  *  Source_dV
    end do

    do iA = 1, nA
      MyMR ( iA, iRS )  =  MyMR ( iA, iRS )  +  R_S ( iA ) 
      MyMI ( iA, iRS )  =  MyMI ( iA, iRS )  +  I_S ( iA )
    end do

  end subroutine ComputeMomentContributionsKernel


end module LaplacianMultipole_Template
