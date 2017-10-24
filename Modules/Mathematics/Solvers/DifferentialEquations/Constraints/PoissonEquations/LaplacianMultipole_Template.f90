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
      SolidHarmonic_RC, SolidHarmonic_IC, &  !-- Regular, Irregular Cos
      SolidHarmonic_RS, SolidHarmonic_IS, &  !-- Regular, Irregular Sin
      Delta
    real ( KDR ), dimension ( : ), pointer :: &
      MyMoment_RC_1D => null ( ), &
      MyMoment_IC_1D => null ( ), &
      MyMoment_RS_1D => null ( ), &
      MyMoment_IS_1D => null ( ), &
      Moment_RC_1D   => null ( ), &    
      Moment_IC_1D   => null ( ), &
      Moment_RS_1D   => null ( ), &
      Moment_IS_1D   => null ( )
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
      FinalizeTemplate
    procedure, public, pass :: &
      SetMomentStorage
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

  end interface

    private :: &
      AssignPointers

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

    if ( allocated ( LM % Delta ) ) &
      deallocate ( LM % Delta )
    if ( allocated ( LM % SolidHarmonic_IS ) ) &
      deallocate ( LM % SolidHarmonic_IS )
    if ( allocated ( LM % SolidHarmonic_IC ) ) &
      deallocate ( LM % SolidHarmonic_IC )
    if ( allocated ( LM % SolidHarmonic_RS ) ) &
      deallocate ( LM % SolidHarmonic_RS )
    if ( allocated ( LM % SolidHarmonic_RC ) ) &
      deallocate ( LM % SolidHarmonic_RC )
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


end module LaplacianMultipole_Template
