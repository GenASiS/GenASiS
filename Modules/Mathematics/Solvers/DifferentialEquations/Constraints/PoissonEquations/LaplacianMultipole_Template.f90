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
      SolidHarmonic_RC, SolidHarmonic_RS, &  !-- Regular Cos, Sin
      SolidHarmonic_IC, SolidHarmonic_IS, &  !-- Irregular Cos, Sin
      Delta
    character ( LDF ) :: &
      Type = '', &
      Name = ''
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure ( SP ), private, pass, deferred :: &
      SetParameters   
    procedure, public, pass :: &
      FinalizeTemplate
    procedure, private, pass :: &
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
    allocate ( LM % SolidHarmonic_RS ( LM % nAngularMomentCells ) )
    allocate ( LM % SolidHarmonic_IC ( LM % nAngularMomentCells ) )
    allocate ( LM % SolidHarmonic_IS ( LM % nAngularMomentCells ) )
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

    ! if ( associated ( LM % MomentIrregularSin_1D ) ) &
    !    deallocate ( LM % MomentIrregularSin_1D )
    ! if ( associated ( LM % MomentIrregularCos_1D ) ) &
    !    deallocate ( LM % MomentIrregularCos_1D )
    ! if ( associated ( LM % MomentRegularSin_1D ) ) &
    !    deallocate ( LM % MomentRegularSin_1D )
    ! if ( associated ( LM % MomentRegularCos_1D ) ) &
    !    deallocate ( LM % MomentRegularCos_1D )
    ! if ( associated ( LM % MyMomentIrregularSin_1D ) ) &
    !    deallocate ( LM % MyMomentIrregularSin_1D )
    ! if ( associated ( LM % MyMomentIrregularCos_1D ) ) &
    !    deallocate ( LM % MyMomentIrregularCos_1D )
    ! if ( associated ( LM % MyMomentRegularSin_1D ) ) &
    !    deallocate ( LM % MyMomentRegularSin_1D )
    ! if ( associated ( LM % MyMomentRegularCos_1D ) ) &
    !    deallocate ( LM % MyMomentRegularCos_1D )

    ! associate ( nR_nA => LM % nRadialCells * LM % nAngularMomentCells )
    ! allocate ( LM % MyMomentRegularCos_1D ( nR_nA ) )
    ! allocate ( LM % MyMomentRegularSin_1D ( nR_nA ) )
    ! allocate ( LM % MyMomentIrregularCos_1D ( nR_nA ) )
    ! allocate ( LM % MyMomentIrregularSin_1D ( nR_nA ) )
    ! allocate ( LM % MomentRegularCos_1D ( nR_nA ) )
    ! allocate ( LM % MomentRegularSin_1D ( nR_nA ) )
    ! allocate ( LM % MomentIrregularCos_1D ( nR_nA ) )
    ! allocate ( LM % MomentIrregularSin_1D ( nR_nA ) )
    ! end associate !-- nR_nA

    ! !-- FIXME: NAG has trouble with pointer rank reassignment when the 
    ! !          left-hand side is a member
    ! call AssignPointers &
    !        ( LM, LM % MyMRC, LM % MyMRS, LM % MyMIC, LM % MyMIS, &
    !          LM % MRC, LM % MRS, LM % MIC, LM % MIS )

    ! if ( allocated ( LM % Reduction_RC ) ) deallocate ( LM % Reduction_RC )
    ! if ( allocated ( LM % Reduction_RS ) ) deallocate ( LM % Reduction_RS )
    ! if ( allocated ( LM % Reduction_IC ) ) deallocate ( LM % Reduction_IC )
    ! if ( allocated ( LM % Reduction_IS ) ) deallocate ( LM % Reduction_IS )

    ! associate ( PHC => PROGRAM_HEADER % Communicator )
    ! allocate ( LM % Reduction_RC )
    ! allocate ( LM % Reduction_RS )
    ! allocate ( LM % Reduction_IC )
    ! allocate ( LM % Reduction_IS )
    ! call LM % Reduction_RC % Initialize &
    !        ( PHC, OutgoingValue = LM % MyMomentRegularCos_1D, &
    !          IncomingValue = LM % MomentRegularCos_1D )
    ! call LM % Reduction_RS % Initialize &
    !        ( PHC, OutgoingValue = LM % MyMomentRegularSin_1D, &
    !          IncomingValue = LM % MomentRegularSin_1D )
    ! call LM % Reduction_IC % Initialize &
    !        ( PHC, OutgoingValue = LM % MyMomentIrregularCos_1D, &
    !          IncomingValue = LM % MomentIrregularCos_1D )
    ! call LM % Reduction_IS % Initialize &
    !        ( PHC, OutgoingValue = LM % MyMomentIrregularSin_1D, &
    !          IncomingValue = LM % MomentIrregularSin_1D )
    ! end associate !-- PHC

  end subroutine SetMomentStorage


end module LaplacianMultipole_Template
