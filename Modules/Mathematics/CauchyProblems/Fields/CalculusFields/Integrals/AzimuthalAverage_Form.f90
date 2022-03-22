module AzimuthalAverage_Form

  use Basics
  use Manifolds
  use FieldSets
  use Geometries

  implicit none
  private

  type, public :: AzimuthalAverageForm
    integer ( KDI ) :: &
      nAverages
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaAverage
    class ( FieldSetForm ), allocatable :: &
      FieldSet_AA
    class ( FieldSetForm ), pointer :: &
      FieldSet => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type AzimuthalAverageForm

    private :: &
      ComputeAverage_CGS


contains


  subroutine Initialize ( AA, G, FS, A_AA, iaAverageOption )

    class ( AzimuthalAverageForm ), intent ( inout ) :: &
      AA
    class ( Geometry_F_Form ), intent ( in ) :: &
      G
    class ( FieldSetForm ), intent ( in ), target :: &
      FS
    class ( Atlas_SCG_Form ), intent ( in ) :: &
      A_AA
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      iaAverageOption
    
    AA % FieldSet  =>  FS

    if ( present ( iaAverageOption ) ) then
      allocate ( AA % iaAverage, source = iaAverageOption )
    else
      allocate ( AA % iaAverage, source = FS % iaSelected )
    end if

    AA % nAverages  =  size ( AA % iaAverage )

    if ( .not. allocated ( AA % FieldSet_AA ) ) then
      allocate ( AA % FieldSet_AA )
      associate ( FS_AA  =>  AA % FieldSet_AA )
      call FS_AA % Initialize &
             ( A_AA, &
               FieldOption = FS % Field, &
               VectorOption = FS % Vector, &
               NameOption = trim ( FS % Name ) // '_AA', &
               UnitOption = FS % Unit, &
               VectorIndicesOption = FS % VectorIndices, &
               nFieldsOption = FS % nFields, &
               IgnorabilityOption = FS % Ignorability )
      end associate !-- FS_AA
    end if !-- allocated FS_AA

  end subroutine Initialize


  subroutine Compute ( AA, IgnorabilityOption )

    class ( AzimuthalAverageForm ), intent ( inout ) :: &
      AA      
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    integer ( KDI ) :: &
      Ignorability
    real ( KDR ), dimension ( : ), pointer :: &
      dPh
    real ( KDR ), dimension ( :, :, :, : ), pointer :: &
      F, &
      F_AA

    Ignorability  =  CONSOLE % INFO_5
    if ( present ( IgnorabilityOption ) ) &
      Ignorability  =  IgnorabilityOption

    associate &
      ( FS     =>  AA % FieldSet, &
        FS_AA  =>  AA % FieldSet_AA )

    call Show ( 'Computing a AzimuthalAverage', Ignorability )
    call Show ( FS % Name, 'Integrand', Ignorability )
    call Show ( FS % Atlas % Name, 'Atlas', Ignorability )

    select type ( A  =>  FS % Atlas )
      class is ( Atlas_SCG_Form )
    select type ( C  =>  A % Chart_GS )
      class is ( Chart_GS_C_Form )
    select type ( A_AA  =>  FS_AA % Atlas )
      class is ( Atlas_SCG_Form )
    select type ( C_AA  =>  A_AA % Chart_GS )
      class is ( Chart_GS_C_Form )
    associate &
      ( FV     =>  FS    % Storage_GS % Value, &
        FV_AA  =>  FS_AA % Storage_GS % Value )

    call C    % SetFieldPointer ( FV,    F    )
    call C_AA % SetFieldPointer ( FV_AA, F_AA )

    if ( C % nDimensions  >  2 ) then
      associate ( nPh  =>  C % nCells ( 3 ) )
      dPh  =>  C % Width ( 3 ) % Value ( 1 : nPh )
      call ComputeAverage_CGS &
             ( F_AA, F, dPh, AA % iaAverage, C % nCellsBrick, &
               C % nGhostLayers ) 
      end associate !-- nPh
    end if !-- nDimensions > 2

    end associate !-- FV, etc.
    end select !-- C_AA
    end select !-- A_AA

    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'AzimuthalAverage_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'AzimuthalAverage_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

    end associate !-- FS, etc.

  end subroutine Compute


  impure elemental subroutine Finalize ( AA )

    type ( AzimuthalAverageForm ), intent ( inout ) :: &
      AA

    if ( allocated ( AA % FieldSet_AA ) ) &
      deallocate ( AA % FieldSet_AA )
    if ( allocated ( AA % iaAverage ) ) &
      deallocate ( AA % iaAverage )

    nullify ( AA % FieldSet )

  end subroutine Finalize


  subroutine ComputeAverage_CGS ( F_AA, F, dPhi, iaAvg, nC, oC )

    real ( KDR ), dimension ( :, :, :, : ), intent ( inout ) :: &
      F_AA
    real ( KDR ), dimension ( :, :, :, : ), intent ( in ) :: &
      F
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      dPhi  !-- dAzimuthalAngle
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iaAvg, &
      nC, oC  !-- nCells, oCell

    integer ( KDI ) :: &
      iR, iT, iP, &  !-- iRadius, iTheta, iPhi
      iS, &  !-- iSelected
      iF     !-- iField
    real ( KDR ) :: &
      TwoPi

    TwoPi  =  2.0_KDR  *  CONSTANT % PI

    !-- Clear variables to be averaged
    do iS  =  1, size ( iaAvg )
      iF  =  iaAvg ( iS )
      F_AA ( :, :, 1, iF )  =  0.0_KDR
    end do !-- iS

    !$OMP parallel do collapse ( 4 ) &
    !$OMP schedule ( OMP_SCHEDULE_HOST ) private ( iF ) &
    !$OMP reduction ( + : F_AA )
    do iS  =  1, size ( iaAvg )
      do iP  =  1,  nC ( 3 )
        do iT  =  1,  nC ( 2 )
          do iR  =  1,  nC ( 1 )

            iF  =  iaAvg ( iS )

            F_AA ( oC ( 1 )  +  iR,  oC ( 2 )  +  iT,  1,  iF )  &
              =   F_AA ( oC ( 1 )  +  iR,  oC ( 2 )  +  iT,  1 ,  iF )  &
                  +  dPhi ( iP )  &
                     *  F ( oC ( 1 )  +  iR, &
                            oC ( 2 )  +  iT, &
                            oC ( 3 )  +  iP, &
                            iF )

          end do !-- iR
        end do !-- iT
      end do !-- iP
    end do !-- iS
    !$OMP  end parallel do      

    !-- Normalize variables to be averaged
    do iS  =  1, size ( iaAvg )
      iF  =  iaAvg ( iS )
      F_AA ( :, :, 1, iF )  =  F_AA  ( :, :, 1, iF )  /  TwoPi
    end do !-- iS

  end subroutine ComputeAverage_CGS


end module AzimuthalAverage_Form
