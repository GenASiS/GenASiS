module SphericalAverage_Form

  use Basics
  use Manifolds
  use FieldSets
  use Geometries

  implicit none
  private

  type, public :: SphericalAverageForm
    integer ( KDI ) :: &
      nAverages
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaAverage
    real ( KDR ), dimension ( :, : ), allocatable :: &
      dSolidAngle
    class ( FieldSetForm ), allocatable :: &
      FieldSet_SA
    class ( FieldSetForm ), pointer :: &
      FieldSet => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type SphericalAverageForm

    private :: &
      ComputeSolidAngle, &
      ComputeAverage_CGS


contains


  subroutine Initialize ( SA, G, FS, A_SA, iaAverageOption )

    class ( SphericalAverageForm ), intent ( inout ) :: &
      SA
    class ( Geometry_F_Form ), intent ( in ) :: &
      G
    class ( FieldSetForm ), intent ( in ), target :: &
      FS
    class ( Atlas_SCG_Form ), intent ( in ) :: &
      A_SA
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      iaAverageOption
    
    SA % FieldSet  =>  FS

    if ( present ( iaAverageOption ) ) then
      allocate ( SA % iaAverage, source = iaAverageOption )
    else
      allocate ( SA % iaAverage, source = FS % iaSelected )
    end if

    SA % nAverages  =  size ( SA % iaAverage )

    call ComputeSolidAngle ( SA, G )

    if ( .not. allocated ( SA % FieldSet_SA ) ) then
      allocate ( SA % FieldSet_SA )
      associate ( FS_SA  =>  SA % FieldSet_SA )
      call FS_SA % Initialize &
             ( A_SA, &
               FieldOption = FS % Field, &
               VectorOption = FS % Vector, &
               NameOption = trim ( FS % Name ) // '_SA', &
               UnitOption = FS % Unit, &
               VectorIndicesOption = FS % VectorIndices, &
               nFieldsOption = FS % nFields, &
               IgnorabilityOption = FS % Ignorability )
      end associate !-- FS_SA
    end if !-- allocated FS_SA

  end subroutine Initialize


  subroutine Compute ( SA, IgnorabilityOption )

    class ( SphericalAverageForm ), intent ( inout ) :: &
      SA      
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    integer ( KDI ) :: &
      Ignorability
    real ( KDR ), dimension ( :, :, :, : ), pointer :: &
      F, &
      F_SA

    Ignorability  =  CONSOLE % INFO_5
    if ( present ( IgnorabilityOption ) ) &
      Ignorability  =  IgnorabilityOption

    associate &
      ( FS     =>  SA % FieldSet, &
        FS_SA  =>  SA % FieldSet_SA )

    call Show ( 'Computing a SphericalAverage', Ignorability )
    call Show ( FS % Name, 'Integrand', Ignorability )
    call Show ( FS % Atlas % Name, 'Atlas', Ignorability )

    select type ( A  =>  FS % Atlas )
      class is ( Atlas_SCG_Form )
    select type ( C  =>  A % Chart_GS )
      class is ( Chart_GS_C_Form )
    select type ( A_SA  =>  FS_SA % Atlas )
      class is ( Atlas_SCG_Form )
    select type ( C_SA  =>  A_SA % Chart_GS )
      class is ( Chart_GS_C_Form )
    associate &
      ( FV     =>  FS    % Storage_GS % Value, &
        FV_SA  =>  FS_SA % Storage_GS % Value )

    call C    % SetFieldPointer ( FV,    F    )
    call C_SA % SetFieldPointer ( FV_SA, F_SA )

    call ComputeAverage_CGS &
           ( F_SA, F, SA % dSolidAngle, SA % iaAverage, C % nCellsBrick, &
             C % nGhostLayers ) 

    end associate !-- FV, etc.
    end select !-- C_SA
    end select !-- A_SA

    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'SphericalAverage_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'SphericalAverage_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

    end associate !-- FS, etc.

  end subroutine Compute


  impure elemental subroutine Finalize ( SA )

    type ( SphericalAverageForm ), intent ( inout ) :: &
      SA

    if ( allocated ( SA % FieldSet_SA ) ) &
      deallocate ( SA % FieldSet_SA )
    if ( allocated ( SA % dSolidAngle ) ) &
      deallocate ( SA % dSolidAngle )
    if ( allocated ( SA % iaAverage ) ) &
      deallocate ( SA % iaAverage )

    nullify ( SA % FieldSet )

  end subroutine Finalize


  subroutine ComputeSolidAngle ( SA, G )

    class ( SphericalAverageForm ), intent ( inout ) :: &
      SA
    class ( Geometry_F_Form ), intent ( in ), target :: &
      G

    integer ( KDI ) :: &
      iTh, &  !-- iTheta
      iPh     !-- iPhi
    real ( KDR ) :: &
      TwoPi, FourPi
    real ( KDR ), dimension ( : ), pointer :: &
       Th_I, &
      dPh

     TwoPi  =  2.0_KDR  *  CONSTANT % PI
    FourPi  =  4.0_KDR  *  CONSTANT % PI

    select type ( A  =>  G % Atlas )
      class is ( Atlas_SCG_Form )
    select type ( C  =>  A % Chart_GS )
      class is ( Chart_GS_C_Form )

    associate &
      ( nD   =>  C % nDimensions, &
        nTh  =>  C % nCells ( 2 ), &
        nPh  =>  C % nCells ( 3 ) )

    Th_I  =>  null ( )
    if ( C % nDimensions  >  1 ) &
      Th_I  =>  C % Edge ( 2 ) % Value ( 1 : nTh + 1 )

    dPh  => null ( )
    if ( C % nDimensions  >  2 ) &
      dPh  =>  C % Width ( 3 ) % Value ( 1 : nPh )

    ! if ( associated ( Th_I ) ) &
    !   call Show ( Th_I, 'Theta_I' )
    ! if ( associated ( dPh ) ) &
    !   call Show ( dPh, 'dPhi' )

    allocate ( SA % dSolidAngle ( nTh, nPh ) )

    select case ( nD )
    case ( 1 )
      SA % dSolidAngle ( 1, 1 )  =  FourPi
    case ( 2 )
      do iTh  =  1, nTh
        SA % dSolidAngle ( iTh, 1 )  &
          =  TwoPi  *  ( cos ( Th_I ( iTh ) )  -  cos ( Th_I ( iTh + 1 ) ) )
      end do !-- iTh
    case ( 3 )
      do iPh  =  1, nPh
        do iTh  =  1, nTh
          SA % dSolidAngle ( iTh, iPh )  &
            =  dPh ( iPh )  &
               *  ( cos ( Th_I ( iTh ) )  -  cos ( Th_I ( iTh + 1 ) ) )
        end do !-- iTh
      end do !-- iPh
    end select !-- nD

    end associate !-- nD, etc.

    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'SphericalAverage_Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeSolidAngle', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'SphericalAverage_Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeSolidAngle', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

    nullify ( Th_I, dPh )

  end subroutine ComputeSolidAngle


  subroutine ComputeAverage_CGS ( F_SA, F, dSA, iaAvg, nC, oC )

    real ( KDR ), dimension ( :, :, :, : ), intent ( inout ) :: &
      F_SA
    real ( KDR ), dimension ( :, :, :, : ), intent ( in ) :: &
      F
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      dSA  !-- dSolidAngle
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iaAvg, &
      nC, oC  !-- nCells, oCell

    integer ( KDI ) :: &
      iR, iT, iP, &  !-- iRadius, iTheta, iPhi
      iS, &  !-- iSelected
      iF     !-- iField
    real ( KDR ) :: &
      FourPi

    FourPi  =  4.0_KDR  *  CONSTANT % PI

    !-- Clear variables to be averaged
    do iS  =  1, size ( iaAvg )
      iF  =  iaAvg ( iS )
      F_SA ( :, 1, 1, iF )  =  0.0_KDR
    end do !-- iS

    !$OMP parallel do collapse ( 4 ) &
    !$OMP schedule ( OMP_SCHEDULE_HOST ) private ( iF ) &
    !$OMP reduction ( + : F_SA )
    do iS  =  1, size ( iaAvg )
      do iP  =  1,  nC ( 3 )
        do iT  =  1,  nC ( 2 )
          do iR  =  1,  nC ( 1 )

            iF  =  iaAvg ( iS )

            F_SA ( oC ( 1 )  +  iR,  1,  1,  iF )  &
              =   F_SA ( oC ( 1 )  +  iR,  1,  1 ,  iF )  &
                  +  dSA ( iT, iP )  &
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
      F_SA ( :, 1, 1, iF )  =  F_SA  ( :, 1, 1, iF )  /  FourPi
    end do !-- iS

  end subroutine ComputeAverage_CGS


end module SphericalAverage_Form
