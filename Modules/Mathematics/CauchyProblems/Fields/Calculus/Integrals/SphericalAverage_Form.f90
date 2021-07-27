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
    real ( KDI ), dimension ( :, : ), allocatable :: &
      SolidAngle
    class ( FieldSetForm ), allocatable :: &
      FieldSet_SA
    class ( FieldSetForm ), pointer :: &
      FieldSet => null ( )
  contains
    procedure, private, pass :: &
      Initialize_SA
    generic, public :: &
      Initialize => Initialize_SA
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type SphericalAverageForm

    private :: &
      ComputeSolidAngle


contains


  subroutine Initialize_SA ( SA, G, FS, A_SA, iaAverageOption )

    class ( SphericalAverageForm ), intent ( inout ) :: &
      SA
    class ( Geometry_F_Form ), intent ( in ), target :: &
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
               NameOption = FS % Name, &
               UnitOption = FS % Unit, &
               VectorIndicesOption = FS % VectorIndices, &
               nFieldsOption = FS % nFields, &
               IgnorabilityOption = FS % Ignorability )
      end associate !-- FS_SA
    end if !-- allocated FS_SA

  end subroutine Initialize_SA


  subroutine Compute ( SA )

    class ( SphericalAverageForm ), intent ( inout ) :: &
      SA      

  end subroutine Compute


  impure elemental subroutine Finalize ( SA )

    type ( SphericalAverageForm ), intent ( inout ) :: &
      SA

    if ( allocated ( SA % FieldSet_SA ) ) &
      deallocate ( SA % FieldSet_SA )
    if ( allocated ( SA % SolidAngle ) ) &
      deallocate ( SA % SolidAngle )
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
    associate &
      ( C  =>  A % Chart_GS )
    select case ( trim ( C % CoordinateSystem ) )
      case ( 'SPHERICAL' )

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

    allocate ( SA % SolidAngle ( nTh, nPh ) )

    select case ( nD )
    case ( 1 )
      SA % SolidAngle ( 1, 1 )  =  FourPi
    case ( 2 )
      do iTh  =  1, nTh
        SA % SolidAngle ( iTh, 1 )  &
          =  TwoPi  *  ( cos ( Th_I ( iTh ) )  -  cos ( Th_I ( iTh + 1 ) ) )
      end do !-- iTh
    case ( 3 )
      do iPh  =  1, nPh
        do iTh  =  1, nTh
          SA % SolidAngle ( iTh, iPh )  &
            =  dPh ( iPh )  &
               *  ( cos ( Th_I ( iTh ) )  -  cos ( Th_I ( iTh + 1 ) ) )
        end do !-- iTh
      end do !-- iPh
    end select !-- nD

    end associate !-- nD, etc.

    case default
      call Show ( 'CoordinateSystem not recognized', CONSOLE % ERROR )
      call Show ( C % CoordinateSystem, 'CoordinateSystem', CONSOLE % ERROR )
      call Show ( 'SphericalAverage_Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeSolidAngle', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- CoordinateSystem

    end associate !-- C

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'SphericalAverage_Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeSolidAngle', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

    nullify ( Th_I, dPh )

  end subroutine ComputeSolidAngle


end module SphericalAverage_Form
