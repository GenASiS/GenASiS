module SphericalAverage_Form

  use Basics
  use Manifolds
  use FieldSets
  use Geometries

  implicit none
  private

  type, public, extends ( FieldSetForm ) :: SphericalAverageForm
    integer ( KDI ) :: &
      nAverages
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaAverage
    real ( KDI ), dimension ( :, : ), allocatable :: &
      SolidAngle
    class ( FieldSetForm ), pointer :: &
      Integrand => null ( )
    class ( Geometry_F_Form ), pointer :: &
      Geometry => null ( )
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


  subroutine Initialize_SA ( SA, G, I, A_SA, iaAverageOption )

    class ( SphericalAverageForm ), intent ( inout ) :: &
      SA
    class ( Geometry_F_Form ), intent ( in ), target :: &
      G
    class ( FieldSetForm ), intent ( in ), target :: &
      I
    class ( Atlas_SCG_Form ), intent ( in ) :: &
      A_SA
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      iaAverageOption
    
    SA % Integrand  =>  I
    SA % Geometry   =>  G

    if ( present ( iaAverageOption ) ) then
      allocate ( SA % iaAverage, source = iaAverageOption )
    else
      allocate ( SA % iaAverage, source = I % iaSelected )
    end if

    SA % nAverages  =  size ( SA % iaAverage )

    call ComputeSolidAngle ( SA )

    call SA % Initialize &
           ( A_SA, &
             FieldOption = I % Field, &
             VectorOption = I % Vector, &
             NameOption = I % Name, &
             UnitOption = I % Unit, &
             VectorIndicesOption = I % VectorIndices, &
             nFieldsOption = I % nFields, &
             IgnorabilityOption = I % Ignorability )

  end subroutine Initialize_SA


  subroutine Compute ( SA )

    class ( SphericalAverageForm ), intent ( inout ) :: &
      SA      

  end subroutine Compute


  impure elemental subroutine Finalize ( SA )

    type ( SphericalAverageForm ), intent ( inout ) :: &
      SA

    if ( allocated ( SA % SolidAngle ) ) &
      deallocate ( SA % SolidAngle )
    if ( allocated ( SA % iaAverage ) ) &
      deallocate ( SA % iaAverage )

    nullify ( SA % Geometry )
    nullify ( SA % Integrand )

  end subroutine Finalize


  subroutine ComputeSolidAngle ( SA )

    class ( SphericalAverageForm ), intent ( inout ) :: &
      SA

    integer ( KDI ) :: &
      iT, &  !-- iTheta
      iP     !-- iPhi

    select type ( A  =>  SA % Geometry % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS )
    select case ( trim ( C % CoordinateSystem ) )
      case ( 'SPHERICAL' )

    if ( C % nDimensions  >  1 ) &
      call Show ( C % Edge ( 2 ) % Value ( 1 : C % nCellsBrick ( 2 )  +  1 ), &
                  'Theta' )
    if ( C % nDimensions  >  2 ) &
      call Show ( C % Edge ( 3 ) % Value ( 1 : C % nCellsBrick ( 3 )  +  1 ), &
                  'Phi' )

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

  end subroutine ComputeSolidAngle


end module SphericalAverage_Form
