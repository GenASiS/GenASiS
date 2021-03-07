!-- Manifold_Form is a skeletal representation of a Manifold.

module Manifold_Form

  use Basics
  use MANIFOLD_Singleton
  use ManifoldHeader_Form
  use FieldsHeader_M__Form

  implicit none
  private

  type, public, extends ( ManifoldHeaderForm ) :: ManifoldForm
    integer ( KDI ) :: &
      nFieldSets = 0
    type ( FieldsHeader_M_Pointer ), dimension ( : ), allocatable :: &
      Fields
  contains
    procedure, private, pass :: &
      InitializeBasic
    procedure, public, pass :: &
      AddFields
    procedure, private, pass :: &
      Show_M
    final :: &
      Finalize
  end type ManifoldForm

    integer ( KDI ), private, parameter :: &
      MAX_FIELDS  =  MANIFOLD % MAX_FIELDS

contains


  subroutine InitializeBasic &
               ( M, Name, CommunicatorOption, nDimensionsOption, &
                 iDimensionalityOption )

    class ( ManifoldForm ), intent ( inout ) :: &
      M
    character ( * ), intent ( in )  :: &
      Name
    type ( CommunicatorForm ), intent ( in ), target, optional :: &
      CommunicatorOption
    integer ( KDI ), intent ( in ), optional :: &
      nDimensionsOption, &
      iDimensionalityOption

    call M % ManifoldHeaderForm % Initialize &
           ( Name, CommunicatorOption, nDimensionsOption, &
             iDimensionalityOption )

    allocate ( M % Fields ( MAX_FIELDS ) )

  end subroutine InitializeBasic


  subroutine AddFields ( M, Fields )

    class ( ManifoldForm ), intent ( inout ) :: &
      M
    class ( FieldsHeader_M_Form ), intent ( in ), target :: &
      Fields

    associate ( nF  =>  M % nFieldSets )

    nF  =  nF + 1

    M % Fields ( nF ) % Pointer  =>  Fields

    end associate !-- nF

  end subroutine AddFields


  subroutine Show_M ( M )

    class ( ManifoldForm ), intent ( in ) :: &
      M

    integer ( KDI ) :: &
      iF

    call M % ManifoldHeaderForm % Show ( )

    call Show ( M % nFieldSets, 'nFieldSets', M % IGNORABILITY )
    call Show ( [ ( M % Fields ( iF ) % Pointer % Name, &
                    iF = 1, M % nFieldSets) ], &
                'FieldSets', M % IGNORABILITY )

  end subroutine Show_M


  subroutine Finalize ( M )

    type ( ManifoldForm ), intent ( inout ) :: &
      M

    if ( allocated ( M % Fields ) ) &
      deallocate ( M % Fields )

  end subroutine Finalize


end module Manifold_Form
