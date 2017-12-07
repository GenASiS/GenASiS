!-- Geometry_ASC represents the geometry on a position space represented
!   by a single chart.

module Geometry_ASC__Form

  !-- Geometry_AtlasSingleChart_Form

  use Basics
  use Mathematics
  use Geometry_G__Form
  use Geometry_N__Form
  use Geometry_CSL__Form

  implicit none
  private

  type, public, extends ( GeometryFlat_ASC_Form ) :: Geometry_ASC_Form
    real ( KDR ) :: &
      GravitationalConstant
    type ( Storage_ASC_Form ), allocatable :: &
      Storage_ASC
    type ( GradientForm ), allocatable :: &
      Gradient
    type ( Poisson_ASC_Form ), allocatable :: &
      Poisson_ASC
  contains
    procedure, public, pass :: &
      Initialize
    procedure, private, pass :: &
      Geometry_G_CSL
    generic, public :: &
      Geometry_G => Geometry_G_CSL
    procedure, private, pass :: &
      Geometry_N_CSL
    generic, public :: &
      Geometry_N => Geometry_N_CSL
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type Geometry_ASC_Form

contains


  subroutine Initialize &
               ( GA, A, GeometryType, NameShortOption, &
                 GravitySolverTypeOption, GravitationalConstantOption, &
                 IgnorabilityOption )

    class ( Geometry_ASC_Form ), intent ( inout ) :: &
      GA
    class ( Atlas_SC_Form ), intent ( in ), target :: &
      A
    character ( * ), intent ( in ) :: &
      GeometryType
    character ( * ), intent ( in ), optional :: &
      NameShortOption, &
      GravitySolverTypeOption
    real ( KDR ), intent ( in ), optional :: &
      GravitationalConstantOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    integer ( KDI ) :: &
      MaxDegree
    character ( LDL ) :: &
      NameShort
    class ( Geometry_N_Form ), pointer :: &
      G

    if ( GA % Type == '' ) &
      GA % Type = 'a Geometry_ASC'

    GA % GeometryType = GeometryType    

    call GA % InitializeFlat ( A, NameShortOption, IgnorabilityOption )

    select case ( trim ( GA % GeometryType ) )
    case ( 'NEWTONIAN' )

      GA % GravitationalConstant  =  CONSTANT % GRAVITATIONAL
      if ( present ( GravitationalConstantOption ) ) &
        GA % GravitationalConstant  =  GravitationalConstantOption

      G => GA % Geometry_N ( )

      allocate ( GA % Storage_ASC )
      associate ( SA => GA % Storage_ASC )
      call SA % Initialize &
             ( GA, NameShort = 'PoissonStorage', &
               iaSelectedOption = [ G % GRAVITATIONAL_POTENTIAL ] )
      end associate !-- SA

      allocate ( GA % Gradient )
      associate ( Grad => GA % Gradient )
      call Grad % Initialize &
             ( 'GeometryGradient', [ G % nValues, 1 ] )
      end associate !-- Grad

      if ( present ( GravitySolverTypeOption ) ) then

        MaxDegree = 10
        call PROGRAM_HEADER % GetParameter ( MaxDegree, 'MaxDegree' )

        allocate ( GA % Poisson_ASC )
        associate ( PA => GA % Poisson_ASC )
        call PA % Initialize &
               ( A, SolverType = GravitySolverTypeOption, &
                 MaxDegreeOption = MaxDegree, &
                 nEquationsOption = 1 )
        end associate !-- PA

      else
        call Show ( 'NEWTONIAN geometry requires GravitySolverType', &
                    CONSOLE % ERROR )
        call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
        call Show ( 'Geometry_ASC__Form', 'module', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if

      nullify ( G )

    end select !-- GeometryType

  end subroutine Initialize


  function Geometry_G_CSL ( GA ) result ( G )

    class ( Geometry_ASC_Form ), intent ( in ) :: &
      GA
    class ( Geometry_G_Form ), pointer :: &
      G

    select type ( GC => GA % Chart )
    class is ( Geometry_CSL_Form )
      G => GC % Geometry_G ( )
    class default
      call Show ( 'Geometry Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Geometry_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Geometry_G_CSL', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- GC

  end function Geometry_G_CSL


  function Geometry_N_CSL ( GA ) result ( G )

    class ( Geometry_ASC_Form ), intent ( in ) :: &
      GA
    class ( Geometry_N_Form ), pointer :: &
      G

    select type ( GC => GA % Chart )
    class is ( Geometry_CSL_Form )
      G => GC % Geometry_N ( )
    class default
      call Show ( 'Geometry Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Geometry_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Geometry_N_CSL', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- GC

  end function Geometry_N_CSL


  impure elemental subroutine Finalize ( GA )

    type ( Geometry_ASC_Form ), intent ( inout ) :: &
      GA

    if ( allocated ( GA % Poisson_ASC ) ) &
      deallocate ( GA % Poisson_ASC )
    if ( allocated ( GA % Gradient ) ) &
      deallocate ( GA % Gradient )
    if ( allocated ( GA % Storage_ASC ) ) &
      deallocate ( GA % Storage_ASC )

  end subroutine Finalize


  subroutine SetField ( FA )

    class ( Geometry_ASC_Form ), intent ( inout ) :: &
      FA

    select type ( A => FA % Atlas )
    class is ( Atlas_SC_Template )

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
    associate ( nValues => C % nProperCells + C % nGhostCells )

    allocate ( Geometry_CSL_Form :: FA % Chart )

    select type ( FC => FA % Chart )
    class is ( Geometry_CSL_Form )
      call FC % Initialize &
             ( C, FA % NameShort, FA % GeometryType, nValues, &
               IgnorabilityOption = FA % IGNORABILITY )
    end select !-- FC

    end associate !-- nValues
    end select !-- C
    end select !-- A

  end subroutine SetField


end module Geometry_ASC__Form
