module RayleighTaylor_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( UniverseTemplate ) :: RayleighTaylorForm
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type RayleighTaylorForm

contains


  subroutine Initialize ( RT, Name )

    class ( RayleighTaylorForm ), intent ( inout ) :: &
      RT
    character ( * ), intent ( in )  :: &
      Name

    integer ( KDI ) :: &
      iD  !-- iDimension
    real ( KDR ) :: &
      Acceleration
    type ( Character_1D_Form ), dimension ( 3 ) :: &
      BoundaryConditionsFace

    if ( RT % Type == '' ) &
      RT % Type = 'a RayleighTaylor'

    call RT % InitializeTemplate ( Name )

    Acceleration = 0.1_KDR
    call PROGRAM_HEADER % GetParameter ( Acceleration, 'Acceleration' )

    !-- Integrator

    allocate ( FluidBoxForm :: RT % Integrator )
    select type ( FB => RT % Integrator )
    type is ( FluidBoxForm )

    associate ( BCF => BoundaryConditionsFace )
    do iD = 1, 3
      call BCF ( iD ) % Initialize ( [ 'REFLECTING', 'REFLECTING' ] )     
    end do
    call FB % Initialize &
           ( Name, FluidType = 'IDEAL', GeometryType = 'NEWTONIAN', &
             BoundaryConditionsFaceOption = BCF, &
             GravitySolverTypeOption = 'UNIFORM', &
             UniformAccelerationOption = Acceleration )
    end associate !-- BCF

    select type ( PS => FB % PositionSpace )
    class is ( Atlas_SC_Form )


    !-- Cleanup

    end select !-- PS
    end select !-- FB

  end subroutine Initialize


  subroutine Finalize ( RT )

    type ( RayleighTaylorForm ), intent ( inout ) :: &
      RT

    call RT % FinalizeTemplate ( )

  end subroutine Finalize


end module RayleighTaylor_Form
