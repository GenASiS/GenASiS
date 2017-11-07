module Poisson_ASC__Form_Test__Form

  use Basics
  use Manifolds
  use Poisson_ASC__Form
  use CreateProportionalChart_Command
  use SetHomogeneousSphere_Command

  implicit none
  private

  integer ( KDI ), private, parameter :: &
    O_HS                 = 0, &  !-- OFFSET_HOMOGENEOUS_SPHERE
    N_HOMOGENEOUS_SPHERE = 3, &
    N_EQUATIONS          = N_HOMOGENEOUS_SPHERE
  character ( LDL ), dimension ( N_EQUATIONS ), private, parameter :: &
    VARIABLE = [ 'HomegeneousSphere_1', &
                 'HomogeneousSphere_2', &
                 'HomogeneousSphere_3' ]
  character ( LDF ), public, parameter :: &
    ProgramName = 'Poisson_ASC__Form_Test'

  type, public :: Poisson_ASC__Form_Test_Form
    type ( GridImageStreamForm ), allocatable :: &
      Stream
    type ( Atlas_SC_Form ), allocatable :: &
      Atlas
    type ( Storage_ASC_Form ), allocatable :: &
      Source, &
      Solution, &
      Reference
    type ( Poisson_ASC_Form ), allocatable :: &
      Poisson
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type Poisson_ASC__Form_Test_Form

contains


  subroutine Initialize ( PFT, Name )

    class ( Poisson_ASC__Form_Test_Form ) :: &
      PFT
    character ( * ), intent ( in ) :: &
      Name

    integer ( KDI ) :: &
      iE, &  !-- iEquation
      MaxDegree
    real ( KDR ), dimension ( N_EQUATIONS ) :: &
      RadiusDensity, &
      Density

    !-- Atlas

    allocate ( PFT % Atlas )
    associate ( A => PFT % Atlas )
    call A % Initialize ( Name, PROGRAM_HEADER % Communicator )

    call CreateProportionalChart ( A )


    !-- Poisson

    MaxDegree = 2
    call PROGRAM_HEADER % GetParameter ( MaxDegree, 'MaxDegree' )

    allocate ( PFT % Poisson )
    associate ( P => PFT % Poisson )
    call P % Initialize &
           ( A, SolverType = 'MULTIPOLE', MaxDegreeOption = MaxDegree, &
             nEquationsOption = N_EQUATIONS )


    !-- Source, Reference

    allocate ( PFT % Source )
    associate ( SA => PFT % Source )
    call SA % Initialize &
           ( A, 'Source', N_EQUATIONS, &
             VariableOption = VARIABLE, &
             WriteOption = .true. )

    allocate ( PFT % Reference )
    associate ( RA => PFT % Reference )
    call RA % Initialize &
           ( A, 'Reference', N_EQUATIONS, &
             VariableOption = VARIABLE, &
             WriteOption = .true. )


    !-- Homogeneous sphere

    RadiusDensity  &
      =  A % Chart % MaxCoordinate ( 1 ) / [ 1.1_KDR, 2.0_KDR, 10.0_KDR ] 
    call PROGRAM_HEADER % GetParameter ( RadiusDensity, 'RadiusDensity' )

    Density = 1.0_KDR / ( 4.0_KDR  *  CONSTANT % PI  *  RadiusDensity ** 3 )
    call PROGRAM_HEADER % GetParameter ( Density, 'Density' )

    do iE = O_HS + 1, O_HS + N_HOMOGENEOUS_SPHERE
      call SetHomogeneousSphere &
             ( SA, RA, A, Density ( iE ), RadiusDensity ( iE ), iE )
    end do


    !-- Solution

    allocate ( PFT % Solution )
    associate ( SA => PFT % Solution )
    call SA % Initialize &
           ( A, 'Solution', N_EQUATIONS, &
             VariableOption = VARIABLE, &
             WriteOption = .true. )
    end associate !-- SA

    call P % Solve ( PFT % Solution, PFT % Source )


    !-- Write

    allocate ( PFT % Stream )
    call PFT % Stream % Initialize &
           ( Name, CommunicatorOption = PROGRAM_HEADER % Communicator )    
    associate ( GIS => PFT % Stream )

    call A % OpenStream ( GIS, 'Stream', iStream = 1 )

    call GIS % Open ( GIS % ACCESS_CREATE )
    call A % Write ( iStream = 1 )
    call GIS % Close ( )

    end associate !-- GIS


    !-- Cleanup

    end associate !-- RA
    end associate !-- SA
    end associate !-- P
    end associate !-- A

  end subroutine Initialize


  subroutine Finalize ( PFT )

    type ( Poisson_ASC__Form_Test_Form ) :: &
      PFT

    if ( allocated ( PFT % Stream ) ) &
      deallocate ( PFT % Stream )
    if ( allocated ( PFT % Reference ) ) &
      deallocate ( PFT % Reference )
    if ( allocated ( PFT % Solution ) ) &
      deallocate ( PFT % Solution )
    if ( allocated ( PFT % Source ) ) &
      deallocate ( PFT % Source )
    if ( allocated ( PFT % Poisson ) ) &
      deallocate ( PFT % Poisson )
    if ( allocated ( PFT % Atlas ) ) &
      deallocate ( PFT % Atlas )

  end subroutine Finalize


end module Poisson_ASC__Form_Test__Form



program Poisson_ASC__Form_Test

  use Basics
  use Poisson_ASC__Form_Test__Form

  implicit none

  type ( Poisson_ASC__Form_Test_Form ), allocatable :: &
    PFT

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( ProgramName )

  allocate ( PFT )
  call PFT % Initialize ( PROGRAM_HEADER % Name )
  deallocate ( PFT )

  deallocate ( PROGRAM_HEADER )

end program Poisson_ASC__Form_Test
