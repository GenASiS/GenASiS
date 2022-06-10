module Poisson_ASC__Form_Test__Form

  use Basics
  use Manifolds
  use Operations
  use Poisson_ASC__Form
  use SetHomogeneousSphere_Command

  implicit none
  private

  integer ( KDI ), private, parameter :: &
    N_EQUATIONS            = 3
  character ( LDL ), dimension ( N_EQUATIONS ), private, parameter :: &
    VARIABLE = [ 'HomogeneousSphere_1  ', &
                 'HomogeneousSphere_2  ', &
                 'HomogeneousSphere_3  ' ]
  character ( LDF ), public, parameter :: &
    ProgramName = 'Poisson_ASC__Form_Test'

  type, public :: Poisson_ASC__Form_Test_Form
    type ( GridImageStreamForm ), allocatable :: &
      Stream
    type ( Atlas_SC_CC_Form ), allocatable :: &
      Atlas
    type ( Storage_ASC_Form ), allocatable :: &
      Source, &
      Solution, &
      Reference, &
      Difference
    type ( Poisson_ASC_Form ), allocatable :: &
      Poisson
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      ComputeError
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
    character ( LDL ) :: &
      SolverType

    !-- Atlas

    allocate ( PFT % Atlas )
    associate ( A => PFT % Atlas )
    call A % Initialize ( Name, PROGRAM_HEADER % Communicator )
    call A % CreateChart_CC ( )
    call A % SetGeometry ( )


    !-- Poisson

    SolverType = 'MULTIPOLE'
    call PROGRAM_HEADER % GetParameter ( SolverType, 'SolverType' )

    MaxDegree = 2
    call PROGRAM_HEADER % GetParameter ( MaxDegree, 'MaxDegree' )

    allocate ( PFT % Poisson )
    associate ( P => PFT % Poisson )
    call P % Initialize &
           ( A, SolverType = SolverType, MaxDegreeOption = MaxDegree, &
             nEquationsOption = N_EQUATIONS )
    call P % InitializeTimers ( BaseLevel = 1 )


    !-- Source, Reference

    allocate ( PFT % Source )
    associate ( SA => PFT % Source )
    call SA % Initialize &
           ( A, 'Source', N_EQUATIONS, &
             VariableOption = VARIABLE, &
             WriteOption = .true. )
    end associate !-- SA

    allocate ( PFT % Reference )
    associate ( RA => PFT % Reference )
    call RA % Initialize &
           ( A, 'Reference', N_EQUATIONS, &
             VariableOption = VARIABLE, &
             WriteOption = .true. )
    end associate !-- RA


    !-- Homogeneous sphere

    RadiusDensity  &
      =  A % Chart % MaxCoordinate ( 1 ) / [ 1.1_KDR, 2.0_KDR, 10.0_KDR ] 
    call PROGRAM_HEADER % GetParameter ( RadiusDensity, 'RadiusDensity' )

    Density ( 1 : N_EQUATIONS ) &
      = 1.0_KDR / ( 4.0_KDR  *  CONSTANT % PI  *  RadiusDensity ** 3 )
    call PROGRAM_HEADER % GetParameter ( Density, 'Density' )

    do iE = 1, N_EQUATIONS
      call SetHomogeneousSphere &
             ( PFT % Source, PFT % Reference, A, Density ( iE ), &
               RadiusDensity ( iE ), iE )
    end do


    !-- Solution, Difference

    allocate ( PFT % Solution )
    associate ( SA => PFT % Solution )
    call SA % Initialize &
           ( A, 'Solution', N_EQUATIONS, &
             VariableOption = VARIABLE, &
             WriteOption = .true. )
    end associate !-- SA

    allocate ( PFT % Difference)
    associate ( DA => PFT % Difference )
    call DA % Initialize &
           ( A, 'Difference', N_EQUATIONS, &
             VariableOption = VARIABLE, &
             WriteOption = .true. )
    end associate !-- DA

    call P % Solve ( PFT % Solution, PFT % Source )
    
    call PFT % ComputeError ( ) 


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

    end associate !-- P
    end associate !-- A

  end subroutine Initialize


  subroutine ComputeError ( PFT )

    class ( Poisson_ASC__Form_Test_Form ), intent ( in ) :: &
      PFT
    
    real ( KDR ) :: &
      L1_HS_1, &
      L1_HS_2, &
      L1_HS_3
    class ( StorageForm ), pointer :: &
      Solution, &
      Reference, &
      Difference         
    type ( CollectiveOperation_R_Form ) :: &
      CO
    
    Solution   => PFT % Solution   % Storage ( )
    Reference  => PFT % Reference  % Storage ( )
    Difference => PFT % Difference % Storage ( )

    call MultiplyAdd &
           ( Solution % Value, Reference % Value, -1.0_KDR, Difference % Value )

    associate ( A => PFT % Atlas ) 
    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
    
    associate ( HS_1   => Difference % Value ( :, 1 ), &
                HS_2   => Difference % Value ( :, 2 ), &
                HS_3   => Difference % Value ( :, 3 ) )

    call CO % Initialize ( A % Communicator, [ 4 ], [ 4 ] )
    CO % Outgoing % Value ( 1 ) = sum ( abs ( HS_1 ), &
                                        mask = C % IsProperCell )
    CO % Outgoing % Value ( 2 ) = sum ( abs ( HS_2 ), &
                                       mask = C % IsProperCell )
    CO % Outgoing % Value ( 3 ) = sum ( abs ( HS_3 ), &
                                        mask = C % IsProperCell )
    CO % Outgoing % Value ( 4 ) = C % nProperCells
    call CO % Reduce ( REDUCTION % SUM )

    end associate !-- HS_1, etc.
    
    associate &
      ( HS_1    => CO % Incoming % Value ( 1 ), &
        HS_2    => CO % Incoming % Value ( 2 ), &
        HS_3    => CO % Incoming % Value ( 3 ), &
        nValues => CO % Incoming % Value ( 4 ) )

    L1_HS_1    = HS_1  / nValues
    L1_HS_2    = HS_2  / nValues
    L1_HS_3    = HS_3  / nValues

    end associate

    call Show ( L1_HS_1, '*** L1_HS_1 error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )
    call Show ( L1_HS_2, '*** L1_HS_2 error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )
    call Show ( L1_HS_3, '*** L1_HS_3 error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )

    end select !-- C
    end associate !-- A

    Difference % Value = abs ( Difference % Value / Reference % Value )

    nullify ( Solution, Reference, Difference )

  end subroutine ComputeError


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
