module DensityDistribution_Template

  use Basics
  use Mathematics
  
  implicit none

  type, public :: DensityDistributionTemplate
    integer ( KDI ) :: &
      N_Equations
    character ( LDL ), dimension ( : ), allocatable :: &
      Variable
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
      Compute
    procedure, public, pass :: &
      ComputeError
    final :: &
      Finalize
  end type DensityDistributionTemplate

  private :: &
    ComputeError, &
    ShiftSolution

contains


  subroutine Initialize ( DD, Name, N_Eq, Variable, RadiusMaxOption )
    class ( DensityDistributionTemplate ) :: &
      DD
    character ( * ), intent ( in ) :: &
      Name
    integer ( KDI ), intent ( in ) :: &
      N_Eq
    character ( LDL ), dimension ( N_Eq ), intent ( in ) :: &
      Variable
    real ( KDR ), optional :: &
      RadiusMaxOption
    
    integer ( KDI ) :: &
      iE, & !iEquation
      MaxDegree
    real ( KDR ) :: &
      Density, &
      TotalMass, &
      RadiusMax
    real ( KDR ), dimension ( N_Eq ) :: &
      Eccentricity, &
      SemiMajor, &
      SemiMinor

    DD % N_Equations = N_eq
    DD % Variable = Variable
    
    RadiusMax = 10.0_KDR
    if ( present ( RadiusMaxOption ) ) RadiusMax = RadiusMaxOption

    !-- Atlas
    
    allocate ( DD % Atlas )
    associate ( A => DD % Atlas )
    call A % Initialize ( Name, PROGRAM_HEADER % Communicator )
    call A % CreateChart_CC ( RadiusMaxOption = RadiusMax )
    call A % SetGeometry ( )

    !-- Poisson
    MaxDegree = 10
    call PROGRAM_HEADER % GetParameter ( MaxDegree, 'MaxDegree' )

    allocate ( DD % Poisson )
    associate ( P => DD % Poisson )
    call P % Initialize &
           ( A, SolverType = 'MULTIPOLE', MaxDegreeOption = MaxDegree, &
             nEquationsOption = DD % N_Equations )
    end associate !-- P

    !-- Source, Reference

    allocate ( DD % Source )
    associate ( SA => DD % Source )
    call SA % Initialize &
           ( A, 'Source', DD % N_Equations, &
             VariableOption = DD % Variable, &
             WriteOption = .true. )

    allocate ( DD % Reference )
    associate ( RA => DD % Reference )
    call RA % Initialize &
           ( A, 'Reference', DD % N_Equations, &
             VariableOption = DD % Variable, &
             WriteOption = .true. )

    end associate !-- SA
    end associate !-- RA

    !-- Solution, Difference

    allocate ( DD % Solution )
    associate ( SA => DD % Solution )
    call SA % Initialize &
           ( A, 'Solution', DD % N_Equations, &
             VariableOption = DD % Variable, &
             WriteOption = .true. )
    end associate !-- SA

    allocate ( DD % Difference)
    associate ( DA => DD % Difference )
    call DA % Initialize &
           ( A, 'Difference', DD % N_Equations, &
             VariableOption = DD % Variable, &
             WriteOption = .true. )
    end associate !-- DA

    end associate !-- A
  
  end subroutine Initialize 


  subroutine Compute ( DD, ShiftSolutionOption )
    class ( DensityDistributionTemplate ), intent ( inout ) :: &
      DD
    logical ( KDL ), optional :: &
      ShiftSolutionOption

    associate ( A => DD % Atlas )

    associate ( P => DD % Poisson )

    call P % Solve ( DD % Solution, DD % Source )

    if ( present ( ShiftSolutionOption ) ) call ShiftSolution ( DD )
    
    call ComputeError ( DD ) 

    !-- Write

    allocate ( DD % Stream )
    call DD % Stream % Initialize &
           ( A % Name, CommunicatorOption = PROGRAM_HEADER % Communicator )    
    associate ( GIS => DD % Stream )

    call A % OpenStream ( GIS, 'Stream', iStream = 1 )

    call GIS % Open ( GIS % ACCESS_CREATE )
    call A % Write ( iStream = 1 )
    call GIS % Close ( )

    end associate !-- GIS

    !-- Cleanup

    end associate !-- P
    end associate !-- A

  end subroutine Compute


  subroutine Finalize ( DD )
    type ( DensityDistributionTemplate ) :: &
      DD

    if ( allocated ( DD % Stream ) ) &
      deallocate ( DD % Stream )
    if ( allocated ( DD % Reference ) ) &
      deallocate ( DD % Reference )
    if ( allocated ( DD % Solution ) ) &
      deallocate ( DD % Solution )
    if ( allocated ( DD % Source ) ) &
      deallocate ( DD % Source )
    if ( allocated ( DD % Poisson ) ) &
      deallocate ( DD % Poisson )
    if ( allocated ( DD % Atlas ) ) &
      deallocate ( DD % Atlas )

  end subroutine Finalize


  subroutine ShiftSolution ( DD )
    type ( DensityDistributionTemplate ), intent ( inout ) :: &
      DD

    integer ( KDI ) :: &
      iE
    real ( KDR ) :: &
      D
    class ( StorageForm ), pointer :: &
      Solution, &
      Reference, &
      Difference
    type ( CollectiveOperation_R_Form ) :: &
      CO
    
    Solution   => DD % Solution   % Storage ( )
    Reference  => DD % Reference  % Storage ( )

    associate ( A => DD % Atlas ) 
    select type ( C => A % Chart )
    class is ( Chart_SL_Template )

    call CO % Initialize &
                ( A % Communicator, [ DD % N_Equations * 2 ], &
                    [ DD % N_Equations * 2 ] )
    do iE = 1, DD % N_Equations
      CO % Outgoing % Value ( 2 * iE - 1 ) &
        = minval ( Solution % Value ( :, iE ), mask = C % IsProperCell )
      CO % Outgoing % Value ( 2 * iE ) &
        = minval ( Reference % Value ( :, iE ), mask = C % IsProperCell )
    end do
    call CO % Reduce ( REDUCTION % MIN )

    do iE = 1, DD % N_Equations
      D = abs ( CO % Incoming % Value ( 2 * iE - 1 ) &
                 - CO % Incoming % Value ( 2 * iE ) )
      Solution % Value ( :, iE ) = Solution % Value ( :, iE ) + D
    end do

    end select !-- C
    end associate !-- A

    nullify ( Solution, Reference )

  end subroutine ShiftSolution


  subroutine ComputeError ( DD )
    class ( DensityDistributionTemplate ), intent ( in ) :: &
      DD

    integer ( KDI ) :: &
      nM, & 
      iM !-- iMessage
    real ( KDR ) :: &
      L1
    class ( StorageForm ), pointer :: &
      Solution, &
      Reference, &
      Difference         
    type ( CollectiveOperation_R_Form ) :: &
      CO
    
    Solution   => DD % Solution   % Storage ( )
    Reference  => DD % Reference  % Storage ( )
    Difference => DD % Difference % Storage ( )

    call MultiplyAdd &
           ( Solution % Value, Reference % Value, -1.0_KDR, Difference % Value )

    associate ( A => DD % Atlas ) 
    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
    
    associate ( D   => Difference % Value )

    nM = DD % N_Equations + 1

    call CO % Initialize ( A % Communicator, [ nM ], [ nM ] )

    do iM = 1, nM - 1
      CO % Outgoing % Value ( iM ) = sum ( abs ( D ( :, iM ) ), &
                                            mask = C % IsProperCell )
    end do

    CO % Outgoing % Value ( nM ) = C % nProperCells

    call CO % Reduce ( REDUCTION % SUM )

    end associate !-- D
    
    !-- L1 error

    associate ( IN   => CO % Incoming % Value )

    do iM = 1, nM - 1
      L1 = IN ( iM ) / IN ( nM )
      call Show ( L1, '*** L1 error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )
    end do

    end associate !-- IN 

    end select !-- C
    end associate !-- A

    !-- Relative difference

    Difference % Value = abs ( Difference % Value / Reference % Value )

    nullify ( Solution, Reference, Difference )

  end subroutine ComputeError  


end module DensityDistribution_Template
