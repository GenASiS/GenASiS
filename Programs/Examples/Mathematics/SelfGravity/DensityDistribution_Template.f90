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
    final :: &
      Finalize
  end type DensityDistributionTemplate

  private :: &
    ComputeError, &
    ShiftSolutionKernel
    
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
      iE, &  !-- iEquation
      MaxDegree
    real ( KDR ) :: &
      Density, &
      TotalMass, &
      RadiusMax
    real ( KDR ), dimension ( N_Eq ) :: &
      Eccentricity, &
      SemiMajor, &
      SemiMinor
    character ( LDL ) :: &
      PoissonSolverType
    class ( GeometryFlatForm ), pointer :: &
      G

    DD % N_Equations = N_eq
    DD % Variable = Variable
    
    RadiusMax = 10.0_KDR
    if ( present ( RadiusMaxOption ) ) RadiusMax = RadiusMaxOption

    !-- Atlas
    
    allocate ( DD % Atlas )
    associate ( A => DD % Atlas )
    call A % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )
    call A % CreateChart_CC ( RadiusMaxOption = RadiusMax )
    call A % SetGeometry ( UsePinnedMemoryOption = .true. )
    call A % Geometry_ASC % AllocateDevice ( )
    G => A % Geometry ( )
    call G % UpdateDevice ( )

    !-- Poisson

    MaxDegree = 10
    PoissonSolverType = 'MULTIPOLE_OLD'
    !-- FIXME: XL 16.1.1-5 does not work without association.
    associate ( PH => PROGRAM_HEADER )
!    call PROGRAM_HEADER % GetParameter ( MaxDegree, 'MaxDegree' )
    call PH % GetParameter ( MaxDegree, 'MaxDegree' )
    call PH % GetParameter ( PoissonSolverType, 'PoissonSolverType' )
    end associate !-- PH

    allocate ( DD % Poisson )
    associate ( P => DD % Poisson )
    call P % Initialize &
           ( A, SolverType = PoissonSolverType, MaxDegreeOption = MaxDegree, &
             nEquationsOption = DD % N_Equations )
    call P % InitializeTimers ( BaseLevel = 0 )
    end associate !-- P

    !-- Source, Reference

    allocate ( DD % Source )
    associate ( SA => DD % Source )
    call SA % Initialize &
           ( A, 'Source', DD % N_Equations, &
             VariableOption = DD % Variable, &
             WriteOption = .true., &
             UsePinnedMemoryOption = .true. )
    call SA % AllocateDevice ( AssociateVariablesOption = .false. )
    end associate !-- SA

    allocate ( DD % Reference )
    associate ( RA => DD % Reference )
    call RA % Initialize &
           ( A, 'Reference', DD % N_Equations, &
             VariableOption = DD % Variable, &
             WriteOption = .true. )
    end associate !-- RA

    !-- Solution, Difference

    allocate ( DD % Solution )
    associate ( SA => DD % Solution )
    call SA % Initialize &
           ( A, 'Solution', DD % N_Equations, &
             VariableOption = DD % Variable, &
             WriteOption = .true., &
             UsePinnedMemoryOption = .true. )
    call SA % AllocateDevice ( AssociateVariablesOption = .false. )
    end associate !-- SA

    allocate ( DD % Difference)
    associate ( DA => DD % Difference )
    call DA % Initialize &
           ( A, 'Difference', DD % N_Equations, &
             VariableOption = DD % Variable, &
             WriteOption = .true. )
    end associate !-- DA

    end associate !-- A
  
    nullify ( G )

  end subroutine Initialize 


  subroutine Compute ( DD, ShiftSolutionOption )

    class ( DensityDistributionTemplate ), intent ( inout ) :: &
      DD
    logical ( KDL ), optional :: &
      ShiftSolutionOption

    integer ( KDI ) :: &
      iS, &  !-- iSolve
      nSolve
    logical ( KDL ) :: &
      ShiftSolution, &
      NormalizeSolution

    ShiftSolution = .false.
    NormalizeSolution = .false. 

    associate ( A => DD % Atlas )

    associate ( P => DD % Poisson )

    nSolve = 1
    !-- FIXME: XL 16.1.1-5 does not work without association.
    associate ( PH => PROGRAM_HEADER )
!    call PROGRAM_HEADER % GetParameter ( nSolve, 'nSolve' )
    call PH % GetParameter ( nSolve, 'nSolve' )
    end associate !-- PH

    call Show ( 'Solving Poisson equation' )
    call Show ( nSolve, 'nSolve' )
    do iS = 1, nSolve
      call Show ( iS, 'iS' )
      call P % Solve ( DD % Solution, DD % Source )
    end do !-- iS

    if ( present ( ShiftSolutionOption )  ) &
      ShiftSolution = ShiftSolutionOption
    if ( ShiftSolution ) call ShiftSolutionKernel ( DD )

    call Show ( 'Computing error' )
    call ComputeError ( DD, NormalizeSolution ) 

    call PROGRAM_HEADER % ShowStatistics &
           ( CONSOLE % INFO_1, &
             CommunicatorOption = PROGRAM_HEADER % Communicator )

    !-- Write

    call Show ( 'Writing results' )
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


  subroutine ComputeError ( DD, Normalized )
    class ( DensityDistributionTemplate ), intent ( inout ) :: &
      DD
    logical ( KDL ), intent ( in ) :: &
      Normalized

    integer ( KDI ) :: &
      nM, &
      iV, &
      iE, &
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

    if ( .not. Normalized ) &
      call MultiplyAdd &
             ( Solution % Value, Reference % Value, &
                 -1.0_KDR, Difference % Value )

    associate ( A => DD % Atlas ) 
    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
    
    associate ( D   => Difference % Value, &
                R   => Reference  % Value )

    nM = DD % N_Equations

    call CO % Initialize ( A % Communicator, [ 2 * nM ], [ 2 * nM ] )

    do iM = 1, nM
      CO % Outgoing % Value ( 2 * iM ) &
        = sum ( abs ( R ( :, iM ) ), mask = C % IsProperCell )
      CO % Outgoing % Value ( 2 * iM - 1 )&
        = sum ( abs ( D ( :, iM ) ), mask = C % IsProperCell )
    end do

    call CO % Reduce ( REDUCTION % SUM )

    end associate !-- D
    
    !-- L1 error

    associate ( IN   => CO % Incoming % Value )

    do iM = 1, nM
      L1 = IN ( 2 * iM - 1 ) / IN ( 2 * iM )
      call Show ( L1, '*** L1 error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )
    end do

    end associate !-- IN 

    end select !-- C
    end associate !-- A

    !-- Relative difference

    do iE = 1, nM
      do iV = 1, size ( Difference % Value ( :, 1 ) )
      Difference % Value ( iV, iE) &
        = abs ( Difference % Value ( iV, iE) ) &
                  / max ( abs ( Reference % Value ( iV, iE ) ), &
                          sqrt ( tiny ( 0.0_KDR ) ) )
      end do
    end do

    nullify ( Solution, Reference, Difference )

  end subroutine ComputeError


  subroutine ShiftSolutionKernel ( DD )
    type ( DensityDistributionTemplate ), intent ( inout ) :: &
      DD

    integer ( KDI ) :: &
      iE
    real ( KDR ) :: &
      D
    class ( StorageForm ), pointer :: &
      Solution, &
      Reference
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
      D = CO % Incoming % Value ( 2 * iE - 1 ) &
                 - CO % Incoming % Value ( 2 * iE )
      call Show ( D, '>>> C_shift=' )
      Reference % Value ( :, iE ) = Reference % Value ( :, iE ) + D
    end do

    end select !-- C
    end associate !-- A

    nullify ( Solution, Reference )

  end subroutine ShiftSolutionKernel


end module DensityDistribution_Template
