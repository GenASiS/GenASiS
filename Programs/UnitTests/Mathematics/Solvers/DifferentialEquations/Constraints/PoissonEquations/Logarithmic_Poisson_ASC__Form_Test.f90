module Logarithmic_Poisson_ASC__Form_Test__Form

  use Basics
  use Manifolds
  use Operations
  use Poisson_ASC__Form
  use SetLogarithmic_Command

  implicit none
  private

  integer ( KDI ), private, parameter :: &
    O_L                   = 0, &  !-- OFFSET_LOGARITHMIC
    N_LOGARITHMIC         = 3, &
    N_EQUATIONS           = N_LOGARITHMIC
  character ( LDL ), dimension ( N_EQUATIONS ), private, parameter :: &
    VARIABLE = [ 'Logarithmic_1        ', &
                 'Logarithmic_2        ', &
                 'Logarithmic_3        ' ]
  character ( LDF ), public, parameter :: &
    ProgramName = 'Logarithmic_Poisson_ASC__Form_Test'

  type, public :: Logarithmic_Poisson_ASC__Form_Test_Form
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
    procedure, public, pass :: &
      ShiftSolution
    final :: &
      Finalize
  end type Logarithmic_Poisson_ASC__Form_Test_Form

contains


  subroutine Initialize ( LPFT, Name )

    class ( Logarithmic_Poisson_ASC__Form_Test_Form ) :: &
      LPFT
    character ( * ), intent ( in ) :: &
      Name

    integer ( KDI ) :: &
      iE, &  !-- iEquation
      MaxDegree
    real ( KDR ), dimension ( N_LOGARITHMIC ) :: &
      q
    real ( KDR ), dimension ( N_EQUATIONS ) :: &
      Density
    real ( KDR ) :: &
      Total_Mass, &
      v0, &
      rho_c

    !-- Atlas

    allocate ( LPFT % Atlas )
    associate ( A => LPFT % Atlas )
    call A % Initialize &
           ( Name, PROGRAM_HEADER % Communicator )
    call A % CreateChart_CC ( RadiusMaxOption = 1.0_KDR )
    call A % SetGeometry ( )


    !-- Poisson

    MaxDegree = 2
    call PROGRAM_HEADER % GetParameter ( MaxDegree, 'MaxDegree' )

    allocate ( LPFT % Poisson )
    associate ( P => LPFT % Poisson )
    call P % Initialize &
           ( A, SolverType = 'MULTIPOLE', MaxDegreeOption = MaxDegree, &
             nEquationsOption = N_EQUATIONS )


    !-- Source, Reference

    allocate ( LPFT % Source )
    associate ( SA => LPFT % Source )
    call SA % Initialize &
           ( A, 'Source', N_EQUATIONS, &
             VariableOption = VARIABLE, &
             WriteOption = .true. )

    allocate ( LPFT % Reference )
    associate ( RA => LPFT % Reference )
    call RA % Initialize &
           ( A, 'Reference', N_EQUATIONS, &
             VariableOption = VARIABLE, &
             WriteOption = .true. )

   !-- Logarithmic Potential 

    v0 = 1.0_KDR
    call PROGRAM_HEADER % GetParameter ( v0, 'v0' )

    rho_c = 0.1_KDR
    call PROGRAM_HEADER % GetParameter ( rho_c, 'rho_c' )

    q = [ 0.52_KDR, 0.8_KDR, 0.95_KDR ]
    call PROGRAM_HEADER % GetParameter ( q, 'q' )

    do iE = O_L + 1, O_L + N_LOGARITHMIC
      call SetLogarithmic &
             ( SA, RA, A, v0, rho_c, q ( iE - O_L ) , iE )
    end do


    !-- Solution, Difference

    allocate ( LPFT % Solution )
    associate ( SA => LPFT % Solution )
    call SA % Initialize &
           ( A, 'Solution', N_EQUATIONS, &
             VariableOption = VARIABLE, &
             WriteOption = .true. )
    end associate !-- SA

    allocate ( LPFT % Difference)
    associate ( DA => LPFT % Difference )
    call DA % Initialize &
           ( A, 'Difference', N_EQUATIONS, &
             VariableOption = VARIABLE, &
             WriteOption = .true. )
    end associate !-- DA

    call P % Solve ( LPFT % Solution, LPFT % Source )

    call LPFT % ShiftSolution ( )

    call LPFT % ComputeError ( ) 


    !-- Write

    allocate ( LPFT % Stream )
    call LPFT % Stream % Initialize &
           ( Name, CommunicatorOption = PROGRAM_HEADER % Communicator )    
    associate ( GIS => LPFT % Stream )

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

  
  subroutine ShiftSolution ( LPFT )

    class ( Logarithmic_Poisson_ASC__Form_Test_Form ), intent ( in ) :: &
      LPFT

    real ( KDR ) :: &
      D_1, &
      D_2, &
      D_3
    class ( StorageForm ), pointer :: &
      Solution, &
      Reference, &
      Difference
    type ( CollectiveOperation_R_Form ) :: &
      CO
    
    Solution   => LPFT % Solution   % Storage ( )
    Reference  => LPFT % Reference  % Storage ( )

    associate ( A => LPFT % Atlas ) 
    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
    
    associate ( S_1   => Solution % Value ( :, O_L + 1 ), &
                S_2   => Solution % Value ( :, O_L + 2 ), &
                S_3   => Solution % Value ( :, O_L + 3 ), &
                R_1   => Reference % Value ( :, O_L + 1 ), &
                R_2   => Reference % Value ( :, O_L + 2 ), &
                R_3   => Reference % Value ( :, O_L + 3 ) )

    call CO % Initialize ( A % Communicator, [ 6 ], [ 6 ] )
      CO % Outgoing % Value ( 1 ) = minval ( S_1, mask = C % IsProperCell )
      CO % Outgoing % Value ( 2 ) = minval ( S_2, mask = C % IsProperCell )
      CO % Outgoing % Value ( 3 ) = minval ( S_3, mask = C % IsProperCell )
      CO % Outgoing % Value ( 4 ) = minval ( R_1, mask = C % IsProperCell )
      CO % Outgoing % Value ( 5 ) = minval ( R_2, mask = C % IsProperCell )
      CO % Outgoing % Value ( 6 ) = minval ( R_3, mask = C % IsProperCell )
    call CO % Reduce ( REDUCTION % MIN )

    end associate !-- S_1, etc.
    
    associate &
      ( S_1 => CO % Incoming % Value ( 1 ), &
        S_2 => CO % Incoming % Value ( 2 ), &
        S_3 => CO % Incoming % Value ( 3 ), &
        R_1 => CO % Incoming % Value ( 4 ), &
        R_2 => CO % Incoming % Value ( 5 ), &
        R_3 => CO % Incoming % Value ( 6 ) )

    D_1   = abs ( S_1 - R_1 )
    D_2   = abs ( S_2 - R_2 )
    D_3   = abs ( S_3 - R_3 )

    end associate

    Solution % Value ( :, O_L + 1 ) = Solution % Value ( :, O_L + 1 ) + D_1
    Solution % Value ( :, O_L + 2 ) = Solution % Value ( :, O_L + 2 ) + D_2
    Solution % Value ( :, O_L + 3 ) = Solution % Value ( :, O_L + 3 ) + D_3
   
    end select !-- C
    end associate !-- A

    nullify ( Solution, Reference )

  end subroutine ShiftSolution


  subroutine ComputeError ( LPFT )

    class ( Logarithmic_Poisson_ASC__Form_Test_Form ), intent ( in ) :: &
      LPFT
    
    real ( KDR ) :: &
      L1_L_1, &
      L1_L_2, &
      L1_L_3
    class ( StorageForm ), pointer :: &
      Solution, &
      Reference, &
      Difference         
    type ( CollectiveOperation_R_Form ) :: &
      CO
    
    Solution   => LPFT % Solution   % Storage ( )
    Reference  => LPFT % Reference  % Storage ( )
    Difference => LPFT % Difference % Storage ( )

    call MultiplyAdd &
           ( Solution % Value, Reference % Value, -1.0_KDR, Difference % Value )

    associate ( A => LPFT % Atlas ) 
    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
    
    associate ( L_1   => Difference % Value ( :, O_L + 1 ), &
                L_2   => Difference % Value ( :, O_L + 2 ), &
                L_3   => Difference % Value ( :, O_L + 3 ) )

    call CO % Initialize ( A % Communicator, [ 4 ], [ 4 ] )
    CO % Outgoing % Value ( 1 ) = sum ( abs ( L_1 ), &
                                        mask = C % IsProperCell )
    CO % Outgoing % Value ( 2 ) = sum ( abs ( L_2 ), &
                                       mask = C % IsProperCell )
    CO % Outgoing % Value ( 3 ) = sum ( abs ( L_3 ), &
                                        mask = C % IsProperCell )
    CO % Outgoing % Value ( 4 ) = C % nProperCells
    call CO % Reduce ( REDUCTION % SUM )

    end associate !-- DS_1, etc.
    
    associate &
      ( L_1   => CO % Incoming % Value ( 1 ), &
        L_2   => CO % Incoming % Value ( 2 ), &
        L_3   => CO % Incoming % Value ( 3 ), &
        nValues => CO % Incoming % Value ( 4 ) )

    L1_L_1   = L_1   / nValues
    L1_L_2   = L_2   / nValues
    L1_L_3   = L_3   / nValues

    end associate

    call Show ( L1_L_1, '*** L1_L_1 error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )
    call Show ( L1_L_2, '*** L1_L_2 error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )
    call Show ( L1_L_3, '*** L1_L_3 error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )

    end select !-- C
    end associate !-- A

    Difference % Value = abs ( Difference % Value / Reference % Value )

    nullify ( Solution, Reference, Difference )

  end subroutine ComputeError


  subroutine Finalize ( LPFT )

    type ( Logarithmic_Poisson_ASC__Form_Test_Form ) :: &
      LPFT

    if ( allocated ( LPFT % Stream ) ) &
      deallocate ( LPFT % Stream )
    if ( allocated ( LPFT % Reference ) ) &
      deallocate ( LPFT % Reference )
    if ( allocated ( LPFT % Solution ) ) &
      deallocate ( LPFT % Solution )
    if ( allocated ( LPFT % Source ) ) &
      deallocate ( LPFT % Source )
    if ( allocated ( LPFT % Poisson ) ) &
      deallocate ( LPFT % Poisson )
    if ( allocated ( LPFT % Atlas ) ) &
      deallocate ( LPFT % Atlas )

  end subroutine Finalize


end module Logarithmic_Poisson_ASC__Form_Test__Form



program Logarithmic_Poisson_ASC__Form_Test

  use Basics
  use Logarithmic_Poisson_ASC__Form_Test__Form

  implicit none

  type ( Logarithmic_Poisson_ASC__Form_Test_Form ), allocatable :: &
    LPFT

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( ProgramName )

  allocate ( LPFT )
  call LPFT % Initialize ( PROGRAM_HEADER % Name )
  deallocate ( LPFT )

  deallocate ( PROGRAM_HEADER )

end program Logarithmic_Poisson_ASC__Form_Test
