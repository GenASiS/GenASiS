module Poisson_ASC__Form_Test__Form

  use Basics
  use Manifolds
  use Operations
  use Poisson_ASC__Form
  use SetHomogeneousSphere_Command
  use SetHomogeneousSpheroid_Command
  use SetKuzminPlummerDisk_Command

  implicit none
  private

  integer ( KDI ), private, parameter :: &
    O_HS                   = 0, &  !-- OFFSET_HOMOGENEOUS_SPHERE
    N_HOMOGENEOUS_SPHERE   = 3, &
    O_S                    = O_HS + N_HOMOGENEOUS_SPHERE, & !-- OFFSET_SPHEROID
    N_HOMOGENEOUS_SPHEROID = 3, &
    O_KPD                  = O_S + N_HOMOGENEOUS_SPHEROID, &
    N_KUZMIN_PLUMMER_DISK  = 4, &
    N_EQUATIONS            &
      = N_HOMOGENEOUS_SPHERE + N_HOMOGENEOUS_SPHEROID + N_KUZMIN_PLUMMER_DISK
  character ( LDL ), dimension ( N_EQUATIONS ), private, parameter :: &
    VARIABLE = [ 'HomogeneousSphere_1  ', &
                 'HomogeneousSphere_2  ', &
                 'HomogeneousSphere_3  ', &
                 'HomogeneousSpheroid_1', &
                 'HomogeneousSpheroid_2', &
                 'HomogeneousSpheroid_3', &
                 'KuzminPlummerDisk_1  ', &
                 'KuzminPlummerDisk_2  ', &
                 'KuzminPlummerDisk_3  ', &
                 'KuzminPlummerDisk_4  ']
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
    real ( KDR ), dimension ( N_HOMOGENEOUS_SPHERE ) :: &
      RadiusDensity
    real ( KDR ), dimension ( N_HOMOGENEOUS_SPHEROID ) :: &
      Eccentricity, &
      SemiMajor, &
      SemiMinor
    real ( KDR ), dimension ( N_KUZMIN_PLUMMER_DISK ) :: &
      KPD_ratio, &
      KPD_a, &
      KPD_b
    real ( KDR ), dimension ( N_EQUATIONS ) :: &
      Density
    real ( KDR ) :: &
      Total_Mass

    !-- Atlas

    allocate ( PFT % Atlas )
    associate ( A => PFT % Atlas )
    call A % Initialize ( Name, PROGRAM_HEADER % Communicator )
    call A % CreateChart_CC ( )
    call A % SetGeometry ( )


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

    Density ( O_HS + 1 : O_HS + N_HOMOGENEOUS_SPHERE ) &
      = 1.0_KDR / ( 4.0_KDR  *  CONSTANT % PI  *  RadiusDensity ** 3 )
    call PROGRAM_HEADER % GetParameter ( Density, 'Density' )

    do iE = O_HS + 1, O_HS + N_HOMOGENEOUS_SPHERE
      call SetHomogeneousSphere &
             ( SA, RA, A, Density ( iE ), RadiusDensity ( iE ), iE )
    end do

    !-- Homogeneous spheroid

    SemiMajor = A % Chart % MaxCoordinate ( 1 ) * 0.5_KDR
    call PROGRAM_HEADER % GetParameter ( SemiMajor, 'SemiMajor' )

    Eccentricity = sqrt ( 1.0_KDR &
                         - [ 0.2_KDR ** 2, 0.7_KDR ** 2, 0.2_KDR ** 2 ] )
    call PROGRAM_HEADER % GetParameter ( Eccentricity, 'Eccentricity' )

    SemiMinor = sqrt ( 1.0_KDR - Eccentricity ** 2 ) * SemiMajor

    !-- Prolate Spheroid 
    SemiMajor ( N_HOMOGENEOUS_SPHEROID ) &
      = SemiMinor ( N_HOMOGENEOUS_SPHEROID )
    SemiMinor ( N_HOMOGENEOUS_SPHEROID ) &
      = A % Chart % MaxCoordinate ( 1 ) * 0.5_KDR

    Density = 1.0_KDR / ( 4.0_KDR  *  CONSTANT % PI  )
    call PROGRAM_HEADER % GetParameter ( Density, 'Density' )

    do iE = O_S + 1, O_S + N_HOMOGENEOUS_SPHEROID
      call SetHomogeneousSpheroid &
             ( SA, RA, A, Density ( iE ), &
               SemiMajor ( iE - O_S ), SemiMinor ( iE - O_S ), iE )
    end do
    
    !-- Kuzmin-Plummer Disk

    KPD_a = 1.0_KDR
    call PROGRAM_HEADER % GetParameter ( KPD_a, 'KPD_a' )
    KPD_a ( N_KUZMIN_PLUMMER_DISK ) = 0.4_KDR

    KPD_ratio = [ 0.2_KDR, 1.0_KDR, 5.0_KDR, 0.25_KDR ]
    call PROGRAM_HEADER % GetParameter ( KPD_ratio, 'KPD_ratio' )

    KPD_b = KPD_ratio * KPD_a
    KPD_b ( N_KUZMIN_PLUMMER_DISK ) = 0.1_KDR

    Total_Mass = 1.0_KDR 
    call PROGRAM_HEADER % GetParameter ( Total_Mass, 'TotalMass' )

    do iE = O_KPD + 1, O_KPD + N_KUZMIN_PLUMMER_DISK
      call SetKuzminPlummerDisk &
             ( SA, RA, A, Total_Mass, &
               KPD_a ( iE - O_KPD ), KPD_b ( iE - O_KPD ) , iE )
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

    end associate !-- RA
    end associate !-- SA
    end associate !-- P
    end associate !-- A

  end subroutine Initialize


  subroutine ComputeError ( PFT )

    class ( Poisson_ASC__Form_Test_Form ), intent ( in ) :: &
      PFT
    
    real ( KDR ) :: &
      L1_DS_1, &
      L1_DS_2, &
      L1_DS_3, &
      L1_KPD_1, &
      L1_KPD_2, &
      L1_KPD_3, &
      L1_KPD_4
    class ( VariableGroupForm ), pointer :: &
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
    
    associate ( DS_1   => Difference % Value ( :, O_S + 1 ), &
                DS_2   => Difference % Value ( :, O_S + 2 ), &
                DS_3   => Difference % Value ( :, O_S + 3 ), &
                KPD_1  => Difference % Value ( :, O_KPD + 1 ), &
                KPD_2  => Difference % Value ( :, O_KPD + 2 ), &
                KPD_3  => Difference % Value ( :, O_KPD + 3 ), &
                KPD_4  => Difference % Value ( :, O_KPD + 4 ) )

    call CO % Initialize ( A % Communicator, [ 8 ], [ 8 ] )
    CO % Outgoing % Value ( 1 ) = sum ( abs ( DS_1 ), &
                                        mask = C % IsProperCell )
    CO % Outgoing % Value ( 2 ) = sum ( abs ( DS_2 ), &
                                       mask = C % IsProperCell )
    CO % Outgoing % Value ( 3 ) = sum ( abs ( DS_3 ), &
                                        mask = C % IsProperCell )
    CO % Outgoing % Value ( 4 ) = sum ( abs ( KPD_1 ), &
                                        mask = C % IsProperCell )
    CO % Outgoing % Value ( 5 ) = sum ( abs ( KPD_2 ), &
                                        mask = C % IsProperCell )
    CO % Outgoing % Value ( 6 ) = sum ( abs ( KPD_3 ), &
                                        mask = C % IsProperCell )
    CO % Outgoing % Value ( 7 ) = sum ( abs ( KPD_4 ), &
                                        mask = C % IsProperCell )
    CO % Outgoing % Value ( 8 ) = C % nProperCells
    call CO % Reduce ( REDUCTION % SUM )

    end associate !-- DS_1, etc.
    
    associate &
      ( DS_1   => CO % Incoming % Value ( 1 ), &
        DS_2   => CO % Incoming % Value ( 2 ), &
        DS_3   => CO % Incoming % Value ( 3 ), &
        KPD_1  => CO % Incoming % Value ( 4 ), &
        KPD_2  => CO % Incoming % Value ( 5 ), &
        KPD_3  => CO % Incoming % Value ( 6 ), &
        KPD_4  => CO % Incoming % Value ( 7 ), &
        nValues => CO % Incoming % Value ( 8 ) )

    L1_DS_1   = DS_1  / nValues
    L1_DS_2   = DS_2  / nValues
    L1_DS_3   = DS_3  / nValues
    L1_KPD_1   = KPD_1  / nValues
    L1_KPD_2   = KPD_2  / nValues
    L1_KPD_3   = KPD_3  / nValues
    L1_KPD_4   = KPD_4  / nValues

    end associate

    call Show ( L1_DS_1, '*** L1_DS_1 error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )
    call Show ( L1_DS_2, '*** L1_DS_2 error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )
    call Show ( L1_DS_3, '*** L1_DS_3 error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )
    call Show ( L1_KPD_1, '*** L1_KPD_1 error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )
    call Show ( L1_KPD_2, '*** L1_KPD_2 error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )
    call Show ( L1_KPD_3, '*** L1_KPD_3 error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )
    call Show ( L1_KPD_4, '*** L1_KPD_4 error', nLeadingLinesOption = 2, &
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
