module LaplacianMultipole_ASC__Form_Test__Form

  use Basics
  use Manifolds
  use LaplacianMultipole_ASC__Form
  use SetHomogeneousSphere_Command

  implicit none
  private

  character ( LDF ), public, parameter :: &
    ProgramName = 'LaplacianMultipole_ASC__Form_Test'
    
  type, public :: LaplacianMultipole_ASC__Form_Test_Form
    type ( GridImageStreamForm ), allocatable :: &
      StreamAngular
    type ( Atlas_SC_CC_Form ), allocatable :: &
      Atlas
    type ( Atlas_SC_Form ), allocatable :: &
      AtlasAngular
    type ( LaplacianMultipole_ASC_Form ), allocatable :: &
      Laplacian
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type LaplacianMultipole_ASC__Form_Test_Form

contains


  subroutine Initialize ( LMFT, Name )

    class ( LaplacianMultipole_ASC__Form_Test_Form ) :: &
      LMFT
    character ( * ), intent ( in ) :: &
      Name

    integer ( KDI ) :: &
      iR, &  !-- iRank
      nEquations, &
      MaxDegree
    integer ( KDI ), dimension ( : ), allocatable :: &
      RankAngular
    real ( KDR ) :: &
      X_Random, &
      Cos_X, Sin_X, &
      Pi
    character ( LDF ) :: &
      NameAngular
    type ( CommunicatorForm ), allocatable :: &
      CommunicatorAngular
    class ( GeometryFlatForm ), pointer :: &
      G


    !-- Atlas

    allocate ( LMFT % Atlas )
    associate ( A => LMFT % Atlas )
    call A % Initialize ( Name, PROGRAM_HEADER % Communicator )
    call A % CreateChart_CC ( )
    call A % SetGeometry ( )

    G => A % Geometry ( )

    select type ( C => A % Chart )
    class is ( Chart_SLD_Form )


    !-- Laplacian

    nEquations = 1

    MaxDegree = 3
    call PROGRAM_HEADER % GetParameter ( MaxDegree, 'MaxDegree' )

    allocate ( LMFT % Laplacian )
    associate ( L => LMFT % Laplacian )
    call L % Initialize ( A, MaxDegree, nEquations )
    call L % InitializeTimers ( BaseLevel = 1 )

    call Show ( L % RadialEdges % Value ( :, 1 ), 'RadialEdge', &
                CONSOLE % INFO_2 )

    !-- Associated Legendre polynomials

    call Show ( 'Testing normalized Associated Legendre polynomials' )

    call random_number ( X_Random )
    X_Random  =  -1.0_KDR  +  2.0 * X_Random
       Cos_X  =  X_Random
       Sin_X  =  sqrt ( 1.0_KDR  -  X_Random ** 2 )
          Pi  =  CONSTANT % PI

    call Show ( X_Random, 'X_Random' )

    call Show ( L % AssociatedLegendre ( X_Random, 0, 0 ), 'P_0_0 computed' )
    call Show ( sqrt ( 1.0_KDR / ( 4.0_KDR * Pi ) ), 'P_0_0 expected' )

    call Show ( L % AssociatedLegendre ( X_Random, 1, 0 ), 'P_1_0 computed' )
    call Show ( sqrt ( 3.0_KDR / ( 4.0_KDR * Pi ) ) * Cos_X, &
                'P_1_0 expected' )

    call Show ( L % AssociatedLegendre ( X_Random, 1, 1 ), 'P_1_1 computed' )
    call Show ( - sqrt ( 3.0_KDR / ( 8.0_KDR * Pi ) ) * Sin_X, &
                'P_1_1 expected' )

    call Show ( L % AssociatedLegendre ( X_Random, 2, 0 ), 'P_2_0 computed' )
    call Show ( sqrt ( 5.0_KDR / ( 4.0_KDR * Pi ) ) &
                *  (    ( 3.0_KDR / 2.0_KDR )  *  Cos_X ** 2  &
                     -  ( 1.0_KDR / 2.0_KDR ) ), &
                'P_2_0 expected' )

    call Show ( L % AssociatedLegendre ( X_Random, 2, 1 ), 'P_2_1 computed' )
    call Show ( - sqrt ( 15.0_KDR / ( 8.0_KDR * Pi ) ) * Sin_X * Cos_X, &
                'P_2_1 expected' )

    call Show ( L % AssociatedLegendre ( X_Random, 2, 2 ), 'P_2_2 computed' )
    call Show ( ( 1.0_KDR / 4.0_KDR ) * sqrt ( 15.0_KDR / ( 2.0_KDR * Pi ) ) &
                *  Sin_X ** 2, &
                'P_2_2 expected' )


    !-- Angular functions

    NameAngular  =  trim ( Name ) // '_Angular'

    allocate ( RankAngular ( A % Communicator % Size  /  C % nBricks ( 1 ) ) )
    RankAngular  =  [ ( iR  *  C % nBricks ( 1 ), &
                        iR = 0, size ( RankAngular ) - 1 ) ]
 
    allocate ( CommunicatorAngular )
    call CommunicatorAngular % Initialize &
           ( A % Communicator, RankAngular, NameAngular ) 

    if ( any ( A % Communicator % Rank == RankAngular ) ) then

      allocate ( LMFT % AtlasAngular )
      associate ( AA => LMFT % AtlasAngular )
      call AA % Initialize &
             ( NameAngular, &
               CommunicatorOption = CommunicatorAngular, &
               nDimensionsOption = 3 )
      call AA % CreateChart &
             ( CoordinateLabelOption = [ 'r    ', 'Theta', 'Phi  ' ], &
               MinCoordinateOption = [ 0.0_KDR, 0.0_KDR, 0.0_KDR ], &
               MaxCoordinateOption = [ 1.0_KDR, Pi, 2.0_KDR * Pi ], &
               nCellsOption = [ 1, C % nCells ( 2 : 3 ) ], &
               nBricksOption = [ 1, C % nBricks ( 2 : 3 ) ], &
               nGhostLayersOption = [ 0, 0, 0 ] )  
      call AA % SetGeometry ( AddTimersOption = .false. )
    
      select type ( CA => AA % Chart )
      class is ( Chart_SLD_Form )

      allocate ( LMFT % StreamAngular )
      call LMFT % StreamAngular % Initialize &
             ( NameAngular, CommunicatorOption = CommunicatorAngular )    
      associate ( GIS => LMFT % StreamAngular )

      call AA % OpenStream ( GIS, 'Stream', iStream = 1 )

      call CA % AddFieldImage ( L % AngularFunctions, iStream = 1 )

      call GIS % Open ( GIS % ACCESS_CREATE )
      call AA % Write ( iStream = 1 )
      call GIS % Close ( )

      end associate !-- GIS
      end select !-- CA
      end associate !-- AA

    end if !-- Rank in RankAngular
    deallocate ( CommunicatorAngular )

    !-- Cleanup

    end associate !-- L
    end select !-- C
    end associate !-- A

    nullify ( G )!, Source, SolidHarmonics_RC, SolidHarmonics_IC, &
!              SolidHarmonics_RS, SolidHarmonics_IS, Solution )


  end subroutine Initialize


  subroutine Finalize ( LMFT )

    type ( LaplacianMultipole_ASC__Form_Test_Form ) :: &
      LMFT

   if ( allocated ( LMFT % Laplacian ) ) &
      deallocate ( LMFT % Laplacian )
    if ( allocated ( LMFT % AtlasAngular ) ) &
      deallocate ( LMFT % AtlasAngular )
    if ( allocated ( LMFT % Atlas ) ) &
      deallocate ( LMFT % Atlas )
    if ( allocated ( LMFT % StreamAngular ) ) &
      deallocate ( LMFT % StreamAngular )

  end subroutine Finalize


end module LaplacianMultipole_ASC__Form_Test__Form


program LaplacianMultipole_ASC__Form_Test

  use Basics
  use LaplacianMultipole_ASC__Form_Test__Form

  implicit none

  type ( LaplacianMultipole_ASC__Form_Test_Form ), allocatable :: &
    LMFT

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( ProgramName )
  call CONSOLE % SetVerbosity ( 'INFO_2' )

  allocate ( LMFT )
  call LMFT % Initialize ( PROGRAM_HEADER % Name )
  deallocate ( LMFT )

  deallocate ( PROGRAM_HEADER )

end program LaplacianMultipole_ASC__Form_Test
