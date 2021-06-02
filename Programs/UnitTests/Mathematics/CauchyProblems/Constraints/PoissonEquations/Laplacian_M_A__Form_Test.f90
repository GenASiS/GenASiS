program Laplacian_M_A__Form_Test

  !-- Laplacian_Multipole_Atlas__Form_Test

  use Basics
  use Manifolds
  use Fields
  use PoissonEquations

  implicit none

  integer ( KDI ) :: &
!    iR, &  !-- iRank
    nEquations, &
    MaxDegree
  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Atlas_SCG_CC_Form ), allocatable :: &
    A
  type ( Stream_A_Form ), allocatable :: &
    SA
  type ( Geometry_F_A_Form ), allocatable :: &
    GA
  type ( Laplacian_M_A_Form ), allocatable :: &
    LA

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Laplacian_M_A__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( A )
  call A % Initialize &
         ( RadiusMax = 10.0_KDR, &
           RadiusCore = 10.0_KDR / 8.0_KDR, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( SA )
  call SA % Initialize ( A, GIS )

  allocate ( GA )
  call GA % Initialize ( A )
  call GA % SetStream ( SA )

  nEquations = 1

  MaxDegree = 3
  call PROGRAM_HEADER % GetParameter ( MaxDegree, 'MaxDegree' )

  allocate ( LA )
  call LA % Initialize ( GA, MaxDegree, nEquations )

  call  A % Show ( )
  call SA % Show ( )
  call GA % Show ( )
  call LA % Show ( )

  call TestAssociatedLegendre ( LA )

  deallocate ( LA )
  deallocate ( GA )
  deallocate ( SA )
  deallocate ( A )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )


contains


  subroutine TestAssociatedLegendre ( LA )

    class ( Laplacian_M_A_Form ), intent ( inout ) :: &
      LA

    real ( KDR ) :: &
      X_Random, &
      Cos_X, Sin_X, &
      Pi

    call Show ( 'Testing normalized Associated Legendre polynomials' )

    call InitializeRandomSeed ( PROGRAM_HEADER % Communicator )
    call random_number ( X_Random )
    X_Random  =  -1.0_KDR  +  2.0 * X_Random
       Cos_X  =  X_Random
       Sin_X  =  sqrt ( 1.0_KDR  -  X_Random ** 2 )
          Pi  =  CONSTANT % PI

    call Show ( X_Random, 'X_Random' )

    call Show ( LA % AssociatedLegendre ( X_Random, 0, 0 ), 'P_0_0 computed' )
    call Show ( sqrt ( 1.0_KDR / ( 4.0_KDR * Pi ) ), 'P_0_0 expected' )

    call Show ( LA % AssociatedLegendre ( X_Random, 1, 0 ), 'P_1_0 computed' )
    call Show ( sqrt ( 3.0_KDR / ( 4.0_KDR * Pi ) ) * Cos_X, &
                'P_1_0 expected' )

    call Show ( LA % AssociatedLegendre ( X_Random, 1, 1 ), 'P_1_1 computed' )
    call Show ( - sqrt ( 3.0_KDR / ( 8.0_KDR * Pi ) ) * Sin_X, &
                'P_1_1 expected' )

    call Show ( LA % AssociatedLegendre ( X_Random, 2, 0 ), 'P_2_0 computed' )
    call Show ( sqrt ( 5.0_KDR / ( 4.0_KDR * Pi ) ) &
                *  (    ( 3.0_KDR / 2.0_KDR )  *  Cos_X ** 2  &
                     -  ( 1.0_KDR / 2.0_KDR ) ), &
                'P_2_0 expected' )

    call Show ( LA % AssociatedLegendre ( X_Random, 2, 1 ), 'P_2_1 computed' )
    call Show ( - sqrt ( 15.0_KDR / ( 8.0_KDR * Pi ) ) * Sin_X * Cos_X, &
                'P_2_1 expected' )

    call Show ( LA % AssociatedLegendre ( X_Random, 2, 2 ), 'P_2_2 computed' )
    call Show ( ( 1.0_KDR / 4.0_KDR ) * sqrt ( 15.0_KDR / ( 2.0_KDR * Pi ) ) &
                *  Sin_X ** 2, &
                'P_2_2 expected' )

  end subroutine TestAssociatedLegendre


end program Laplacian_M_A__Form_Test
