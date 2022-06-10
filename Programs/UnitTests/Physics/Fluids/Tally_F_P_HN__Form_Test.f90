program Tally_F_P_HN__Form_Test

  use Basics
  use Mathematics
  use Gravitations
  use Fluids

  implicit none

  character ( LDF ) :: &
    ProgramName = 'Tally_F_P_HN__Form_Test'

  integer ( KDI ) :: &
    iS  !-- iSelected
  real ( KDR ) :: &
    RadiusMax, &
    RadiusCore, &
    RadialRatio
  type ( Atlas_SCG_Form ), allocatable :: &
    A_G
  type ( Atlas_SCG_CC_Form ), allocatable :: &
    A_N
  type ( Gravitation_G_Form ), allocatable :: &
    G_G
  type ( Gravitation_N_SG_Form ), allocatable :: &
    G_N
  type ( Units_F_Form ), allocatable :: &
    U_G, U_N
  type ( Tally_F_P_HN_Form ), allocatable :: &
    T_G, T_N

  allocate ( PROGRAM_HEADER )  
  call PROGRAM_HEADER % Initialize ( ProgramName, DimensionalityOption = '2D' )

  allocate ( U_G, U_N )
  call U_G % Initialize ( )
  call U_N % Initialize ( TypeOption = 'ASTROPHYSICS' )

  RadiusCore   =   16.0_KDR  *  UNIT % KILOMETER
  RadiusMax    =  1.0e4_KDR  *  UNIT % KILOMETER
  RadialRatio  =  5.9_KDR

  allocate ( A_G, A_N )
  call A_G % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator, &
           NameOption = 'Atlas_G' )
  call A_N % Initialize &
         ( RadiusMax = 10.0_KDR, RadiusCore = 10.0_KDR / 8.0_KDR, &
           CommunicatorOption = PROGRAM_HEADER % Communicator, &
           NameOption = 'Atlas_N', &
           CoordinateUnitOption = U_N % Coordinate_PS, &
           RadialRatioOption = RadialRatio )

  allocate ( G_G, G_N )
  call G_G % Initialize ( A_G )
  call G_N % Initialize ( A_N, GravitationalConstant = 1.0_KDR )

  allocate ( T_G, T_N )
  call T_G % Initialize ( G_G, U_G )
  call T_N % Initialize ( G_N, U_N )

  call Show ( 'Tally_G' )
  call Show ( T_G % nSelected, 'nSelected' )
  call Show ( pack ( T_G % Variable, &
                     mask = [ ( any ( iS == T_G % iaSelected ), &
                                iS = 1, T_G % nIntegrals ) ] ),  &
              'Variable' )
  call Show ( pack ( T_G % Unit, &
                     mask = [ ( any ( iS == T_G % iaSelected ), &
                                iS = 1, T_G % nIntegrals ) ] ),  &
              'Unit' )

  call Show ( 'Tally_N' )
  call Show ( T_N % nSelected, 'nSelected' )
  call Show ( pack ( T_N % Variable, &
                     mask = [ ( any ( iS == T_N % iaSelected ), &
                                iS = 1, T_N % nIntegrals ) ] ),  &
              'Variable' )
  call Show ( pack ( T_N % Unit, &
                     mask = [ ( any ( iS == T_N % iaSelected ), &
                                iS = 1, T_N % nIntegrals ) ] ),  &
              'Unit' )

  deallocate ( T_N, T_G )
  deallocate ( G_N, G_G )
  deallocate ( A_N, A_G )
  deallocate ( U_N, U_G )
  deallocate ( PROGRAM_HEADER )

end program Tally_F_P_HN__Form_Test
