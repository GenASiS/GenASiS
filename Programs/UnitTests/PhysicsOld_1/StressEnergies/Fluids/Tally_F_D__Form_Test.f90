program Tally_F_D__Form_Test

  use Basics
  use Mathematics
  use Spaces
  use Tally_F_D__Form

  implicit none

  character ( LDF ) :: &
    ProgramName = 'Tally_F_D__Form_Test'

  integer ( KDI ) :: &
    iS  !-- iSelected
  type ( Atlas_SC_Form ), allocatable :: &
    A_G, A_N
  type ( Geometry_ASC_Form ), allocatable :: &
    GA_G, GA_N
  type ( Tally_F_D_Form ), allocatable :: &
    T_G, T_N

  allocate ( PROGRAM_HEADER )  
  call PROGRAM_HEADER % Initialize ( ProgramName )

  allocate ( A_G, A_N )
  call A_G % Initialize ( 'Atlas_SC_G', PROGRAM_HEADER % Communicator )
  call A_N % Initialize ( 'Atlas_SC_N', PROGRAM_HEADER % Communicator )
  call A_G % CreateChart ( )
  call A_N % CreateChart ( )

  allocate ( GA_G, GA_N )
  call GA_G % Initialize ( A_G, 'GALILEAN' )
  call GA_N % Initialize ( A_N, 'NEWTONIAN', &
                           GravitySolverTypeOption = 'MULTIPOLE' )
  call A_G % SetGeometry ( GA_G )
  call A_N % SetGeometry ( GA_N )

  allocate ( T_G, T_N )
  call T_G % Initialize ( A_G )
  call T_N % Initialize ( A_N )

  call Show ( 'Tally_G' )
  call Show ( T_G % nSelected, 'nSelected' )
  do iS = 1, T_G % nSelected
    call Show ( T_G % Variable ( T_G % iaSelected ( iS ) ), 'Variable' ) 
  end do !-- iI

  call Show ( 'Tally_N' )
  call Show ( T_N % nSelected, 'nSelected' )
  do iS = 1, T_N % nSelected
    call Show ( T_N % Variable ( T_N % iaSelected ( iS ) ), 'Variable' ) 
  end do !-- iI

  deallocate ( T_N, T_G )
  deallocate ( GA_N, GA_G )
  deallocate ( A_N, A_G )
  deallocate ( PROGRAM_HEADER )

end program Tally_F_D__Form_Test
