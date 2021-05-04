module RiemannSolver_HLL_C__Form

  !-- RiemannSolver_HartenLaxVanLeer_Chart_Form

  use Basics
  use Manifolds
  use Fields
  use Reconstruction_C__Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_SOLVER_SPEEDS_HLL  =  2

  type, public, extends ( FieldSet_C_Form ) :: RiemannSolver_HLL_C_Form
    integer ( KDI ) :: &
      iTimer = 0
    integer ( KDI ) :: &
      N_SOLVER_SPEEDS_HLL = N_SOLVER_SPEEDS_HLL
    integer ( KDI ) :: &
      ALPHA_PLUS_U    = 0, &
      ALPHA_MINUS_U   = 0, &
      N_SOLVER_SPEEDS = 0
    class ( CurrentSet_C_Form ), pointer :: &
      CurrentSet_C => null ( )
    class ( FluxSet_C_Form ), pointer :: &
      FluxSet_C => null ( )
    class ( Eigenspeeds_F_C_Form ), pointer :: &
      Eigenspeeds_C => null ( )
    class ( Reconstruction_C_Form ), pointer :: &
      Reconstruction_B_C => null ( ), &
      Reconstruction_F_C => null ( ), &
      Reconstruction_E_C => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_RS
    generic, public :: &
      Initialize => InitializeAllocate_RS
    procedure, private, pass :: &
      Show_FSC
    final :: &
      Finalize
  end type RiemannSolver_HLL_C_Form


contains


  subroutine InitializeAllocate_RS &
               ( RSC, RBC, RFC, REC, EC, FSC, CSC, FieldOption, NameOption, &
                 nFieldsOption )

    class ( RiemannSolver_HLL_C_Form ), intent ( inout ) :: &
      RSC
    class ( Reconstruction_C_Form ), intent ( in ), target :: &
      RFC, &
      REC, &
      RBC
    class ( Eigenspeeds_F_C_Form ), intent ( in ), target :: &
      EC
    class ( FluxSet_C_Form ), intent ( in ), target :: &
      FSC
    class ( CurrentSet_C_Form ), intent ( in ), target :: &
      CSC
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption

    integer ( KDI ) :: &
      iB, &  !-- iBalanced
      iF, &  !-- iField
      nFields
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    if ( RSC % Type  ==  '' ) &
      RSC % Type  =  'a RiemannSolver_HLL_C' 
    
    Name  =  'RS_' // trim ( CSC % Name )
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    RSC % CurrentSet_C        =>  CSC
    RSC % FluxSet_C           =>  FSC
    RSC % Eigenspeeds_C       =>  EC
    RSC % Reconstruction_B_C  =>  RBC
    RSC % Reconstruction_F_C  =>  RFC
    RSC % Reconstruction_E_C  =>  REC

    associate &
      ( nB  =>  CSC % nBalanced, &
        DeviceMemory  =>  CSC % Storage_FSC % DeviceMemory, &
        PinnedMemory  =>  CSC % Storage_FSC % PinnedMemory, &
        DevicesCommunicate  =>  CSC % GhostExchange_FSC % DevicesCommunicate ) 

    !-- Field indices

    if ( present ( nFieldsOption ) ) then
      nFields  =  nFieldsOption
    else
      RSC % N_SOLVER_SPEEDS  =  RSC % N_SOLVER_SPEEDS_HLL
      nFields  =  nB  +  RSC % N_SOLVER_SPEEDS
    end if

    RSC % ALPHA_PLUS_U   =  nB  +  1
    RSC % ALPHA_MINUS_U  =  nB  +  2

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    do iB  =  1,  nB
      iF  =  CSC % iaBalanced ( iB )
      Field ( iB )  =  CSC % Field ( iF )
    end do !-- iS

    Field ( nB + 1 : nB + RSC % N_SOLVER_SPEEDS ) &
      =  [ 'AlphaPlus_U ', &
           'AlphaMinus_U' ]
          
    !-- FieldSet

    call RSC % FieldSet_C_Form % Initialize &
           ( CSC % Chart, &
             FieldOption = Field, &
             NameOption = Name, &
             DeviceMemoryOption = DeviceMemory, &
             PinnedMemoryOption = PinnedMemory, &
             DevicesCommunicateOption = DevicesCommunicate, &
             nFieldsOption = nFields, &
             IgnorabilityOption = CSC % IGNORABILITY )

    end associate !-- nB, etc.

  end subroutine InitializeAllocate_RS


  subroutine Show_FSC ( FSC )

    class ( RiemannSolver_HLL_C_Form ), intent ( in ) :: &
      FSC

    call FSC % FieldSet_C_Form % Show ( )
    call FSC % FluxSet_C % Show ( )
    call FSC % Eigenspeeds_C % Show ( )
    call FSC % Reconstruction_B_C % Show ( )
    call FSC % Reconstruction_F_C % Show ( )
    call FSC % Reconstruction_E_C % Show ( )

  end subroutine Show_FSC


  impure elemental subroutine Finalize ( RSC )

    type ( RiemannSolver_HLL_C_Form ), intent ( inout ) :: &
      RSC

    nullify ( RSC % Reconstruction_E_C )
    nullify ( RSC % Reconstruction_F_C )
    nullify ( RSC % Reconstruction_B_C )
    nullify ( RSC % Eigenspeeds_C )
    nullify ( RSC % FluxSet_C )
    nullify ( RSC % CurrentSet_C )

  end subroutine Finalize


end module RiemannSolver_HLL_C__Form
