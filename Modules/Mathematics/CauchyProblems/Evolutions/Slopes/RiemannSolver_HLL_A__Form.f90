module RiemannSolver_HLL_A__Form

  !-- RiemannSolver_HartenLaxVanLeer_Atlas_Form

  use Basics
  use Fields
  use Reconstruction_A__Form
  use RiemannSolver_HLL_C__Form

  implicit none
  private

  type, public, extends ( FieldSet_A_Form ) :: RiemannSolver_HLL_A_Form
    class ( FieldSet_A_Form ), allocatable :: &
      Balanced_A
    class ( CurrentSet_A_Form ), pointer :: &
      CurrentSet_A => null ( )
    class ( FluxSet_A_Form ), allocatable :: &
      FluxSet_A
    class ( Eigenspeeds_F_A_Form ), allocatable :: &
      Eigenspeeds_A
    class ( Reconstruction_A_Form ), allocatable :: &
      Reconstruction_B_A, &
      Reconstruction_F_A, &
      Reconstruction_E_A
  contains
    procedure, private, pass :: &
      InitializeAllocate_RS
    generic, public :: &
      Initialize => InitializeAllocate_RS
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type RiemannSolver_HLL_A_Form


contains


  subroutine InitializeAllocate_RS &
               ( RSA, CSA, FieldOption, NameOption, nFieldsOption )

    class ( RiemannSolver_HLL_A_Form ), intent ( inout ) :: &
      RSA
    class ( CurrentSet_A_Form ), intent ( in ), target :: &
      CSA
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption

    integer ( KDI ) :: &
      iC !-- iChart
    integer ( KDI ), dimension ( : ), pointer :: &
      iaBalanced
    logical :: &
      PreviouslyAllocated
    character ( LDL ) :: &
      Name

    if ( RSA % Type  ==  '' ) &
      RSA % Type  =  'a RiemannSolver_HLL_A'

    Name  =  'RS_' // trim ( CSA % Name )
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    RSA % CurrentSet_A  =>  CSA

    select type ( CSC  =>  CSA % FieldSet_C ( 1 ) % Element )
    class is ( CurrentSet_C_Form )
    iaBalanced  =>  CSC % iaBalanced
    end select !-- CSC

    allocate ( RSA % Balanced_A )
    associate ( BA  =>  RSA % Balanced_A )
    call BA % Initialize &
           ( CSA, iaBalanced, &
             NameOption = trim ( CSA % Name ) // '_Balanced' )

    allocate ( RSA % FluxSet_A )
    associate ( FSA  =>  RSA % FluxSet_A )
    call FSA % Initialize ( CSA )

    allocate ( RSA % Eigenspeeds_A )
    associate ( EA  =>  RSA % Eigenspeeds_A )
    call EA % Initialize ( CSA )

    allocate ( RSA % Reconstruction_B_A )
    associate ( RBA  =>  RSA % Reconstruction_B_A )
    call RBA % Initialize ( CSA % Geometry_A, BA )

    allocate ( RSA % Reconstruction_F_A )
    associate ( RFA  =>  RSA % Reconstruction_F_A )
    call RFA % Initialize ( CSA % Geometry_A, FSA )

    allocate ( RSA % Reconstruction_E_A )
    associate ( REA  =>  RSA % Reconstruction_E_A )
    call REA % Initialize ( CSA % Geometry_A, EA )

    associate ( nC  =>  CSA % Atlas % nCharts )

    if ( allocated ( RSA % FieldSet_C ) ) then
      PreviouslyAllocated  =  .true.
    else
      PreviouslyAllocated  =  .false.
      allocate ( RSA % FieldSet_C ( nC ) )
    end if

    call RSA % FieldSet_A_Form % Initialize &
           ( CSA % Atlas, &
             NameOption = Name, &
             IgnorabilityOption = CSA % IGNORABILITY )

    if ( .not. PreviouslyAllocated ) then
      do iC  =  1,  nC

        allocate &
          ( RiemannSolver_HLL_C_Form :: RSA % FieldSet_C ( iC ) % Element ) 
        select type ( RSC  =>  RSA % FieldSet_C ( iC ) % Element )
        class is ( RiemannSolver_HLL_C_Form )

        associate &
          ( RBC  =>  RBA % Reconstruction_C ( iC ) % Element, &
            RFC  =>  RFA % Reconstruction_C ( iC ) % Element, &
            REC  =>  REA % Reconstruction_C ( iC ) % Element )

        select type ( EC  =>  EA % FieldSet_C ( iC ) % Element )
        class is ( Eigenspeeds_F_C_Form )

        select type ( FSC  =>  FSA % FieldSet_C ( iC ) % Element )
        class is ( FluxSet_C_Form )

        select type ( CSC  =>  CSA % FieldSet_C ( iC ) % Element )
        class is ( CurrentSet_C_Form )

        call RSC % Initialize &
               ( RBC, RFC, REC, EC, FSC, CSC, FieldOption, NameOption, &
                 nFieldsOption )

        end select !-- EC
        end select !-- FSC
        end select !-- CSC
        end associate !-- RBC, etc.
        end select !-- RSC

      end do !-- iC
    end if !-- PreviouslyAllocated

    end associate !-- nC
    end associate !-- REA
    end associate !-- RFA
    end associate !-- RBA
    end associate !-- EA
    end associate !-- FSA
    end associate !-- BA

  end subroutine InitializeAllocate_RS


  subroutine Compute ( RSA, iD, TimerLevelOption )

    class ( RiemannSolver_HLL_A_Form ), intent ( inout ) :: &
      RSA
    integer ( KDI ), intent ( in ) :: &
      iD  !-- iDimensions
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    integer ( KDI ) :: &
      iC  !-- iChart

    do iC  =  1, size ( RSA % FieldSet_C )
      select type ( RSC  =>  RSA % FieldSet_C ( iC ) % Element )
      class is ( RiemannSolver_HLL_C_Form )
      call RSC % Compute ( iD, TimerLevelOption )
      end select !-- EC
    end do !-- iC

  end subroutine Compute


  impure elemental subroutine Finalize ( RSA )

    type ( RiemannSolver_HLL_A_Form ), intent ( inout ) :: &
      RSA

    if ( allocated ( RSA % Reconstruction_E_A ) ) &
      deallocate ( RSA % Reconstruction_E_A )
    if ( allocated ( RSA % Reconstruction_F_A ) ) &
      deallocate ( RSA % Reconstruction_F_A )
    if ( allocated ( RSA % Reconstruction_B_A ) ) &
      deallocate ( RSA % Reconstruction_B_A )
    if ( allocated ( RSA % Eigenspeeds_A ) ) &
      deallocate ( RSA % Eigenspeeds_A )
    if ( allocated ( RSA % FluxSet_A ) ) &
      deallocate ( RSA % FluxSet_A )

    nullify ( RSA % CurrentSet_A )

    if ( allocated ( RSA % Balanced_A ) ) &
      deallocate ( RSA % Balanced_A )

  end subroutine Finalize


end module RiemannSolver_HLL_A__Form
