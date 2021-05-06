module Slope_DFV_A__Form

  !-- Slope_DivergenceFiniteVolume_Atlas_Form

  use Basics
  use Fields
  use RiemannSolver_HLL_C__Form
  use RiemannSolver_HLL_A__Form
  use Slope_DFV_C__Form

  implicit none
  private

  type, public, extends ( FieldSet_A_Form ) :: Slope_DFV_A_Form
    class ( CurrentSet_A_Form ), pointer :: &
      CurrentSet_A => null ( )
    class ( RiemannSolver_HLL_A_Form ), allocatable :: &
      RiemannSolver_A
  contains
    procedure, private, pass :: &
      InitializeAllocate_S
    generic, public :: &
      Initialize => InitializeAllocate_S
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Slope_DFV_A_Form


contains


  subroutine InitializeAllocate_S ( SA, CSA, NameOption )

    class ( Slope_DFV_A_Form ), intent ( inout ) :: &
      SA
    class ( CurrentSet_A_Form ), intent ( in ), target :: &
      CSA
    character ( * ), intent ( in ), optional :: &
      NameOption

    integer ( KDI ) :: &
      iC !-- iChart
    logical :: &
      PreviouslyAllocated
    character ( LDL ) :: &
      Name

    if ( SA % Type  ==  '' ) &
      SA % Type  =  'a Slope_DFV_A'

    Name  =  'SDFV_' // trim ( CSA % Name )
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    SA % CurrentSet_A  =>  CSA

    associate ( nC  =>  CSA % Atlas % nCharts )

    if ( allocated ( SA % FieldSet_C ) ) then
      PreviouslyAllocated  =  .true.
    else
      PreviouslyAllocated  =  .false.
      allocate ( SA % FieldSet_C ( nC ) )
    end if

    call SA % FieldSet_A_Form % Initialize &
           ( CSA % Atlas, &
             NameOption = Name, &
             IgnorabilityOption = CSA % IGNORABILITY )

    allocate ( SA % RiemannSolver_A )
    associate ( RSA  =>  SA % RiemannSolver_A )
    call RSA % Initialize ( CSA )

    if ( .not. PreviouslyAllocated ) then
      do iC  =  1,  nC

        allocate ( Slope_DFV_C_Form :: SA % FieldSet_C ( iC ) % Element ) 
        select type ( SC  =>  SA % FieldSet_C ( iC ) % Element )
        class is ( Slope_DFV_C_Form )

        select type ( RSC  =>  RSA % FieldSet_C ( iC ) % Element )
        class is ( RiemannSolver_HLL_C_Form )

        select type ( CSC  =>  CSA % FieldSet_C ( iC ) % Element )
        class is ( CurrentSet_C_Form )

        call SC % Initialize ( RSC, CSC, NameOption )

        end select !-- CSC
        end select !-- RSC
        end select !-- SC

      end do !-- iC
    end if !-- PreviouslyAllocated

    end associate !-- RSA
    end associate !-- nC

  end subroutine InitializeAllocate_S


  subroutine Compute ( SA, TimerLevelOption )

    class ( Slope_DFV_A_Form ), intent ( inout ) :: &
      SA
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    integer ( KDI ) :: &
      iC  !-- iChart

    do iC  =  1, size ( SA % FieldSet_C )
      select type ( SC  =>  SA % FieldSet_C ( iC ) % Element )
      class is ( Slope_DFV_C_Form )
      call SC % Compute ( TimerLevelOption )
      end select !-- SC
    end do !-- iC

  end subroutine Compute


  impure elemental subroutine Finalize ( SA )

    type ( Slope_DFV_A_Form ), intent ( inout ) :: &
      SA

    if ( allocated ( SA % RiemannSolver_A ) ) &
      deallocate ( SA % RiemannSolver_A )

    nullify ( SA % CurrentSet_A )

  end subroutine Finalize


end module Slope_DFV_A__Form
