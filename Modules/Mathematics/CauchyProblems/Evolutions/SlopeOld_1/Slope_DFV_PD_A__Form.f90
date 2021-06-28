module Slope_DFV_PD_A__Form

  !-- Slope_DivergenceFiniteVolume_PartialDerivative_Atlas_Form

  use Basics
  use Fields
  use RiemannSolver_HLL_C__Form
  use RiemannSolver_HLL_A__Form
  use Slope_H_A__Form
  use Slope_DFV_PD_C__Form

  implicit none
  private

  type, public, extends ( Slope_H_A_Form ) :: Slope_DFV_PD_A_Form
    class ( CurrentSet_A_Form ), pointer :: &
      CurrentSet_A => null ( )
    class ( RiemannSolver_HLL_A_Form ), pointer :: &
      RiemannSolver_A => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_PD
    generic, public :: &
      Initialize => InitializeAllocate_PD
    final :: &
      Finalize
  end type Slope_DFV_PD_A_Form


contains


  subroutine InitializeAllocate_PD ( SA, RSA, SuffixOption )

    class ( Slope_DFV_PD_A_Form ), intent ( inout ) :: &
      SA
    class ( RiemannSolver_HLL_A_Form ), intent ( in ), target :: &
      RSA
    character ( * ), intent ( in ), optional :: &
      SuffixOption

    integer ( KDI ) :: &
      iC !-- iChart
    logical :: &
      PreviouslyAllocated
    character ( LDL ) :: &
      Name

    if ( SA % Type  ==  '' ) &
      SA % Type  =  'a Slope_DFV_PD_A'

    Name  =  'S_DFV_PD_' // trim ( RSA % CurrentSet_A % Name )
    if ( present ( SuffixOption ) ) &
      Name  =  trim ( Name ) // '_' // trim ( SuffixOption )

    SA % CurrentSet_A     =>  RSA % CurrentSet_A
    SA % RiemannSolver_A  =>  RSA

    associate ( nC  =>  RSA % Atlas % nCharts )

    if ( allocated ( SA % FieldSet_C ) ) then
      PreviouslyAllocated  =  .true.
    else
      PreviouslyAllocated  =  .false.
      allocate ( SA % FieldSet_C ( nC ) )
    end if

    call SA % Slope_H_A_Form % Initialize &
           ( RSA % Atlas, &
             NameOption = Name, &
             IgnorabilityOption = RSA % IGNORABILITY )

    if ( .not. PreviouslyAllocated ) then
      do iC  =  1,  nC

        allocate ( Slope_DFV_PD_C_Form :: SA % FieldSet_C ( iC ) % Element ) 
        select type ( SC  =>  SA % FieldSet_C ( iC ) % Element )
          class is ( Slope_DFV_PD_C_Form )
        select type ( RSC  =>  RSA % FieldSet_C ( iC ) % Element )
          class is ( RiemannSolver_HLL_C_Form )

        call SC % Initialize ( RSC, SuffixOption )

        end select !-- RSC
        end select !-- SC

      end do !-- iC
    end if !-- PreviouslyAllocated

    end associate !-- nC

  end subroutine InitializeAllocate_PD


  impure elemental subroutine Finalize ( SA )

    type ( Slope_DFV_PD_A_Form ), intent ( inout ) :: &
      SA

    nullify ( SA % RiemannSolver_A )
    nullify ( SA % CurrentSet_A )

  end subroutine Finalize


end module Slope_DFV_PD_A__Form
