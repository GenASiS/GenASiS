module Step_RK_CSA__Form

  !-- Step_RungeKutta_FieldSetAtlas_Form

  use Basics
  use Fields
  use Slopes
  use Step_RK_H__Form

  implicit none
  private

  type, public, extends ( Step_RK_H_Form ) :: Step_RK_CSA_Form
    class ( CurrentSet_A_Form ), pointer :: &
      CurrentSet_A
    class ( FieldSet_A_Form ), allocatable :: &
      Slope_A
  contains
    procedure, private, pass :: &
      Initialize_CSA
    generic, public :: &
      Initialize => Initialize_CSA
    procedure, private, pass :: &
      Show_S
    final :: &
      Finalize
  end type Step_RK_CSA_Form


contains


  subroutine Initialize_CSA ( S, CSA, NameOption )

    class ( Step_RK_CSA_Form ), intent ( inout ) :: &
      S
    class ( CurrentSet_A_Form ), intent ( in ), target :: &
      CSA
    character ( * ), intent ( in ), optional :: &
      NameOption

    real ( KDR ), dimension ( 2 : 2, 1 : 1 ) :: &
      A
    real ( KDR ), dimension ( 2 : 2 ) :: &
      C
    real ( KDR ), dimension ( 1 : 2 ) :: &
      B

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Step_RK_CSA'

    S % CurrentSet_A  =>  CSA

    call Clear ( A )
    A ( 2, 1 ) = 1.0_KDR

    B ( 1 ) = 0.5_KDR
    B ( 2 ) = 0.5_KDR

    C ( 2 ) = 1.0_KDR
    
    call S % Step_RK_H_Form % Initialize ( A, B, C, NameOption )

    if ( .not. allocated ( S % Slope_A ) ) then
      allocate ( Slope_DFV_A_Form :: S % Slope_A )
      select type ( SA  =>  S % Slope_A )
      class is ( Slope_DFV_A_Form )
      call SA % Initialize ( CSA )
      end select !-- SA
    end if !-- allocated Slope

  end subroutine Initialize_CSA


  subroutine Show_S ( S )

    class ( Step_RK_CSA_Form ), intent ( in ) :: &
      S

    call S % Step_RK_H_Form % Show ( )
    call S % Slope_A % Show ( )

  end subroutine Show_S


  impure elemental subroutine Finalize ( S )

    type ( Step_RK_CSA_Form ), intent ( inout ) :: &
      S

    if ( allocated ( S % Slope_A ) ) &
      deallocate ( S % Slope_A )

    nullify ( S % CurrentSet_A )

  end subroutine Finalize


end module Step_RK_CSA__Form
