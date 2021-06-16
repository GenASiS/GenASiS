module Slope_DFV_F_A__Form

  !-- Slope_DivergenceFiniteVolume_Flat_Atlas_Form

  use Basics
  use RiemannSolver_HLL_C__Form
  use RiemannSolver_HLL_A__Form
  use Slope_H_C__Form
  use Slope_H_A__Form
  use Slope_DFV_PD_A__Form

  implicit none
  private

  type, public, extends ( Slope_H_A_Form ) :: Slope_DFV_F_A_Form
  contains
    procedure, private, pass :: &
      InitializeAllocate_F
    generic, public :: &
      Initialize => InitializeAllocate_F
    final :: &
      Finalize
  end type Slope_DFV_F_A_Form


contains


  subroutine InitializeAllocate_F ( SA, RSA, NameOption )

    class ( Slope_DFV_F_A_Form ), intent ( inout ) :: &
      SA
    class ( RiemannSolver_HLL_A_Form ), intent ( in ), target :: &
      RSA
    character ( * ), intent ( in ), optional :: &
      NameOption

    integer ( KDI ) :: &
      iC !-- iChart
    logical :: &
      PreviouslyAllocated
    character ( LDL ) :: &
      Name

    if ( SA % Type  ==  '' ) &
      SA % Type  =  'a Slope_DFV_F_A'

    Name  =  'S_DFV_' // trim ( RSA % CurrentSet_A % Name )
    if ( present ( NameOption ) ) &
      Name  =  NameOption

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

        allocate ( Slope_H_C_Form :: SA % FieldSet_C ( iC ) % Element ) 
        select type ( SC  =>  SA % FieldSet_C ( iC ) % Element )
          class is ( Slope_H_C_Form )
        select type ( RSC  =>  RSA % FieldSet_C ( iC ) % Element )
          class is ( RiemannSolver_HLL_C_Form )
        associate &
          ( CSC  =>  RSC % CurrentSet_C )
        associate &
          ( nB  =>  CSC % nBalanced, &
            DeviceMemory  =>  CSC % Storage_FSC % DeviceMemory, &
            PinnedMemory  =>  CSC % Storage_FSC % PinnedMemory, &
            DevicesCommunicate  =>  CSC % GhostExchange_FSC &
                                        % DevicesCommunicate ) 

        call SC % Initialize &
               ( CSC % Chart, &
                 FieldOption = CSC % Balanced, &
                 NameOption = Name, &
                 DeviceMemoryOption = DeviceMemory, &
                 PinnedMemoryOption = PinnedMemory, &
                 DevicesCommunicateOption = DevicesCommunicate, &
                 nFieldsOption = nB, &
                 IgnorabilityOption = CSC % IGNORABILITY )

        end associate !-- nB, etc.
        end associate !-- CSC
        end select !-- RSC
        end select !-- SC

      end do !-- iC
    end if !-- PreviouslyAllocated

    end associate !-- nC

    !-- Slope component: Partial derivative

    associate ( nC  =>  SA % nComponents )
    nC  =  nC + 1
    allocate ( Slope_DFV_PD_A_Form :: SA % Component_A ( iC ) % Element )
    select type ( SPDA  =>  SA % Component_A ( iC ) % Element )
      class is ( Slope_DFV_PD_A_Form )
    call SPDA % Initialize ( RSA, NameOption )
    end select !-- SPDA
    end associate !-- nC

  end subroutine InitializeAllocate_F


  impure elemental subroutine Finalize ( SA )

    type ( Slope_DFV_F_A_Form ), intent ( inout ) :: &
      SA

  end subroutine Finalize


end module Slope_DFV_F_A__Form
