module Slope_DFV_F_DP__Form

  !-- Slope_DivergenceFiniteVolume_Flat_DivergenceParts__Form

  use Basics
  use Fields
  use RiemannSolver_HLL__Form
  use Slope_H__Form
  use Slope_DFV_DP__Form
  use Slope_DFV_DD__Form

  implicit none
  private

  type, public, extends ( Slope_H_Form ) :: Slope_DFV_F_DP_Form
    class ( RiemannSolver_HLL_Form ), pointer :: &
      RiemannSolver => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_F
    generic, public :: &
      Initialize => InitializeAllocate_F
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Slope_DFV_F_DP_Form


contains


  subroutine InitializeAllocate_F ( S, RS, DP_1D, IgnorabilityOption )

    class ( Slope_DFV_F_DP_Form ), intent ( inout ) :: &
      S
    class ( RiemannSolver_HLL_Form ), intent ( in ), target :: &
      RS
    type ( DivergencePartElement ), dimension ( : ), intent ( in ) :: &
      DP_1D
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    integer ( KDI ) :: &
      iDP
    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Slope_DFV_F_DP'

    Name  =  trim ( RS % CurrentSet % Name ) // '_Slp_DFV_F_DP'

    S % RiemannSolver  =>  RS

    associate ( CS  =>  RS % CurrentSet )
    call S % Slope_H_Form % Initialize &
           ( CS % Atlas, &
             FieldOption = CS % Balanced, &
             NameOption = Name, &
             DeviceMemoryOption = CS % DeviceMemory, &
             PinnedMemoryOption = CS % PinnedMemory, &
             DevicesCommunicateOption = CS % DevicesCommunicate, &
             nFieldsOption = CS % nBalanced, &
             IgnorabilityOption = IgnorabilityOption )
    end associate !-- CS

    !-- Slope components: divergence contributions

    associate ( nSC  =>  S % nComponents )
    do iDP  =  1, size ( DP_1D )
      nSC  =  nSC + 1
      allocate ( Slope_DFV_DP_Form :: S % Component ( nSC ) % Element )
      select type ( SDP  =>  S % Component ( nSC ) % Element )
      class is ( Slope_DFV_DP_Form )
        call SDP % Initialize ( RS, DP_1D ( iDP ) % Element )
      end select !-- SDP
    end do !-- iDP
    end associate !-- nSC

    !-- Slope component: divergence diffusion

    associate ( nSC  =>  S % nComponents )
    nSC  =  nSC + 1
    allocate ( Slope_DFV_DD_Form :: S % Component ( nSC ) % Element )
    select type ( SDD  =>  S % Component ( nSC ) % Element )
    class is ( Slope_DFV_DD_Form )
      call SDD % Initialize ( RS )
    end select !-- SDD
    end associate !-- nSC

  end subroutine InitializeAllocate_F


  subroutine Compute ( S, T_Option )

    class ( Slope_DFV_F_DP_Form ), intent ( inout ) :: &
      S
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iC, &  !-- iChart
      iD, &  !-- iDimension
      iDP    !-- iDivergencePart
    type ( TimerForm ), pointer :: &
      T_RPP

    associate ( nC  =>  S % Atlas % nCharts )
    do iC  =  1,  nC
      associate ( C  =>  S % Atlas % Chart ( 1 ) % Element )

      !-- Partial derivative contributions
      associate ( RS  =>  S % RiemannSolver )
      
      if ( present ( T_Option ) ) then
        T_RPP  =>  RS % Timer_P ( Level = T_Option % Level + 1 )
      else
        T_RPP  =>  null ( )
      end if

      do iD  =  1,  C % nDimensions

        !-- Reconstruction and eigenspeeds
        if ( associated ( T_RPP ) ) call T_RPP % Start ( )
        call RS % Prepare ( iC = iC, iD = iD, T_Option = T_RPP )
        if ( associated ( T_RPP ) ) call T_RPP % Stop ( )

        !-- Flux contributions
        do iDP  =  1,  S % nComponents - 1
          select type ( SDP  =>  S % Component ( iDP ) % Element )
          class is ( Slope_DFV_DP_Form )
            call SDP % ComputePartialDerivative ( iC, iD, T_Option = T_Option )
          end select !-- SDP
        end do !-- iDP

        !-- Diffusive term
        select type ( SDD  =>  S % Component ( S % nComponents ) % Element )
        class is ( Slope_DFV_DD_Form )
          call SDD % ComputeDimension ( iC, iD, T_Option = T_Option )
        end select !-- SDD

      end do !-- iD
      
      end associate !-- RS

      !-- Connection contributions and sums
      do iDP  =  1,  S % nComponents - 1
        select type ( SDP  =>  S % Component ( iDP ) % Element )
        class is ( Slope_DFV_DP_Form )
          call SDP % ComputeConnectionFlat ( iC, T_Option = T_Option )
          call SDP % AddComponents ( T_Option = T_Option )
        end select !-- SDP
      end do !-- iDP

      !-- Diffusive term

      end associate !-- C
    end do !-- iC

    call S % AddComponents ( T_Option = T_Option )

    end associate !-- nC

  end subroutine Compute


  impure elemental subroutine Finalize ( S )

    type ( Slope_DFV_F_DP_Form ), intent ( inout ) :: &
      S

  end subroutine Finalize


end module Slope_DFV_F_DP__Form
