module Slope_DFV_DP__Form

  !-- Slope_DivergenceFiniteVolume_DivergencePart__Form

  use Basics
  use Fields
  use RiemannSolver_HLL__Form
  use Slope_H__Form
  use Slope_DFV_PD__Form
  use Slope_DFV_C_F__Form

  implicit none
  private

  type, public, extends ( Slope_H_Form ) :: Slope_DFV_DP_Form
    logical ( KDL ) :: &
      StreamFluxes
    type ( FieldSetElement ), dimension ( : ), allocatable :: &
      FluxSetDimension, &
      FluxSet_IL_Dimension, &
      FluxSet_IR_Dimension, &
      FluxSet_RS_Dimension
  contains
    procedure, private, pass :: &
      InitializeAllocate_DP
    generic, public :: &
      Initialize => InitializeAllocate_DP
    procedure, public, pass :: &
      SetStream
    procedure, public, pass :: &
      ComputePartialDerivative
    procedure, public, pass :: &
      ComputeConnectionFlat
    procedure, public, pass :: &
      ClearRecursive
    procedure, public, pass :: &
      MultiplyAddRecursive
    final :: &
      Finalize
  end type Slope_DFV_DP_Form


contains


  subroutine InitializeAllocate_DP ( S, RS, DP, IgnorabilityOption )

    class ( Slope_DFV_DP_Form ), intent ( inout ) :: &
      S
    class ( RiemannSolver_HLL_Form ), intent ( in ) :: &
      RS
    class ( DivergencePart_CS_Form ), intent ( in ) :: &
      DP
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    integer ( KDI ) :: &
      iC, &  !-- iChart
      iD, &  !-- iDimension
      nD     !-- nDimensions
    character ( 1 ) :: &
      DimensionNumber
    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Slope_DFV_DP'

    associate ( CS  =>  RS % CurrentSet )

    Name  =  trim ( CS % Name ) // '_Slp_DFV_DP_' // trim ( DP % Name )

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

    !-- Slope components: Partial derivative and flat connection

    associate ( nSC  =>  S % nComponents )

    nSC  =  nSC + 1
    allocate ( Slope_DFV_PD_Form :: S % Component ( nSC ) % Element )
    select type ( SPD  =>  S % Component ( nSC ) % Element )
      class is ( Slope_DFV_PD_Form )
    call SPD % Initialize ( RS, DP )
    end select !-- SPD

    nSC  =  nSC + 1
    allocate ( Slope_DFV_C_F_Form :: S % Component ( nSC ) % Element )
    select type ( SCF  =>  S % Component ( nSC ) % Element )
      class is ( Slope_DFV_C_F_Form )
    call SCF % Initialize ( DP )
    end select !-- SCF

    end associate !-- nSC

    !-- Stream parameters

    S % StreamComponents  =  .false.

    S % StreamFluxes  =  .false.
    call PROGRAM_HEADER % GetParameter ( S % StreamFluxes, 'StreamFluxes' )

    if ( S % StreamFluxes ) then

      associate ( A  =>  S % Atlas )
      nD  =  A % Chart ( 1 ) % Element % nDimensions
      do iC  =  2, A % nCharts
        nD  =  max ( nD, A % Chart ( iC ) % Element % nDimensions )
      end do
      end associate !-- A

      allocate &
        ( S % FluxSetDimension ( nD ), &
          S % FluxSet_IL_Dimension ( nD ), &
          S % FluxSet_IR_Dimension ( nD ), &
          S % FluxSet_RS_Dimension ( nD ) )
      do iD  =  1, nD
        write ( DimensionNumber, fmt = '(i1.1)' ) iD

        allocate ( S % FluxSetDimension ( iD ) % Element )
        associate ( FSD  =>  S % FluxSetDimension ( iD ) % Element )
        call FSD % Initialize &
               ( S % Atlas, &
                 FieldOption = S % Field, &
                 NameOption = 'FS_' // trim ( S % Name ) // '_' &
                              // DimensionNumber, &
                 DeviceMemoryOption = S % DeviceMemory, &
                 DevicesCommunicateOption = S % DevicesCommunicate, &
                 nFieldsOption = S % nFields, &
                 IgnorabilityOption = S % IGNORABILITY + 1 )
        end associate !-- FSD

        allocate ( S % FluxSet_IL_Dimension ( iD ) % Element )
        associate ( FSD  =>  S % FluxSet_IL_Dimension ( iD ) % Element )
        call FSD % Initialize &
               ( S % Atlas, &
                 FieldOption = S % Field, &
                 NameOption = 'FS_IL_' // trim ( S % Name ) // '_' &
                              // DimensionNumber, &
                 DeviceMemoryOption = S % DeviceMemory, &
                 DevicesCommunicateOption = S % DevicesCommunicate, &
                 nFieldsOption = S % nFields, &
                 IgnorabilityOption = S % IGNORABILITY + 1 )
        end associate !-- FSD

        allocate ( S % FluxSet_IR_Dimension ( iD ) % Element )
        associate ( FSD  =>  S % FluxSet_IR_Dimension ( iD ) % Element )
        call FSD % Initialize &
               ( S % Atlas, &
                 FieldOption = S % Field, &
                 NameOption = 'FS_IR_' // trim ( S % Name ) // '_' &
                              // DimensionNumber, &
                 DeviceMemoryOption = S % DeviceMemory, &
                 DevicesCommunicateOption = S % DevicesCommunicate, &
                 nFieldsOption = S % nFields, &
                 IgnorabilityOption = S % IGNORABILITY + 1 )
        end associate !-- FSD

        allocate ( S % FluxSet_RS_Dimension ( iD ) % Element )
        associate ( FSD  =>  S % FluxSet_RS_Dimension ( iD ) % Element )
        call FSD % Initialize &
               ( S % Atlas, &
                 FieldOption = RS % Field, &
                 NameOption = 'FS_RS_' // trim ( S % Name ) // '_' &
                              // DimensionNumber, &
                 DeviceMemoryOption = RS % DeviceMemory, &
                 DevicesCommunicateOption = RS % DevicesCommunicate, &
                 nFieldsOption = RS % nFields, &
                 IgnorabilityOption = S % IGNORABILITY + 1 )
        end associate !-- FSD

      end do !-- iD

    end if

  end subroutine InitializeAllocate_DP


  subroutine SetStream ( S, Sm )

    class ( Slope_DFV_DP_Form ), intent ( inout ) :: &
      S
    class ( StreamForm ), intent ( inout ) :: &
      Sm

    integer ( KDI ) :: &
      iD  !-- iDimension

    call S % Slope_H_Form % SetStream ( Sm )

    if ( S % StreamFluxes ) then
      do iD  =  1, size ( S % FluxSetDimension )
        associate &
          ( FSD     =>  S  % FluxSetDimension ( iD ) % Element, &
            FSD_IL  =>  S  % FluxSet_IL_Dimension ( iD ) % Element, &
            FSD_IR  =>  S  % FluxSet_IR_Dimension ( iD ) % Element, &
            FSD_RS  =>  S  % FluxSet_RS_Dimension ( iD ) % Element )
        call Sm % AddFieldSet ( FSD )
        call Sm % AddFieldSet ( FSD_IL )
        call Sm % AddFieldSet ( FSD_IR )
        call Sm % AddFieldSet ( FSD_RS )
        end associate !-- FSD, etc.
      end do !-- iD
    end if

  end subroutine SetStream


  subroutine ComputePartialDerivative ( S, iC, iD, T_Option )

    class ( Slope_DFV_DP_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iC, &  !-- iChart
      iD     !-- iDimension
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    select type ( S_PD  =>  S % Component ( 1 ) % Element )
      class is ( Slope_DFV_PD_Form )

    call S_PD % ComputeDimension ( iC, iD, T_Option )

    if ( allocated ( S % FluxSetDimension ) ) then
      associate &
        ( FS_RS     =>  S_PD % RiemannSolver % FluxSet, &
          FS_IL_RS  =>  S_PD % RiemannSolver % FluxSet_IL, &
          FS_IR_RS  =>  S_PD % RiemannSolver % FluxSet_IR, &
          FS_RS_RS  =>  S_PD % RiemannSolver, &
          FS_D      =>  S % FluxSetDimension     ( iD ) % Element, &
          FS_IL_D   =>  S % FluxSet_IL_Dimension ( iD ) % Element, &
          FS_IR_D   =>  S % FluxSet_IR_Dimension ( iD ) % Element, &
          FS_RS_D   =>  S % FluxSet_RS_Dimension ( iD ) % Element )
      call FS_RS    % Copy ( FS_D )
      call FS_IL_RS % Copy ( FS_IL_D )
      call FS_IR_RS % Copy ( FS_IR_D )
      call FS_RS_RS % Copy ( FS_RS_D )
      end associate !-- FS_RS, etc.
    end if !-- StreamFluxes

    end select !-- S_PD

  end subroutine ComputePartialDerivative


  subroutine ComputeConnectionFlat ( S, iC, T_Option )

    class ( Slope_DFV_DP_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iC  !-- iChart
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    select type ( S_C_F  =>  S % Component ( 2 ) % Element )
      class is ( Slope_DFV_C_F_Form )

    call S_C_F % ComputeChart ( iC, T_Option )

    end select !-- S_C_F

  end subroutine ComputeConnectionFlat


  subroutine ClearRecursive ( S )

    class ( Slope_DFV_DP_Form ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iD  !-- iDimension

    call S % Slope_H_Form % ClearRecursive ( )

    if ( S % StreamFluxes ) then
      do iD  =  1, size ( S % FluxSetDimension )
        associate &
          ( FSD     =>  S  % FluxSetDimension ( iD ) % Element, &
            FSD_IL  =>  S  % FluxSet_IL_Dimension ( iD ) % Element, &
            FSD_IR  =>  S  % FluxSet_IR_Dimension ( iD ) % Element, &
            FSD_RS  =>  S  % FluxSet_RS_Dimension ( iD ) % Element )
        call FSD    % Clear ( )
        call FSD_IL % Clear ( )
        call FSD_IR % Clear ( )
        call FSD_RS % Clear ( )
        end associate !-- FSD, etc.
      end do !-- iD
    end if

  end subroutine ClearRecursive


  subroutine MultiplyAddRecursive ( S, SS, B )

    class ( Slope_DFV_DP_Form ), intent ( inout ) :: &
      S  !-- Slope
    class ( Slope_H_Form ), intent ( in ) :: &
      SS  !-- SlopeStage
    real ( KDR ) :: &
      B  !-- RungeKutta weight

    integer ( KDI ) :: &
      iD  !-- iDimension

    call S % Slope_H_Form % MultiplyAddRecursive ( SS, B )

    if ( S % StreamFluxes ) then
      do iD  =  1, size ( S % FluxSetDimension )
        select type ( SS )
          class is ( Slope_DFV_DP_Form )
        associate &
          ( FSD        =>  S  % FluxSetDimension ( iD ) % Element, &
            FSD_SS     =>  SS % FluxSetDimension ( iD ) % Element, &
            FSD_IL     =>  S  % FluxSet_IL_Dimension ( iD ) % Element, &
            FSD_IL_SS  =>  SS % FluxSet_IL_Dimension ( iD ) % Element, &
            FSD_IR     =>  S  % FluxSet_IR_Dimension ( iD ) % Element, &
            FSD_IR_SS  =>  SS % FluxSet_IR_Dimension ( iD ) % Element, &
            FSD_RS     =>  S  % FluxSet_RS_Dimension ( iD ) % Element, &
            FSD_RS_SS  =>  SS % FluxSet_RS_Dimension ( iD ) % Element )
        call FSD    % MultiplyAdd ( FSD_SS,    B )
        call FSD_IL % MultiplyAdd ( FSD_IL_SS, B )
        call FSD_IR % MultiplyAdd ( FSD_IR_SS, B )
        call FSD_RS % MultiplyAdd ( FSD_RS_SS, B )
        end associate !-- FSD, etc.
        end select !-- SS
      end do !-- iD
    end if

    !-- This slope

  end subroutine MultiplyAddRecursive


  impure elemental subroutine Finalize ( S )

    type ( Slope_DFV_DP_Form ), intent ( inout ) :: &
      S

    if ( allocated ( S % FluxSet_RS_Dimension ) ) &
      deallocate ( S % FluxSet_RS_Dimension )
    if ( allocated ( S % FluxSet_IR_Dimension ) ) &
      deallocate ( S % FluxSet_IR_Dimension )
    if ( allocated ( S % FluxSet_IL_Dimension ) ) &
      deallocate ( S % FluxSet_IL_Dimension )
    if ( allocated ( S % FluxSetDimension ) ) &
      deallocate ( S % FluxSetDimension )

  end subroutine Finalize


end module Slope_DFV_DP__Form
