module Slope_H_C__Form

  !-- Slope_Header_Chart_Form

  use Basics
  use Algebra
  use Manifolds
  use Fields

  implicit none
  private

    integer, private, parameter :: &
      MAX_COMPONENTS = 16

  type, public, extends ( FieldSet_C_Form ) :: Slope_H_C_Form
    integer ( KDI ) :: &
      nComponents = 0
    class ( Slope_H_C_Pointer ), dimension ( : ), pointer :: &
      Component_C => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_FS
    procedure, public, pass :: &
      AddComponent
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Slope_H_C_Form

  type, public :: Slope_H_C_Pointer
    class ( Slope_H_C_Form ), pointer :: &
      Pointer => null ( )
  contains
    final :: &
      Finalize_P
  end type Slope_H_C_Pointer


contains


  subroutine InitializeAllocate_FS &
               ( FSC, C, FieldOption, VectorOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, UnitOption, VectorIndicesOption, &
                 nFieldsOption, IgnorabilityOption )

    class ( Slope_H_C_Form ), intent ( inout ) :: &
      FSC
    class ( Chart_H_Form ), intent ( in ), target :: &
      C
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      DeviceMemoryOption, &
      PinnedMemoryOption, &
      DevicesCommunicateOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption, &
      IgnorabilityOption

    if ( FSC % Type  ==  '' ) &
      FSC % Type  =  'a Slope_H_C' 
    
    call FSC % FieldSet_C_Form % Initialize &
           ( C, FieldOption, VectorOption, NameOption, &
             DeviceMemoryOption, PinnedMemoryOption, &
             DevicesCommunicateOption, UnitOption, VectorIndicesOption, &
             nFieldsOption, IgnorabilityOption )

    allocate ( FSC % Component_C ( MAX_COMPONENTS ) )

  end subroutine InitializeAllocate_FS


  subroutine AddComponent ( SC, SCC )

    class ( Slope_H_C_Form ), intent ( inout ) :: &
      SC
    class ( Slope_H_C_Form ), intent ( in ), target :: &
      SCC

    associate ( nC  =>  SC % nComponents )
    nC  =  nC + 1
    SC % Component_C ( nC ) % Pointer  =>  SCC
    end associate !-- nC

  end subroutine AddComponent


  subroutine Compute ( SC, TimerLevelOption, iS_Option )

    class ( Slope_H_C_Form ), intent ( inout ) :: &
      SC
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption, &
      iS_Option

    integer ( KDI ) :: &
      iC  !-- iComponent

call Show ( SC % Name, '>>> Name' )
    if ( SC % nComponents  >  0 ) then

      call SC % Clear ( )

      do iC  =  1, SC % nComponents
        associate &
          ( SCC  =>  SC % Component_C ( iC ) % Pointer )
        associate &
          (  SV  =>   SC % Storage_FSC % Storage % Value, &
            SCV  =>  SCC % Storage_FSC % Storage % Value )

        call SCC % Compute ( TimerLevelOption, iS_Option )
        call MultiplyAdd &
               ( SV, SCV, 1.0_KDR, &
                 UseDeviceOption = SC % Storage_FSC % DeviceMemory )

        end associate !-- SV, etc.
        end associate !-- SCC
      end do !-- iC

    else
      call Show ( 'Slope has no components', CONSOLE % ERROR )
      call Show ( 'Compute must be overridden', CONSOLE % ERROR )
      call Show ( SC % Name, 'Name', CONSOLE % ERROR )
      call Show ( 'Slope_H_C_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

  end subroutine Compute


  subroutine Finalize ( SC )

    type ( Slope_H_C_Form ), intent ( inout ) :: &
      SC

    if ( associated ( SC % Component_C ) ) &
      deallocate ( SC % Component_C )

  end subroutine Finalize


  impure elemental subroutine Finalize_P ( SP )
    
    type ( Slope_H_C_Pointer ), intent ( inout ) :: &
      SP

    nullify ( SP % Pointer )

  end subroutine Finalize_P


end module Slope_H_C__Form
