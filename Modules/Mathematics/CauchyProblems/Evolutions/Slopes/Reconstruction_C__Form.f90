module Reconstruction_C__Form

  !-- Reconstruction_Chart_Form

  use Basics
  use Manifolds
  use Fields

  implicit none
  private

  type, public :: Reconstruction_C_Form
    integer ( KDI ) :: &
      IGNORABILITY, &
      iTimer = 0, &
      Order
    character ( LDL ) :: &
      Name
    class ( FieldSet_C_Form ), pointer :: &
      Output_IL_C => null ( ), &
      Output_IR_C => null ( ), & 
      FieldSet_C => null ( )
    class ( Geometry_F_C_Form ), pointer :: &
      Geometry_C => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Compute
    procedure, public, pass :: &
      Show => Show_RC
    final :: &
      Finalize
  end type Reconstruction_C_Form

  type, public :: Reconstruction_C_Element
    !-- Reconstruction_Chart_Element
    class ( Reconstruction_C_Form ), allocatable :: &
      Element
  contains
    final :: &
      Finalize_E
  end type Reconstruction_C_Element

    private :: &
      ComputeConstant_CGS_Kernel, &
      ComputeLinear_CGS_Kernel, &
      ComputeParabolic_CGS_Kernel

    interface
  
      module subroutine ComputeConstant_CGS_Kernel &
               ( F, iaSlctd, iD, oV, F_IL, F_IR, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, :, : ), intent ( in ) :: &
          F
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          iaSlctd
        integer ( KDI ), intent ( in ) :: &
          iD, &
          oV   
        real ( KDR ), dimension ( :, :, :, : ), intent ( out ) :: &
          F_IL, F_IR
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeConstant_CGS_Kernel

      module subroutine ComputeLinear_CGS_Kernel &
               ( F, X, dX, XA, iaSlctd, iD, oV, F_IL, F_IR, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, :, : ), intent ( in ) :: &
          F
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
           X, &
          dX, &
           XA
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          iaSlctd
        integer ( KDI ), intent ( in ) :: &
          iD, &
          oV   
        real ( KDR ), dimension ( :, :, :, : ), intent ( out ) :: &
          F_IL, F_IR
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeLinear_CGS_Kernel

      module subroutine ComputeParabolic_CGS_Kernel &
               ( F, X, dX, XA, X2A, iaSlctd, iD, oV, F_IL, F_IR, &
                 UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, :, : ), intent ( in ) :: &
          F
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
           X, &
          dX, &
           XA, &
           X2A
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          iaSlctd
        integer ( KDI ), intent ( in ) :: &
          iD, &
          oV   
        real ( KDR ), dimension ( :, :, :, : ), intent ( out ) :: &
          F_IL, F_IR
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeParabolic_CGS_Kernel

    end interface

contains


  subroutine Initialize &
               ( RC, GC, FSC, O_IL_C, O_IR_C, NameOption, OrderOption )

    class ( Reconstruction_C_Form ), intent ( inout ) :: &
      RC
    class ( Geometry_F_C_Form ), intent ( in ), target :: &
      GC
    class ( FieldSet_C_Form ), intent ( in ), target :: &
      FSC, &
      O_IL_C, O_IR_C
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), intent ( in ), optional :: &
      OrderOption

    RC % IGNORABILITY  =  FSC % IGNORABILITY

    RC % Name  =  'R_' // trim ( FSC % Name )
    if ( present ( NameOption ) ) &
      RC % Name  =  trim ( NameOption ) // '_' // trim ( FSC % Name )

    call Show ( 'Initializing a Reconstruction_C', RC % IGNORABILITY )
    call Show ( RC % Name, 'Name', RC % IGNORABILITY )
   
    RC % FieldSet_C   =>  FSC
    RC % Output_IL_C  =>  O_IL_C
    RC % Output_IR_C  =>  O_IR_C

    RC % Geometry_C  =>   GC

    RC % Order  =  2
    if ( present ( OrderOption ) ) &
      RC % Order  =  OrderOption
    call PROGRAM_HEADER % GetParameter ( RC % Order, 'Order' )

  end subroutine Initialize


  subroutine Compute ( RC, iD, TimerLevelOption )

    class ( Reconstruction_C_Form ), intent ( inout ) :: &
      RC
    integer ( KDI ), intent ( in ) :: &
      iD  !-- iDimensions
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    real ( KDR ), dimension ( :, :, : ), pointer :: &
       X, &
      dX
    real ( KDR ), dimension ( :, :, :, : ), pointer :: &
      F, &
      F_IL, &
      F_IR
    character ( LDL ) :: &
      TimerName
    type ( TimerForm ), pointer :: &
      T 

    associate ( iT  =>  RC % iTimer )
    if ( iT == 0 ) then
      TimerName  =  RC % Name
      if ( present ( TimerLevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, TimerLevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if
    end associate !-- iT

    T  =>  PROGRAM_HEADER % TimerPointer ( RC % iTimer )
    call T % Start ( )

    call Show ( 'Computing a Reconstruction_C', RC % IGNORABILITY + 4 )
    call Show ( RC % Name, 'Name', RC % IGNORABILITY + 4 )

    associate &
      ( GC     =>  RC % Geometry_C, &
        FC     =>  RC % FieldSet_C, &
        FC_IL  =>  RC % Output_IL_C, &
        FC_IR  =>  RC % Output_IR_C )
    associate &
      ( GV     =>  GC    % Storage_FSC % Storage % Value, &
        FV     =>  FC    % Storage_FSC % Storage % Value, &
        FV_IL  =>  FC_IL % Storage_FSC % Storage % Value, &
        FV_IR  =>  FC_IR % Storage_FSC % Storage % Value, &
        DeviceMemory  =>  FC % Storage_FSC % DeviceMemory )

    select type ( C  =>  RC % FieldSet_C % Chart )
    class is ( Chart_GS_Form )

      call C % SetFieldPointer ( GV ( :, GC % CENTER_U ( iD ) ),  X )
      call C % SetFieldPointer ( GV ( :, GC % WIDTH_U  ( iD ) ), dX )
      call C % SetFieldPointer ( FV,    F    )
      call C % SetFieldPointer ( FV_IL, F_IL )
      call C % SetFieldPointer ( FV_IR, F_IR )

      select case ( RC % Order )
      case ( 0 )
        call ComputeConstant_CGS_Kernel &
               ( F, FC % iaSelected, iD, C % nGhostLayers ( iD ), F_IL, F_IR, &
                 UseDeviceOption = DeviceMemory )
      case ( 1 )
        call ComputeLinear_CGS_Kernel &
               ( F, X, dX, X, FC % iaSelected, iD, C % nGhostLayers ( iD ), &
                 F_IL, F_IR, UseDeviceOption = DeviceMemory )
      case ( 2 )
        call ComputeParabolic_CGS_Kernel &
               ( F, X, dX, X, X ** 2, FC % iaSelected, iD, & 
                 C % nGhostLayers ( iD ), F_IL, F_IR, &
                 UseDeviceOption = DeviceMemory )
      case default
        call Show ( 'Order not implemented', CONSOLE % ERROR )
        call Show ( 'Reconstruction_C__Form', 'module', CONSOLE % ERROR )
        call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- Order

    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Reconstruction_C__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

    end associate !-- FV, etc.
    end associate !-- FC, etc.

    call T % Stop ( )

  end subroutine Compute


  subroutine Show_RC ( RC )

    class ( Reconstruction_C_Form ), intent ( in ) :: &
      RC

    call Show ( 'Reconstruction_C Parameters', RC % IGNORABILITY )

    call Show ( RC % Name, 'Name',  RC % IGNORABILITY )
    call Show ( RC % Order, 'Order', RC % IGNORABILITY )
    call RC % Output_IL_C % Show ( )
    call RC % Output_IR_C % Show ( )

  end subroutine Show_RC


  impure elemental subroutine Finalize ( RC )

    type ( Reconstruction_C_Form ), intent ( inout ) :: &
      RC

    nullify ( RC % Geometry_C )
    nullify ( RC % FieldSet_C )
    nullify ( RC % Output_IR_C )
    nullify ( RC % Output_IL_C )

    call Show ( 'Finalizing a Reconstruction_C', RC % IGNORABILITY )
    call Show ( RC % Name, 'Name', RC % IGNORABILITY )
   
  end subroutine Finalize

  
  impure elemental subroutine Finalize_E ( RE )
    
    type ( Reconstruction_C_Element ), intent ( inout ) :: &
      RE

    if ( allocated ( RE % Element ) ) &
      deallocate ( RE % Element )

  end subroutine Finalize_E


end module Reconstruction_C__Form
