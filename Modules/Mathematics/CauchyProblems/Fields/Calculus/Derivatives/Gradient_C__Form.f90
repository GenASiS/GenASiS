module Gradient_C__Form

  !-- Gradient_Chart_Form

  use Basics
  use Manifolds
  use FieldSets
  use Streams
  use Geometries

  implicit none
  private

  type, public, extends ( FieldSet_C_Form ) :: Gradient_C_Form
    class ( FieldSet_C_Form ), pointer :: &
      FieldSet_C  => null ( )
    class ( Geometry_F_C_Form ), pointer :: &
      Geometry_C => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_G
    generic, public :: &
      Initialize => InitializeAllocate_G
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Gradient_C_Form

    private :: &
      Compute_CGS_Kernel

    interface

      module subroutine Compute_CGS_Kernel &
               ( F, XA, iaSlctd, iD, oV, dFdX, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, :, : ), intent ( in ) :: &
          F
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          XA
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          iaSlctd
        integer ( KDI ), intent ( in ) :: &
          iD, &
          oV   
        real ( KDR ), dimension ( :, :, :, : ), intent ( out ) :: &
          dFdX
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_CGS_Kernel

    end interface


contains


  subroutine InitializeAllocate_G ( GC, GyC, FSC, NameOption )

    class ( Gradient_C_Form ), intent ( inout ) :: &
      GC
    class ( Geometry_F_C_Form ), intent ( in ), target :: &
      GyC
    class ( FieldSet_C_Form ), intent ( in ), target :: &
      FSC
    character ( * ), intent ( in ), optional :: &
      NameOption
 
    integer ( KDI ) :: &
      iS, &    !-- iSelected
      iF, &    !-- iField
      nFields
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    if ( GC % Type  ==  '' ) &
      GC % Type  =  'a Gradient_C' 
    
    Name  =  'Grad_' // trim ( FSC % Name )
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    GC % FieldSet_C  =>  FSC
    GC % Geometry_C  =>  GyC

    associate &
      ( DeviceMemory  =>  FSC % Storage_FSC % DeviceMemory, &
        PinnedMemory  =>  FSC % Storage_FSC % PinnedMemory, &
        DevicesCommunicate  =>  FSC % GhostExchange_FSC % DevicesCommunicate ) 

    nFields  =  FSC % nFields
    
    allocate ( Field ( nFields ) )
    do iS  =  1, nFields
      iF  =  FSC % iaSelected ( iS )
      Field ( iS )  =  FSC % Field ( iF )
    end do !-- iS

    call GC % FieldSet_C_Form % Initialize &
           ( FSC % Chart, &
             FieldOption = Field, &
             NameOption = Name, &
             DeviceMemoryOption = DeviceMemory, &
             PinnedMemoryOption = PinnedMemory, &
             DevicesCommunicateOption = DevicesCommunicate, &
             nFieldsOption = nFields, &
             IgnorabilityOption = FSC % IGNORABILITY )

    end associate !-- DeviceMemory, etc.

  end subroutine InitializeAllocate_G


  subroutine Compute ( GC, iD, TimerLevelOption )

    class ( Gradient_C_Form ), intent ( inout ) :: &
      GC
    integer ( KDI ), intent ( in ) :: &
      iD  !-- iDimensions
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    real ( KDR ), dimension ( :, :, : ), pointer :: &
       X
    real ( KDR ), dimension ( :, :, :, : ), pointer :: &
       F, &
      dFdX
    character ( LDL ) :: &
      TimerName
    type ( TimerForm ), pointer :: &
      T 

    ! associate ( iT  =>  GC % iTimer )
    ! if ( iT == 0 ) then
    !   TimerName  =  GC % Name
    !   if ( present ( TimerLevelOption ) ) then
    !     call PROGRAM_HEADER % AddTimer ( TimerName, iT, TimerLevelOption )
    !   else
    !     call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
    !   end if
    ! end if
    ! end associate !-- iT

    ! T  =>  PROGRAM_HEADER % TimerPointer ( GC % iTimer )
    ! call T % Start ( )

    call Show ( 'Computing a Reconstruction_C', GC % IGNORABILITY + 4 )
    call Show ( GC % Name, 'Name', GC % IGNORABILITY + 4 )

    associate &
      ( GyC  =>  GC % Geometry_C, &
         FC  =>  GC % FieldSet_C )
    associate &
      (  GV  =>   GC % Storage_FSC % Storage % Value, &
        GyV  =>  GyC % Storage_FSC % Storage % Value, &
         FV  =>   FC % Storage_FSC % Storage % Value, &
        DeviceMemory  =>  FC % Storage_FSC % DeviceMemory )

    select type ( C  =>  GC % FieldSet_C % Chart )
    class is ( Chart_GS_Form )

      call C % SetFieldPointer ( GyV ( :, GyC % CENTER_U ( iD ) ),  X )
      call C % SetFieldPointer (  FV,  F   )
      call C % SetFieldPointer (  GV, dFdX )

      call Compute_CGS_Kernel &
             ( F, X, FC % iaSelected, iD, C % nGhostLayers ( iD ), dFdX, &
               UseDeviceOption = DeviceMemory )

    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Gradient_C__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

    end associate !-- GV, etc.
    end associate !-- GyC, etc.

    ! if ( allocated ( GC % StageDimension_C ) .and. present ( iS_Option ) ) &
    ! then
    !   associate &
    !     ( SDC  =>  GC % StageDimension_C ( iS_Option, iD ) % Element, &
    !       FSC  =>  GC % FieldSet_C )
    !   call FSC % Copy ( SDC )
    !   end associate !-- SDC, etc.
    ! end if

    ! if ( allocated ( GC % StageDimension_IL_C ) .and. present ( iS_Option ) ) &
    ! then
    !   associate &
    !     ( SDC  =>  GC % StageDimension_IL_C ( iS_Option, iD ) % Element, &
    !        OC  =>  GC % Output_IL_C )
    !   call OC % Copy ( SDC )
    !   end associate !-- SDC, etc.
    ! end if

    ! if ( allocated ( GC % StageDimension_IR_C ) .and. present ( iS_Option ) ) &
    ! then
    !   associate &
    !     ( SDC  =>  GC % StageDimension_IR_C ( iS_Option, iD ) % Element, &
    !        OC  =>  GC % Output_IR_C )
    !   call OC % Copy ( SDC )
    !   end associate !-- SDC, etc.
    ! end if

    ! call T % Stop ( )

  end subroutine Compute


  impure elemental subroutine Finalize ( GC )

    type ( Gradient_C_Form ), intent ( inout ) :: &
      GC

    nullify ( GC % Geometry_C )
    nullify ( GC % FieldSet_C )

  end subroutine Finalize


end module Gradient_C__Form
