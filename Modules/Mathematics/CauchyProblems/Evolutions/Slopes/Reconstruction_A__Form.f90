module Reconstruction_A__Form

  !-- Reconstruction_Atlas_Form

  use Basics
  use Manifolds
  use Fields
  use Reconstruction_C__Form

  implicit none
  private

  type, public :: Reconstruction_A_Form
    integer ( KDI ) :: &
      IGNORABILITY
    character ( LDL ) :: &
      Name
    class ( FieldSet_A_Form ), allocatable :: &
      Output_IL_A, &
      Output_IR_A       
    class ( FieldSet_A_Form ), pointer :: &
      FieldSet_A => null ( )
    class ( Geometry_F_A_Form ), pointer :: &
      Geometry_A => null ( )
    class ( Reconstruction_C_Element ), dimension ( : ), allocatable :: &
      Reconstruction_C
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Compute
    procedure, public, pass :: &
      Show => Show_RA
    final :: &
      Finalize
  end type Reconstruction_A_Form


contains


  subroutine Initialize ( RA, GA, FSA, NameOption, OrderOption )

    class ( Reconstruction_A_Form ), intent ( inout ) :: &
      RA
    class ( Geometry_F_A_Form ), intent ( in ), target :: &
      GA
    class ( FieldSet_A_Form ), intent ( in ), target :: &
      FSA
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), intent ( in ), optional :: &
      OrderOption

    integer ( KDI ) :: &
      iC, &  !-- iChart
      iS, &  !-- iSelected
      iF     !-- iField
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      Unit
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    RA % IGNORABILITY  =  FSA % IGNORABILITY

    RA % Name  =  'Reconstruction_' // trim ( FSA % Name )
    if ( present ( NameOption ) ) &
      RA % Name  =  trim ( NameOption ) // '_' // trim ( FSA % Name )

    call Show ( 'Initializing a Reconstruction_C', RA % IGNORABILITY )
    call Show ( RA % Name, 'Name', RA % IGNORABILITY )
   
    RA % FieldSet_A  =>  FSA
    RA % Geometry_A  =>   GA

    associate &
      ( FSC  =>  FSA % FieldSet_C ( 1 ) % Element )
    associate &
      ( nF  =>  FSC % nFields, &
        DeviceMemory  =>  FSC % Storage_FSC % DeviceMemory, &
        PinnedMemory  =>  FSC % Storage_FSC % PinnedMemory, &
        DevicesCommunicate  =>  FSC % GhostExchange_FSC % DevicesCommunicate ) 

    allocate ( Field ( nF ) )
    allocate ( Unit ( nF ) )
    do iS  =  1,  nF
      iF  =  FSC % iaSelected ( iS )
      Field ( iS )  =  FSC % Field ( iF )
      Unit  ( iS )  =  FSC % Unit  ( iF )
    end do !-- iS

    allocate ( RA % Output_IL_A )
    allocate ( RA % Output_IR_A )
    call RA % Output_IL_A % Initialize &
           ( FSA % Atlas, &
             FieldOption = Field, &
             NameOption = trim ( RA % Name ) // '_IL', &
             DeviceMemoryOption = DeviceMemory, &
             PinnedMemoryOption = PinnedMemory, &
             DevicesCommunicateOption = DevicesCommunicate, &
             UnitOption = Unit, &
             nFieldsOption = nF, &
             IgnorabilityOption = FSC % Ignorability )
    call RA % Output_IR_A % Initialize &
           ( FSA % Atlas, &
             FieldOption = Field, &
             NameOption = trim ( RA % Name ) // '_IR', &
             DeviceMemoryOption = DeviceMemory, &
             PinnedMemoryOption = PinnedMemory, &
             DevicesCommunicateOption = DevicesCommunicate, &
             UnitOption = Unit, &
             nFieldsOption = nF, &
             IgnorabilityOption = FSC % IGNORABILITY )

    end associate !-- nF, etc.
    end associate !-- FSC

    associate ( nC  =>  FSA % Atlas % nCharts )
    allocate ( RA % Reconstruction_C ( nC ) )
    do iC  =  1, nC
      allocate ( RA % Reconstruction_C ( iC ) % Element )
      associate &
        (     RC  =>  RA % Reconstruction_C ( iC ) % Element, &
             FSC  =>  RA % FieldSet_A % FieldSet_C ( iC ) % Element, &
          O_IL_C  =>  RA % Output_IL_A % FieldSet_C ( iC ) % Element, &
          O_IR_C  =>  RA % Output_IR_A % FieldSet_C ( iC ) % Element )
      select type ( GC  =>  GA % FieldSet_C ( iC ) % Element )
      class is ( Geometry_F_C_Form )
        call RC % Initialize &
               ( GC, FSC, O_IL_C, O_IR_C, NameOption, OrderOption )
      end select !-- GC
      end associate !-- RC, etc.
    end do !-- iC
    end associate !-- nC

  end subroutine Initialize


  subroutine Compute ( RA, iD, TimerLevelOption )

    class ( Reconstruction_A_Form ), intent ( inout ) :: &
      RA
    integer ( KDI ), intent ( in ) :: &
      iD  !-- iDimensions
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

   integer ( KDI ) :: &
     iC  !-- iChart

    do iC  =  1, size ( RA % Reconstruction_C )
      associate ( RC  =>  RA % Reconstruction_C ( iC ) % Element )
      call RC % Compute ( iD, TimerLevelOption )
      end associate !-- RC
    end do !-- iC

  end subroutine Compute


  subroutine Show_RA ( RA )

    class ( Reconstruction_A_Form ), intent ( in ) :: &
      RA

   integer ( KDI ) :: &
     iC  !-- iChart

    call Show ( 'Reconstruction_A Parameters', RA % IGNORABILITY )

    call Show ( RA % Name, 'Name',  RA % IGNORABILITY )
    call Show ( RA % FieldSet_A % Atlas % Name, 'Atlas', RA % IGNORABILITY )
    do iC  =  1, size ( RA % Reconstruction_C )
      associate ( RC  =>  RA % Reconstruction_C ( iC ) % Element )
      call RC % Show ( )
      end associate !-- RC
    end do !-- iC

  end subroutine Show_RA


  impure elemental subroutine Finalize ( RA )

    type ( Reconstruction_A_Form ), intent ( inout ) :: &
      RA

    if ( allocated ( RA % Reconstruction_C ) ) &
      deallocate ( RA % Reconstruction_C )

    nullify ( RA % Geometry_A )
    nullify ( RA % FieldSet_A )

    if ( allocated ( RA % Output_IR_A ) ) &
      deallocate ( RA % Output_IR_A )
    if ( allocated ( RA % Output_IL_A ) ) &
      deallocate ( RA % Output_IL_A )

    call Show ( 'Finalizing a Reconstruction_A', RA % IGNORABILITY )
    call Show ( RA % Name, 'Name', RA % IGNORABILITY )
   
  end subroutine Finalize

  
end module Reconstruction_A__Form
