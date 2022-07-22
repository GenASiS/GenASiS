module Reconstruction_Form

  use Basics
  use Manifolds
  use Fields

  implicit none
  private

  type, public :: ReconstructionForm
    integer ( KDI ) :: &
      IGNORABILITY, &
      Order
    integer ( KDI ) :: &
      iTimer = 0
    integer ( KDL ), dimension ( : ), allocatable :: &
      iaSelected
    logical ( KDL ) :: &
      AllocatedOutput
    character ( LDL ) :: &
      Name
    class ( FieldSetForm ), pointer :: &
      FieldSet  => null ( )
    class ( FieldSetForm ), pointer :: &
      Output_IL => null ( ), &
      Output_IR => null ( )
    class ( Geometry_F_Form ), pointer :: &
      Geometry => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate
    procedure, private, pass :: &
      InitializeAssociate
    generic, public :: &
      Initialize => InitializeAllocate, InitializeAssociate
    procedure, public, pass :: &
      Show => Show_R
    procedure, public, pass :: &
      Timer
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type ReconstructionForm

    private :: &
      ComputeConstant_CGS_Kernel, &
      ComputeLinear_CGS_Kernel, &
      ComputeParabolic_CGS_Kernel

    interface
  
      module subroutine ComputeConstant_CGS_Kernel &
               ( F, iaSlctd, iaSlctd_R, iD, oV, F_IL, F_IR, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, :, : ), intent ( in ) :: &
          F
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          iaSlctd, &
          iaSlctd_R
        integer ( KDI ), intent ( in ) :: &
          iD, &
          oV   
        real ( KDR ), dimension ( :, :, :, : ), intent ( out ) :: &
          F_IL, F_IR
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeConstant_CGS_Kernel

      module subroutine ComputeLinear_CGS_Kernel &
               ( F, X, dX, XA, iaSlctd, iaSlctd_R, iD, oV, F_IL, F_IR, &
                 UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, :, : ), intent ( in ) :: &
          F
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
           X, &
          dX, &
           XA
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          iaSlctd, &
          iaSlctd_R
        integer ( KDI ), intent ( in ) :: &
          iD, &
          oV   
        real ( KDR ), dimension ( :, :, :, : ), intent ( out ) :: &
          F_IL, F_IR
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeLinear_CGS_Kernel

      module subroutine ComputeParabolic_CGS_Kernel &
               ( F, X, dX, XA, X2A, iaSlctd, iaSlctd_R, iD, oV, F_IL, F_IR, &
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
          iaSlctd, &
          iaSlctd_R
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


  subroutine InitializeAllocate ( R, G, FS, OrderOption )

    class ( ReconstructionForm ), intent ( inout ) :: &
      R
    class ( Geometry_F_Form ), intent ( in ), target :: &
      G
    class ( FieldSetForm ), intent ( in ), target :: &
      FS
    integer ( KDI ), intent ( in ), optional :: &
      OrderOption

    integer ( KDI ) :: &
      iC, &  !-- iChart
      iS, &  !-- iSelected
      iF     !-- iField
    type ( QuantityForm ), dimension ( :, : ), allocatable :: &
      Unit
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    R % IGNORABILITY  =  FS % Atlas % IGNORABILITY

    R % Name  =  trim ( FS % Name ) // '_Rcnstrctn'

    call Show ( 'InitializeAllocating a Reconstruction', R % IGNORABILITY )
    call Show ( R % Name, 'Name', R % IGNORABILITY )
   
    R % FieldSet  =>  FS
    R % Geometry  =>   G

    R % AllocatedOutput  =  .true.

    associate &
      ( nC  =>  FS % Atlas % nCharts, &
        nF  =>  FS % nFields )

    allocate ( Field ( nF ) )
    allocate ( Unit ( nF, nC ) )
    do iS  =  1,  nF
      iF  =  FS % iaSelected ( iS )
      Field ( iS )  =  FS % Field ( iF )
      do iC  =  1,  nC
        Unit  ( iS, iC )  =  FS % Unit ( iF, iC )
      end do !-- iC
    end do !-- iS

    allocate ( R % Output_IL )
    allocate ( R % Output_IR )
    call R % Output_IL % Initialize &
           ( FS % Atlas, &
             FieldOption = Field, &
             NameOption = trim ( R % Name ) // '_IL', &
             DeviceMemoryOption = FS % DeviceMemory, &
             DevicesCommunicateOption = FS % DevicesCommunicate, &
             UnitOption = Unit, &
             nFieldsOption = nF, &
             IgnorabilityOption = R % IGNORABILITY + 1 )
    call R % Output_IR % Initialize &
           ( FS % Atlas, &
             FieldOption = Field, &
             NameOption = trim ( R % Name ) // '_IR', &
             DeviceMemoryOption = FS % DeviceMemory, &
             DevicesCommunicateOption = FS % DevicesCommunicate, &
             UnitOption = Unit, &
             nFieldsOption = nF, &
             IgnorabilityOption = R % IGNORABILITY + 1 )

    R % Order  =  2
    if ( present ( OrderOption ) ) &
      R % Order  =  OrderOption
    call PROGRAM_HEADER % GetParameter ( R % Order, 'ReconstructionOrder' )

    allocate ( R % iaSelected, source = R % Output_IL % iaSelected )
 
    end associate !-- nC, etc.
  
  end subroutine InitializeAllocate


  subroutine InitializeAssociate ( R, G, FS, O_IL, O_IR, iaS, OrderOption )

    class ( ReconstructionForm ), intent ( inout ) :: &
      R
    class ( Geometry_F_Form ), intent ( in ), target :: &
      G
    class ( FieldSetForm ), intent ( in ), target :: &
      FS, &
      O_IL, O_IR
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iaS
    integer ( KDI ), intent ( in ), optional :: &
      OrderOption

    R % IGNORABILITY  =  FS % Atlas % IGNORABILITY

    R % Name  =  trim ( FS % Name ) // '_Rcnstrctn' 

    call Show ( 'InitializeAssociating a Reconstruction', R % IGNORABILITY )
    call Show ( R % Name, 'Name', R % IGNORABILITY )
   
    R % FieldSet  =>  FS
    R % Geometry  =>   G

    R % AllocatedOutput  =  .false.

    R % Output_IL  =>  O_IL
    R % Output_IR  =>  O_IR

    R % Order  =  2
    if ( present ( OrderOption ) ) &
      R % Order  =  OrderOption
    call PROGRAM_HEADER % GetParameter ( R % Order, 'ReconstructionOrder' )

    allocate ( R % iaSelected, source = iaS )
 
  end subroutine InitializeAssociate


  subroutine Show_R ( R )

    class ( ReconstructionForm ), intent ( in ) :: &
      R

    call Show ( 'Reconstruction Parameters', R % IGNORABILITY )

    call Show ( R % Name,       'Name',  R % IGNORABILITY )
    call Show ( R % Order,      'Order', R % IGNORABILITY )
    call Show ( R % iaSelected, 'iaSelected', R % IGNORABILITY + 1 )

    call R % Output_IL % Show ( )
    call R % Output_IR % Show ( )

  end subroutine Show_R


  function Timer ( R, Level ) result ( T )

    class ( ReconstructionForm ), intent ( inout ) :: &
      R
    integer ( KDI ), intent ( in ) :: &
      Level
    type ( TimerForm ), pointer :: &
      T

    T  =>  PROGRAM_HEADER % Timer &
             ( Handle = R % iTimer, &
               Name = trim ( R % Name ), &
               Level = Level )

  end function Timer


  subroutine Compute ( R, iC, iD )

    class ( ReconstructionForm ), intent ( inout ) :: &
      R
    integer ( KDI ), intent ( in ) :: &
      iC, &  !-- iChart
      iD     !-- iDimensions

    real ( KDR ), dimension ( :, :, : ), pointer :: &
       X, &
      dX, &
       XA, &
       X2A
    real ( KDR ), dimension ( :, :, :, : ), pointer :: &
      F, &
      F_IL, &
      F_IR

    call Show ( 'Computing a Reconstruction', R % IGNORABILITY + 3 )
    call Show ( R % Name, 'Name', R % IGNORABILITY + 3 )

    associate &
      ( G      =>  R % Geometry, &
        FS     =>  R % FieldSet, &
        FS_IL  =>  R % Output_IL, &
        FS_IR  =>  R % Output_IR )
    associate &
      ( GV     =>  G     % Storage ( iC ) % Value, &
        FV     =>  FS    % Storage ( iC ) % Value, &
        FV_IL  =>  FS_IL % Storage ( iC ) % Value, &
        FV_IR  =>  FS_IR % Storage ( iC ) % Value )
        
    call FS % Storage ( iC ) % ReassociateHost &
            ( AssociateVariablesOption = .false. )
    call FS_IL % Storage ( iC ) % ReassociateHost &
            ( AssociateVariablesOption = .false. )
    call FS_IR % Storage ( iC ) % ReassociateHost &
            ( AssociateVariablesOption = .false. )

    select type ( C  =>  FS % Atlas % Chart ( iC ) % Element )
    class is ( Chart_GS_Form )

      call C % SetFieldPointer ( GV ( :, G % CENTER_U ( iD ) ),    X )
      call C % SetFieldPointer ( GV ( :, G % WIDTH_U  ( iD ) ),   dX )
      call C % SetFieldPointer ( GV ( :, G % AVERAGE_1_U ( iD ) ), XA )
      call C % SetFieldPointer ( GV ( :, G % AVERAGE_2_U ( iD ) ), X2A )
      call C % SetFieldPointer ( FV,    F    )
      call C % SetFieldPointer ( FV_IL, F_IL )
      call C % SetFieldPointer ( FV_IR, F_IR )

      select case ( R % Order )
      case ( 0 )
        call ComputeConstant_CGS_Kernel &
               ( F, FS % iaSelected, R % iaSelected, iD, &
                 C % nGhostLayers ( iD ), F_IL, F_IR, &
                 UseDeviceOption = FS % DeviceMemory )
      case ( 1 )
        call ComputeLinear_CGS_Kernel &
               ( F, X, dX, XA, FS % iaSelected, R % iaSelected, iD, &
                 C % nGhostLayers ( iD ), F_IL, F_IR, &
                 UseDeviceOption = FS % DeviceMemory )
      case ( 2 )
        call ComputeParabolic_CGS_Kernel &
               ( F, X, dX, XA, X2A, FS % iaSelected, R % iaSelected, iD, & 
                 C % nGhostLayers ( iD ), F_IL, F_IR, &
                 UseDeviceOption = FS % DeviceMemory )
      case default
        call Show ( 'Order not implemented', CONSOLE % ERROR )
        call Show ( 'Reconstruction_Form', 'module', CONSOLE % ERROR )
        call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- Order

    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Reconstruction_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C
    
    call FS % Storage ( iC ) % ReassociateHost &
            ( AssociateVariablesOption = .true. )
    call FS_IL % Storage ( iC ) % ReassociateHost &
            ( AssociateVariablesOption = .true. )
    call FS_IR % Storage ( iC ) % ReassociateHost &
            ( AssociateVariablesOption = .true. )

    end associate !-- GV, etc.
    end associate !-- G, etc.

  end subroutine Compute


  impure elemental subroutine Finalize ( R )

    type ( ReconstructionForm ), intent ( inout ) :: &
      R

    if ( R % AllocatedOutput ) then
      deallocate ( R % Output_IR )
      deallocate ( R % Output_IL )
    else
      nullify ( R % Output_IR )
      nullify ( R % Output_IL )
    end if

    nullify ( R % Geometry )
    nullify ( R % FieldSet )

    if ( allocated ( R % iaSelected ) ) &
      deallocate ( R % iaSelected )

    call Show ( 'Finalizing a Reconstruction', R % IGNORABILITY )
    call Show ( R % Name, 'Name', R % IGNORABILITY )
   
  end subroutine Finalize

  
end module Reconstruction_Form
