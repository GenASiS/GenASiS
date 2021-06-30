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
    character ( LDL ) :: &
      Name
    class ( FieldSetForm ), pointer :: &
      FieldSet  => null ( )
    class ( FieldSetForm ), allocatable :: &
      Output_IL, &
      Output_IR 
    class ( Geometry_F_Form ), pointer :: &
      Geometry => null ( )
    type ( FieldSetElement ), dimension ( :, : ), allocatable :: &
      StageDimension_IL, &
      StageDimension_IR, &
      StageDimension
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      SetStream
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
               ( R, G, FS, PrefixOption, OrderOption )

    class ( ReconstructionForm ), intent ( inout ) :: &
      R
    class ( Geometry_F_Form ), intent ( in ), target :: &
      G
    class ( FieldSetForm ), intent ( in ), target :: &
      FS
    character ( * ), intent ( in ), optional :: &
      PrefixOption
    integer ( KDI ), intent ( in ), optional :: &
      OrderOption

    integer ( KDI ) :: &
      iC, &  !-- iChart
      iS, &  !-- iSelected
      iF     !-- iField
    type ( MeasuredValueForm ), dimension ( :, : ), allocatable :: &
      Unit
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    R % IGNORABILITY  =  FS % IGNORABILITY + 1

    R % Name  =  'R_' // trim ( FS % Name )
    if ( present ( PrefixOption ) ) &
      R % Name  =  trim ( PrefixOption ) // '_' // trim ( FS % Name )

    call Show ( 'Initializing a Reconstruction', R % IGNORABILITY )
    call Show ( R % Name, 'Name', R % IGNORABILITY )
   
    R % FieldSet  =>  FS
    R % Geometry  =>   G

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
             IgnorabilityOption = R % IGNORABILITY )
    call R % Output_IR % Initialize &
           ( FS % Atlas, &
             FieldOption = Field, &
             NameOption = trim ( R % Name ) // '_IR', &
             DeviceMemoryOption = FS % DeviceMemory, &
             DevicesCommunicateOption = FS % DevicesCommunicate, &
             UnitOption = Unit, &
             nFieldsOption = nF, &
             IgnorabilityOption = R % IGNORABILITY )

    R % Order  =  2
    if ( present ( OrderOption ) ) &
      R % Order  =  OrderOption
    call PROGRAM_HEADER % GetParameter ( R % Order, 'ReconstructionOrder' )

  end associate !-- nC, etc.
  
  end subroutine Initialize


  subroutine SetStream ( R, S, nS )

    class ( ReconstructionForm ), intent ( inout ) :: &
      R
    class ( StreamForm ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      nS  !-- nStages

    integer ( KDI ) :: &
      iS, &  !-- iStage
      iD     !-- iDimension
    character ( 1 ) :: &
      StageNumber, &
      DimensionNumber

    associate ( FS  =>  R % FieldSet )
    associate ( nD  =>  3 )

    allocate ( R % StageDimension ( nS, nD ) )
    allocate ( R % StageDimension_IL ( nS, nD ) )
    allocate ( R % StageDimension_IR ( nS, nD ) )
    do iS  =  1, nS
      do iD  =  1, nD

        write ( StageNumber, fmt = '(i1.1)' ) iS
        write ( DimensionNumber, fmt = '(i1.1)' ) iD

        allocate ( R % StageDimension ( iS, iD ) % Element )
        allocate ( R % StageDimension_IL ( iS, iD ) % Element )
        allocate ( R % StageDimension_IR ( iS, iD ) % Element )

        associate &
          ( SD  =>  R % StageDimension ( iS, iD ) % Element, &
            FS  =>  R % FieldSet )
        call SD % Initialize &
               ( FS % Atlas, &
                 FieldOption = FS % Field, &
                 NameOption = 'R_' // trim ( FS % Name ) // '_' &
                              // StageNumber // '_' // DimensionNumber, &
                 DeviceMemoryOption = FS % DeviceMemory, &
                 DevicesCommunicateOption = FS % DevicesCommunicate, &
                 nFieldsOption = size ( FS % Field ) )
        call S % AddFieldSet ( SD, iaSelectedOption = FS % iaSelected )
        end associate !-- SD, etc.

        associate &
          ( SD  =>  R % StageDimension_IL ( iS, iD ) % Element, &
             O  =>  R % Output_IL )
        call SD % Initialize &
               ( O % Atlas, &
                 FieldOption = O % Field, &
                 NameOption = trim ( O % Name ) // '_' // StageNumber // '_' &
                              // DimensionNumber, &
                 DeviceMemoryOption = FS % DeviceMemory, &
                 DevicesCommunicateOption = FS % DevicesCommunicate, &
                 nFieldsOption = O % nFields )
        call S % AddFieldSet ( SD )
        end associate !-- SD, etc.

        associate &
          ( SD  =>  R % StageDimension_IR ( iS, iD ) % Element, &
             O  =>  R % Output_IR )
        call SD % Initialize &
               ( O % Atlas, &
                 FieldOption = O % Field, &
                 NameOption = trim ( O % Name ) // '_' // StageNumber // '_' &
                              // DimensionNumber, &
                 DeviceMemoryOption = FS % DeviceMemory, &
                 DevicesCommunicateOption = FS % DevicesCommunicate, &
                 nFieldsOption = O % nFields )
        call S % AddFieldSet ( SD )
        end associate !-- SD, etc.

      end do !-- iD
    end do !-- iS

    end associate !-- nD
    end associate !-- FS

  end subroutine SetStream


  subroutine Show_R ( R )

    class ( ReconstructionForm ), intent ( in ) :: &
      R

    call Show ( 'Reconstruction Parameters', R % IGNORABILITY )

    call Show ( R % Name, 'Name',  R % IGNORABILITY )
    call Show ( R % Order, 'Order', R % IGNORABILITY )
!    call R % FieldSet % Show ( )
    call R % Output_IL % Show ( )
    call R % Output_IR % Show ( )

  end subroutine Show_R


  function Timer ( R, LevelOption ) result ( T )

    class ( ReconstructionForm ), intent ( inout ) :: &
      R
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDF ) :: &
      TimerName

    associate ( iT  =>  R % iTimer )

    if ( iT == 0 ) then
      TimerName  =  R % Name
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function Timer


  subroutine Compute ( R, iD, iS_Option )

    class ( ReconstructionForm ), intent ( inout ) :: &
      R
    integer ( KDI ), intent ( in ) :: &
      iD  !-- iDimensions
    integer ( KDI ), intent ( in ), optional :: &
      iS_Option

    integer ( KDI ) :: &
      iC  !-- iChart
    real ( KDR ), dimension ( :, :, : ), pointer :: &
       X, &
      dX
    real ( KDR ), dimension ( :, :, :, : ), pointer :: &
      F, &
      F_IL, &
      F_IR

    call Show ( 'Computing a Reconstruction', R % IGNORABILITY + 4 )
    call Show ( R % Name, 'Name', R % IGNORABILITY + 4 )

    associate &
      ( G      =>  R % Geometry, &
        FS     =>  R % FieldSet, &
        FS_IL  =>  R % Output_IL, &
        FS_IR  =>  R % Output_IR )
    do iC  =  1,  FS % Atlas % nCharts
      associate &
        ( GV     =>  G     % Storage ( iC ) % Value, &
          FV     =>  FS    % Storage ( iC ) % Value, &
          FV_IL  =>  FS_IL % Storage ( iC ) % Value, &
          FV_IR  =>  FS_IR % Storage ( iC ) % Value )

      select type ( C  =>  FS % Atlas % Chart ( iC ) % Element )
      class is ( Chart_GS_Form )

        call C % SetFieldPointer ( GV ( :, G % CENTER_U ( iD ) ),  X )
        call C % SetFieldPointer ( GV ( :, G % WIDTH_U  ( iD ) ), dX )
        call C % SetFieldPointer ( FV,    F    )
        call C % SetFieldPointer ( FV_IL, F_IL )
        call C % SetFieldPointer ( FV_IR, F_IR )

        select case ( R % Order )
        case ( 0 )
          call ComputeConstant_CGS_Kernel &
                 ( F, FS % iaSelected, iD, C % nGhostLayers ( iD ), &
                   F_IL, F_IR, UseDeviceOption = FS % DeviceMemory )
        case ( 1 )
          call ComputeLinear_CGS_Kernel &
                 ( F, X, dX, X, FS % iaSelected, iD, C % nGhostLayers ( iD ), &
                   F_IL, F_IR, UseDeviceOption = FS % DeviceMemory )
        case ( 2 )
          call ComputeParabolic_CGS_Kernel &
                 ( F, X, dX, X, X ** 2, FS % iaSelected, iD, & 
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
    
      end associate !-- GV, etc.
    end do !-- iC
    end associate !-- G, etc.

    if ( allocated ( R % StageDimension ) .and. present ( iS_Option ) ) &
    then
      associate &
        ( SD  =>  R % StageDimension ( iS_Option, iD ) % Element, &
          FS  =>  R % FieldSet )
      call FS % Copy ( SD )
      end associate !-- SD, etc.
    end if

    if ( allocated ( R % StageDimension_IL ) .and. present ( iS_Option ) ) &
    then
      associate &
        ( SD  =>  R % StageDimension_IL ( iS_Option, iD ) % Element, &
           O  =>  R % Output_IL )
      call O % Copy ( SD )
      end associate !-- SD, etc.
    end if

    if ( allocated ( R % StageDimension_IR ) .and. present ( iS_Option ) ) &
    then
      associate &
        ( SD  =>  R % StageDimension_IR ( iS_Option, iD ) % Element, &
           O  =>  R % Output_IR )
      call O % Copy ( SD )
      end associate !-- SD, etc.
    end if

  end subroutine Compute


  impure elemental subroutine Finalize ( R )

    type ( ReconstructionForm ), intent ( inout ) :: &
      R

    if ( allocated ( R % StageDimension ) ) &
      deallocate ( R % StageDimension )
    if ( allocated ( R % StageDimension_IR ) ) &
      deallocate ( R % StageDimension_IR )
    if ( allocated ( R % StageDimension_IL ) ) &
      deallocate ( R % StageDimension_IL )

    if ( allocated ( R % Output_IR ) ) &
      deallocate ( R % Output_IR )
    if ( allocated ( R % Output_IL ) ) &
      deallocate ( R % Output_IL )

    nullify ( R % Geometry )
    nullify ( R % FieldSet )

    call Show ( 'Finalizing a Reconstruction', R % IGNORABILITY )
    call Show ( R % Name, 'Name', R % IGNORABILITY )
   
  end subroutine Finalize

  
end module Reconstruction_Form
