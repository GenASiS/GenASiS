module Gradient_Form

  use Basics
  use Manifolds
  use FieldSets
  use Geometries

  implicit none
  private

  type, public, extends ( FieldSet_BM_Form ) :: GradientForm
    logical ( KDL ) :: &
      FromDivergence
    class ( FieldSet_BM_Form ), pointer :: &
      FieldSet  => null ( )
    class ( Geometry_F_Form ), pointer :: &
      Geometry => null ( )
    type ( FieldSet_BM_Form ), dimension ( : ), allocatable :: &
      Dimension
  contains
    procedure, private, pass :: &
      InitializeAllocate_G
    generic, public :: &
      Initialize => InitializeAllocate_G
    procedure, public, pass :: &
      SetStream
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type GradientForm

    private :: &
      Compute_CGS_Kernel, &
      ComputeFromDivergence_CGS_Kernel

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

      module subroutine ComputeFromDivergence_CGS_Kernel &
               ( F, A_I, V, iaSlctd, iD, oV, dFdX, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, :, : ), intent ( in ) :: &
          F
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          A_I, &
          V
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          iaSlctd
        integer ( KDI ), intent ( in ) :: &
          iD, &
          oV   
        real ( KDR ), dimension ( :, :, :, : ), intent ( out ) :: &
          dFdX
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeFromDivergence_CGS_Kernel

    end interface


contains


  subroutine InitializeAllocate_G &
               ( G, Gy, FS, NameOption, FromDivergenceOption, &
                 IgnorabilityOption )

    class ( GradientForm ), intent ( inout ) :: &
      G
    class ( Geometry_F_Form ), intent ( in ), target :: &
      Gy
    class ( FieldSet_BM_Form ), intent ( in ), target :: &
      FS
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      FromDivergenceOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
 
    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iF, &  !-- iField
      Ignorability
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    if ( G % Type  ==  '' ) &
      G % Type  =  'a Gradient' 
    
    Name  =  'Grdnt_' // trim ( FS % Name )
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    Ignorability  =  FS % IGNORABILITY  +  1
    if ( present ( IgnorabilityOption ) ) &
      Ignorability  =  IgnorabilityOption

    G % FieldSet  =>  FS
    G % Geometry  =>  Gy

    allocate ( Field ( FS % nFields ) )
    do iS  =  1,  FS % nFields
      iF  =  FS % iaSelected ( iS )
      Field ( iS )  =  FS % Field ( iF )
    end do !-- iS

    call G % FieldSet_BM_Form % Initialize &
           ( FS % Atlas, &
             FieldOption = Field, &
             NameOption = Name, &
             DeviceMemoryOption = FS % DeviceMemory, &
             DevicesCommunicateOption = FS % DevicesCommunicate, &
             nFieldsOption = FS % nFields, &
             IgnorabilityOption = Ignorability )

    G % FromDivergence  =  .false.
    if ( present ( FromDivergenceOption ) ) &
      G % FromDivergence  =  FromDivergenceOption

    call Show ( G % FromDivergence, 'FromDivergence', G % IGNORABILITY )

  end subroutine InitializeAllocate_G


  subroutine SetStream ( G, S )

    class ( GradientForm ), intent ( inout ) :: &
      G
    class ( Stream_BM_Form ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iD     !-- iDimension
    character ( 1 ) :: &
      DimensionNumber

    allocate ( G % Dimension ( 3 ) )
    do iD  =  1, 3

      write ( DimensionNumber, fmt = '(i1.1)' ) iD

      associate ( GD  =>  G % Dimension ( iD ) )
      call GD % Initialize &
             ( G % Atlas, &
               FieldOption = G % Field, &
               NameOption = trim ( G % Name ) // '_' // DimensionNumber, &
               DeviceMemoryOption = G % DeviceMemory, &
               DevicesCommunicateOption = G % DevicesCommunicate, &
               nFieldsOption = G % nFields )
      call S % AddFieldSet ( GD )
      end associate !-- GD, etc.

    end do !-- iD

  end subroutine SetStream


  subroutine Compute ( G, iD )

    class ( GradientForm ), intent ( inout ) :: &
      G
    integer ( KDI ), intent ( in ) :: &
      iD  !-- iDimensions

    integer ( KDI ) :: &
      iC  !-- iChart
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      X, &
      A_I, &
      V
    real ( KDR ), dimension ( :, :, :, : ), pointer :: &
       F, &
      dFdX

    call Show ( 'Computing a Gradient', G % IGNORABILITY + 3 )
    call Show ( G % Name, 'Name', G % IGNORABILITY + 3 )

    associate ( A  =>  G % Atlas )
    do iC  =  1,  A % nCharts

    select type ( C  =>  A % Chart ( iC ) % Element )
      class is ( Chart_GS_Form )
    associate &
      ( Gy  =>  G % Geometry, &
        FS  =>  G % FieldSet )
    associate &
      (  GV  =>  G  % Storage ( iC ) % Value, &
        GyV  =>  Gy % Storage ( iC ) % Value, &
         FV  =>  FS % Storage ( iC ) % Value )

      call C % SetFieldPointer (  FV,  F   )
      call C % SetFieldPointer (  GV, dFdX )

      call G % Storage ( iC ) % ReassociateHost &
            ( AssociateVariablesOption = .false. )
      call FS % Storage ( iC ) % ReassociateHost &
            ( AssociateVariablesOption = .false. )
      
      if ( G % FromDivergence ) then
        call C % SetFieldPointer (  GyV ( :, Gy % AREA_I_D ( iD ) ), A_I )
        call C % SetFieldPointer (  GyV ( :, Gy % VOLUME ), V )
        call ComputeFromDivergence_CGS_Kernel &
               ( F, A_I, V, FS % iaSelected, iD, C % nGhostLayers ( iD ), &
                 dFdX, UseDeviceOption = FS % DeviceMemory )
      else !-- standard central difference
        call C % SetFieldPointer ( GyV ( :, Gy % AVERAGE_1_U ( iD ) ), X )
        call Compute_CGS_Kernel &
               ( F, X, FS % iaSelected, iD, C % nGhostLayers ( iD ), dFdX, &
                 UseDeviceOption = FS % DeviceMemory )
      end if !-- FromDivergence

      call FS % Storage ( iC ) % ReassociateHost &
            ( AssociateVariablesOption = .true. )
      call G % Storage ( iC ) % ReassociateHost &
            ( AssociateVariablesOption = .true. )
      
    end associate !-- GV, etc.
    end associate !-- Gy, etc.

    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Gradient_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

    end do !-- iC
    end associate !-- A

    if ( allocated ( G % Dimension ) ) then
      associate ( GD  =>  G % Dimension ( iD ) )
      call G % Copy ( GD )
      end associate !-- GD, etc.
    end if

  end subroutine Compute


  impure elemental subroutine Finalize ( G )

    type ( GradientForm ), intent ( inout ) :: &
      G

    if ( allocated ( G % Dimension ) ) &
      deallocate ( G % Dimension )

    nullify ( G % Geometry )
    nullify ( G % FieldSet )

  end subroutine Finalize


end module Gradient_Form
