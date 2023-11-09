module DiffusionFactor_RM__Form

  !-- DiffusionFactor_RadiationMoments__Form

  use Basics
  use Mathematics
  use Gravitations
  use Interactions_BM__Form
  use RadiationMoments_BM__Form

  implicit none
  private

  type, public, extends ( DiffusionFactor_CS_Form ) :: DiffusionFactor_RM_Form
  contains
    procedure, private, pass :: &
      InitializeAllocate_DF
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type DiffusionFactor_RM_Form

    private :: &
      ComputeKernel

    interface

      module subroutine ComputeKernel ( SF, TO, M_DD, dX, DF, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          SF, &
          TO, &
          M_DD, &
          dX
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          DF
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeKernel

    end interface

contains


  subroutine InitializeAllocate_DF &
               ( DF, CS, FieldOption, nFieldsOption )

    class ( DiffusionFactor_RM_Form ), intent ( inout ) :: &
      DF
    class ( CurrentSetForm ), intent ( in ), target :: &
      CS
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption

    if ( DF % Type  ==  '' ) &
      DF % Type  =  'a DiffusionFactor_RM' 
    
    call DF % DiffusionFactor_CS_Form % Initialize &
           ( CS, FieldOption, nFieldsOption )

  end subroutine InitializeAllocate_DF


  subroutine Compute ( DF, iC, iD )

    class ( DiffusionFactor_RM_Form ), intent ( inout ) :: &
      DF
    integer ( KDI ), intent ( in ) :: &
      iC, &  !-- iChart
      iD     !-- iDimensions

    call Show ( 'Computing ' // trim ( DF % Type ), DF % IGNORABILITY + 3 )
    call Show ( DF % Name, 'Name', DF % IGNORABILITY + 3 )

    select type ( RM  =>  DF % CurrentSet )
      class is ( RadiationMoments_BM_Form )
    select type ( I  =>  RM % Interactions )
      class is ( Interactions_BM_Form )
    associate &
      ( DFV  =>  DF % Storage ( iC ) % Value, &
        RMV  =>  RM % Storage ( iC ) % Value, &
         IV  =>   I % Storage ( iC ) % Value )
    associate &
      ( DFR  =>  DFV ( :, DF % DIFFUSION_FACTOR ), &
         SF  =>  RMV ( :, RM % STRESS_FACTOR ), &
         TO  =>   IV ( :,  I % OPACITY_H ) )

    select type ( G  =>  RM % Geometry )
    class is ( Gravitation_G_Form )

        associate &
          ( GSV  =>  G % Storage ( iC ) % Value )
        associate &
          ( dX     =>  GSV ( :, G % WIDTH_U ( iD ) ), &
             M_DD  =>  GSV ( :, G % METRIC_F_DD ( iD ) ) )

       call ComputeKernel &
              ( SF, TO, M_DD, dX, DFR, UseDeviceOption = DF % DeviceMemory )

        end associate !-- dX, etc.
        end associate !-- GSV

    class default
      call Show ( 'Gravitation type not recognized', CONSOLE % ERROR )
      call Show ( 'DiffusionFactor_RM__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- G

    end associate !-- DFR, etc.
    end associate !-- DFV, etc.
    end select !-- I
    end select !-- RM

  end subroutine Compute


  impure elemental subroutine Finalize ( DF )

    type ( DiffusionFactor_RM_Form ), intent ( inout ) :: &
      DF

  end subroutine Finalize


end module DiffusionFactor_RM__Form
