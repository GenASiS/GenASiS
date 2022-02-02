module FluxSet_CS__Form

  use Basics
  use FieldSets
  use FluxSet_Form
  use CurrentSet_Form

  implicit none
  private

  type, public, extends ( FluxSetForm ) :: FluxSet_CS_Form
    class ( CurrentSetForm ), pointer :: &
      CurrentSet => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass ( FS ) :: &
      Compute
    final :: &
      Finalize
  end type FluxSet_CS_Form

    interface

      module subroutine ComputeKernel ( D, V_Dim, F_D, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          D, &
          V_Dim
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          F_D
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeKernel

    end interface


contains


  subroutine Initialize ( FS, CS, NameOption, IgnorabilityOption )

    class ( FluxSet_CS_Form ), intent ( inout ) :: &
      FS
    class ( CurrentSetForm ), intent ( in ), target :: &
      CS
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    call FS % Initialize_H ( NameOption, IgnorabilityOption )
 
    FS % CurrentSet  =>  CS

  end subroutine Initialize


  subroutine Compute ( FS_FS, FS_SS, FS, FS_CS, iC, iD )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS_FS, &  !-- FluxSet
      FS_SS     !-- StressSet
    class ( FluxSet_CS_Form ), intent ( in ) :: &
      FS
    class ( FieldSetForm ), intent ( in ) :: &
      FS_CS
    integer ( KDI ), intent ( in ) :: &
      iC, &  !-- iChart
      iD     !-- iDimension

    integer ( KDI ) :: &
      iDensity

    associate ( CS  =>  FS % CurrentSet )

    if ( CS % DENSITY_CS > 0 ) then

      call Search ( CS % iaBalanced, CS % DENSITY_CS, iDensity )

      associate &
        ( FSS  =>  FS_FS % Storage ( iC ), &
          CSS  =>  FS_CS % Storage ( iC ) )
      associate &
        ( F_D      =>  FSS % Value ( :, iDensity ), &
            D      =>  CSS % Value ( :, CS % DENSITY_CS ), & 
            V_Dim  =>  CSS % Value ( :, CS % VELOCITY_CS_U ( iD ) ) ) 
 
      call ComputeKernel &
             ( D, V_Dim, F_D, UseDeviceOption = CS % DeviceMemory )
  
      end associate !-- F_D, etc.
      end associate !-- FSS, etc.
  
    end if !-- Density default

    end associate !-- CS

  end subroutine Compute


  impure elemental subroutine Finalize ( FS )

    type ( FluxSet_CS_Form ), intent ( inout ) :: &
      FS

  end subroutine Finalize


end module FluxSet_CS__Form
