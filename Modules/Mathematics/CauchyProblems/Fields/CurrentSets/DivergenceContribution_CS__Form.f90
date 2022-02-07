module DivergenceContribution_CS__Form

  use Basics
  use FieldSets
  use DivergenceContribution_Form
  use CurrentSet_Form

  implicit none
  private

  type, public, extends ( DivergenceContributionForm ) :: &
    DivergenceContribution_CS_Form
      class ( CurrentSetForm ), pointer :: &
        CurrentSet => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass ( DC ) :: &
      ComputeFluxes
    final :: &
      Finalize
  end type DivergenceContribution_CS_Form

    interface

      module subroutine ComputeFluxesKernel ( D, V_Dim, F_D, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          D, &
          V_Dim
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          F_D
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeFluxesKernel

    end interface


contains


  subroutine Initialize ( DC, CS, NameOption, IgnorabilityOption )

    class ( DivergenceContribution_CS_Form ), intent ( inout ) :: &
      DC
    class ( CurrentSetForm ), intent ( in ), target :: &
      CS
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( DC % Type  ==  '' ) &
      DC % Type  =  'a DivergenceContribution_CS' 

    call DC % Initialize_H ( NameOption, IgnorabilityOption )
 
    DC % CurrentSet  =>  CS

  end subroutine Initialize


  subroutine ComputeFluxes ( FS_F, DC, FS_CS, iC, iD )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS_F  !-- Fluxes
    class ( DivergenceContribution_CS_Form ), intent ( in ) :: &
      DC
    class ( FieldSetForm ), intent ( in ) :: &
      FS_CS
    integer ( KDI ), intent ( in ) :: &
      iC, &  !-- iChart
      iD     !-- iDimension

    integer ( KDI ) :: &
      iDensity

    associate ( CS  =>  DC % CurrentSet )

    if ( CS % DENSITY_CS > 0 ) then

      call Search ( CS % iaBalanced, CS % DENSITY_CS, iDensity )

      associate &
        ( FSS  =>  FS_F  % Storage ( iC ), &
          CSS  =>  FS_CS % Storage ( iC ) )
      associate &
        ( F_D      =>  FSS % Value ( :, iDensity ), &
            D      =>  CSS % Value ( :, CS % DENSITY_CS ), & 
            V_Dim  =>  CSS % Value ( :, CS % VELOCITY_CS_U ( iD ) ) ) 
 
      call ComputeFluxesKernel &
             ( D, V_Dim, F_D, UseDeviceOption = CS % DeviceMemory )
  
      end associate !-- F_D, etc.
      end associate !-- FSS, etc.
  
    end if !-- Density default

    end associate !-- CS

  end subroutine ComputeFluxes


  impure elemental subroutine Finalize ( DC )

    type ( DivergenceContribution_CS_Form ), intent ( inout ) :: &
      DC

  end subroutine Finalize


end module DivergenceContribution_CS__Form
