module DivergenceContribution_CS__Form

  use Basics
  use FieldSets
  use CurrentSet_Form

  implicit none
  private

  type, public :: DivergenceContribution_CS_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      iTimer = 0
    character ( LDL ) :: &
      Type = '', &
      Name
    class ( CurrentSetForm ), pointer :: &
      CurrentSet => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Show => Show_DC
    procedure, public, pass :: &
      Timer
    procedure, public, pass ( DC ) :: &
      ComputeFluxes
    procedure, public, pass ( DC ) :: &
      ComputeStresses
    final :: &
      Finalize
  end type DivergenceContribution_CS_Form

  type, public :: DivergenceContributionElement
    class ( DivergenceContribution_CS_Form ), allocatable :: &
      Element
  contains
    final :: &
      Finalize_E
  end type DivergenceContributionElement

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

    DC % IGNORABILITY  =  CS % IGNORABILITY
    if ( present ( IgnorabilityOption ) ) &
      DC % IGNORABILITY  =  IgnorabilityOption

    if ( DC % Type  ==  '' ) &
      DC % Type  =  'a DivergenceContribution_CS' 

    DC % Name  =  'DivergenceContribution'
    if ( present ( NameOption ) ) &
      DC % Name  =  NameOption

    call Show ( 'Initializing ' // trim ( DC % Type ), DC % IGNORABILITY )
    call Show ( DC % Name, 'Name', DC % IGNORABILITY )

    DC % CurrentSet  =>  CS

  end subroutine Initialize


  subroutine Show_DC ( DC )

    class ( DivergenceContribution_CS_Form ), intent ( in ) :: &
      DC

    character ( LDL ), dimension ( : ), allocatable :: &
      TypeWord

    call Split ( DC % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', DC % IGNORABILITY )

    call Show ( DC % Name, 'Name', DC % IGNORABILITY )
    call Show ( DC % CurrentSet % Name, 'CurrentSet', DC % IGNORABILITY )

  end subroutine Show_DC


  function Timer ( DC, LevelOption ) result ( T )

    class ( DivergenceContribution_CS_Form ), intent ( inout ) :: &
      DC
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDL ) :: &
      TimerName

    associate ( iT  =>  DC % iTimer )

    if ( iT == 0 ) then
      TimerName  =  DC % Name
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function Timer


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


  subroutine ComputeStresses ( S_UD, DC, iC, iMomentum_1, iMomentum_2 )

    class ( FieldSetForm ), intent ( inout ) :: &
      S_UD
    class ( DivergenceContribution_CS_Form ), intent ( in ) :: &
      DC
    integer ( KDI ), intent ( in ) :: &
      iC  !-- iChart
    integer ( KDI ), intent ( out ) :: &
      iMomentum_1, iMomentum_2

    call Show ( 'ComputeStresses should be overridden', CONSOLE % WARNING )
    call Show ( 'DivergenceContribution_CS__Form', 'module', CONSOLE % WARNING )

    call S_UD % Clear ( )

  end subroutine ComputeStresses


  impure elemental subroutine Finalize ( DC )

    type ( DivergenceContribution_CS_Form ), intent ( inout ) :: &
      DC

    call Show ( 'Finalizing ' // trim ( DC % Type ), DC % IGNORABILITY )
    call Show ( DC % Name, 'Name', DC % IGNORABILITY )
   
  end subroutine Finalize


  impure elemental subroutine Finalize_E ( DCE )
    
    type ( DivergenceContributionElement ), intent ( inout ) :: &
      DCE

    if ( allocated ( DCE % Element ) ) &
      deallocate ( DCE % Element )

  end subroutine Finalize_E


end module DivergenceContribution_CS__Form
