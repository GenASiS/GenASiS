module DivergencePart_CS__Form

  use Basics
  use FieldSets
  use CurrentSet_Form

  implicit none
  private

  type, public :: DivergencePart_CS_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      iTimer_F     = 0, &
      iTimer_S     = 0
    character ( LDL ) :: &
      Type = '', &
      Name
    class ( CurrentSetForm ), pointer :: &
      CurrentSet => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Show => Show_DP
    procedure, public, pass :: &
      Timer_F
    procedure, public, pass :: &
      Timer_S
    procedure, public, pass ( DP ) :: &
      ComputeFluxes
    procedure, public, pass ( DP ) :: &
      ComputeStresses
    final :: &
      Finalize
  end type DivergencePart_CS_Form

  type, public :: DivergencePartElement
    class ( DivergencePart_CS_Form ), allocatable :: &
      Element
  contains
    final :: &
      Finalize_E
  end type DivergencePartElement

    private :: &
      ComputeFluxesKernel

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


  subroutine Initialize ( DP, CS, NameOption, IgnorabilityOption )

    class ( DivergencePart_CS_Form ), intent ( inout ) :: &
      DP
    class ( CurrentSetForm ), intent ( in ), target :: &
      CS
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    DP % IGNORABILITY  =  CS % IGNORABILITY
    if ( present ( IgnorabilityOption ) ) &
      DP % IGNORABILITY  =  IgnorabilityOption

    if ( DP % Type  ==  '' ) &
      DP % Type  =  'a DivergencePart_CS' 

    DP % Name  =  'V'
    if ( present ( NameOption ) ) &
      DP % Name  =  NameOption

    call Show ( 'Initializing ' // trim ( DP % Type ), DP % IGNORABILITY )
    call Show ( DP % Name, 'Name', DP % IGNORABILITY )

    DP % CurrentSet  =>  CS

  end subroutine Initialize


  subroutine Show_DP ( DP )

    class ( DivergencePart_CS_Form ), intent ( in ) :: &
      DP

    character ( LDL ), dimension ( : ), allocatable :: &
      TypeWord

    call Split ( DP % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', DP % IGNORABILITY )

    call Show ( DP % Name, 'Name', DP % IGNORABILITY )
    call Show ( DP % CurrentSet % Name, 'CurrentSet', DP % IGNORABILITY )

  end subroutine Show_DP


  function Timer_F ( DP, Level ) result ( T )

    class ( DivergencePart_CS_Form ), intent ( inout ) :: &
      DP
    integer ( KDI ), intent ( in ) :: &
      Level
    type ( TimerForm ), pointer :: &
      T

    T  =>  PROGRAM_HEADER % Timer &
             ( Handle = DP % iTimer_F, &
               Name = 'Flx_' // trim ( DP % Name ) &
                      // '_' // trim ( DP % CurrentSet % Name ), &
               Level = Level )

  end function Timer_F


  function Timer_S ( DP, Level ) result ( T )

    class ( DivergencePart_CS_Form ), intent ( inout ) :: &
      DP
    integer ( KDI ), intent ( in ) :: &
      Level
    type ( TimerForm ), pointer :: &
      T

    T  =>  PROGRAM_HEADER % Timer &
             ( Handle = DP % iTimer_S, &
               Name = 'Strss_' // trim ( DP % Name ) &
                      // '_' // trim ( DP % CurrentSet % Name ), &
               Level = Level )

  end function Timer_S


  subroutine ComputeFluxes ( FS_F, DP, FS_CS, iC, iD )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS_F  !-- Fluxes
    class ( DivergencePart_CS_Form ), intent ( in ) :: &
      DP
    class ( FieldSetForm ), intent ( in ) :: &
      FS_CS
    integer ( KDI ), intent ( in ) :: &
      iC, &  !-- iChart
      iD     !-- iDimension

    integer ( KDI ) :: &
      iDensity

    associate ( CS  =>  DP % CurrentSet )

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


  subroutine ComputeStresses ( S_UD, DP, iC, iMomentum_1, iMomentum_2 )

    class ( FieldSetForm ), intent ( inout ) :: &
      S_UD
    class ( DivergencePart_CS_Form ), intent ( in ) :: &
      DP
    integer ( KDI ), intent ( in ) :: &
      iC  !-- iChart
    integer ( KDI ), intent ( out ) :: &
      iMomentum_1, iMomentum_2

    call S_UD % Clear ( )

  end subroutine ComputeStresses


  impure elemental subroutine Finalize ( DP )

    type ( DivergencePart_CS_Form ), intent ( inout ) :: &
      DP

    call Show ( 'Finalizing ' // trim ( DP % Type ), DP % IGNORABILITY )
    call Show ( DP % Name, 'Name', DP % IGNORABILITY )
   
  end subroutine Finalize


  impure elemental subroutine Finalize_E ( DPE )
    
    type ( DivergencePartElement ), intent ( inout ) :: &
      DPE

    if ( allocated ( DPE % Element ) ) &
      deallocate ( DPE % Element )

  end subroutine Finalize_E


end module DivergencePart_CS__Form
