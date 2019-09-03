!-- Differences computes variable differentials on a chart: cell centered
!   values on input yield differences on the inner face on output.

module Difference_Form

  use Basics
  use Manifolds

  implicit none
  private

  type, public :: DifferenceForm
    integer ( KDI ) :: &
      IGNORABILITY = 0
    character ( LDF ) :: &
      Name = ''
    type ( StorageForm ), allocatable :: &
      OutputInner
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      AllocateDevice => AllocateDevice_D
    procedure, private, pass :: &
      ComputeChart_SL
    generic, public :: &
      Compute => ComputeChart_SL
    final :: &
      Finalize
  end type DifferenceForm

    private :: &
      ComputeChart_SL_Kernel
      
    interface
      
      module subroutine ComputeChart_SL_Kernel &
                          ( V, iD, oV, dV_I, UseDeviceOption )
        use Basics
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          V
        integer ( KDI ), intent ( in ) :: &
          iD, &
          oV   
        real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
          dV_I
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeChart_SL_Kernel
    
    end interface

contains


  subroutine Initialize ( D, Name, ValueShape )

    class ( DifferenceForm ), intent ( inout ) :: &
      D
    character ( * ), intent ( in ) :: &
      Name
    integer ( KDI ), dimension ( 2 ), intent ( in ) :: &
      ValueShape

    D % IGNORABILITY = CONSOLE % INFO_5
    D % Name = Name

    call Show ( 'Initializing a Difference', D % IGNORABILITY )
    call Show ( trim ( D % Name ), 'Name', D % IGNORABILITY )

    allocate ( D % OutputInner )
    call D % OutputInner % Initialize ( ValueShape, ClearOption = .true. )

  end subroutine Initialize
  
  
  subroutine AllocateDevice_D ( D )
  
    class ( DifferenceForm ), intent ( inout ) :: &
      D
      
    call D % OutputInner % AllocateDevice ( )
    
  end subroutine AllocateDevice_D


  subroutine ComputeChart_SL ( D, CSL, Input, iDimension )

    class ( DifferenceForm ), intent ( inout ) :: &
      D
    class ( Chart_SL_Template ), intent ( in ) :: &
      CSL
    class ( StorageForm ), intent ( in ) :: &
      Input
    integer ( KDI ), intent ( in ) :: &
      iDimension

    integer ( KDI ) :: &
      iS  !-- iSelected
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      V, &
      dV_I

    call Show ( 'Computing Difference', D % IGNORABILITY + 1 )
    call Show ( trim ( D % Name ), 'Name', D % IGNORABILITY + 1 )
    call Show ( iDimension, 'iDimension', D % IGNORABILITY + 1 )

    associate &
      ( I  => Input, &
        OI => D % OutputInner )
    associate &
      ( iaS  => I % iaSelected )
    
    do iS = 1, I % nVariables
      call CSL % SetVariablePointer ( I % Value ( :, iaS ( iS ) ), V )
      call CSL % SetVariablePointer ( OI % Value ( :, iS ), dV_I )
      call ComputeChart_SL_Kernel &
             ( V, iDimension, CSL % nGhostLayers ( iDimension ), dV_I, &
               UseDeviceOption = OI % AllocatedDevice )
    end do !-- iS
    
    end associate !-- D_I, etc
    end associate !-- I, etc.

    nullify ( V, dV_I )

  end subroutine ComputeChart_SL


  impure elemental subroutine Finalize ( D )

    type ( DifferenceForm ), intent ( inout ) :: &
      D

    if ( allocated ( D % OutputInner ) ) &
      deallocate ( D % OutputInner )

    call Show ( 'Finalizing a Difference', D % IGNORABILITY )
    call Show ( trim ( D % Name ), 'Name', D % IGNORABILITY )

  end subroutine Finalize


end module Difference_Form
