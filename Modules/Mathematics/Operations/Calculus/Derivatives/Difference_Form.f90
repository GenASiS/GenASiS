!-- Differences computes variable differentials on a chart: cell centered
!   values on input yield differences on the inner face on output.

#include "Preprocessor"

module Difference_Form

  use iso_c_binding
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
      ComputeChart_SL_KernelHost, &
      ComputeChart_SL_KernelDevice

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
      ( D_I  => I % D_Selected, &
        D_OI => OI % D_Selected, &
        iaS  => I % iaSelected )
    
    do iS = 1, I % nVariables
      call CSL % SetVariablePointer ( I % Value ( :, iaS ( iS ) ), V )
      call CSL % SetVariablePointer ( OI % Value ( :, iS ), dV_I )
      if ( OI % AllocatedDevice ) then
        call ComputeChart_SL_KernelDevice &
               ( V, iDimension, CSL % nGhostLayers ( iDimension ), &
                 D_I ( iaS ( iS ) ), D_OI ( iS ), dV_I )
      else
        call ComputeChart_SL_KernelHost &
               ( V, iDimension, CSL % nGhostLayers ( iDimension ), dV_I )
      end if
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


  subroutine ComputeChart_SL_KernelHost ( V, iD, oV, dV_I )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      V
    integer ( KDI ), intent ( in ) :: &
      iD, &
      oV   
    real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
      dV_I

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVS, &   
      lV, uV
    
!    dV_I = V - cshift ( V, shift = -1, dim = iD )

    lV = 1
    where ( shape ( V ) > 1 )
      lV = oV + 1
    end where
    lV ( iD ) = oV

    uV = 1
    where ( shape ( V ) > 1 )
      uV = shape ( V ) - oV
    end where
    uV ( iD ) = size ( V, dim = iD )
    
    iaS = 0
    iaS ( iD ) = -1

    !$OMP  parallel do collapse ( 3 ) &
    !$OMP& schedule ( OMP_SCHEDULE ) private ( iV, jV, kV, iaVS )
    do kV = lV ( 3 ), uV ( 3 ) 
      do jV = lV ( 2 ), uV ( 2 )
        do iV = lV ( 1 ), uV ( 1 )

          iaVS = [ iV, jV, kV ] + iaS

          dV_I ( iV, jV, kV )  &
            =  V ( iV, jV, kV )  -  V ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )

        end do !-- iV
      end do !-- jV
    end do !-- kV
    !$OMP  end parallel do
     
  end subroutine ComputeChart_SL_KernelHost


  subroutine ComputeChart_SL_KernelDevice ( V, iD, oV, D_V, D_dV_I, dV_I )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      V
    integer ( KDI ), intent ( in ) :: &
      iD, &
      oV
    type ( c_ptr ), intent ( in ) :: &
      D_V, &
      D_dV_I   
    real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
      dV_I

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVS, &   
      lV, uV
      
    call AssociateHost ( D_V, V )
    call AssociateHost ( D_dV_I, dV_I )
    
!    dV_I = V - cshift ( V, shift = -1, dim = iD )

    lV = 1
    where ( shape ( V ) > 1 )
      lV = oV + 1
    end where
    lV ( iD ) = oV

    uV = 1
    where ( shape ( V ) > 1 )
      uV = shape ( V ) - oV
    end where
    uV ( iD ) = size ( V, dim = iD )
    
    iaS = 0
    iaS ( iD ) = -1

    !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
    !$OMP& schedule ( OMP_SCHEDULE ) private ( iV, jV, kV, iaVS )
    do kV = lV ( 3 ), uV ( 3 ) 
      do jV = lV ( 2 ), uV ( 2 )
        do iV = lV ( 1 ), uV ( 1 )

          iaVS = [ iV, jV, kV ] + iaS

          dV_I ( iV, jV, kV )  &
            =  V ( iV, jV, kV )  -  V ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )

        end do !-- iV
      end do !-- jV
    end do !-- kV
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    call DisassociateHost ( dV_I )
    call DisassociateHost ( V )
    
     
  end subroutine ComputeChart_SL_KernelDevice


end module Difference_Form
