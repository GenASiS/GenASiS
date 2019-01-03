!-- Differences computes variable differentials on a chart: cell centered
!   values on input yield differences on the inner face on output.

module Difference_Form

  use Basics

  implicit none
  private

  type, public :: DifferenceForm
    integer ( KDI ) :: &
      IGNORABILITY = 0
    character ( LDF ) :: &
      Name = ''
    type ( VariableGroupForm ), allocatable :: &
      OutputInner
  contains
    procedure, public, pass :: &
      Initialize
    procedure, private, pass :: &
      ComputeChart_SL
    generic, public :: &
      Compute => ComputeChart_SL
    final :: &
      Finalize
  end type DifferenceForm

    private :: &
      ComputeChart_SL_Kernel

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


  subroutine ComputeChart_SL ( D, Input, nValues, iDimension )

    class ( DifferenceForm ), intent ( inout ) :: &
      D
    class ( VariableGroupForm ), intent ( in ) :: &
      Input
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nValues
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
        OI => D % OutputInner, &
        nV => nValues )

    do iS = 1, I % nVariables
      
      V ( 1 : nV ( 1 ), 1 : nV ( 2 ), 1 : nV ( 3 ) ) &
        => I % Value ( :, I % iaSelected ( iS ) )
      dV_I ( 1 : nV ( 1 ), 1 : nV ( 2 ), 1 : nV ( 3 ) ) &
        => OI % Value ( :,  iS )
      
      call ComputeChart_SL_Kernel ( V, iDimension, 0, dV_I )
    end do !-- iS

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


  subroutine ComputeChart_SL_Kernel ( V, iD, oV, dV_I )

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
!;    lV ( iD ) = oV

    uV = 1
    where ( shape ( V ) > 1 )
      uV = shape ( V ) - oV
    end where
    uV ( iD ) = size ( V, dim = iD )
    
    iaS = 0
!;    iaS ( iD ) = -1

    !$OMP OMP_TARGET_DIRECTIVE 
    do kV = lV ( 3 ), uV ( 3 ) 
      !$OMP parallel 
      !$OMP do private ( iaVS ) collapse ( 2 ) schedule (OMP_SCHEDULE)
      do jV = lV ( 2 ), uV ( 2 )
        do iV = lV ( 1 ), uV ( 1 )

          iaVS = [ iV, jV, kV ] + iaS

          dV_I ( iV, jV, kV )  &
            =  V ( iV, jV, kV )  -  V ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )

        end do !-- iV
      end do !-- jV
      !$OMP end do
      !$OMP end parallel
    end do !-- kV
    !$OMP end OMP_TARGET_DIRECTIVE
    
     
  end subroutine ComputeChart_SL_Kernel


end module Difference_Form
