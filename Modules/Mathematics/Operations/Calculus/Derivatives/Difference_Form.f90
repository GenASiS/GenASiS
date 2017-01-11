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
    type ( VariableGroupForm ), allocatable :: &
      Input, &
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


  subroutine Initialize ( D, Input, Name )

    class ( DifferenceForm ), intent ( inout ) :: &
      D
    class ( VariableGroupForm ), intent ( in ) :: &
      Input
    character ( * ), intent ( in ) :: &
      Name

    integer ( KDI ) :: &
      iS  !-- iSelected
    character ( LDL ) :: &
      dName
    character ( LDL ), dimension ( : ), allocatable :: &
      Variable

    D % IGNORABILITY = CONSOLE % INFO_5
    D % Name = Name

    call Show ( 'Initializing a Difference', D % IGNORABILITY )
    call Show ( trim ( D % Name ), 'Name', D % IGNORABILITY )

    allocate &
      ( D % Input, &
        D % OutputInner )
    associate &
      ( I  => D % Input, &
        OI => D % OutputInner )

    call I % Initialize ( Input )

    allocate ( Variable ( I % nVariables ) )
    Variable = [ ( I % Variable ( I % iaSelected ( iS ) ), &
                   iS = 1, I % nVariables ) ]

    dName = 'd_' // trim ( I % Name )
    call OI % Initialize &
           ( [ I % nValues, I % nVariables ], VariableOption = Variable, &
             NameOption = dName )!, ClearOption = .true. )

    end associate !-- I, etc.

  end subroutine Initialize


  subroutine ComputeChart_SL ( D, CSL, iDimension )

    class ( DifferenceForm ), intent ( inout ) :: &
      D
    class ( Chart_SL_Template ), intent ( in ) :: &
      CSL
    integer ( KDI ), intent ( in ) :: &
      iDimension

    integer ( KDI ) :: &
      iS  !-- iSelected
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      V, &
      dV_I

    call Show ( 'Computing a Difference', D % IGNORABILITY )
    call Show ( trim ( D % Name ), 'Name', D % IGNORABILITY )
    call Show ( iDimension, 'iDimension', D % IGNORABILITY )

    associate &
      ( I  => D % Input, &
        OI => D % OutputInner )

    do iS = 1, I % nVariables
      call CSL % SetVariablePointer &
             ( I % Value ( :, I % iaSelected ( iS ) ), V )
      call CSL % SetVariablePointer &
             ( OI % Value ( :, iS ), dV_I )
      call ComputeChart_SL_Kernel &
             ( V, iDimension, CSL % nGhostLayers ( iDimension ), dV_I )
    end do !-- iS

    end associate !-- I, etc.

    nullify ( V, dV_I )

  end subroutine ComputeChart_SL


  impure elemental subroutine Finalize ( D )

    type ( DifferenceForm ), intent ( inout ) :: &
      D

    if ( allocated ( D % OutputInner ) ) deallocate ( D % OutputInner )
    if ( allocated ( D % Input ) ) deallocate ( D % Input )

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
    lV ( iD ) = oV

    uV = 1
    where ( shape ( V ) > 1 )
      uV = shape ( V ) - oV
    end where
    uV ( iD ) = size ( V, dim = iD )
    
    iaS = 0
    iaS ( iD ) = -1

    !$OMP parallel do private ( iV, jV, kV, iaVS )
    do kV = lV ( 3 ), uV ( 3 ) 
      do jV = lV ( 2 ), uV ( 2 )
        do iV = lV ( 1 ), uV ( 1 )

          iaVS = [ iV, jV, kV ] + iaS

          dV_I ( iV, jV, kV )  &
            =  V ( iV, jV, kV )  -  V ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )

        end do !-- iV
      end do !-- jV
    end do !-- kV
    !$OMP end parallel do
     
  end subroutine ComputeChart_SL_Kernel


end module Difference_Form
