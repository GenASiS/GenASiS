module IncrementDamping_Form

  use Basics
  use Operations
  use Fields

  implicit none
  private

  type, public :: IncrementDampingForm
    integer ( KDI ) :: &
      IGNORABILITY = 0
    character ( LDF ) :: &
      Name = ''
  contains   
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type IncrementDampingForm

    private :: &
      ComputeKernel

contains


  subroutine Initialize ( I, NameSuffix, IgnorabilityOption )

    class ( IncrementDampingForm ), intent ( inout ) :: &
      I
    character ( * ), intent ( in ) :: &
      NameSuffix
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    I % IGNORABILITY = CONSOLE % INFO_1
    if ( present ( IgnorabilityOption ) ) &
      I % IGNORABILITY = IgnorabilityOption
    
    I % Name = 'IncrementDamping_' // trim ( NameSuffix )

    call Show ( 'Initializing an Increment_C_R', I % IGNORABILITY )
    call Show ( I % Name, 'Name', I % IGNORABILITY )

  end subroutine Initialize


  subroutine Compute &
               ( I, Increment, C, IncrementExplicit, DampingCoefficient, &
                 TimeStep )

    class ( IncrementDampingForm ), intent ( inout ) :: &
      I
    type ( VariableGroupForm ), intent ( inout ) :: &
      Increment
    class ( CurrentTemplate ), intent ( in ) :: &
      C
    type ( VariableGroupForm ), intent ( in ) :: &
      IncrementExplicit, &
      DampingCoefficient
    real ( KDR ), intent ( in ) :: &
      TimeStep

    integer ( KDI ) :: &
      iF  !-- iField

    call Show ( 'Computing an IncrementDamping', I % IGNORABILITY + 3 )
    call Show ( I % Name, 'Name', I % IGNORABILITY + 3 )

    do iF = 1, C % N_CONSERVED
      call ComputeKernel &
             ( Increment % Value ( :, iF ), &
               C % Value ( :, C % iaConserved ( iF ) ), &
               IncrementExplicit % Value ( :, iF ), &
               DampingCoefficient % Value ( :, iF ), TimeStep )
    end do

  end subroutine Compute


  impure elemental subroutine Finalize ( I )

    type ( IncrementDampingForm ), intent ( inout ) :: &
      I

!    if ( allocated ( I % GridImageStream ) ) &
!      deallocate ( I % GridImageStream )
!    if ( allocated ( I % Output ) ) &
!      deallocate ( I % Output )

    if ( I % Name == '' ) return

    call Show ( 'Finalizing an IncrementDamping', I % IGNORABILITY )
    call Show ( I % Name, 'Name', I % IGNORABILITY )

  end subroutine Finalize


  subroutine ComputeKernel ( dU, U, IE, DC, dT )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      dU
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      U, &
      IE, &
      DC
    real ( KDR ), intent ( in ) :: &
      dT

    integer ( KDI ) :: &
      iV, &
      nV

    nV  =  size ( dU )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      dU ( iV )  =  ( IE ( iV )  -  DC ( iV ) * U ( iV ) * dT ) &
                    /  ( 1.0_KDR  +  DC ( iV ) * dT )
    end do
    !$OMP end parallel do

  end subroutine ComputeKernel


end module IncrementDamping_Form
