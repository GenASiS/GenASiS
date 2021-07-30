#include "Preprocessor"

module OppenheimerSnyder_Form

  !-- For example, Misner, Thorne, Wheeler p. 663, Eqs. (25.28)-(25.29)

  use GenASiS

  implicit none
  private

  type, public, extends ( Universe_F_CC_Form ) :: OppenheimerSnyderForm
  !   real ( KDR ) :: &
  !     DensityInitial, &
  !     RadiusInitial, &
  !     TimeScale
  !   type ( RootFinderForm ), allocatable :: &
  !     RootFinder
    type ( Fluid_D_Form ), allocatable :: &
      Reference, &
      Difference
  contains
    procedure, private, pass :: &
      Initialize_H
    procedure, public, pass :: &
      ComputeError
    final :: &
      Finalize
  end type OppenheimerSnyderForm

    private :: &
      InitializeUniverse!, &
      ! InitializeDiagnostics, &
      ! SetInitial, &
      ! SetReference

      ! private :: &
      !   SetFluid

      !   private :: &
      !     SetFluidKernel

contains

 
  subroutine Initialize_H ( U, NameOption )

    class ( OppenheimerSnyderForm ), intent ( inout ), target :: &
      U
    character ( * ), intent ( in ), optional  :: &
      NameOption

    character ( LDL ) :: &
      Name

    if ( U % Type  ==  '' ) &
      U % Type  =  'a OppenheimerSnyder'

    Name  =  'OppenheimerSnyder'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    call InitializeUniverse ( U, Name )
!    call InitializeDiagnostics ( U )

!    if ( .not. associated ( U % Integrator % SetInitial ) ) &
!      U % Integrator % SetInitial  =>  SetInitial

  end subroutine Initialize_H


  subroutine ComputeError ( OS )

    class ( OppenheimerSnyderForm ), intent ( inout ) :: &
      OS

    real ( KDR ) :: &
      L1
    type ( CollectiveOperation_R_Form ) :: &
      CO
    
    select type ( A  =>  OS % Reference % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C    =>  A % Chart_GS, &
        F_R  =>  OS % Reference, &
        F_D  =>  OS % Difference )
    associate &
      ( FV_R  =>  F_R % Storage_GS % Value, &
        FV_D  =>  F_D % Storage_GS % Value )

    call CO % Initialize ( C % Communicator, [ 2 ], [ 2 ] )

    associate &
      ( D  =>  FV_D ( :, F_D % BARYON_DENSITY_C ), &
        R  =>  FV_R ( :, F_R % BARYON_DENSITY_C ), &
        Norm_D  =>  CO % Incoming % Value ( 1 ), &
        Norm_R  =>  CO % Incoming % Value ( 2 ) )

    CO % Outgoing % Value ( 1 ) &
      =  sum ( abs ( D ), mask = C % ProperCell )
    CO % Outgoing % Value ( 2 ) &
      =  sum ( abs ( R ), mask = C % ProperCell )

    call CO % Reduce ( REDUCTION % SUM )

    L1  =  Norm_D / Norm_R
    call Show ( L1, '*** L1 error', &
                nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )

    end associate !-- D, etc.
    end associate !-- FV_R, etc.
    end associate !-- C, etc.
    end select !-- A

  end subroutine ComputeError


  impure elemental subroutine Finalize ( OS )
    
    type ( OppenheimerSnyderForm ), intent ( inout ) :: &
      OS

    if ( allocated ( OS % Difference ) ) &
      deallocate ( OS % Difference )
    if ( allocated ( OS % Reference ) ) &
      deallocate ( OS % Reference )
    ! if ( allocated ( OS % RootFinder ) ) &
    !   deallocate ( OS % RootFinder )

  end subroutine Finalize


  subroutine InitializeUniverse ( OS, Name )

    class ( OppenheimerSnyderForm ), intent ( inout ) :: &
      OS
    character ( * ), intent ( in )  :: &
      Name

    ! integer ( KDI ) :: &
    !   iD

    call OS % Initialize &
           ( FluidType = 'DUST', &
             GravitationType = 'NEWTON_SG', &
             NameOption = Name, &
             nCellsPolarOption = 128 )

    ! select type ( I  =>  OS % Integrator )
    !   class is ( Integrator_CS_Form )
    ! associate &
    !   ( F  =>  I % CurrentSet_X )
    ! do iD  =  1, 3
    !   call F % SetBoundaryConditionsFace &
    !          ( [ 'PERIODIC', 'PERIODIC' ], iC = 1, iD = iD )
    ! end do !-- iD
    ! end associate !-- F
    ! end select !-- I
             
    ! OS % Integrator % SetReference  =>  SetReference

  end subroutine InitializeUniverse


end module OppenheimerSnyder_Form
