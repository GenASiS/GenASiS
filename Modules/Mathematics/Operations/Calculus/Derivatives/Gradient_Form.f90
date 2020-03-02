!-- Gradient computes gradients on a chart: cell centered
!   values on input yield cell-centered gradients on output.

module Gradient_Form

  use Basics
  use Manifolds
  use Difference_Form

  implicit none
  private

  type, public :: GradientForm
    integer ( KDI ) :: &
      IGNORABILITY = 0
    integer ( KDI ) :: &
      iTimerComputeGradient
    character ( LDF ) :: &
      Name = ''
    type ( StorageForm ), allocatable :: &
      Output
    type ( DifferenceForm ), allocatable :: &
      CoordinateDifference, &
      VariableDifference
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      AllocateDevice => AllocateDevice_G
    procedure, private, pass :: &
      ComputeChart_SL
    generic, public :: &
      Compute => ComputeChart_SL
    final :: &
      Finalize
  end type GradientForm

    private :: &
      ComputeChart_SL_Kernel
      
    interface
      
      module subroutine ComputeChart_SL_Kernel &
                   ( dV_I, dX_I, iD, oV, dVdX, UseDeviceOption, &
                     UseLimiterOption, ThetaOption )
        use Basics
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          dV_I, &
          dX_I
        integer ( KDI ), intent ( in ) :: &
          iD, &
          oV
        real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
          dVdX
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
        real ( KDR ), dimension ( :, :, : ), intent ( in ), optional, target :: &
          UseLimiterOption
        real ( KDR ), intent ( in ), optional, target :: &
          ThetaOption
      end subroutine ComputeChart_SL_Kernel
    
    end interface

contains


  subroutine Initialize ( G, Name, ValueShape )

    class ( GradientForm ), intent ( inout ) :: &
      G
    character ( * ), intent ( in ) :: &
      Name
    integer ( KDI ), dimension ( 2 ), intent ( in ) :: &
      ValueShape

    G % IGNORABILITY = CONSOLE % INFO_4
    G % Name = Name

    call Show ( 'Initializing a Gradient', G % IGNORABILITY )
    call Show ( trim ( G % Name ), 'Name', G % IGNORABILITY )

    allocate ( G % Output )
    call G % Output % Initialize ( ValueShape, ClearOption = .true. )

    allocate ( G % CoordinateDifference )
    allocate ( G % VariableDifference )
    call G % CoordinateDifference % Initialize &
           ( 'CoordinateDifference', [ ValueShape ( 1 ), 1 ] )
    call G % VariableDifference % Initialize &
           ( 'VariableDifference', ValueShape )
           
    call PROGRAM_HEADER % AddTimer &
           ( 'ComputeGradient', G % iTimerComputeGradient, Level = 6 )

  end subroutine Initialize
  
  
  subroutine AllocateDevice_G ( G )
    
    class ( GradientForm ), intent ( inout ) :: &
      G
    
    call G % Output % AllocateDevice ( )
    call G % CoordinateDifference % AllocateDevice ( )
    call G % VariableDifference % AllocateDevice ( )
  
  end subroutine AllocateDevice_G


  subroutine ComputeChart_SL &
               ( G, CSL, Input, iDimension, UseLimiterOption, &
                 LimiterParameterOption )

    class ( GradientForm ), intent ( inout ) :: &
      G
    class ( Chart_SL_Template ), intent ( in ) :: &
      CSL
    class ( StorageForm ), intent ( in ) :: &
      Input
    integer ( KDI ), intent ( in ) :: &
      iDimension
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      UseLimiterOption
    real ( KDR ), intent ( in ), optional :: &
      LimiterParameterOption

    integer ( KDI ) :: &
      iV  !-- iVariable
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      dV_I, &
      dX_I, &
      dVdX, &
      UL_Option
    type ( StorageForm ), allocatable :: &
      Coordinate
    class ( GeometryFlatForm ), pointer :: &
      Gmtry

    call Show ( 'Computing Gradient', G % IGNORABILITY + 1 )
    call Show ( trim ( G % Name ), 'Name', G % IGNORABILITY + 1 )
    call Show ( iDimension, 'iDimension', G % IGNORABILITY + 1 )

    associate &
      ( CD => G % CoordinateDifference, &
        VD => G % VariableDifference )

    !-- Coordinate difference

    Gmtry => CSL % Geometry ( )
    allocate ( Coordinate )
    call Coordinate % Initialize &
           ( Gmtry, iaSelectedOption = [ Gmtry % CENTER_U ( iDimension ) ] )
    
    call CD % Compute ( CSL, Coordinate, iDimension )

    !-- Variable difference

    call VD % Compute ( CSL, Input, iDimension )

    !-- Compute gradient

    associate &
      ( O    => G  % Output, &
        V_OI => VD % OutputInner, &
        C_OI => CD % OutputInner, &
        T_G  => PROGRAM_HEADER % Timer ( G % iTimerComputeGradient ) )
        
    call CSL % SetVariablePointer ( C_OI % Value ( :, 1 ), dX_I )
    
    do iV = 1, Input % nVariables
      call CSL % SetVariablePointer ( V_OI % Value ( :, iV ), dV_I )
      call CSL % SetVariablePointer ( O % Value ( :, iV ), dVdX )
      
      call T_G % Start ()
      
      if ( present ( UseLimiterOption ) ) then
        call CSL % SetVariablePointer ( UseLimiterOption, UL_Option )
        call ComputeChart_SL_Kernel &
               ( dV_I, dX_I, iDimension, CSL % nGhostLayers ( iDimension ), &
                 dVdX, UseLimiterOption = UL_Option, &
                 ThetaOption = LimiterParameterOption, &
                 UseDeviceOption = G % Output % AllocatedDevice )
      else
        call ComputeChart_SL_Kernel &
               ( dV_I, dX_I, iDimension, CSL % nGhostLayers ( iDimension ), &
                 dVdX, ThetaOption = LimiterParameterOption, &
                 UseDeviceOption = G % Output % AllocatedDevice )
      end if
      
      call T_G % Stop ()
      
    end do !-- iV
    
    end associate !-- O, V_OI, etc.

    end associate !-- VD, etc.
    nullify ( dV_I, dX_I, dVdX, Gmtry )

  end subroutine ComputeChart_SL


  impure elemental subroutine Finalize ( G )

    type ( GradientForm ), intent ( inout ) :: &
      G

    if ( allocated ( G % VariableDifference ) ) &
      deallocate ( G % VariableDifference )
    if ( allocated ( G % CoordinateDifference ) ) &
      deallocate ( G % CoordinateDifference )
    if ( allocated ( G % Output ) ) &
      deallocate ( G % Output )

    call Show ( 'Finalizing a Gradient', G % IGNORABILITY )
    call Show ( trim ( G % Name ), 'Name', G % IGNORABILITY )

  end subroutine Finalize


end module Gradient_Form
