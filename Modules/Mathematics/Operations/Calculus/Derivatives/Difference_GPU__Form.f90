!-- Differences computes variable differentials on a chart: cell centered
!   values on input yield differences on the inner face on output.

module Difference_GPU__Form
  
  use ISO_C_BINDING
  use OMP_LIB
  use Basics
  use Manifolds
  use AllocateMemory_Command

  implicit none
  private

  type, public :: Difference_GPU_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0
    character ( LDF ) :: &
      Name = ''
    type ( VariableGroupForm ), allocatable :: &
      Input, &
      OutputInner
    type ( TimerForm ) :: &
      T_Compute, &
      T_Communication
    class ( Chart_SL_Template ), pointer :: &
      Chart
    type ( c_ptr ) :: &
      D_OutputInner
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      SetInput
    procedure, private, pass :: &
      ComputeChart_SL
    procedure, private, pass :: &
      ComputeChart_SL_Kernel
    generic, public :: &
      Compute => ComputeChart_SL
    procedure, public, pass :: &
      GetOutput
    final :: &
      Finalize
  end type Difference_GPU_Form

    private :: &
      ComputeChart_SL_Kernel
      
contains


  subroutine Initialize ( D, Name, ValueShape )

    class ( Difference_GPU_Form ), intent ( inout ) :: &
      D
    character ( * ), intent ( in ) :: &
      Name
    integer ( KDI ), dimension ( 2 ), intent ( in ) :: &
      ValueShape

    D % IGNORABILITY = CONSOLE % INFO_1
    D % Name = Name

    call Show ( 'Initializing a Difference', D % IGNORABILITY )
    call Show ( trim ( D % Name ), 'Name', D % IGNORABILITY )
    
    allocate ( D % Input )

    allocate ( D % OutputInner )
    call D % OutputInner % Initialize ( ValueShape, ClearOption = .true. )
    
    call D % T_Compute % Initialize ( 'D_GPU Compute', Level = 1 )
    call D % T_Communication % Initialize ( 'D_GPU Communication', Level = 1 )

  end subroutine Initialize
  
  
  subroutine SetInput ( D, CSL, Input )
  
    class ( Difference_GPU_Form ), intent ( inout ) :: &
      D
    class ( Chart_SL_Template ), intent ( in ), target :: &
      CSL
    class ( VariableGroupForm ), intent ( in ), target :: &
      Input
    
    
    integer ( KDI ) :: &
      nValues
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      dV_I
    real ( KDR ), dimension ( :, :, :, : ), pointer :: &
      V
    !type ( c_ptr ) :: &
    !  dV_I
    
    call D % T_Communication % Start ( ) 
    
    call D % Input % Initialize ( Input ) 
    D % Chart => CSL
    call CSL % SetVariablePointer_Any_4D &
           ( D % Input % Value, size ( D % Input % Value ),  V )
    call CSL % SetVariablePointer ( D % OutputInner % Value ( :, 1 ), dV_I )
    
    !$OMP target enter data map ( to: V )
    !$OMP target enter data map ( alloc : dV_I ) 
    
    !
    !nValues = size ( D % OutputInner % Value )
    !call Show ( nValues, 'nValues F' )
    !D % D_OutputInner = Allocate_D ( nValues )
    !
    !dV_I = D % D_OutputInner
    !!$OMP target is_device_ptr ( dV_I )
    !dV_I = D % D_OutputInner
    !!$OMP end target 
    
    call D % T_Communication % Stop ( ) 
    
  end subroutine SetInput 


  subroutine ComputeChart_SL ( D, CSL, iDimension )

    class ( Difference_GPU_Form ), intent ( inout ) :: &
      D
    class ( Chart_SL_Template ), intent ( in ) :: &
      CSL
    integer ( KDI ), intent ( in ) :: &
      iDimension

    integer ( KDI ) :: &
      iS  !-- iSelected
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      dV_I
    real ( KDR ), dimension ( :, :, :, : ), pointer :: &
      V

    call Show ( 'Computing Difference', D % IGNORABILITY + 1 )
    call Show ( trim ( D % Name ), 'Name', D % IGNORABILITY + 1 )
    call Show ( iDimension, 'iDimension', D % IGNORABILITY + 1 )

    associate &
      ( I  => D % Input, &
        OI => D % OutputInner )

    do iS = 1, I % nVariables
      call CSL % SetVariablePointer_Any_4D ( I % Value,  size ( I % Value ), V )
      call CSL % SetVariablePointer ( OI % Value ( :, iS ), dV_I )
      call D % ComputeChart_SL_Kernel &
             ( V, iDimension, I % iaSelected ( iS ), &
               CSL % nGhostLayers ( iDimension ), dV_I )
      call D % GetOutput ( iS )
    end do !-- iS

    end associate !-- I, etc.

    nullify ( V, dV_I )

  end subroutine ComputeChart_SL
  
  
  subroutine GetOutput ( D, iS )
  
    class ( Difference_GPU_Form ), intent ( inout ) :: &
      D
    integer ( KDI ), intent ( in ) :: &
      iS  !-- iSelected
    
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      dV_I
    type ( c_ptr ) :: &
      dV_I_D
    
    associate &
      ( I  => D % Input, &
        OI => D % OutputInner )
    
    dV_I_D = D % D_OutputInner
    print*, 'dV_I_', D % D_OutputInner
    call D % Chart % SetVariablePointer ( OI % Value ( :, iS ), dV_I )
    !-- !$OMP target data map ( dV_I ) is_device_ptr ( dV_I_D )
    !-- !$OMP target update from ( dV_I )
    !-- !$OMP end target data
             
    end associate 
    
  end subroutine GetOutput 


  impure elemental subroutine Finalize ( D )

    type ( Difference_GPU_Form ), intent ( inout ) :: &
      D
    
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      dV_I
    real ( KDR ), dimension ( :, :, :, : ), pointer :: &
      V
      
    call D % Chart % SetVariablePointer_Any_4D &
           ( D % Input % Value, size ( D % Input % Value ),  V )
    call D % Chart % SetVariablePointer &
           ( D % OutputInner % Value ( :, 1 ), dV_I )
    
    !-- !$OMP target exit data map ( delete: V )
    !-- !$OMP target exit data map ( delete: dV_I )
    
    !nullify ( V, dV_I )
    !nullify ( D % Chart )
    
    call D % T_Compute % ShowTotal ( CONSOLE % INFO_1 )
    call D % T_Communication % ShowTotal ( CONSOLE % INFO_1 )
    
    if ( allocated ( D % OutputInner ) ) &
      deallocate ( D % OutputInner )

    call Show ( 'Finalizing a Difference', D % IGNORABILITY )
    call Show ( trim ( D % Name ), 'Name', D % IGNORABILITY )

  end subroutine Finalize


  subroutine ComputeChart_SL_Kernel ( D, V, iD, iS, oV, dV_I )
  
    class ( Difference_GPU_Form ), intent ( inout ) :: &
      D 
    real ( KDR ), dimension ( :, :, :, : ), intent ( in ) :: &
      V
    integer ( KDI ), intent ( in ) :: &
      iD, &
      iS, &
      oV   
    real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
      dV_I

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVS, &   
      lV, uV
    type ( c_ptr ) :: &
      dV_I_D
    
!    dV_I = V - cshift ( V, shift = -1, dim = iD )

    lV = 1
    where ( shape ( V ( :, :, :, iS ) ) > 1 )
      lV = oV + 1
    end where
    lV ( iD ) = oV

    uV = 1
    where ( shape ( V ( :, :, :, iS ) ) > 1 )
      uV = shape ( V ( :, :, :, iS ) ) - oV
    end where
    uV ( iD ) = size ( V, dim = iD )
    
    iaS = 0
    iaS ( iD ) = -1
    
    call D % T_Compute % Start ( )
    
    !$OMP OMP_TARGET_DIRECTIVE 
    do kV = lV ( 3 ), uV ( 3 ) 
      !$OMP parallel do 
      do jV = lV ( 2 ), uV ( 2 )
        !$OMP simd
        do iV = lV ( 1 ), uV ( 1 )

          !-- iaVS = [ iV, jV, kV ] + iaS

          dV_I ( iV, jV, kV )  &
            =  V ( iV, jV, kV, iS )  &
               -  V ( iV + iaS ( 1 ), jV + iaS ( 2 ), kV + iaS ( 3 ), iS )

        end do !-- iV
        !$OMP end simd
      end do !-- jV
      !$OMP end parallel do
    end do !-- kV
    !$OMP end OMP_TARGET_DIRECTIVE
    
    call D % T_Compute % Stop ( )
    
    call D % T_Communication % Start ( )
    
    !$OMP target data map ( dV_I_D ) use_device_ptr ( dV_I )
    D % D_OutputInner = dV_I_D
    !$OMP end target data
    
    !--$OMP target update from ( dV_I )
    
    call D % T_Communication % Stop ( )
    
  end subroutine ComputeChart_SL_Kernel
  

end module Difference_GPU__Form
