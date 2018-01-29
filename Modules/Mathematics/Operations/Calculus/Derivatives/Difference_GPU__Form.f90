!-- Differences computes variable differentials on a chart: cell centered
!   values on input yield differences on the inner face on output.

module Difference_GPU__Form

  use Basics
  use Manifolds

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
    class ( Chart_SL_Template ), pointer :: &
      Chart
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      SetInput
    procedure, private, pass :: &
      ComputeChart_SL
    generic, public :: &
      Compute => ComputeChart_SL
    final :: &
      Finalize
  end type Difference_GPU_Form

    private :: &
      ComputeChart_SL_Kernel, &
      ComputeChart_SL_Kernel_2

contains


  subroutine Initialize ( D, Name, ValueShape )

    class ( Difference_GPU_Form ), intent ( inout ) :: &
      D
    character ( * ), intent ( in ) :: &
      Name
    integer ( KDI ), dimension ( 2 ), intent ( in ) :: &
      ValueShape

    D % IGNORABILITY = CONSOLE % INFO_5
    D % Name = Name

    call Show ( 'Initializing a Difference', D % IGNORABILITY )
    call Show ( trim ( D % Name ), 'Name', D % IGNORABILITY )
    
    allocate ( D % Input )

    allocate ( D % OutputInner )
    call D % OutputInner % Initialize ( ValueShape, ClearOption = .true. )

  end subroutine Initialize
  
  
  subroutine SetInput ( D, CSL, Input )
  
    class ( Difference_GPU_Form ), intent ( inout ) :: &
      D
    class ( Chart_SL_Template ), intent ( in ), target :: &
      CSL
    class ( VariableGroupForm ), intent ( in ), target :: &
      Input
    
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      dV_I
    real ( KDR ), dimension ( :, :, :, : ), pointer :: &
      V
   
    call D % Input % Initialize ( Input ) 
    D % Chart => CSL
    call CSL % SetVariablePointer_Any_4D &
           ( D % Input % Value, size ( D % Input % Value ),  V )
    call CSL % SetVariablePointer ( D % OutputInner % Value ( :, 1 ), dV_I )
    
    !$OMP target enter data map ( to: V )
    !$OMP target enter data map ( alloc : dV_I )
    
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
      call ComputeChart_SL_Kernel &
             ( V, iDimension, I % iaSelected ( iS ), &
               CSL % nGhostLayers ( iDimension ), dV_I )
    end do !-- iS

    end associate !-- I, etc.

    nullify ( V, dV_I )

  end subroutine ComputeChart_SL


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
    
    if ( allocated ( D % OutputInner ) ) &
      deallocate ( D % OutputInner )

    call Show ( 'Finalizing a Difference', D % IGNORABILITY )
    call Show ( trim ( D % Name ), 'Name', D % IGNORABILITY )

  end subroutine Finalize


  subroutine ComputeChart_SL_Kernel ( V, iD, iS, oV, dV_I )

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
    
    !$OMP target data map ( from : dV_I )
    
    !$OMP target teams distribute 
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
    !$OMP end target teams distribute
    
    !$OMP end target data
    
  end subroutine ComputeChart_SL_Kernel
  
  
  subroutine ComputeChart_SL_Kernel_2 ( V, iD, oV, dV_I )

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
    
    if ( iD == 1 ) then
    
      !$OMP target teams distribute 
      do kV = lV ( 3 ), uV ( 3 ) 
        !$OMP parallel do 
        do jV = lV ( 2 ), uV ( 2 )
          !$OMP simd
          do iV = lV ( 1 ), uV ( 1 )

            dV_I ( iV, jV, kV )  =  V ( iV, jV, kV )  -  V ( iV - 1, jV, kV )

          end do !-- iV
          !$OMP end simd
        end do !-- jV
        !$OMP end parallel do
      end do !-- kV
      !$OMP end target teams distribute
    
    else if ( iD == 2 ) then
    
      !$OMP target teams distribute 
      do kV = lV ( 3 ), uV ( 3 ) 
        !$OMP parallel do 
        do jV = lV ( 2 ), uV ( 2 )
          !$OMP simd
          do iV = lV ( 1 ), uV ( 1 )

            dV_I ( iV, jV, kV )  =  V ( iV, jV, kV )  -  V ( iV, jV - 1, kV )

          end do !-- iV
          !$OMP end simd
        end do !-- jV
        !$OMP end parallel do
      end do !-- kV
      !$OMP end target teams distribute
      
    else if ( iD == 3 ) then

      !$OMP target teams distribute 
      do kV = lV ( 3 ), uV ( 3 ) 
        !$OMP parallel do 
        do jV = lV ( 2 ), uV ( 2 )
          !$OMP simd
          do iV = lV ( 1 ), uV ( 1 )

            dV_I ( iV, jV, kV )  =  V ( iV, jV, kV )  -  V ( iV, jV, kV - 1 )

          end do !-- iV
          !$OMP end simd
        end do !-- jV
        !$OMP end parallel do
      end do !-- kV
      !$OMP end target teams distribute
      
    end if 
         
  end subroutine ComputeChart_SL_Kernel_2


end module Difference_GPU__Form
