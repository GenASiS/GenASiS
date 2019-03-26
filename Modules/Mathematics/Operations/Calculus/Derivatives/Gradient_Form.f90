!-- Gradient computes gradients on a chart: cell centered
!   values on input yield cell-centered gradients on output.

#include "Preprocessor"

module Gradient_Form

  use iso_c_binding
  use Basics
  use Manifolds
  use Difference_Form

  implicit none
  private

  type, public :: GradientForm
    integer ( KDI ) :: &
      IGNORABILITY = 0
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
      ComputeChart_SL_KernelHost, &
      ComputeChart_SL_KernelDevice

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
        C_OI => CD % OutputInner )
    associate &
      ( D_O    => O    % D_Selected, &
        D_V_OI => V_OI % D_Selected, &
        D_C_OI => C_OI % D_Selected )
        
    call CSL % SetVariablePointer ( C_OI % Value ( :, 1 ), dX_I )
    
    do iV = 1, Input % nVariables
      call CSL % SetVariablePointer ( V_OI % Value ( :, iV ), dV_I )
      call CSL % SetVariablePointer ( O % Value ( :, iV ), dVdX )
      
      if ( present ( UseLimiterOption ) ) then
      
        call CSL % SetVariablePointer ( UseLimiterOption, UL_Option )
        if ( G % Output % AllocatedDevice ) then
          call ComputeChart_SL_KernelDevice &
                 ( dV_I, dX_I, iDimension, CSL % nGhostLayers ( iDimension ), &
                   D_V_OI ( iV ), D_C_OI ( iV ), D_O ( iV ), &
                   dVdX, UseLimiterOption = UL_Option, &
                   ThetaOption = LimiterParameterOption )
        else
          call ComputeChart_SL_KernelHost &
                 ( dV_I, dX_I, iDimension, CSL % nGhostLayers ( iDimension ), &
                   dVdX, UseLimiterOption = UL_Option, &
                   ThetaOption = LimiterParameterOption )
        end if
      
      else
        
        if ( G % Output % AllocatedDevice ) then
          
          call ComputeChart_SL_KernelDevice &
                 ( dV_I, dX_I, iDimension, CSL % nGhostLayers ( iDimension ), &
                   D_V_OI ( iV ), D_C_OI ( iV ), D_O ( iV ), &
                   dVdX, ThetaOption = LimiterParameterOption )
        
        else
        
          call ComputeChart_SL_KernelHost &
                 ( dV_I, dX_I, iDimension, CSL % nGhostLayers ( iDimension ), &
                   dVdX, ThetaOption = LimiterParameterOption )
        
        end if
      
      end if
    end do !-- iV
    
    end associate !-- D_O, D_V_OI, etc
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


  subroutine ComputeChart_SL_KernelHost &
               ( dV_I, dX_I, iD, oV, dVdX, UseLimiterOption, ThetaOption )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      dV_I, &
      dX_I
    integer ( KDI ), intent ( in ) :: &
      iD, &
      oV
    real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
      dVdX
    real ( KDR ), dimension ( :, :, : ), intent ( in ), optional :: &
      UseLimiterOption
    real ( KDR ), intent ( in ), optional :: &
      ThetaOption

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVS, &
      lV, uV
    real ( KDR ) :: &
      dX_L, dX_R, &
      dV_L, dV_R
    
    lV = 1
    where ( shape ( dV_I ) > 1 )
      lV = oV + 1
    end where
    lV ( iD ) = oV
    
    uV = 1
    where ( shape ( dV_I ) > 1 )
      uV = shape ( dV_I ) - oV
    end where
    uV ( iD ) = size ( dV_I, dim = iD ) - oV + 1 
      
    iaS = 0
    iaS ( iD ) = +1

    if ( present ( ThetaOption ) .and. present ( UseLimiterOption ) ) then
      associate &
        ( Theta => ThetaOption, &
          UseLimiter => UseLimiterOption )

      !$OMP  parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE ) &
      !$OMP& private ( iV, jV, kV, iaVS, dV_L, dV_R, dX_L, dX_R )
      do kV = lV ( 3 ), uV ( 3 ) 
        do jV = lV ( 2 ), uV ( 2 )
          do iV = lV ( 1 ), uV ( 1 )

            iaVS = [ iV, jV, kV ] + iaS

            dV_L = dV_I ( iV, jV, kV )
            dV_R = dV_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
            dX_L = dX_I ( iV, jV, kV )
            dX_R = dX_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
              
            if ( UseLimiter ( iV, jV, kV ) > 0.0_KDR ) then
              dVdX ( iV, jV, kV ) &
                = ( sign ( 0.5_KDR, dV_L ) + sign ( 0.5_KDR, dV_R ) ) &
                  * min ( abs ( Theta * dV_L / dX_L ), &
                          abs ( Theta * dV_R / dX_R ), &
                          abs ( ( dX_R ** 2  *  dV_L  +  dX_L ** 2  *  dV_R ) &
                                / ( dX_L * dX_R * ( dX_L + dX_R ) ) ) )
            else
              dVdX ( iV, jV, kV )&
                =  ( dX_R ** 2  *  dV_L  +  dX_L ** 2  *  dV_R ) &
                   / ( dX_L * dX_R * ( dX_L + dX_R ) )
            end if

          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end parallel do

      end associate !-- Theta    
    else if ( present ( ThetaOption ) ) then
      associate ( Theta => ThetaOption )

      !$OMP parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE ) &
      !$OMP& private ( iV, jV, kV, iaVS, dV_L, dV_R, dX_L, dX_R )
      do kV = lV ( 3 ), uV ( 3 ) 
        do jV = lV ( 2 ), uV ( 2 )
          do iV = lV ( 1 ), uV ( 1 )

            iaVS = [ iV, jV, kV ] + iaS

            dV_L = dV_I ( iV, jV, kV )
            dV_R = dV_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
            dX_L = dX_I ( iV, jV, kV )
            dX_R = dX_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
              
            dVdX ( iV, jV, kV ) &
              = ( sign ( 0.5_KDR, dV_L ) + sign ( 0.5_KDR, dV_R ) ) &
                * min ( abs ( Theta * dV_L / dX_L ), &
                        abs ( Theta * dV_R / dX_R ), &
                        abs ( ( dX_R ** 2  *  dV_L  +  dX_L ** 2  *  dV_R ) &
                              / ( dX_L * dX_R * ( dX_L + dX_R ) ) ) )

          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end parallel do

      end associate !-- Theta
    else

      !$OMP parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE ) &
      !$OMP& private ( iV, jV, kV, iaVS, dV_L, dV_R, dX_L, dX_R )
      do kV = lV ( 3 ), uV ( 3 ) 
        do jV = lV ( 2 ), uV ( 2 )
          do iV = lV ( 1 ), uV ( 1 )

            iaVS = [ iV, jV, kV ] + iaS

            dV_L = dV_I ( iV, jV, kV )
            dV_R = dV_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
            dX_L = dX_I ( iV, jV, kV )
            dX_R = dX_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
              
            dVdX ( iV, jV, kV )&
              =  ( dX_R ** 2  *  dV_L  +  dX_L ** 2  *  dV_R ) &
                 / ( dX_L * dX_R * ( dX_L + dX_R ) )
            
          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end parallel do

    end if
    
  end subroutine ComputeChart_SL_KernelHost


  subroutine ComputeChart_SL_KernelDevice &
               ( dV_I, dX_I, iD, oV, D_dV_I, D_dX_I, D_dVdX, dVdX, &
                 UseLimiterOption, ThetaOption )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      dV_I, &
      dX_I
    integer ( KDI ), intent ( in ) :: &
      iD, &
      oV
    type ( c_ptr ), intent ( in ) :: &
      D_dV_I, &
      D_dX_I, &
      D_dVdX
    real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
      dVdX
    real ( KDR ), dimension ( :, :, : ), intent ( in ), optional :: &
      UseLimiterOption
    real ( KDR ), intent ( in ), optional :: &
      ThetaOption

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVS, &
      lV, uV
    real ( KDR ) :: &
      dX_L, dX_R, &
      dV_L, dV_R
      
    call AssociateHost ( D_dV_I, dV_I )
    call AssociateHost ( D_dX_I, dX_I )
    call AssociateHost ( D_dVdX, dVdX )
    
    lV = 1
    where ( shape ( dV_I ) > 1 )
      lV = oV + 1
    end where
    lV ( iD ) = oV
    
    uV = 1
    where ( shape ( dV_I ) > 1 )
      uV = shape ( dV_I ) - oV
    end where
    uV ( iD ) = size ( dV_I, dim = iD ) - oV + 1 
      
    iaS = 0
    iaS ( iD ) = +1

    if ( present ( ThetaOption ) .and. present ( UseLimiterOption ) ) then
      associate &
        ( Theta => ThetaOption, &
          UseLimiter => UseLimiterOption )

      !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE ) &
      !$OMP& private ( iV, jV, kV, iaVS, dV_L, dV_R, dX_L, dX_R )
      do kV = lV ( 3 ), uV ( 3 ) 
        do jV = lV ( 2 ), uV ( 2 )
          do iV = lV ( 1 ), uV ( 1 )

            iaVS = [ iV, jV, kV ] + iaS

            dV_L = dV_I ( iV, jV, kV )
            dV_R = dV_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
            dX_L = dX_I ( iV, jV, kV )
            dX_R = dX_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
              
            if ( UseLimiter ( iV, jV, kV ) > 0.0_KDR ) then
              dVdX ( iV, jV, kV ) &
                = ( sign ( 0.5_KDR, dV_L ) + sign ( 0.5_KDR, dV_R ) ) &
                  * min ( abs ( Theta * dV_L / dX_L ), &
                          abs ( Theta * dV_R / dX_R ), &
                          abs ( ( dX_R ** 2  *  dV_L  +  dX_L ** 2  *  dV_R ) &
                                / ( dX_L * dX_R * ( dX_L + dX_R ) ) ) )
            else
              dVdX ( iV, jV, kV )&
                =  ( dX_R ** 2  *  dV_L  +  dX_L ** 2  *  dV_R ) &
                   / ( dX_L * dX_R * ( dX_L + dX_R ) )
            end if

          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do

      end associate !-- Theta    
    else if ( present ( ThetaOption ) ) then
      associate ( Theta => ThetaOption )

      !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE ) &
      !$OMP& private ( iV, jV, kV, iaVS, dV_L, dV_R, dX_L, dX_R )
      do kV = lV ( 3 ), uV ( 3 ) 
        do jV = lV ( 2 ), uV ( 2 )
          do iV = lV ( 1 ), uV ( 1 )

            iaVS = [ iV, jV, kV ] + iaS

            dV_L = dV_I ( iV, jV, kV )
            dV_R = dV_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
            dX_L = dX_I ( iV, jV, kV )
            dX_R = dX_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
              
            dVdX ( iV, jV, kV ) &
              = ( sign ( 0.5_KDR, dV_L ) + sign ( 0.5_KDR, dV_R ) ) &
                * min ( abs ( Theta * dV_L / dX_L ), &
                        abs ( Theta * dV_R / dX_R ), &
                        abs ( ( dX_R ** 2  *  dV_L  +  dX_L ** 2  *  dV_R ) &
                              / ( dX_L * dX_R * ( dX_L + dX_R ) ) ) )

          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do

      end associate !-- Theta
    else

      !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE ) &
      !$OMP& private ( iV, jV, kV, iaVS, dV_L, dV_R, dX_L, dX_R )
      do kV = lV ( 3 ), uV ( 3 ) 
        do jV = lV ( 2 ), uV ( 2 )
          do iV = lV ( 1 ), uV ( 1 )

            iaVS = [ iV, jV, kV ] + iaS

            dV_L = dV_I ( iV, jV, kV )
            dV_R = dV_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
            dX_L = dX_I ( iV, jV, kV )
            dX_R = dX_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
              
            dVdX ( iV, jV, kV )&
              =  ( dX_R ** 2  *  dV_L  +  dX_L ** 2  *  dV_R ) &
                 / ( dX_L * dX_R * ( dX_L + dX_R ) )
            
          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do

    end if
    
    call DisassociateHost ( dVdX )
    call DisassociateHost ( dX_I )
    call DisassociateHost ( dV_I )
    
  end subroutine ComputeChart_SL_KernelDevice


end module Gradient_Form
