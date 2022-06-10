module Interpolation_Form

  use Basics
  
  implicit none
  private
  
  type, public :: InterpolationForm
    integer ( KDI ), private :: &
      IGNORABILITY = 0
    integer ( KDI ) :: &
      nValues
    real ( KDR ), dimension ( : ), allocatable :: &
      Input, &
      Value
    real ( KDR ), dimension ( : ), allocatable, private :: &
      DD_Interpolant
  contains 
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Evaluate
    final :: &
      Finalize
  end type InterpolationForm
    
    !   Subroutines Spline() and Splint() are adapted from "Numerical Recipes
    !   in Fortran 77", 2nd ed., pg. 109 - 110, for curve fitting with 
    !   cubic spline interpolation.
    
    private :: &
      Spline, &
      Splint

contains

  
  subroutine Initialize &
               ( I, Input, Value, ValueDerivativeEndPointOption, &
                 VerbosityOption )
  
    class ( InterpolationForm ), intent ( inout ) :: &
      I
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Input, &
      Value
    real ( KDR ), dimension ( 2 ), intent ( in ), optional :: &
      ValueDerivativeEndPointOption
    integer ( KDI ), intent ( in ), optional :: & 
      VerbosityOption 
    
    integer ( KDI ) :: &
      nValues
    
    I % IGNORABILITY = CONSOLE % INFO_4
    if ( present ( VerbosityOption ) ) &
      I % IGNORABILITY = VerbosityOption
      
    I % nValues = size ( Input )
    
    allocate ( I % Input, source = Input )   
    allocate ( I % Value, source = Value )
    allocate ( I % DD_Interpolant ( I % nValues ) )
    
    call Spline ( I % Input, I % Value, I % DD_Interpolant, &
                  ValueDerivativeEndPointOption )
    
  end subroutine Initialize
  
  
  subroutine Evaluate ( I, Input, Value, SuccessOption )
    
    class ( InterpolationForm ), intent ( in ) :: &
      I
    real ( KDR ), intent ( in ) :: &
      Input
    real ( KDR ), intent ( out ) :: &
      Value
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
    
    integer ( KDI ) :: &
      iV
    
    call Search ( I % Input, Input, iV )
    
    if ( iV == 0 .or. iV >= I % nValues ) then
      if ( present ( SuccessOption ) ) &
        SuccessOption = .false.
      Value = 0.0_KDR
      call Show ( 'Cannot extrapolate value', &
                  'SplineInterpolationForm % Evaluate', CONSOLE % WARNING )
      call Show ( [ I % Input ( 1 ), I % Input ( I % nValues ) ], &
                  'Valid input range', I % IGNORABILITY )
      call Show ( Input, 'Out-of-range input', I % IGNORABILITY )
      call Show ( Value, 'Zeroing value', I % IGNORABILITY )
      return
    end if
    
    call Splint ( I % Input, I % Value, I % DD_Interpolant, Input, iV, Value )
    
    if ( present ( SuccessOption ) ) &
      SuccessOption = .true.
    
  end subroutine Evaluate
  
  
  impure elemental subroutine Finalize ( I )
    
    type ( InterpolationForm ), intent ( inout ) :: &
      I

    if ( allocated ( I % DD_Interpolant ) ) &
      deallocate ( I % DD_Interpolant )
    if ( allocated ( I % Value ) ) &
      deallocate ( I % Value )
    if ( allocated ( I % Input ) ) &
      deallocate ( I % Input )
    
  end subroutine Finalize 
  
  
  subroutine Spline ( X, Y, Y2, YP_Option )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      X, &
      Y
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      Y2
    real ( KDR ), dimension ( 2 ), intent ( in ), optional :: &
      YP_Option
    
    integer ( KDI ) :: &
      iV, &
      kV, &
      nValues
    real ( KDR )  :: &
      P, &
      Qn, &
      Sig, &
      Un
    real ( KDR ) , dimension ( size ( X ) )  :: &
      U
      
    nValues = size ( X ) 

    if ( present ( YP_Option ) ) then
      Y2 ( 1 )  = -0.5_KDR
      U  ( 1 )  =  ( 3. / ( X ( 2 ) -X ( 1 ) ) ) &
                   * ( ( Y ( 2 ) - Y ( 1 ) ) / ( X ( 2 ) - X ( 1 ) ) &
                       - YP_Option ( 1 ) ) 
    else
      Y2 ( 1 )  = 0.0_KDR
      U  ( 1 )  = 0.0_KDR
    end if
    
    do  iV = 2, nValues - 1
      Sig = ( X ( iV ) - X ( iV - 1 ) ) / ( X ( iV + 1 ) - X ( iV - 1 ) ) 
      P   = Sig * Y2 ( iV - 1 )  + 2.0_KDR
      Y2 ( iV ) = ( Sig - 1.0_KDR )  / P
      U  ( iV ) &
        =  ( 6.0_KDR * ( ( Y ( iV + 1 )  - Y ( iV ) ) &
                           /  ( X ( iV+1 )  - X ( iV ) ) &
                         -  ( Y ( iV )  - Y ( iV - 1 ) ) &
                              /  ( X ( iV )  - X ( iV - 1 ) )  )  &
             / ( X ( iV + 1 )  - X ( iV - 1 ) ) - Sig * U ( iV - 1 ) )  &
           / P
    end do
    
    if ( present ( YP_Option ) ) then
      Qn = 0.5_KDR
      Un = ( 3.0_KDR /  ( X ( nValues )  - X ( nValues - 1 ) ) )  &
           * ( YP_Option ( nValues )  &
               - ( Y ( nValues )  - Y ( nValues - 1 ) ) &
                   /  ( X ( nValues ) - X ( nValues - 1 ) ) ) 
    else
      Qn = 0.0_KDR
      Un = 0.0_KDR
    end if
      
    Y2 ( nValues ) &
      = ( Un - Qn * U ( nValues - 1 ) ) / ( Qn * Y2 ( nValues - 1 ) + 1.0_KDR ) 
    
    do kV = nValues - 1, 1, -1
      Y2 ( kV )  =  Y2 ( kV ) * Y2 ( kV + 1 ) + U ( kV ) 
    end do
    
  end subroutine Spline
  
  
  subroutine Splint ( X_Known, Y_Known, Y2_Known, X, iV, Y ) 
    
    real ( KDR ) , dimension ( : ) , intent ( in )  :: &
      X_Known,  &
      Y_Known, &
      Y2_Known
    real ( KDR ) , intent ( in )  :: &
      X
    integer ( KDI ) , intent ( in )  :: &
      iV
    real ( KDR ) , intent ( out )  :: &
      Y
    
    real ( KDR )  :: &
      A, &
      B, &
      H
    
    H = X_Known ( iV + 1 ) - X_Known ( iV ) 
    A = ( X_Known ( iV + 1 )  - X ) / H
    B = ( X - X_Known ( iV ) ) / H
    Y = A * Y_Known ( iV ) + B * Y_Known ( iV + 1 )  &
        + ( ( A ** 3 - A ) * Y2_Known ( iV ) + ( B ** 3 - B ) &
            * Y2_Known ( iV + 1 ) ) * ( h ** 2 ) / 6.0_KDR
  
  end subroutine Splint
  
  
end module Interpolation_Form
