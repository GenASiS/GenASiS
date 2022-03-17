module Integral_Form

  !-- https://en.wikipedia.org/wiki/Romberg%27s_method

  use Basics

  implicit none
  private

  type, public :: IntegralForm
    integer ( KDI ), private :: &
      IGNORABILITY
    integer ( KDI ) :: &
      MaxIterations, &
      nIterations
    real ( KDR ) :: &
      RequestedPrecision, &
      AbsolutePrecision, &
      RelativePrecision
    real ( KDR ), dimension ( :, : ), allocatable :: &
      Romberg
    logical ( KDL ) :: &
      Success = .false.
    procedure ( IF ), public, pointer, nopass :: &
      Integrand => null ( )
    class ( * ), private, pointer :: &
      Parameters => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Compute
    procedure, public, pass :: &
      Show => Show_I
    final :: &
      Finalize
  end type IntegralForm

  abstract interface
    function IF ( Parameters, X ) result ( F )
      use Basics
      implicit none
      class ( * ), intent ( in ) :: &
        Parameters
      real ( KDR ), intent ( in ) :: &
        X
      real ( KDR ) :: &
        F
    end function IF
  end interface
  

contains


  subroutine Initialize &
               ( I, Parameters, PrecisionOption, MaxIterationsOption, &
                 VerbosityOption )

    class ( IntegralForm ), intent ( inout ) :: &
      I
    class ( * ), intent ( in ), target :: &
      Parameters
    real ( KDR ), intent ( in ), optional :: &
      PrecisionOption
    integer ( KDI ), intent ( in ), optional :: &
      MaxIterationsOption, &
      VerbosityOption

    I % IGNORABILITY  =  CONSOLE % INFO_3 
    if ( present ( VerbosityOption ) ) &
      I % IGNORABILITY  =  VerbosityOption 

    I % MaxIterations  =  20
    if ( present ( MaxIterationsOption ) ) &
      I % MaxIterations  =  MaxIterationsOption

    allocate ( I % Romberg ( 0 : I % MaxIterations - 1, 2 ) )

    I % RequestedPrecision  =  epsilon ( 1.0_KDR ) * 10.0_KDR 
    if ( present ( PrecisionOption ) ) &
      I % RequestedPrecision  =  PrecisionOption
    
    I % Parameters  =>  Parameters
  
  end subroutine Initialize


  subroutine Compute ( I, A, B, Integral )

    class ( IntegralForm ), intent ( inout ) :: &
      I
    real ( KDR ), intent ( in ) :: &
      A, B
    real ( KDR ), intent ( out ) :: &
      Integral

    integer ( KDI ) :: &
      iI, &   !-- iIteration
      iS, &   !-- iSum
      iR, &   !-- iRecursion
      iP, iC  !-- iPrevious, iCurrent
    real ( KDR ) :: &
      H, &  !-- Step size
      S, &  !-- Sum
      P4    !-- Power of 4
    procedure ( IF ), pointer :: &
      F

    associate &
      ( R      =>  I % Romberg, &
        P      =>  I % Parameters, &
        P_Abs  =>  I % AbsolutePrecision, &
        P_Rel  =>  I % RelativePrecision, &
        P_Req  =>  I % RequestedPrecision )

    I % nIterations  =  1

    F  =>  I % Integrand

    iC  =  1
    iP  =  2

    H  =  B - A
    
    R ( 0, iC )  =  0.5_KDR * H * ( F ( P, A )  +  F ( P, B ) )
    call Show ( R ( : 0, iC ), 'R', CONSOLE % INFO_7 )

    I % Success  =  .false.

    do iI  =  1,  I % MaxIterations - 1

      I % nIterations  =  I % nIterations  +  1

      iC  =  mod ( iC, 2 )  +  1
      iP  =  mod ( iP, 2 )  +  1

      H  =  H  /  2.0_KDR

      S  =  0.0_KDR
      do iS  =  1,  2 ** ( iI - 1 )
        S  =  S  +  F ( P, A + ( 2*iS - 1 ) * H )
      end do !-- iS

      R ( 0, iC )  =  0.5_KDR * R ( 0, iP )  +  H * S

      do iR  =  1,  iI

        P4  =  4.0_KDR ** iR

        R ( iR, iC )  =  ( P4 * R ( iR - 1, iC )  -  R ( iR - 1, iP ) )  &
                         /  ( P4 - 1.0_KDR )

      end do !-- iR

      Integral  =  R ( iI, iC )

      P_Abs  =  abs ( R ( iI, iC )  -  R ( iI - 1, iP ) )
      P_Rel  =  P_Abs  /  max ( abs ( R ( iI, iC ) ), tiny ( 0.0_KDR ) )

      call Show ( R ( : iI, iC ), 'R', CONSOLE % INFO_7 )
      call Show ( I % AbsolutePrecision, 'AbsolutePrecision', CONSOLE % INFO_7 )
      call Show ( I % RelativePrecision, 'RelativePrecision', CONSOLE % INFO_7 )

      if ( any ( [ P_Abs, P_Rel ]  <=  P_Req ) ) then
        I % Success  =  .true.
        exit
      end if

    end do !-- iI

    end associate !-- R, etc.

  end subroutine Compute


  subroutine Show_I ( I, Name )

    class ( IntegralForm ), intent ( in ) :: &
      I
    character ( * ), intent ( in ) :: &
      Name

    call Show ( Name )
    call Show ( I % Success, 'Success' )
    call Show ( I % nIterations, 'nIterations' )
    call Show ( I % RequestedPrecision, 'RequestedPrecision' )
    call Show ( I % AbsolutePrecision, 'AbsolutePrecision' )
    call Show ( I % RelativePrecision, 'RelativePrecision' )

  end subroutine Show_I


  subroutine Finalize ( I )
  
    type ( IntegralForm ), intent ( inout ) :: &
      I
    
    nullify ( I % Parameters )
    nullify ( I % Integrand )

    if ( allocated ( I % Romberg ) )  &
      deallocate ( I % Romberg )

  end subroutine Finalize  


end module Integral_Form
