module NormalRandomNumber_Function

  !-- Joseph L. Leva, ACM Transactions on Mathematical Software, Vol. 18, No. 4, 
  !   December 1992,449-453 
  !-- Or see Numerical Recipes 3rd ed. Section 7.3.9

  use Basics
  
  implicit none
  private

  public :: &
    NormalRandomNumber

   interface NormalRandomNumber
     module procedure NormalRandomNumberZeroMean
     module procedure NormalRandomNumberNonzeroMean
   end interface NormalRandomNumber
contains


  function NormalRandomNumberZeroMean ( Variance ) result ( NRN )

    real ( KDR ), intent ( in ) :: &
      Variance !-- Sigma ( standard deviation ) squared
    real ( KDR ) :: NRN

    NRN = NormalRandomNumberNonzeroMean ( 0.0_KDR, Variance )

  end function NormalRandomNumberZeroMean


  function NormalRandomNumberNonzeroMean ( Mean, Variance ) result ( NRN )

    real ( KDR ), intent ( in ) :: &
      Mean, &
      Variance !-- Sigma ( standard deviation ) squared
    real ( KDR ) :: NRN

    real ( KDR ) :: &
      U, V, &  !-- Uniform random numbers
      X, Y, Q  !-- Quadratic form
    logical ( KDL ) :: &
      TryAgain

    TryAgain = .true.
    do while ( TryAgain )

      call random_number ( U )
      call random_number ( V )
      V = 1.7156_KDR * ( V - 0.5_KDR )
      X = U - 0.449871_KDR
      Y = abs ( V ) + 0.386595_KDR
      Q = X ** 2  +  Y * ( 0.19600_KDR * Y  -  0.25472 * X )

      !-- Avoid log evaluation if possible

      if ( Q < 0.27597 ) exit

      if ( Q > 0.27846 ) then
        TryAgain = .true.
      else if ( V ** 2 > -4.0_KDR * log ( U ) * ( U ** 2 ) ) then
        TryAgain = .true.
      else
        TryAgain = .false.
      end if

    end do

    NRN = Mean + sqrt ( Variance ) * V / U

  end function NormalRandomNumberNonzeroMean


end module NormalRandomNumber_Function
