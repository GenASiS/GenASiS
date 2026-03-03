#include "Preprocessor"

submodule ( Gradient_Form ) Gradient_Kernel

  use Basics
  
  implicit none
  
contains


  module procedure Compute_CGS_Kernel

    integer ( KDI ) :: &
      iS, &
      iF, &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVP, iaVM, &
      lV, uV
    real ( KDR ) :: &
       fM,  fP, &  !-- f_Minus, f_Plus
      xAM, xAP
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice  =  .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice  =  UseDeviceOption

    lV  =  1
    where ( shape ( XA )  >  1 )
      lV  =  oV
    end where
    
    uV  =  1
    where ( shape ( XA )  >  1 )
      uV  =  shape ( XA )  -  oV
    end where
    uV ( iD )  =  size ( XA, dim = iD )  -  oV  +  1 
      
    iaS  =  0
    iaS ( iD )  =  1
    
    if ( UseDevice ) then
    
      !$OMP OMP_TARGET_DIRECTIVE parallel do collapse ( 4 ) &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP private ( iF, iaVP, iaVM, fM, fP, xAM, xAP )
      do iS  =  1,  size ( iaSlctd )
        do kV  =  lV ( 3 ),  uV ( 3 ) 
          do jV  =  lV ( 2 ),  uV ( 2 )
            do iV  =  lV ( 1 ),  uV ( 1 )

              iF  =  iaSlctd ( iS )

              iaVM  =  [ iV, jV, kV ]  -  iaS
              iaVP  =  [ iV, jV, kV ]  +  iaS

              fM  =  F ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ), iF )
              fP  =  F ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ), iF )

              xAM  =  XA ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ) )
              xAP  =  XA ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )

              dFdX ( iV, jV, kV, iS )  =  ( fP - fM ) / ( xAP - xAM )

            end do !-- iV
          end do !-- jV
        end do !-- kV
      end do !-- iS
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    else !-- use host
              
      !$OMP parallel do collapse ( 4 ) &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP private ( iF, iaVP, iaVM, fM, fP, xAM, xAP )
      do iS  =  1,  size ( iaSlctd )
        do kV  =  lV ( 3 ),  uV ( 3 ) 
          do jV  =  lV ( 2 ),  uV ( 2 )
            do iV  =  lV ( 1 ),  uV ( 1 )

              iF  =  iaSlctd ( iS )

              iaVM  =  [ iV, jV, kV ]  -  iaS
              iaVP  =  [ iV, jV, kV ]  +  iaS

              fM  =  F ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ), iF )
              fP  =  F ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ), iF )

              xAM  =  XA ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ) )
              xAP  =  XA ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )

              dFdX ( iV, jV, kV, iS )  =  ( fP - fM ) / ( xAP - xAM )

            end do !-- iV
          end do !-- jV
        end do !-- kV
      end do !-- iS
      !$OMP end parallel do
      
    end if !-- UseDevice
        
  end procedure Compute_CGS_Kernel


  module procedure ComputeFromDivergence_CGS_Kernel

    integer ( KDI ) :: &
      iS, &
      iF, &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVP, iaVM, &
      lV, uV
    real ( KDR ) :: &
       fM, fC, fP,       &  !-- f_Minus, f_Center, f_Plus
       fI,   fO,         &  !-- f_Inner, f_Outer
       xAM,  xAC,  xAP,  & 
      x2AM, x2AC, x2AP,  & 
        xI,   xO,        &  !-- X_Inner, X_Outer
         d,              &  !-- Determinant / Denominator
        c0,   c1,   c2      !-- Parabola coefficients
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice  =  .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice  =  UseDeviceOption

    lV  =  1
    where ( shape ( V )  >  1 )
      lV  =  oV
    end where
    
    uV  =  1
    where ( shape ( V )  >  1 )
      uV  =  shape ( V )  -  oV
    end where
    uV ( iD )  =  size ( V, dim = iD )  -  oV  +  1 
      
    iaS  =  0
    iaS ( iD )  =  1
    
    if ( UseDevice ) then
    
      !$OMP OMP_TARGET_DIRECTIVE parallel do collapse ( 4 ) &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP private ( iF, iaVP, iaVM, fM, fC, fP, fI, fO ) &
      !$OMP private ( xAM, xAC, xAP, x2AM, x2AC, x2AP ) &
      !$OMP private ( xI, xO, c0, c1, c2, d )
      do iS  =  1,  size ( iaSlctd )
        do kV  =  lV ( 3 ),  uV ( 3 ) 
          do jV  =  lV ( 2 ),  uV ( 2 )
            do iV  =  lV ( 1 ),  uV ( 1 )

              iF  =  iaSlctd ( iS )

              iaVM  =  [ iV, jV, kV ]  -  iaS
              iaVP  =  [ iV, jV, kV ]  +  iaS

              fM  =  F ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ), iF )
              fC  =  F ( iV        , jV        , kV        , iF )
              fP  =  F ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ), iF )

!              fI  =  0.5_KDR * ( fM + fC )
!              fO  =  0.5_KDR * ( fC + fP )

              xAM  =  XA ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ) )
              xAC  =  XA ( iV, jV, kV )
              xAP  =  XA ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )

              x2AM  =  X2A ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ) )
              x2AC  =  X2A ( iV, jV, kV )
              x2AP  =  X2A ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )

              xI  =  X ( iV, jV, kV )  -  0.5 * dX ( iV, jV, kV )
              xO  =  X ( iV, jV, kV )  +  0.5 * dX ( iV, jV, kV )

              d  =     x2AP * ( -xAC + xAM ) &
                    +  x2AM * (  xAC - xAP ) &
                    +  x2AC * ( -xAM + xAP )

              c0  =  (    fP * (   x2AM * xAC   -  x2AC * xAM ) &
                       +  fM * ( -(x2AP * xAC)  +  x2AC * xAP ) &
                       +  fC * (   x2AP * xAM   -  x2AM * xAP ) )  /  d

              c1  =  (    fP * (  x2AC - x2AM ) &
                       +  fC * (  x2AM - x2AP ) & 
                       +  fM * ( -x2AC + x2AP ) )  /  d

              c2  =  (    fP * ( -xAC + xAM ) &
                       +  fM * (  xAC - xAP ) &
                       +  fC * ( -xAM + xAP ) )  /  d

              fI  =  c0  +  c1 * xI  +  c2 * xI**2
              fO  =  c0  +  c1 * xO  +  c2 * xO**2

              dFdX ( iV, jV, kV, iS )  &
                =  (    A_I ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )  *  fO  &
                     -  A_I ( iV        , jV        , kV         )  *  fI )  &
                   /  V ( iV, jV, kV )  &
                   -  fC  *  (    A_I ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )  &
                               -  A_I ( iV        , jV        , kV         ) ) &
                      /  V ( iV, jV, kV )

            end do !-- iV
          end do !-- jV
        end do !-- kV
      end do !-- iS
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    else !-- use host
              
      !$OMP parallel do collapse ( 4 ) &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP private ( iF, iaVP, iaVM, fM, fC, fP, fI, fO ) &
      !$OMP private ( xAM, xAC, xAP, x2AM, x2AC, x2AP ) &
      !$OMP private ( xI, xO, c0, c1, c2, d )
      do iS  =  1,  size ( iaSlctd )
        do kV  =  lV ( 3 ),  uV ( 3 ) 
          do jV  =  lV ( 2 ),  uV ( 2 )
            do iV  =  lV ( 1 ),  uV ( 1 )

              iF  =  iaSlctd ( iS )

              iaVM  =  [ iV, jV, kV ]  -  iaS
              iaVP  =  [ iV, jV, kV ]  +  iaS

              fM  =  F ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ), iF )
              fC  =  F ( iV        , jV        , kV        , iF )
              fP  =  F ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ), iF )

!              fI  =  0.5_KDR * ( fM + fC )
!              fO  =  0.5_KDR * ( fC + fP )

              xAM  =  XA ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ) )
              xAC  =  XA ( iV, jV, kV )
              xAP  =  XA ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )

              x2AM  =  X2A ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ) )
              x2AC  =  X2A ( iV, jV, kV )
              x2AP  =  X2A ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )

              xI  =  X ( iV, jV, kV )  -  0.5 * dX ( iV, jV, kV )
              xO  =  X ( iV, jV, kV )  +  0.5 * dX ( iV, jV, kV )

              d  =     x2AP * ( -xAC + xAM ) &
                    +  x2AM * (  xAC - xAP ) &
                    +  x2AC * ( -xAM + xAP )

              c0  =  (    fP * (   x2AM * xAC   -  x2AC * xAM ) &
                       +  fM * ( -(x2AP * xAC)  +  x2AC * xAP ) &
                       +  fC * (   x2AP * xAM   -  x2AM * xAP ) )  /  d

              c1  =  (    fP * (  x2AC - x2AM ) &
                       +  fC * (  x2AM - x2AP ) & 
                       +  fM * ( -x2AC + x2AP ) )  /  d

              c2  =  (    fP * ( -xAC + xAM ) &
                       +  fM * (  xAC - xAP ) &
                       +  fC * ( -xAM + xAP ) )  /  d

              fI  =  c0  +  c1 * xI  +  c2 * xI**2
              fO  =  c0  +  c1 * xO  +  c2 * xO**2

              dFdX ( iV, jV, kV, iS )  &
                =  (    A_I ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )  *  fO  &
                     -  A_I ( iV        , jV        , kV         )  *  fI )  &
                   /  V ( iV, jV, kV )  &
                   -  fC  *  (    A_I ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )  &
                               -  A_I ( iV        , jV        , kV         ) ) &
                      /  V ( iV, jV, kV )

            end do !-- iV
          end do !-- jV
        end do !-- kV
      end do !-- iS
      !$OMP end parallel do
      
    end if !-- UseDevice
        
  end procedure ComputeFromDivergence_CGS_Kernel


end submodule Gradient_Kernel
