#include "Preprocessor"

submodule ( Reconstruction_Form ) Reconstruction_Kernel

  use Basics
  
  implicit none
  
contains


  module procedure ComputeConstant_CGS_Kernel

    integer ( KDI ) :: &
      iS, &
      iF, &
      iF_R, &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVP, &
      lV, uV
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice  =  .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice  =  UseDeviceOption

    lV  =  1
    where ( shape ( F ( :, :, :, 1 ) )  >  1 )
      lV  =  oV
    end where
    
    uV  =  1
    where ( shape ( F ( :, :, :, 1 ) )  >  1 )
      uV  =  shape ( F ( :, :, :, 1 ) )  -  oV
    end where
    uV ( iD )  =  size ( F, dim = iD )  -  oV  +  1 
      
    iaS  =  0
    iaS ( iD )  =  1
    
    if ( UseDevice ) then
    
      !$OMP OMP_TARGET_DIRECTIVE parallel do simd collapse ( 4 ) &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) private ( iF, iF_R, iaVP )
      do iS  =  1,  size ( iaSlctd )
        do kV  =  lV ( 3 ),  uV ( 3 ) 
          do jV  =  lV ( 2 ),  uV ( 2 )
            do iV  =  lV ( 1 ),  uV ( 1 )

              iF    =  iaSlctd   ( iS )
              iF_R  =  iaSlctd_R ( iS )

              iaVP  =  [ iV, jV, kV ]  +  iaS

              F_IR ( iV, jV, kV, iF_R )  &
                =  F ( iV, jV, kV, iF )

              F_IL ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ), iF_R )  &
                =  F ( iV, jV, kV, iF )                    

            end do !-- iV
          end do !-- jV
        end do !-- kV
      end do !-- iS
      !$OMP end OMP_TARGET_DIRECTIVE parallel do simd
      
    else !-- use host
              
      !$OMP parallel do collapse ( 4 ) &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) private ( iF, iF_R, iaVP )
      do iS  =  1,  size ( iaSlctd )
        do kV  =  lV ( 3 ),  uV ( 3 ) 
          do jV  =  lV ( 2 ),  uV ( 2 )
            do iV  =  lV ( 1 ),  uV ( 1 )

              iF    =  iaSlctd   ( iS )
              iF_R  =  iaSlctd_R ( iS )

              iaVP  =  [ iV, jV, kV ]  +  iaS

              F_IR ( iV, jV, kV, iF_R )  &
                =  F ( iV, jV, kV, iF )

              F_IL ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ), iF_R )  &
                =  F ( iV, jV, kV, iF )                    

            end do !-- iV
          end do !-- jV
        end do !-- kV
      end do !-- iS
      !$OMP end parallel do
      
    end if !-- UseDevice
        
  end procedure ComputeConstant_CGS_Kernel


  module procedure ComputeLinear_CGS_Kernel

    integer ( KDI ) :: &
      iS, &
      iF, &
      iF_R, &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVP, iaVM, &
      lV, uV
    real ( KDR ) :: &
       fM,  fC,  fP, &  !-- f_Minus, f_Center, f_Plus
      xAM, xAC, xAP, & 
       fI,  fO, &       !-- F_Inner, F_Outer
       xI,  xO, &       !-- X_Inner, X_Outer
       c0,  c1          !-- Line coefficients
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice  =  .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice  =  UseDeviceOption

    lV  =  1
    where ( shape ( X )  >  1 )
      lV  =  oV
    end where
    
    uV  =  1
    where ( shape ( X )  >  1 )
      uV  =  shape ( X )  -  oV
    end where
    uV ( iD )  =  size ( X, dim = iD )  -  oV  +  1 
      
    iaS  =  0
    iaS ( iD )  =  1
    
    if ( UseDevice ) then

      !$OMP OMP_TARGET_DIRECTIVE parallel do simd collapse ( 4 ) &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP private ( iF, iF_R, iaVP, iaVM, fM, fC, fP, fI, fO ) &
      !$OMP private ( xAM, xAC, xAP, xI, xO, c0, c1 )
      do iS  =  1,  size ( iaSlctd )
        do kV  =  lV ( 3 ),  uV ( 3 ) 
          do jV  =  lV ( 2 ),  uV ( 2 )
            do iV  =  lV ( 1 ),  uV ( 1 )

              iF    =  iaSlctd   ( iS )
              iF_R  =  iaSlctd_R ( iS )

              iaVP  =  [ iV, jV, kV ]  +  iaS
              iaVM  =  [ iV, jV, kV ]  -  iaS

              fM  =  F ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ), iF )
              fC  =  F ( iV, jV, kV, iF )
              fP  =  F ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ), iF )

              xAM  =  XA ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ) )
              xAC  =  XA ( iV, jV, kV )
              xAP  =  XA ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )

              xI  =  X ( iV, jV, kV )  -  0.5 * dX ( iV, jV, kV )
              xO  =  X ( iV, jV, kV )  +  0.5 * dX ( iV, jV, kV )

!call Show ( iV, '>>> iV' )
!call Show ( [ fM, fC, fP ], '>>> fM, fC, fP' )
              !-- Local extremum of cell average values? 
              !   Then reconstruction is constant.
              if ( ( fC - fM ) * ( fP - fC )  <  0.0_KDR ) then
!call Show ( '>>> Local extremum' )

                c1  =  0.0_KDR
                c0  =  fC

              else  !-- Linear reconstruction

                c1  =  ( fP - fM ) / ( xAP - xAM )
                c0  =  fC  -  c1 * xAC

                fI  =  c0  +  c1 * xI
                fO  =  c0  +  c1 * xO

                !-- Overshoot at inner face?
                !   Reduce slope.
                if ( c1 * ( fI - fM )  <  0.0_KDR ) then
!call Show ( '>>> Overshoot inner' )
!call Show ( [ fM, fI, fC, fO, fP ], '>>> fM, fI, fC, fO, fP' )

                  c1  =  ( fC - fM ) / ( xAC - xAM )
                  c0  =  fC  -  c1 * xAC

                end if  !-- Overshoot inner

                !-- Overshoot at outer face?
                !   Reduce slope.
                if ( c1 * ( fP - fO )  <  0.0_KDR ) then
!call Show ( '>>> Overshoot outer' )
!call Show ( [ fM, fI, fC, fO, fP ], '>>> fM, fI, fC, fO, fP' )

                  c1  =  ( fP - fC ) / ( xAP - xAC )
                  c0  =  fC  -  c1 * xAC

                end if  !-- Overshoot outer

              end if  !-- Local extremum

              F_IR ( iV, jV, kV, iF_R )  &
                =  c0  +  c1 * xI

              F_IL ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ), iF_R )  &
                =  c0  +  c1 * xO

!call Show ( '>>> Final values' )
!call Show ( [ fM, F_IR ( iV, jV, kV, iS ), fC, &
!              F_IL ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ), iS ), fP ], &
!            '>>> fM, fI, fC, fO, fP' )
            end do !-- iV
          end do !-- jV
        end do !-- kV
      end do !-- iS
      !$OMP end OMP_TARGET_DIRECTIVE parallel do simd

    else !-- use host
              
      !$OMP parallel do collapse ( 4 ) &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP private ( iF, iF_R, iaVP, iaVM, fM, fC, fP, fI, fO ) &
      !$OMP private ( xAM, xAC, xAP, xI, xO, c0, c1 )
      do iS  =  1,  size ( iaSlctd )
        do kV  =  lV ( 3 ),  uV ( 3 ) 
          do jV  =  lV ( 2 ),  uV ( 2 )
            do iV  =  lV ( 1 ),  uV ( 1 )

              iF    =  iaSlctd   ( iS )
              iF_R  =  iaSlctd_R ( iS )

              iaVP  =  [ iV, jV, kV ]  +  iaS
              iaVM  =  [ iV, jV, kV ]  -  iaS

              fM  =  F ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ), iF )
              fC  =  F ( iV, jV, kV, iF )
              fP  =  F ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ), iF )

              xAM  =  XA ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ) )
              xAC  =  XA ( iV, jV, kV )
              xAP  =  XA ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )

              xI  =  X ( iV, jV, kV )  -  0.5 * dX ( iV, jV, kV )
              xO  =  X ( iV, jV, kV )  +  0.5 * dX ( iV, jV, kV )

!call Show ( iV, '>>> iV' )
!call Show ( [ fM, fC, fP ], '>>> fM, fC, fP' )
              !-- Local extremum of cell average values? 
              !   Then reconstruction is constant.
              if ( ( fC - fM ) * ( fP - fC )  <  0.0_KDR ) then
!call Show ( '>>> Local extremum' )

                c1  =  0.0_KDR
                c0  =  fC

              else  !-- Linear reconstruction

                c1  =  ( fP - fM ) / ( xAP - xAM )
                c0  =  fC  -  c1 * xAC

                fI  =  c0  +  c1 * xI
                fO  =  c0  +  c1 * xO

                !-- Overshoot at inner face?
                !   Reduce slope.
                if ( c1 * ( fI - fM )  <  0.0_KDR ) then
!call Show ( '>>> Overshoot inner' )
!call Show ( [ fM, fI, fC, fO, fP ], '>>> fM, fI, fC, fO, fP' )

                  c1  =  ( fC - fM ) / ( xAC - xAM )
                  c0  =  fC  -  c1 * xAC

                end if  !-- Overshoot inner

                !-- Overshoot at outer face?
                !   Reduce slope.
                if ( c1 * ( fP - fO )  <  0.0_KDR ) then
!call Show ( '>>> Overshoot outer' )
!call Show ( [ fM, fI, fC, fO, fP ], '>>> fM, fI, fC, fO, fP' )

                  c1  =  ( fP - fC ) / ( xAP - xAC )
                  c0  =  fC  -  c1 * xAC

                end if  !-- Overshoot outer

              end if  !-- Local extremum

              F_IR ( iV, jV, kV, iF_R )  &
                =  c0  +  c1 * xI

              F_IL ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ), iF_R )  &
                =  c0  +  c1 * xO

!call Show ( '>>> Final values' )
!call Show ( [ fM, F_IR ( iV, jV, kV, iS ), fC, &
!              F_IL ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ), iS ), fP ], &
!            '>>> fM, fI, fC, fO, fP' )
            end do !-- iV
          end do !-- jV
        end do !-- kV
      end do !-- iS
      !$OMP end parallel do
      
    end if !-- UseDevice
        
  end procedure ComputeLinear_CGS_Kernel


  module procedure ComputeParabolic_CGS_Kernel

    integer ( KDI ) :: &
      iS, &
      iF, &
      iF_R, &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVP, iaVM, &
      lV, uV
    real ( KDR ) :: &
        fM,   fC,   fP,  &  !-- f_Minus, f_Center, f_Plus
       xAM,  xAC,  xAP,  & 
      x2AM, x2AC, x2AP,  & 
        fI,   fO,        &  !-- F_Inner, F_Outer
        xI,   xO,   xE,  &  !-- X_Inner, X_Outer, X_Extremum
       xIM,   xC,  xOP,  &
         d,              &  !-- Determinant / Denominator
        c0,   c1,   c2,  &  !-- Parabola coefficients,
      c2_S,   SqrtTiny      !-- c2_Safe
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice  =  .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice  =  UseDeviceOption

    lV  =  1
    where ( shape ( X )  >  1 )
      lV  =  oV
    end where
    
    uV  =  1
    where ( shape ( X )  >  1 )
      uV  =  shape ( X )  -  oV
    end where
    uV ( iD )  =  size ( X, dim = iD )  -  oV  +  1 
      
    iaS  =  0
    iaS ( iD )  =  1
    
    SqrtTiny  =  tiny ( 0.0_KDR )

    if ( UseDevice ) then
    
      !$OMP OMP_TARGET_DIRECTIVE parallel do simd collapse ( 4 ) &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP private ( iF, iF_R, iaVP, iaVM, fM, fC, fP, fI, fO ) &
      !$OMP private ( xAM, xAC, xAP, x2AM, x2AC, x2AP ) &
      !$OMP private ( xI, xO, xE, xIM, xC, xOP ) &
      !$OMP private ( c0, c1, c2, c2_S, d ) &
      !$OMP firstprivate ( SqrtTiny )
      do iS  =  1,  size ( iaSlctd )
        do kV  =  lV ( 3 ),  uV ( 3 ) 
          do jV  =  lV ( 2 ),  uV ( 2 )
            do iV  =  lV ( 1 ),  uV ( 1 )

              iF    =  iaSlctd   ( iS )
              iF_R  =  iaSlctd_R ( iS )

              iaVP  =  [ iV, jV, kV ]  +  iaS
              iaVM  =  [ iV, jV, kV ]  -  iaS

              fM  =  F ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ), iF )
              fC  =  F ( iV, jV, kV, iF )
              fP  =  F ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ), iF )

              xAM  =  XA ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ) )
              xAC  =  XA ( iV, jV, kV )
              xAP  =  XA ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )

              x2AM  =  X2A ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ) )
              x2AC  =  X2A ( iV, jV, kV )
              x2AP  =  X2A ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )

              xI  =  X ( iV, jV, kV )  -  0.5 * dX ( iV, jV, kV )
              xO  =  X ( iV, jV, kV )  +  0.5 * dX ( iV, jV, kV )

              xC  =  X ( iV, jV, kV )

              xIM  =  X ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ) )  &
                      -  0.5 * dX ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ) )
              xOP  =  X ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )  &
                      +  0.5 * dX ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )

!call Show ( iV, '>>> iV' )
!call Show ( [ fM, fC, fP ], '>>> fM, fC, fP' )
              !-- Local extremum of cell average values? 
              !   Then reconstruction is constant.
              if ( ( fC - fM ) * ( fP - fC )  <=  0.0_KDR ) then
!call Show ( '>>> Local extremum' )

                fI  =  fC
                fO  =  fC

              else  !-- Parabolic reconstruction
                 
                !-- Fortran expressions from Mathematica notebook 
                !   "FORNAX reconstruction 5.nb"

                !-- First parabola

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

                c2_S  =  sign ( max ( abs ( c2 ), SqrtTiny ), c2 )
                  xE  =  - c1 / ( 2.0 * c2_S )

                fI  =  c0  +  c1 * xI  +  c2 * xI**2
                fO  =  c0  +  c1 * xO  +  c2 * xO**2

                !-- Extremum near inner face?
                !   New inner parabola, flat slope at inner face, revise fI
                if ( xE  >  xIM  .and.  xE  <=  xC ) then
!call Show ( '>>> Extremum near inner face' )
!call Show ( [ fM, fI, fC, fO, fP ], '>>> fM, fI, fC, fO, fP' )

                   d  =  -x2AC + x2AM  +  2 * ( xAC - xAM ) * xIM

                  c0  =  (    fM * ( -x2AC  +  2 * xAC * xIM ) &
                           +  fC * (  x2AM  -  2 * xAM * xIM ) )  /  d

                  c1  =  2 * ( fC - fM ) * xIM  /  d

                  c2  =  ( -fC + fM )  /  d

                  fI  =  c0  +  c1 * xI  +  c2 * xI**2

!call Show ( '>>> Revised fI' )
!call Show ( [ fM, fI, fC ], '>>> fM, fI, fC' )

                !-- Extremum outer face?
                !   New outer parabola, flat slope at outer face, revise fO
                else if ( xE  >  xC  .and.  xE  <  xOP ) then
!call Show ( '>>> Extremum near outer face' )
!call Show ( [ fM, fI, fC, fO, fP ], '>>> fM, fI, fC, fO, fP' )

                   d  =  x2AC - x2AP  +  2 * ( -xAC + xAP ) * xOP
  
                  c0  =  (    fP * (  x2AC  -  2 * xAC * xOP ) &
                           +  fC * ( -x2AP  +  2 * xAP * xOP ) )  /  d

                  c1  =  -2 * ( fC - fP ) * xOP  /  d

                  c2  =  ( fC - fP )  /  d

                  fO  =  c0  +  c1 * xO  +  c2 * xO**2

!call Show ( '>>> Revised fO' )
!call Show ( [ fC, fO, fP ], '>>> fC, fO, fP' )

                end if  !-- First parabola extremum

                !-- Second parabola

                d  =   ( xI - xO ) * ( x2AC  +  xI * xO  -  xAC * ( xI + xO ) )

                c0  =  (    fO * xI * (  x2AC - xAC * xI ) &
                         +  fC * xI * ( xI - xO ) * xO &
                         +  fI * xO * ( -x2AC + xAC * xO ) )  /  d

                c1  =  (    fO * ( -x2AC  +  xI**2 ) &
                         +  fI * (  x2AC  -  xO**2 ) &
                         +  fC * ( -xI**2 +  xO**2 ) )  /  d

                c2  =  (    fO * (  xAC - xI ) &
                         +  fC * (  xI  - xO ) &
                         +  fI * ( -xAC + xO ) )  /  d

                c2_S  =  sign ( max ( abs ( c2 ), SqrtTiny ), c2 )
                  xE  =  - c1 / ( 2.0 * c2_S )

                !   Need c1, c2 to check for extremum

                ! !-- Not necessary to reset, just a consistency check
                ! fI  =  c0  +  c1 * xI  +  c2 * xI**2
                ! fO  =  c0  +  c1 * xO  +  c2 * xO**2

                !-- Extremum near inner face?
                !   New parabola, flat slope at inner face, revise fO
                if ( xE  >  xI  .and.  xE  <=  xC ) then
!call Show ( '>>> Extremum near inner face' )
!call Show ( [ fI, fC, fO ], '>>> fI, fC, fO' )

                  d  =  x2AC + xI * ( -2 * xAC  +  xI )

                  c0  =  ( fC * xI**2  +  fI * ( x2AC  -  2 * xAC * xI ) ) &
                         /  d

                  c1  = -2 * ( fC - fI ) * xI  /  d

                  c2  =  ( fC - fI )  /  d

                  fO  =  c0  +  c1 * xO  +  c2 * xO**2

                !-- Extremum near outer face?
                !   New parabola, flat slope at outer face, revise fI
                else if ( xE  >  xC  .and.  xE  <  xO ) then
!call Show ( '>>> Extremum near outer face' )
!call Show ( [ fI, fC, fO ], '>>> fI, fC, fO' )

                  d  =  x2AC + xO * ( -2 * xAC  +  xO )

                  c0  =  ( fC * xO**2  +  fO * ( x2AC  -  2 * xAC * xO ) ) &
                         /  d

                  c1  =  -2 * ( fC - fO ) * xO  /  d

                  c2  =  ( fC - fO )  /  d

                  fI  =  c0  +  c1 * xI  +  c2 * xI**2

                end if  !-- Second parabola extremum

              end if  !-- Local extremum

              F_IR ( iV, jV, kV, iF_R )  &
                =  fI
              F_IL ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ), iF_R )  &
                =  fO

!call Show ( '>>> Final values' )
!call Show ( [ fM, F_IR ( iV, jV, kV, iS ), fC, &
!              F_IL ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ), iS ), fP ], &
!            '>>> fM, fI, fC, fO, fP' )
            end do !-- iV
          end do !-- jV
        end do !-- kV
      end do !-- iS
      !$OMP end OMP_TARGET_DIRECTIVE parallel do simd

    else !-- use host
              
      !$OMP parallel do collapse ( 4 ) &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP private ( iF, iF_R, iaVP, iaVM, fM, fC, fP, fI, fO ) &
      !$OMP private ( xAM, xAC, xAP, x2AM, x2AC, x2AP ) &
      !$OMP private ( xI, xO, xE, xIM, xC, xOP ) &
      !$OMP private ( c0, c1, c2, c2_S, d ) &
      !$OMP firstprivate ( SqrtTiny )
      do iS  =  1,  size ( iaSlctd )
        do kV  =  lV ( 3 ),  uV ( 3 ) 
          do jV  =  lV ( 2 ),  uV ( 2 )
            do iV  =  lV ( 1 ),  uV ( 1 )

              iF    =  iaSlctd   ( iS )
              iF_R  =  iaSlctd_R ( iS )

              iaVP  =  [ iV, jV, kV ]  +  iaS
              iaVM  =  [ iV, jV, kV ]  -  iaS

              fM  =  F ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ), iF )
              fC  =  F ( iV, jV, kV, iF )
              fP  =  F ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ), iF )

              xAM  =  XA ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ) )
              xAC  =  XA ( iV, jV, kV )
              xAP  =  XA ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )

              x2AM  =  X2A ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ) )
              x2AC  =  X2A ( iV, jV, kV )
              x2AP  =  X2A ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )

              xI  =  X ( iV, jV, kV )  -  0.5 * dX ( iV, jV, kV )
              xO  =  X ( iV, jV, kV )  +  0.5 * dX ( iV, jV, kV )

              xC  =  X ( iV, jV, kV )

              xIM  =  X ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ) )  &
                      -  0.5 * dX ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ) )
              xOP  =  X ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )  &
                      +  0.5 * dX ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )

!call Show ( iV, '>>> iV' )
!call Show ( [ fM, fC, fP ], '>>> fM, fC, fP' )
              !-- Local extremum of cell average values? 
              !   Then reconstruction is constant.
              if ( ( fC - fM ) * ( fP - fC )  <=  0.0_KDR ) then
!call Show ( '>>> Local extremum' )

                fI  =  fC
                fO  =  fC

              else  !-- Parabolic reconstruction
                 
                !-- Fortran expressions from Mathematica notebook 
                !   "FORNAX reconstruction 5.nb"

                !-- First parabola

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

                c2_S  =  sign ( max ( abs ( c2 ), SqrtTiny ), c2 )
                  xE  =  - c1 / ( 2.0 * c2_S )

                fI  =  c0  +  c1 * xI  +  c2 * xI**2
                fO  =  c0  +  c1 * xO  +  c2 * xO**2

                !-- Extremum near inner face?
                !   New inner parabola, flat slope at inner face, revise fI
                if ( xE  >  xIM  .and.  xE  <=  xC ) then
!call Show ( '>>> Extremum near inner face' )
!call Show ( [ fM, fI, fC, fO, fP ], '>>> fM, fI, fC, fO, fP' )

                   d  =  -x2AC + x2AM  +  2 * ( xAC - xAM ) * xIM

                  c0  =  (    fM * ( -x2AC  +  2 * xAC * xIM ) &
                           +  fC * (  x2AM  -  2 * xAM * xIM ) )  /  d

                  c1  =  2 * ( fC - fM ) * xIM  /  d

                  c2  =  ( -fC + fM )  /  d

                  fI  =  c0  +  c1 * xI  +  c2 * xI**2

!call Show ( '>>> Revised fI' )
!call Show ( [ fM, fI, fC ], '>>> fM, fI, fC' )

                !-- Extremum outer face?
                !   New outer parabola, flat slope at outer face, revise fO
                else if ( xE  >  xC  .and.  xE  <  xOP ) then
!call Show ( '>>> Extremum near outer face' )
!call Show ( [ fM, fI, fC, fO, fP ], '>>> fM, fI, fC, fO, fP' )

                   d  =  x2AC - x2AP  +  2 * ( -xAC + xAP ) * xOP
  
                  c0  =  (    fP * (  x2AC  -  2 * xAC * xOP ) &
                           +  fC * ( -x2AP  +  2 * xAP * xOP ) )  /  d

                  c1  =  -2 * ( fC - fP ) * xOP  /  d

                  c2  =  ( fC - fP )  /  d

                  fO  =  c0  +  c1 * xO  +  c2 * xO**2

!call Show ( '>>> Revised fO' )
!call Show ( [ fC, fO, fP ], '>>> fC, fO, fP' )

                end if  !-- First parabola extremum

                !-- Second parabola

                d  =   ( xI - xO ) * ( x2AC  +  xI * xO  -  xAC * ( xI + xO ) )

                c0  =  (    fO * xI * (  x2AC - xAC * xI ) &
                         +  fC * xI * ( xI - xO ) * xO &
                         +  fI * xO * ( -x2AC + xAC * xO ) )  /  d

                c1  =  (    fO * ( -x2AC  +  xI**2 ) &
                         +  fI * (  x2AC  -  xO**2 ) &
                         +  fC * ( -xI**2 +  xO**2 ) )  /  d

                c2  =  (    fO * (  xAC - xI ) &
                         +  fC * (  xI  - xO ) &
                         +  fI * ( -xAC + xO ) )  /  d

                c2_S  =  sign ( max ( abs ( c2 ), SqrtTiny ), c2 )
                  xE  =  - c1 / ( 2.0 * c2_S )

                !   Need c1, c2 to check for extremum

                ! !-- Not necessary to reset, just a consistency check
                ! fI  =  c0  +  c1 * xI  +  c2 * xI**2
                ! fO  =  c0  +  c1 * xO  +  c2 * xO**2

                !-- Extremum near inner face?
                !   New parabola, flat slope at inner face, revise fO
                if ( xE  >  xI  .and.  xE  <=  xC ) then
!call Show ( '>>> Extremum near inner face' )
!call Show ( [ fI, fC, fO ], '>>> fI, fC, fO' )

                  d  =  x2AC + xI * ( -2 * xAC  +  xI )

                  c0  =  ( fC * xI**2  +  fI * ( x2AC  -  2 * xAC * xI ) ) &
                         /  d

                  c1  = -2 * ( fC - fI ) * xI  /  d

                  c2  =  ( fC - fI )  /  d

                  fO  =  c0  +  c1 * xO  +  c2 * xO**2

                !-- Extremum near outer face?
                !   New parabola, flat slope at outer face, revise fI
                else if ( xE  >  xC  .and.  xE  <  xO ) then
!call Show ( '>>> Extremum near outer face' )
!call Show ( [ fI, fC, fO ], '>>> fI, fC, fO' )

                  d  =  x2AC + xO * ( -2 * xAC  +  xO )

                  c0  =  ( fC * xO**2  +  fO * ( x2AC  -  2 * xAC * xO ) ) &
                         /  d

                  c1  =  -2 * ( fC - fO ) * xO  /  d

                  c2  =  ( fC - fO )  /  d

                  fI  =  c0  +  c1 * xI  +  c2 * xI**2

                end if  !-- Second parabola extremum

              end if  !-- Local extremum

              F_IR ( iV, jV, kV, iF_R )  &
                =  fI
              F_IL ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ), iF_R )  &
                =  fO

!call Show ( '>>> Final values' )
!call Show ( [ fM, F_IR ( iV, jV, kV, iS ), fC, &
!              F_IL ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ), iS ), fP ], &
!            '>>> fM, fI, fC, fO, fP' )
            end do !-- iV
          end do !-- jV
        end do !-- kV
      end do !-- iS
      !$OMP end parallel do
      
    end if !-- UseDevice
        
  end procedure ComputeParabolic_CGS_Kernel


end submodule Reconstruction_Kernel
