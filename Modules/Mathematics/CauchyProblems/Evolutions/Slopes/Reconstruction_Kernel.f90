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
      !$OMP private ( xAM, xAC, xAP, x2AM, x2AC, x2AP, xI, xO, xE ) &
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

!call Show ( iV, '>>> iV' )
!call Show ( [ fM, fC, fP ], '>>> fM, fC, fP' )
              !-- Local extremum of cell average values? 
              !   Then reconstruction is constant.
              if ( ( fC - fM ) * ( fP - fC )  <  0.0_KDR ) then
!call Show ( '>>> Local extremum' )

                c2  =  0.0_KDR
                c1  =  0.0_KDR
                c0  =  fC

              else  !-- Parabolic reconstruction

                !-- First parabola

                d  =    ( x2AM - x2AP ) * xAC  &
                      + ( x2AP - x2AC ) * xAM  &
                      + ( x2AC - x2AM ) * xAP

                c0  =  (   fP * ( x2AM * xAC  -  x2AC * xAM )  &
                         + fM * ( x2AC * xAP  -  x2AP * xAC )  &
                         + fC * ( x2AP * xAM  -  x2AM * xAP ) )  /  d

                c1  =  (   fP * ( x2AC - x2AM )  &
                         + fC * ( x2AM - x2AP )  &
                         + fM * ( x2AP - x2AC ) )  /  d

                c2  =  (   fP * ( xAM - xAC )  &
                         + fM * ( xAC - xAP )  &
                         + fC * ( xAP - xAM ) )  /  d

                c2_S  =  sign ( max ( abs ( c2 ), SqrtTiny ), c2 )
                  xE  =  - c1 / ( 2.0 * c2_S )

                fI  =  c0  +  c1 * xI  +  c2 * xI**2
                fO  =  c0  +  c1 * xO  +  c2 * xO**2

                !-- Overshoot at inner face?
                !   New inner parabola, revise fI
                if ( xE  >  xAM  .and.  xE  <=  xAC ) then
!call Show ( '>>> Overshoot inner' )
!call Show ( [ fM, fI, fC, fO, fP ], '>>> fM, fI, fC, fO, fP' )

                   d  =  ( x2AC - x2AM )  +  2.0 * ( xAM - xAC ) * xAM

                  c0  =  (    ( fM * x2AC  -  fC * x2AM )  &
                           +  2.0 * ( fC * xAM  -  fM * xAC ) * xAM )  /  d

                  c1  =  -2.0 * ( fC - fM ) * xAM  /  d

                  c2  =  ( fC - fM )  /  d

                  fI  =  c0  +  c1 * xI  +  c2 * xI**2

!call Show ( '>>> Revised fI' )
!call Show ( [ fM, fI, fC ], '>>> fM, fI, fC' )

                !-- Overshoot at outer face?
                !   New outer parabola, revise fO
                else if ( xE  >  xAC  .and.  xE  <  xAP ) then
!call Show ( '>>> Overshoot outer' )
!call Show ( [ fM, fI, fC, fO, fP ], '>>> fM, fI, fC, fO, fP' )

                   d  =  ( x2AC - x2AP )  +  2.0 * ( xAP - xAC ) * xAP
  
                  c0  =  (    ( fP * x2AC  -  fC * x2AP )  &
                           +  2.0 * ( fC * xAP  -  fP * xAC ) * xAP )  /  d

                  c1  =  -2.0 * ( fC - fP ) * xAP  /  d

                  c2  =  ( fC - fP )  /  d

                  fO  =  c0  +  c1 * xO  +  c2 * xO**2

!call Show ( '>>> Revised fO' )
!call Show ( [ fC, fO, fP ], '>>> fC, fO, fP' )

                end if  !-- First parabola extremum

                !-- Second parabola

                d  =  ( xI - xO ) * ( x2AC  +  xI * xO  -  xAC * ( xI + xO ) )

                c0  =  (    fO * xI * ( x2AC  -  xAC * xI )  &
                         +  fC * xI * xO * ( xI  -  xO )  &
                         +  fI * xO * ( xAC * xO  -  x2AC ) )  /  d

                c1  =  (    fO * ( xI**2  -  x2AC )  &
                         +  fI * (  x2AC  -  xO**2 )  &
                         +  fC * ( xO**2  -  xI**2 ) )  /  d

                c2  =  (    fO * ( xAC -  xI )  &
                         +  fC * (  xI -  xO )  &
                         +  fI * (  xO - xAC ) )  /  d

                c2_S  =  sign ( max ( abs ( c2 ), SqrtTiny ), c2 )
                  xE  =  - c1 / ( 2.0 * c2_S )

                fI  =  c0  +  c1 * xI  +  c2 * xI**2
                fO  =  c0  +  c1 * xO  +  c2 * xO**2

                !-- Extremum near inner face?
                !   New parabola, flat slope at inner face
                if ( xE  >  xI  .and.  xE  <=  xAC ) then
!call Show ( '>>> Extremum near inner face' )
!call Show ( [ fI, fC, fO ], '>>> fI, fC, fO' )

                  d  =  x2AC  -  2.0 * xAC * xI  +  xI**2

                  c0  =  ( fI * x2AC  -  2.0 * fI * xAC * xI  +  fC * xI**2 ) &
                         /  d

                  c1  =  -2.0 * ( fC - fI ) * xI  /  d

                  c2  =  ( fC - fI )  /  d

                !-- Extremum near outer face?
                !   New parabola, flat slope at outer face
                else if ( xE  >  xAC  .and.  xE  <  xO ) then
!call Show ( '>>> Extremum near outer face' )
!call Show ( [ fI, fC, fO ], '>>> fI, fC, fO' )

                  d  =  x2AC  -  2.0 * xAC * xO  +  xO**2

                  c0  =  ( fO * x2AC  -  2.0 * fO *xAC * xO  +  fC * xO**2 ) &
                         /  d

                  c1  =  -2.0 * ( fC - fO ) * xO  /  d

                  c2  =  ( fC - fO )  /  d

                end if  !-- Second parabola extremum

              end if  !-- Local extremum

              F_IR ( iV, jV, kV, iF_R )  &
                =  c0  +  c1 * xI  +  c2 * xI**2

              F_IL ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ), iF_R )  &
                =  c0  +  c1 * xO  +  c2 * xO**2

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
      !$OMP private ( xAM, xAC, xAP, x2AM, x2AC, x2AP, xI, xO, xE ) &
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

!call Show ( iV, '>>> iV' )
!call Show ( [ fM, fC, fP ], '>>> fM, fC, fP' )
              !-- Local extremum of cell average values? 
              !   Then reconstruction is constant.
              if ( ( fC - fM ) * ( fP - fC )  <  0.0_KDR ) then
!call Show ( '>>> Local extremum' )

                c2  =  0.0_KDR
                c1  =  0.0_KDR
                c0  =  fC

              else  !-- Parabolic reconstruction

                !-- First parabola

                d  =    ( x2AM - x2AP ) * xAC  &
                      + ( x2AP - x2AC ) * xAM  &
                      + ( x2AC - x2AM ) * xAP

                c0  =  (   fP * ( x2AM * xAC  -  x2AC * xAM )  &
                         + fM * ( x2AC * xAP  -  x2AP * xAC )  &
                         + fC * ( x2AP * xAM  -  x2AM * xAP ) )  /  d

                c1  =  (   fP * ( x2AC - x2AM )  &
                         + fC * ( x2AM - x2AP )  &
                         + fM * ( x2AP - x2AC ) )  /  d

                c2  =  (   fP * ( xAM - xAC )  &
                         + fM * ( xAC - xAP )  &
                         + fC * ( xAP - xAM ) )  /  d

                c2_S  =  sign ( max ( abs ( c2 ), SqrtTiny ), c2 )
                  xE  =  - c1 / ( 2.0 * c2_S )

                fI  =  c0  +  c1 * xI  +  c2 * xI**2
                fO  =  c0  +  c1 * xO  +  c2 * xO**2

                !-- Overshoot at inner face?
                !   New inner parabola, revise fI
                if ( xE  >  xAM  .and.  xE  <=  xAC ) then
!call Show ( '>>> Overshoot inner' )
!call Show ( [ fM, fI, fC, fO, fP ], '>>> fM, fI, fC, fO, fP' )

                   d  =  ( x2AC - x2AM )  +  2.0 * ( xAM - xAC ) * xAM

                  c0  =  (    ( fM * x2AC  -  fC * x2AM )  &
                           +  2.0 * ( fC * xAM  -  fM * xAC ) * xAM )  /  d

                  c1  =  -2.0 * ( fC - fM ) * xAM  /  d

                  c2  =  ( fC - fM )  /  d

                  fI  =  c0  +  c1 * xI  +  c2 * xI**2

!call Show ( '>>> Revised fI' )
!call Show ( [ fM, fI, fC ], '>>> fM, fI, fC' )

                !-- Overshoot at outer face?
                !   New outer parabola, revise fO
                else if ( xE  >  xAC  .and.  xE  <  xAP ) then
!call Show ( '>>> Overshoot outer' )
!call Show ( [ fM, fI, fC, fO, fP ], '>>> fM, fI, fC, fO, fP' )

                   d  =  ( x2AC - x2AP )  +  2.0 * ( xAP - xAC ) * xAP
  
                  c0  =  (    ( fP * x2AC  -  fC * x2AP )  &
                           +  2.0 * ( fC * xAP  -  fP * xAC ) * xAP )  /  d

                  c1  =  -2.0 * ( fC - fP ) * xAP  /  d

                  c2  =  ( fC - fP )  /  d

                  fO  =  c0  +  c1 * xO  +  c2 * xO**2

!call Show ( '>>> Revised fO' )
!call Show ( [ fC, fO, fP ], '>>> fC, fO, fP' )

                end if  !-- First parabola extremum

                !-- Second parabola

                d  =  ( xI - xO ) * ( x2AC  +  xI * xO  -  xAC * ( xI + xO ) )

                c0  =  (    fO * xI * ( x2AC  -  xAC * xI )  &
                         +  fC * xI * xO * ( xI  -  xO )  &
                         +  fI * xO * ( xAC * xO  -  x2AC ) )  /  d

                c1  =  (    fO * ( xI**2  -  x2AC )  &
                         +  fI * (  x2AC  -  xO**2 )  &
                         +  fC * ( xO**2  -  xI**2 ) )  /  d

                c2  =  (    fO * ( xAC -  xI )  &
                         +  fC * (  xI -  xO )  &
                         +  fI * (  xO - xAC ) )  /  d

                c2_S  =  sign ( max ( abs ( c2 ), SqrtTiny ), c2 )
                  xE  =  - c1 / ( 2.0 * c2_S )

                fI  =  c0  +  c1 * xI  +  c2 * xI**2
                fO  =  c0  +  c1 * xO  +  c2 * xO**2

                !-- Extremum near inner face?
                !   New parabola, flat slope at inner face
                if ( xE  >  xI  .and.  xE  <=  xAC ) then
!call Show ( '>>> Extremum near inner face' )
!call Show ( [ fI, fC, fO ], '>>> fI, fC, fO' )

                  d  =  x2AC  -  2.0 * xAC * xI  +  xI**2

                  c0  =  ( fI * x2AC  -  2.0 * fI * xAC * xI  +  fC * xI**2 ) &
                         /  d

                  c1  =  -2.0 * ( fC - fI ) * xI  /  d

                  c2  =  ( fC - fI )  /  d

                !-- Extremum near outer face?
                !   New parabola, flat slope at outer face
                else if ( xE  >  xAC  .and.  xE  <  xO ) then
!call Show ( '>>> Extremum near outer face' )
!call Show ( [ fI, fC, fO ], '>>> fI, fC, fO' )

                  d  =  x2AC  -  2.0 * xAC * xO  +  xO**2

                  c0  =  ( fO * x2AC  -  2.0 * fO *xAC * xO  +  fC * xO**2 ) &
                         /  d

                  c1  =  -2.0 * ( fC - fO ) * xO  /  d

                  c2  =  ( fC - fO )  /  d

                end if  !-- Second parabola extremum

              end if  !-- Local extremum

              F_IR ( iV, jV, kV, iF_R )  &
                =  c0  +  c1 * xI  +  c2 * xI**2

              F_IL ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ), iF_R )  &
                =  c0  +  c1 * xO  +  c2 * xO**2

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
