#include "Preprocessor"

submodule ( IncrementDivergence_FV__Form ) IncrementDivergence_FV__Kernel

  use Basics
  
  implicit none
  
contains

  module procedure ComputeReconstructionLinear_CSL_Kernel

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVS, &
      lV, uV
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

!    V_IL  =  cshift ( V  +  0.5_KDR * dX * dVdX, shift = -1, dim = iD )

!    V_IR  =  V  -  0.5_KDR * dX * dVdX

    lV = 1
    where ( shape ( V ) > 1 )
      lV = oV + 1
    end where
    
    uV = 1
    where ( shape ( V ) > 1 )
      uV = shape ( V ) - oV
    end where
    uV ( iD ) = size ( V, dim = iD ) - oV + 1 
      
    iaS = 0
    iaS ( iD ) = -1
    
    if ( UseDevice ) then
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iV, jV, kV, iaVS )
      do kV = lV ( 3 ), uV ( 3 ) 
        do jV = lV ( 2 ), uV ( 2 )
          do iV = lV ( 1 ), uV ( 1 )

            iaVS = [ iV, jV, kV ] + iaS

            V_IL ( iV, jV, kV )  &
              =  V ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )  &
                 +  dX_R ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) &
                    *  dVdX ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
            
            V_IR ( iV, jV, kV )  &
              =  V ( iV, jV, kV )  &
                 -  dX_L ( iV, jV, kV )  *  dVdX ( iV, jV, kV )

          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
      
    else
              
      !$OMP  parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iV, jV, kV, iaVS )
      do kV = lV ( 3 ), uV ( 3 ) 
        do jV = lV ( 2 ), uV ( 2 )
          do iV = lV ( 1 ), uV ( 1 )

            iaVS = [ iV, jV, kV ] + iaS

            V_IL ( iV, jV, kV )  &
              =  V ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )  &
                 +  dX_R ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) &
                    *  dVdX ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
                    
            V_IR ( iV, jV, kV )  &
              =  V ( iV, jV, kV )  &
                 -  dX_L ( iV, jV, kV )  *  dVdX ( iV, jV, kV )

          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end parallel do
      
    end if
        
  end procedure ComputeReconstructionLinear_CSL_Kernel


!   module procedure ComputeReconstructionParabolic_CSL_Kernel

!     integer ( KDI ) :: &
!       iV, jV, kV  !-- iValue, etc.
!     integer ( KDI ), dimension ( 3 ) :: &
!       iaS, &         !-- iaShift
!       iaVP, iaVM, &  !-- iaValuePlus, iaValueMinus
!       lV, uV         !-- lowerValue, upperValue
!     real ( KDR ) :: &
!         vM,   vC,   vP, &  !-- V_Minus, V_Center, V_Plus
!         xM,   xC,   xP, &
!        xAM,  xAC,  xAP, & 
!       x2AM, x2AC, x2AP, &
!       xI, xO, &            !-- X_Inner, X_Outer
!       D, &                 !-- Denominator / Determinant
!       c0, c1, c2, &        !-- Parabola coefficients
!       SqrtTiny, c2_Safe, &
!       xE                   !-- X_Extremum
!     logical ( KDL ) :: &
!       UseDevice

!     SqrtTiny  =  tiny ( 0.0_KDR )

!     UseDevice = .false.
!     if ( present ( UseDeviceOption ) ) &
!       UseDevice = UseDeviceOption

!     lV = 1
!     where ( shape ( V ) > 1 )
!       lV = oV
!     end where
    
!     uV = 1
!     where ( shape ( V ) > 1 )
!       uV = shape ( V ) - oV
!     end where
!     uV ( iD ) = size ( V, dim = iD ) - oV + 1 
      
!     iaS = 0
!     iaS ( iD ) = 1
    
!     if ( UseDevice ) then
    
!     else
              
! !call Show ( '>>> New variable' )
!       !$OMP parallel do collapse ( 3 ) &
!       !$OMP schedule ( OMP_SCHEDULE_HOST ) &
!       !$OMP private ( iV, jV, kV, iaVP, iaVM, vM, vC, vP ) &
!       !$OMP private ( xM, xC, xP, xAM, xAC, xAP, x2AM, x2AC, x2AP, xI, xO ) &
!       !$OMP private ( D, c0, c1, c2, c2_Safe, xE ) &
!       !$OMP firstprivate ( SqrtTiny )
!       do kV = lV ( 3 ), uV ( 3 ) 
!         do jV = lV ( 2 ), uV ( 2 )
!           do iV = lV ( 1 ), uV ( 1 )

! !call Show ( '>>> Reconstructing cell' )
! !call Show ( [ iV, jV, kV ], '>>> iV, jV, kV' )
!             iaVP = [ iV, jV, kV ] + iaS
!             iaVM = [ iV, jV, kV ] - iaS

!             vM  =  V ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ) )
!             vC  =  V ( iV, jV, kV )
!             vP  =  V ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )
! !call Show ( [ vM, vC, vP ], '>>> vM, vC, vP' )

!             xM  =  X ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ) )
!             xC  =  X ( iV, jV, kV )
!             xP  =  X ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )
! !call Show ( [ xM, xC, xP ], '>>> xM, xC, xP' )

!             xAM  =  XA ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ) )
!             xAC  =  XA ( iV, jV, kV )
!             xAP  =  XA ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )

!             x2AM  =  X2A ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ) )
!             x2AC  =  X2A ( iV, jV, kV )
!             x2AP  =  X2A ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )

!             xI  =  X ( iV, jV, kV )  -  dX_L ( iV, jV, kV )
!             xO  =  X ( iV, jV, kV )  +  dX_R ( iV, jV, kV )

!             !-- Local extremum of central values? 
!             !   Then reconstruction is constant.
!             if ( ( vC - vM ) * ( vP - vC )  <=  0.0_KDR ) then
! !call Show ( '>>> Local extremum' )

!               V_IR (         iV,         jV,         kV )  =  vC
!               V_IL ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )  =  vC

!             else  !-- Parabolic reconstruction
               
!               !-- First parabola
! !call Show ( '>>> First parabola' )

!               D  =   ( x2AM - x2AP ) * xAC  &
!                    + ( x2AP - x2AC ) * xAM  &
!                    + ( x2AC - x2AM ) * xAP
! !call Show ( D, '>>> D' )

!               c0  =   (   vP * ( x2AM * xAC  -  x2AC * xAM )  &
!                         + vM * ( x2AC * xAP  -  x2AP * xAC )  &
!                         + vC * ( x2AP * xAM  -  x2AM * xAP ) )  /  D
! !call Show ( c0, '>>> c0' )

!               c1  =  (   vP * ( x2AC - x2AM )  &
!                        + vC * ( x2AM - x2AP )  &
!                        + vM * ( x2AP - x2AC ) )  /  D
! !call Show ( c1, '>>> c1' )

!               c2  =  (   vP * ( xAM - xAC )  &
!                        + vM * ( xAC - xAP )  &
!                        + vC * ( xAP - xAM ) )  /  D
! !call Show ( c2, '>>> c2' )

!               c2_Safe =  sign ( max ( abs ( c2 ), SqrtTiny ), c2 )
! !call Show ( c2_Safe, '>>> c2_Safe' )
!                   xE  =  - c1 / ( 2.0_KDR * c2_Safe )

! !call Show ( xE, '>>> xE' ) 
!               !-- Check for and eliminate local extrema in first parabola
!               if ( xE > xM .and. xE <= xC ) then
! !call Show ( '>>> Left extremum 1' )
! ! call Show ( [ vM, vC, vP ], '>>> vM, vC, vP' )
! ! call Show ( [ xM, xC, xP ], '>>> xM, xC, xP' )
! ! call Show ( xE, '>>> xE' )
! ! call Show ( xI, '>>> xI' )
! ! call Show ( [ c0, c1, c2 ], '>>> Old c0, c1, c2' )
! ! call Show (  c0  +  c1 * xI  +  c2 * xI**2, '>>> Old V_IR' )
! ! call Show (  c0  +  c1 * xE  +  c2 * xE**2, '>>> Old V_Extremum' )

!                 D  =  ( x2AC - x2AP ) + 2 * ( xAP - xAC ) * xM

!                 c0  =  (          vP * x2AC  -  vC * x2AP  &
!                          +  2 * ( vC * xAP   -  vP * xAC  ) * xM )  /  D

!                 c1  =  ( 2 * ( vP - vC ) * xM )  /  D

!                 c2  =  ( vC - vP )  /  D

! ! call Show ( [ c0, c1, c2 ], '>>> New c0, c1, c2' )
! ! call Show (  c0  +  c1 * xI  +  c2 * xI**2, '>>> New V_IR' )
! ! call Show (  - c1 / ( 2.0_KDR * c2 ), '>>> New xE' )
! ! call Show (  c0  +  c1 * xM  +  c2 * xM**2, '>>> New V_Extremum' )
! ! call Show ( c1 + 2 * c2 * xM, '>>> c1 + 2 * c2 * xM' )
! ! call Show (  c0  +  c1 * xAC  +  c2 * x2AC, '>>> c0  +  c1 * xAC  +  c2 * x2AC' )
! ! call Show (  c0  +  c1 * xAP  +  c2 * x2AP, '>>> c0  +  c1 * xAP  +  c2 * x2AP' )
! ! call PROGRAM_HEADER % Abort ( )
!               else if ( xE > xC .and. xE < xP ) then
! !call Show ( '>>> Right extremum 1' )
! ! call Show ( [ vM, vC, vP ], '>>> vM, vC, vP' )
! ! call Show ( [ xM, xC, xP ], '>>> xM, xC, xP' )
! ! call Show ( xE, '>>> xE' )
! ! call Show ( xO, '>>> xO' )
! ! call Show ( [ c0, c1, c2 ], '>>> Old c0, c1, c2' )
! ! call Show (  c0  +  c1 * xO  +  c2 * xO**2, '>>> Old V_IL_P' )
! ! call Show (  c0  +  c1 * xE  +  c2 * xE**2, '>>> Old V_Extremum' )

!                 D  =  ( x2AM - x2AC ) + 2 * ( xAC - xAM ) * xP

!                 c0  =  ( ( vC * x2AM - vM * x2AC )  &
!                          + 2 * ( vM * xAC - vC * xAM ) * xP )  /  D

!                 c1  =  ( 2 * ( vC - vM ) * xP ) / D

!                 c2  =  ( vM - vC ) / D

! ! call Show ( [ c0, c1, c2 ], '>>> New c0, c1, c2' )
! ! call Show (  c0  +  c1 * xO  +  c2 * xO**2, '>>> New V_IL_P' )
! ! call Show (  - c1 / ( 2.0_KDR * c2 ), '>>> New xE' )
! ! call Show (  c0  +  c1 * xP  +  c2 * xP**2, '>>> New V_Extremum' )
! ! call Show ( c1 + 2 * c2 * xP, '>>> c1 + 2 * c2 * xP' )
! ! call Show (  c0  +  c1 * xAM  +  c2 * x2AM, '>>> c0  +  c1 * xAM  +  c2 * x2AM' )
! ! call Show (  c0  +  c1 * xAC  +  c2 * x2AC, '>>> c0  +  c1 * xAC  +  c2 * x2AC' )
! ! call PROGRAM_HEADER % Abort ( )
!               end if  !-- Local extrema in first parabola

!               V_IR (         iV,         jV,         kV )  &
!                 =  c0  +  c1 * xI  +  c2 * xI**2

!               V_IL ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )  &
!                 =  c0  +  c1 * xO  +  c2 * xO**2

! !call Show ( '>>> Parabola' )
!             end if

! !call Show ( V_IR (         iV,         jV,         kV ), '>>> V_IR_M' )
! !call Show ( V_IL ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) ), '>>> V_IL_P' ) 

!           end do !-- iV
!         end do !-- jV
!       end do !-- kV
!       !$OMP end parallel do
      
!     end if
        
!   end procedure ComputeReconstructionParabolic_CSL_Kernel


  module procedure ComputeReconstructionParabolic_CSL_Kernel

    integer ( KDI ) :: &
      iV, jV, kV  !-- iValue, etc.
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &         !-- iaShift
      iaVP, iaVM, &  !-- iaValuePlus, iaValueMinus
      lV, uV         !-- lowerValue, upperValue
    real ( KDR ) :: &
        vM,   vC,   vP, &  !-- V_Minus, V_Center, V_Plus
        xM,   xC,   xP, &
       xAM,  xAC,  xAP, & 
      x2AM, x2AC, x2AP, &
      vI, vO, &             !-- V_Inner, V_Outer
      xI, xO, &             !-- X_Inner, X_Outer
      D, &                  !-- Denominator / Determinant
      c0, c1, c2, &         !-- Parabola coefficients
      SqrtTiny, c2_Safe, &
      xE                    !-- X_Extremum
    logical ( KDL ) :: &
      UseDevice

    SqrtTiny  =  tiny ( 0.0_KDR )

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    lV = 1
    where ( shape ( V ) > 1 )
      lV = oV
    end where
    
    uV = 1
    where ( shape ( V ) > 1 )
      uV = shape ( V ) - oV
    end where
    uV ( iD ) = size ( V, dim = iD ) - oV + 1 
      
    iaS = 0
    iaS ( iD ) = 1
    
    if ( UseDevice ) then
    
    else

      !$OMP parallel do collapse ( 3 ) &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP private ( iV, jV, kV, iaVP, iaVM, vM, vC, vP, vI, vO ) &
      !$OMP private ( xM, xC, xP, xAM, xAC, xAP, x2AM, x2AC, x2AP, xI, xO ) &
      !$OMP private ( D, c0, c1, c2, c2_Safe, xE ) &
      !$OMP firstprivate ( SqrtTiny )
      do kV = lV ( 3 ), uV ( 3 ) 
        do jV = lV ( 2 ), uV ( 2 )
          do iV = lV ( 1 ), uV ( 1 )

            iaVP = [ iV, jV, kV ] + iaS
            iaVM = [ iV, jV, kV ] - iaS

            vM  =  V ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ) )
            vC  =  V ( iV, jV, kV )
            vP  =  V ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )

            xM  =  X ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ) )
            xC  =  X ( iV, jV, kV )
            xP  =  X ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )

            xAM  =  XA ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ) )
            xAC  =  XA ( iV, jV, kV )
            xAP  =  XA ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )

            x2AM  =  X2A ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ) )
            x2AC  =  X2A ( iV, jV, kV )
            x2AP  =  X2A ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )

            xI  =  X ( iV, jV, kV )  -  dX_L ( iV, jV, kV )
            xO  =  X ( iV, jV, kV )  +  dX_R ( iV, jV, kV )

            !-- Local extremum of cell average values? 
            !   Then reconstruction is constant.
            if ( ( vC - vM ) * ( vP - vC )  <=  0.0_KDR ) then
!call Show ( '>>> Local extremum' )
!call Show ( [ vM, vC, vP ], '>>> vM, vC, vP' )
!call Show ( [ xM, xC, xP ], '>>> xM, xC, xP' )

              vI  =  vC
              vO  =  vC

!call Show ( [ vI, vO ], '>>> vI, vO' )
!call PROGRAM_HEADER % Abort ( )
            else  !-- Parabolic reconstruction
               
              !-- First parabola
!call Show ( '>>> First parabola' )
!call Show ( [ vM, vC, vP ], '>>> vM, vC, vP' )
!call Show ( [ xM, xC, xP ], '>>> xM, xC, xP' )

              D  =   ( x2AM - x2AP ) * xAC  &
                   + ( x2AP - x2AC ) * xAM  &
                   + ( x2AC - x2AM ) * xAP

              c0  =   (   vP * ( x2AM * xAC  -  x2AC * xAM )  &
                        + vM * ( x2AC * xAP  -  x2AP * xAC )  &
                        + vC * ( x2AP * xAM  -  x2AM * xAP ) )  /  D

              c1  =  (   vP * ( x2AC - x2AM )  &
                       + vC * ( x2AM - x2AP )  &
                       + vM * ( x2AP - x2AC ) )  /  D

              c2  =  (   vP * ( xAM - xAC )  &
                       + vM * ( xAC - xAP )  &
                       + vC * ( xAP - xAM ) )  /  D
!call Show ( [ c0, c1, c2 ], '>>> c0, c1, c2' )

              c2_Safe =  sign ( max ( abs ( c2 ), SqrtTiny ), c2 )
!call Show ( c2_Safe, '>>> c2_Safe' )
                  xE  =  - c1 / ( 2.0_KDR * c2_Safe )
!call Show ( xE, '>>> xE' ) 

              vI  =  c0  +  c1 * xI  +  c2 * xI**2
              vO  =  c0  +  c1 * xO  +  c2 * xO**2

!call Show ( c0  +  c1 * xAM  +  c2 * x2AM, '>>> vM ?' )
!call Show ( vI, '>>> vI' )
!call Show ( c0  +  c1 * xAC  +  c2 * x2AC, '>>> vC ?' )
!call Show ( vO, '>>> vO' )
!call Show ( c0  +  c1 * xAP  +  c2 * x2AP, '>>> vP ?' )
!call PROGRAM_HEADER % Abort ( )

              !-- Check for local extremum in first parabola, 
              !   and adjust vI or vO accordingly
              if ( xE > xM .and. xE <= xC ) then
!call Show ( '>>> Left extremum 1' )
!call Show ( [ vM, vC, vP ], '>>> vM, vC, vP' )
!call Show ( [ xM, xC, xP ], '>>> xM, xC, xP' )
!call Show ( xE, '>>> xE' )
!call Show ( vI, '>>> vI' )
!call Show ( xI, '>>> xI' )

                 D  =  ( x2AC - x2AM )  +  2 * ( xAM - xAC ) * xAM

                c0  =  (   ( vM * x2AC  -  vC * x2AM )  &
                         +  2 * ( vC * xAM  -  vM * xAC ) * xAM )  /  D

                c1  =  -2 * ( vC - vM ) * xAM  /  D

                c2  =  ( vC - vM )  /  D
!call Show ( [ c0, c1, c2 ], '>>> c0, c1, c2' )

                c2_Safe =  sign ( max ( abs ( c2 ), SqrtTiny ), c2 )
!call Show ( c2_Safe, '>>> c2_Safe' )
                    xE  =  - c1 / ( 2.0_KDR * c2_Safe )
!call Show ( xE, '>>> xE New' ) 

                vI  =  c0  +  c1 * xI  +  c2 * xI**2
!call Show (        c1        +  2 * c2 * xAM, '>>> d vM / dx == 0 ?' )
!call Show ( c0  +  c1 * xAM  +      c2 * x2AM, '>>> vM ?' )
!call Show ( vI, '>>> vI New' )
!call Show ( c0  +  c1 * xAC  +  c2 * x2AC, '>>> vC ?' )
!call PROGRAM_HEADER % Abort ( )

              else if ( xE > xC .and. xE < xP ) then
!call Show ( '>>> Right extremum 1' )
!call Show ( [ vM, vC, vP ], '>>> vM, vC, vP' )
!call Show ( [ xM, xC, xP ], '>>> xM, xC, xP' )
!call Show ( xE, '>>> xE' )
!call Show ( vO, '>>> vO' )
!call Show ( xO, '>>> xO' )

                 D  =  ( x2AC - x2AP )  +  2 * ( xAP - xAC ) * xAP

                c0  =  (   ( vP * x2AC  -  vC * x2AP )  &
                         + 2 * ( vC * xAP  -  vP * xAC ) * xAP )  /  D

                c1  =  -2 * ( vC - vP ) * xAP  /  D

                c2  =  ( vC - vP )  /  D
!call Show ( [ c0, c1, c2 ], '>>> c0, c1, c2' )

                c2_Safe =  sign ( max ( abs ( c2 ), SqrtTiny ), c2 )
!call Show ( c2_Safe, '>>> c2_Safe' )
                    xE  =  - c1 / ( 2.0_KDR * c2_Safe )
!call Show ( xE, '>>> xE New' ) 

                vO  =  c0  +  c1 * xO  +  c2 * xO**2
!call Show ( c0  +  c1 * xAC  +  c2 * x2AC, '>>> vC ?' )
!call Show ( vO, '>>> vO New' )
!call Show ( c0  +  c1 * xAP  +      c2 * x2AP, '>>> vP ?' )
!call Show (        c1        +  2 * c2 * xAP, '>>> d vP / dx == 0 ?' )
!call PROGRAM_HEADER % Abort ( )

              end if  !-- Local extrema in first parabola

           end if !-- Local extremum of cell average values

            V_IR (         iV,         jV,         kV )  =  vI
            V_IL ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )  =  vO

          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end parallel do
      
    end if
        
  end procedure ComputeReconstructionParabolic_CSL_Kernel


  module procedure ComputeLogDerivative_CSL_Kernel

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVS, &
      lV, uV
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    lV = 1
    where ( shape ( V ) > 1 )
      lV = oV + 1
    end where
    
    uV = 1
    where ( shape ( V ) > 1 )
      uV = shape ( V ) - oV
    end where
      
    iaS = 0
    iaS ( iD ) = +1
    
    if ( UseDevice ) then
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iV, jV, kV, iaVS )
      do kV = lV ( 3 ), uV ( 3 ) 
        do jV = lV ( 2 ), uV ( 2 )
          do iV = lV ( 1 ), uV ( 1 )

            iaVS = [ iV, jV, kV ] + iaS

            dLVdX ( iV, jV, kV ) &
              =  (    A_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) &
                   -  A_I ( iV, jV, kV )  ) &
                 /  V ( iV, jV, kV )

          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
      
    else
        
      !$OMP  parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iV, jV, kV, iaVS )
      do kV = lV ( 3 ), uV ( 3 ) 
        do jV = lV ( 2 ), uV ( 2 )
          do iV = lV ( 1 ), uV ( 1 )

            iaVS = [ iV, jV, kV ] + iaS

            dLVdX ( iV, jV, kV ) &
              =  (    A_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) &
                   -  A_I ( iV, jV, kV )  ) &
                 /  V ( iV, jV, kV )

          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end parallel do
      
    end if
      
  end procedure ComputeLogDerivative_CSL_Kernel

  
  module procedure ComputeIncrement_CSL_Kernel
    
    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVS, &
      lV, uV
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

!    dU  =  dU  +  dT * ( VJ_I * F_I  &
!                         -  cshift ( VJ_I * F_I, shift = +1, dim = iD ) ) &
!                       / ( VJ * dX )

    lV = 1
    where ( shape ( dU ) > 1 )
      lV = oV + 1
    end where
    
    uV = 1
    where ( shape ( dU ) > 1 )
      uV = shape ( dU ) - oV
    end where
      
    iaS = 0
    iaS ( iD ) = +1
    
    if ( UseDevice ) then
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iV, jV, kV, iaVS )
      do kV = lV ( 3 ), uV ( 3 ) 
        do jV = lV ( 2 ), uV ( 2 )
          do iV = lV ( 1 ), uV ( 1 )

            iaVS = [ iV, jV, kV ] + iaS

            dU ( iV, jV, kV ) &
              = dU ( iV, jV, kV )  &
                +  dT * (    A_I ( iV, jV, kV ) &
                             *  F_I ( iV, jV, kV ) &
                          -  A_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) &
                             *  F_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) ) &
                        /  V ( iV, jV, kV )

          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do

    else
        
      !$OMP  parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iV, jV, kV, iaVS )
      do kV = lV ( 3 ), uV ( 3 ) 
        do jV = lV ( 2 ), uV ( 2 )
          do iV = lV ( 1 ), uV ( 1 )

            iaVS = [ iV, jV, kV ] + iaS

            dU ( iV, jV, kV ) &
              = dU ( iV, jV, kV )  &
                +  dT * (    A_I ( iV, jV, kV ) &
                             *  F_I ( iV, jV, kV ) &
                          -  A_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) &
                             *  F_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) ) &
                        /  V ( iV, jV, kV )

          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP  end parallel do
    
    end if
    
  end procedure ComputeIncrement_CSL_Kernel
  
  
  module procedure RecordBoundaryFluence_CSL_Kernel

    integer ( KDI ) :: &
      iV, jV, kV
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    if ( UseDevice ) then 
      
      !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iV, jV, kV ) 
      do kV = 1, nB ( 3 )
        do jV = 1, nB ( 2 )
          do iV = 1, nB ( 1 )
            BF ( iV, jV, kV ) &
              =  BF ( iV, jV, kV ) &
                 +  Factor &
                    *   F ( oB ( 1 ) + iV, oB ( 2 ) + jV, oB ( 3 ) + kV )
          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do

    else 
    
      !$OMP  parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iV, jV, kV ) 
      do kV = 1, nB ( 3 )
        do jV = 1, nB ( 2 )
          do iV = 1, nB ( 1 )
            BF ( iV, jV, kV ) &
              =  BF ( iV, jV, kV ) &
                 +  Factor &
                    *   F ( oB ( 1 ) + iV, oB ( 2 ) + jV, oB ( 3 ) + kV )
          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP  end parallel do
    
    end if

  end procedure RecordBoundaryFluence_CSL_Kernel


end submodule IncrementDivergence_FV__Kernel
