#include "Preprocessor"

submodule ( Features_F_P__Form ) Features_F_P__Kernel

  use Basics

  implicit none

contains


  module procedure DetectShocksKernel

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS_i, iaS_j, iaS_k, &
      iaV_i, iaV_j, iaV_ij, iaV_k, iaV_ik, &
      lV, uV
    real ( KDR ) :: &
      dP, &
      P_Min, &
      dLnP, &
      dV_iD, &
      SqrtTiny
    logical ( KDL ) :: &
      UseDevice
         
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
    UseDevice = UseDeviceOption

    lV = 1
    where ( shape ( S ) > 1 )
      lV = oV
    end where

    uV = 1
    where ( shape ( S ) > 1 )
      uV = shape ( S ) - oV + 1
    end where
    
    iaS_i = 0
    iaS_i ( iD ) = -1
      
    iaS_j = 0
    if ( size ( S, dim = jD ) > 1 ) &
      iaS_j ( jD ) = +1

    iaS_k = 0
    if ( size ( S, dim = kD ) > 1 ) &
      iaS_k ( kD ) = +1

    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )
    
    if ( UseDevice ) then

      !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP& private ( iV, jV, kV, iaV_i, iaV_j, iaV_ij, iaV_k, iaV_ik ) &
      !$OMP& private ( dP, P_Min, dLnP, dV_iD ) &
      !$OMP& firstprivate ( SqrtTiny, lV, uV )
      do kV = lV ( 3 ), uV ( 3 ) 
        do jV = lV ( 2 ), uV ( 2 )
          do iV = lV ( 1 ), uV ( 1 )

            iaV_i  = [ iV, jV, kV ] + iaS_i
            iaV_j  = [ iV, jV, kV ] + iaS_j
            iaV_ij = [ iV, jV, kV ] + iaS_i + iaS_j
            iaV_k  = [ iV, jV, kV ] + iaS_k
            iaV_ik = [ iV, jV, kV ] + iaS_i + iaS_k

            dP  =  abs ( P ( iV, jV, kV )  &
                         -  P ( iaV_i ( 1 ), iaV_i ( 2 ), iaV_i ( 3 ) ) )
            P_Min  =  max ( min ( P ( iV, jV, kV ), &
                                  P ( iaV_i ( 1 ), iaV_i ( 2 ), iaV_i ( 3 ) )), &
                            SqrtTiny )
            dLnP  =  dP / P_Min

            dV_iD  =  V_iD ( iV, jV, kV ) &
                      -  V_iD ( iaV_i ( 1 ), iaV_i ( 2 ), iaV_i ( 3 ) )
   
            if ( dLnP > ST .and. dV_iD <= 0.0_KDR ) then

              S_I_iD ( iV, jV, kV )  &
                =  1.0_KDR
              S ( iV, jV, kV )  &
                =  1.0_KDR
              S ( iaV_i ( 1 ), iaV_i ( 2 ), iaV_i ( 3 ) )  &
                =  1.0_KDR

              !-- Use diffuse flux in transverse directions, on both sides of 
              !   the shock

              DF_I_jD ( iV, jV, kV ) &
                =  1.0_KDR
              DF_I_jD ( iaV_i ( 1 ), iaV_i ( 2 ), iaV_i ( 3 ) ) &
                =  1.0_KDR
              DF_I_jD ( iaV_j ( 1 ), iaV_j ( 2 ), iaV_j ( 3 ) ) &
                =  1.0_KDR
              DF_I_jD ( iaV_ij ( 1 ), iaV_ij ( 2 ), iaV_ij ( 3 ) ) &
                =  1.0_KDR

              DF_I_kD ( iV, jV, kV ) &
                =  1.0_KDR
              DF_I_kD ( iaV_i ( 1 ), iaV_i ( 2 ), iaV_i ( 3 ) ) &
                =  1.0_KDR
              DF_I_kD ( iaV_k ( 1 ), iaV_k ( 2 ), iaV_k ( 3 ) ) &
                =  1.0_KDR
              DF_I_kD ( iaV_ik ( 1 ), iaV_ik ( 2 ), iaV_ik ( 3 ) ) &
                =  1.0_KDR

            end if

          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
    
    else
    
      !$OMP  parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP& private ( iV, jV, kV, iaV_i, iaV_j, iaV_ij, iaV_k, iaV_ik ) &
      !$OMP& private ( dP, P_Min, dLnP, dV_iD ) &
      !$OMP& firstprivate ( SqrtTiny, lV, uV )
      do kV = lV ( 3 ), uV ( 3 ) 
        do jV = lV ( 2 ), uV ( 2 )
          do iV = lV ( 1 ), uV ( 1 )

            iaV_i  = [ iV, jV, kV ] + iaS_i
            iaV_j  = [ iV, jV, kV ] + iaS_j
            iaV_ij = [ iV, jV, kV ] + iaS_i + iaS_j
            iaV_k  = [ iV, jV, kV ] + iaS_k
            iaV_ik = [ iV, jV, kV ] + iaS_i + iaS_k

            dP  =  abs ( P ( iV, jV, kV )  &
                         -  P ( iaV_i ( 1 ), iaV_i ( 2 ), iaV_i ( 3 ) ) )
            P_Min  =  max ( min ( P ( iV, jV, kV ), &
                                  P ( iaV_i ( 1 ), iaV_i ( 2 ), iaV_i ( 3 ) )), &
                            SqrtTiny )
            dLnP  =  dP / P_Min

            dV_iD  =  V_iD ( iV, jV, kV ) &
                      -  V_iD ( iaV_i ( 1 ), iaV_i ( 2 ), iaV_i ( 3 ) )
   
            if ( dLnP > ST .and. dV_iD <= 0.0_KDR ) then

              S_I_iD ( iV, jV, kV )  &
                =  1.0_KDR
              S ( iV, jV, kV )  &
                =  1.0_KDR
              S ( iaV_i ( 1 ), iaV_i ( 2 ), iaV_i ( 3 ) )  &
                =  1.0_KDR

              !-- Use diffuse flux in transverse directions, on both sides of 
              !   the shock

              DF_I_jD ( iV, jV, kV ) &
                =  1.0_KDR
              DF_I_jD ( iaV_i ( 1 ), iaV_i ( 2 ), iaV_i ( 3 ) ) &
                =  1.0_KDR
              DF_I_jD ( iaV_j ( 1 ), iaV_j ( 2 ), iaV_j ( 3 ) ) &
                =  1.0_KDR
              DF_I_jD ( iaV_ij ( 1 ), iaV_ij ( 2 ), iaV_ij ( 3 ) ) &
                =  1.0_KDR

              DF_I_kD ( iV, jV, kV ) &
                =  1.0_KDR
              DF_I_kD ( iaV_i ( 1 ), iaV_i ( 2 ), iaV_i ( 3 ) ) &
                =  1.0_KDR
              DF_I_kD ( iaV_k ( 1 ), iaV_k ( 2 ), iaV_k ( 3 ) ) &
                =  1.0_KDR
              DF_I_kD ( iaV_ik ( 1 ), iaV_ik ( 2 ), iaV_ik ( 3 ) ) &
                =  1.0_KDR

            end if

          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP  end parallel do
    
    end if
        
  end procedure DetectShocksKernel

  
  module procedure DetectPhaseTransitionKernel

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS_im, iaS_ip, iaS_jp, iaS_kp, &  !-- iaShift
      iaV_im, iaV_ip, iaV_jp, iaV_kp, iaV_im_jp, iaV_im_kp, &  !--iaValue
      lV, uV
    real ( KDR ) :: &
      dGamma, &
      GammaMin, &
      dLnGamma, &
      SqrtTiny
    logical ( KDL ) :: &
      UseDevice
         
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
    UseDevice = UseDeviceOption

    lV = 1
    where ( shape ( PT ) > 1 )
      lV = oV
    end where

    uV = 1
    where ( shape ( PT ) > 1 )
      uV = shape ( PT ) - oV + 1
    end where
    
    iaS_im = 0
    iaS_im ( iD ) = -1
      
    iaS_ip = 0
    iaS_ip ( iD ) = +1
      
    iaS_jp = 0
    if ( size ( PT, dim = jD ) > 1 ) &
      iaS_jp ( jD ) = +1

    iaS_kp = 0
    if ( size ( PT, dim = kD ) > 1 ) &
      iaS_kp ( kD ) = +1

    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )
    
    if ( UseDevice ) then

      !$OMP OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP private ( iV, jV, kV, iaV_im, iaV_ip, iaV_jp, iaV_kp ) &
      !$OMP private ( iaV_im_jp, iaV_im_kp ) &
      !$OMP private ( dGamma, GammaMin, dLnGamma ) &
      !$OMP firstprivate ( SqrtTiny, lV, uV )
      do kV = lV ( 3 ), uV ( 3 ) 
        do jV = lV ( 2 ), uV ( 2 )
          do iV = lV ( 1 ), uV ( 1 )

            iaV_im    = [ iV, jV, kV ] + iaS_im
            iaV_ip    = [ iV, jV, kV ] + iaS_ip
            iaV_jp    = [ iV, jV, kV ] + iaS_jp
            iaV_kp    = [ iV, jV, kV ] + iaS_kp
            iaV_im_jp = [ iV, jV, kV ] + iaS_im + iaS_jp
            iaV_im_kp = [ iV, jV, kV ] + iaS_im + iaS_kp

            dGamma  =  abs ( Gamma ( iV, jV, kV )  &
                         -  Gamma ( iaV_im ( 1 ), iaV_im ( 2 ), iaV_im ( 3 ) ) )
            GammaMin  &
              =  max ( &
                   min ( Gamma ( iV, jV, kV ), &
                         Gamma ( iaV_im ( 1 ), iaV_im ( 2 ), iaV_im ( 3 ) ) ), &
                   SqrtTiny )
            dLnGamma  =  dGamma / GammaMin

            !-- Fractional difference across inner face
            if ( dLnGamma > PTT ) then

              PT ( iV, jV, kV )  &
                =  1.0_KDR
              PT ( iaV_im ( 1 ), iaV_im ( 2 ), iaV_im ( 3 ) )  &
                =  1.0_KDR

              !-- Use diffuse flux in longitudinal direction

              DF_I_iD ( iV, jV, kV ) &
                =  1.0_KDR
              DF_I_iD ( iaV_im ( 1 ), iaV_im ( 2 ), iaV_im ( 3 ) )  &
                =  1.0_KDR
              DF_I_iD ( iaV_ip ( 1 ), iaV_ip ( 2 ), iaV_ip ( 3 ) ) &
                =  1.0_KDR

              !-- Use diffuse flux in transverse directions, on both sides of 
              !   the shock

              DF_I_jD ( iV, jV, kV ) &
                =  1.0_KDR
              DF_I_jD ( iaV_im ( 1 ), iaV_im ( 2 ), iaV_im ( 3 ) ) &
                =  1.0_KDR
              DF_I_jD ( iaV_jp ( 1 ), iaV_jp ( 2 ), iaV_jp ( 3 ) ) &
                =  1.0_KDR
              DF_I_jD ( iaV_im_jp ( 1 ), iaV_im_jp ( 2 ), iaV_im_jp ( 3 ) ) &
                =  1.0_KDR

              DF_I_kD ( iV, jV, kV ) &
                =  1.0_KDR
              DF_I_kD ( iaV_im ( 1 ), iaV_im ( 2 ), iaV_im ( 3 ) ) &
                =  1.0_KDR
              DF_I_kD ( iaV_kp ( 1 ), iaV_kp ( 2 ), iaV_kp ( 3 ) ) &
                =  1.0_KDR
              DF_I_kD ( iaV_im_kp ( 1 ), iaV_im_kp ( 2 ), iaV_im_kp ( 3 ) ) &
                =  1.0_KDR

            end if !-- Fractional difference

            !-- Absolute value
            if ( Gamma ( iV, jV, kV ) < 1.0_KDR ) then

              PT ( iV, jV, kV )  &
                =  1.0_KDR

              !-- Use diffuse flux on all faces

              DF_I_iD ( iV, jV, kV ) &
                =  1.0_KDR
              DF_I_iD ( iaV_ip ( 1 ), iaV_ip ( 2 ), iaV_ip ( 3 ) ) &
                =  1.0_KDR

              DF_I_jD ( iV, jV, kV ) &
                =  1.0_KDR
              DF_I_jD ( iaV_jp ( 1 ), iaV_jp ( 2 ), iaV_jp ( 3 ) ) &
                =  1.0_KDR

              DF_I_kD ( iV, jV, kV ) &
                =  1.0_KDR
              DF_I_kD ( iaV_kp ( 1 ), iaV_kp ( 2 ), iaV_kp ( 3 ) ) &
                =  1.0_KDR

            end if !-- Absolute value

          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
    
    else
    
      !$OMP parallel do collapse ( 3 ) &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP private ( iV, jV, kV, iaV_im, iaV_ip, iaV_jp, iaV_kp ) &
      !$OMP private ( iaV_im_jp, iaV_im_kp ) &
      !$OMP private ( dGamma, GammaMin, dLnGamma ) &
      !$OMP firstprivate ( SqrtTiny, lV, uV )
      do kV = lV ( 3 ), uV ( 3 ) 
        do jV = lV ( 2 ), uV ( 2 )
          do iV = lV ( 1 ), uV ( 1 )

            iaV_im    = [ iV, jV, kV ] + iaS_im
            iaV_ip    = [ iV, jV, kV ] + iaS_ip
            iaV_jp    = [ iV, jV, kV ] + iaS_jp
            iaV_kp    = [ iV, jV, kV ] + iaS_kp
            iaV_im_jp = [ iV, jV, kV ] + iaS_im + iaS_jp
            iaV_im_kp = [ iV, jV, kV ] + iaS_im + iaS_kp

            dGamma  =  abs ( Gamma ( iV, jV, kV )  &
                         -  Gamma ( iaV_im ( 1 ), iaV_im ( 2 ), iaV_im ( 3 ) ) )
            GammaMin  &
              =  max ( &
                   min ( Gamma ( iV, jV, kV ), &
                         Gamma ( iaV_im ( 1 ), iaV_im ( 2 ), iaV_im ( 3 ) ) ), &
                   SqrtTiny )
            dLnGamma  =  dGamma / GammaMin

            !-- Fractional difference across inner face
            if ( dLnGamma > PTT ) then

              PT ( iV, jV, kV )  &
                =  1.0_KDR
              PT ( iaV_im ( 1 ), iaV_im ( 2 ), iaV_im ( 3 ) )  &
                =  1.0_KDR

              !-- Use diffuse flux in longitudinal direction

              DF_I_iD ( iV, jV, kV ) &
                =  1.0_KDR
              DF_I_iD ( iaV_im ( 1 ), iaV_im ( 2 ), iaV_im ( 3 ) )  &
                =  1.0_KDR
              DF_I_iD ( iaV_ip ( 1 ), iaV_ip ( 2 ), iaV_ip ( 3 ) ) &
                =  1.0_KDR

              !-- Use diffuse flux in transverse directions, on both sides of 
              !   the shock

              DF_I_jD ( iV, jV, kV ) &
                =  1.0_KDR
              DF_I_jD ( iaV_im ( 1 ), iaV_im ( 2 ), iaV_im ( 3 ) ) &
                =  1.0_KDR
              DF_I_jD ( iaV_jp ( 1 ), iaV_jp ( 2 ), iaV_jp ( 3 ) ) &
                =  1.0_KDR
              DF_I_jD ( iaV_im_jp ( 1 ), iaV_im_jp ( 2 ), iaV_im_jp ( 3 ) ) &
                =  1.0_KDR

              DF_I_kD ( iV, jV, kV ) &
                =  1.0_KDR
              DF_I_kD ( iaV_im ( 1 ), iaV_im ( 2 ), iaV_im ( 3 ) ) &
                =  1.0_KDR
              DF_I_kD ( iaV_kp ( 1 ), iaV_kp ( 2 ), iaV_kp ( 3 ) ) &
                =  1.0_KDR
              DF_I_kD ( iaV_im_kp ( 1 ), iaV_im_kp ( 2 ), iaV_im_kp ( 3 ) ) &
                =  1.0_KDR

            end if !-- Fractional difference

            !-- Absolute value
            if ( Gamma ( iV, jV, kV ) < 1.0_KDR ) then

              PT ( iV, jV, kV )  &
                =  1.0_KDR

              !-- Use diffuse flux on all faces

              DF_I_iD ( iV, jV, kV ) &
                =  1.0_KDR
              DF_I_iD ( iaV_ip ( 1 ), iaV_ip ( 2 ), iaV_ip ( 3 ) ) &
                =  1.0_KDR

              DF_I_jD ( iV, jV, kV ) &
                =  1.0_KDR
              DF_I_jD ( iaV_jp ( 1 ), iaV_jp ( 2 ), iaV_jp ( 3 ) ) &
                =  1.0_KDR

              DF_I_kD ( iV, jV, kV ) &
                =  1.0_KDR
              DF_I_kD ( iaV_kp ( 1 ), iaV_kp ( 2 ), iaV_kp ( 3 ) ) &
                =  1.0_KDR

            end if !-- Absolute value

          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP  end parallel do
    
    end if
        
  end procedure DetectPhaseTransitionKernel

  
  module procedure ClearBoundaryKernel 

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      lV, &
      uV
    logical ( KDL ) :: &
      UseDevice
         
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
    UseDevice = UseDeviceOption
    
    if ( InnerBoundary ) then

      lV = 1
      where ( shape ( S ) > 1 )
        lV = oV
      end where

      uV = 1
      where ( shape ( S ) > 1 )
        uV = shape ( S ) - oV + 1 
      end where
      uV ( iD ) = lV ( iD ) + 1
      
      if ( UseDevice ) then
        !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
        !$OMP& schedule ( OMP_SCHEDULE_TARGET ) &
        !$OMP& private ( iV, jV, kV ) &
        !$OMP& firstprivate ( lV, uV )
        do kV = lV ( 3 ), uV ( 3 ) 
          do jV = lV ( 2 ), uV ( 2 )
            do iV = lV ( 1 ), uV ( 1 )

                    S ( iV, jV, kV )  =  0.0_KDR
                   PT ( iV, jV, kV )  =  0.0_KDR
               S_I_iD ( iV, jV, kV )  =  0.0_KDR
              DF_I_iD ( iV, jV, kV )  =  0.0_KDR
              DF_I_jD ( iV, jV, kV )  =  0.0_KDR
              DF_I_kD ( iV, jV, kV )  =  0.0_KDR

            end do !-- iV
          end do !-- jV
        end do !-- kV
        !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      else
        !$OMP  parallel do collapse ( 3 ) &
        !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
        !$OMP& private ( iV, jV, kV ) &
        !$OMP& firstprivate ( lV, uV )
        do kV = lV ( 3 ), uV ( 3 ) 
          do jV = lV ( 2 ), uV ( 2 )
            do iV = lV ( 1 ), uV ( 1 )

                    S ( iV, jV, kV )  =  0.0_KDR
                   PT ( iV, jV, kV )  =  0.0_KDR
               S_I_iD ( iV, jV, kV )  =  0.0_KDR
              DF_I_iD ( iV, jV, kV )  =  0.0_KDR
              DF_I_jD ( iV, jV, kV )  =  0.0_KDR
              DF_I_kD ( iV, jV, kV )  =  0.0_KDR

            end do !-- iV
          end do !-- jV
        end do !-- kV
        !$OMP  end parallel do
      end if

    end if !-- InnerBoundary
       
    if ( OuterBoundary ) then

      lV = 1
      where ( shape ( S ) > 1 )
        lV = oV - 1
      end where
      lV ( iD ) = size ( S, dim = iD ) - oV

      uV = 1
      where ( shape ( S ) > 1 )
        uV = shape ( S ) - oV + 1 
      end where
      
      if ( UseDevice ) then 
        !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
        !$OMP& schedule ( OMP_SCHEDULE_TARGET ) &
        !$OMP& private ( iV, jV, kV ) &
        !$OMP& firstprivate ( lV, uV )
        do kV = lV ( 3 ), uV ( 3 ) 
          do jV = lV ( 2 ), uV ( 2 )
            do iV = lV ( 1 ), uV ( 1 )

                    S ( iV, jV, kV )  =  0.0_KDR
                   PT ( iV, jV, kV )  =  0.0_KDR
               S_I_iD ( iV, jV, kV )  =  0.0_KDR
              DF_I_iD ( iV, jV, kV )  =  0.0_KDR
              DF_I_jD ( iV, jV, kV )  =  0.0_KDR
              DF_I_kD ( iV, jV, kV )  =  0.0_KDR

            end do !-- iV
          end do !-- jV
        end do !-- kV
        !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      else
        !$OMP  parallel do collapse ( 3 ) &
        !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
        !$OMP& private ( iV, jV, kV ) &
        !$OMP& firstprivate ( lV, uV )
        do kV = lV ( 3 ), uV ( 3 ) 
          do jV = lV ( 2 ), uV ( 2 )
            do iV = lV ( 1 ), uV ( 1 )

                    S ( iV, jV, kV )  =  0.0_KDR
                   PT ( iV, jV, kV )  =  0.0_KDR
               S_I_iD ( iV, jV, kV )  =  0.0_KDR
              DF_I_iD ( iV, jV, kV )  =  0.0_KDR
              DF_I_jD ( iV, jV, kV )  =  0.0_KDR
              DF_I_kD ( iV, jV, kV )  =  0.0_KDR

            end do !-- iV
          end do !-- jV
        end do !-- kV
        !$OMP  end parallel do
      end if

    end if !-- OuterBoundary
       
  end procedure ClearBoundaryKernel


end submodule Features_F_P__Kernel
