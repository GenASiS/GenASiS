#include "Preprocessor"

submodule ( FluidFeatures_P__Form ) FluidFeatures_P__Kernel

  use Basics

  implicit none

contains


  module procedure DetectShocks_CSL_Kernel

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
                =  S ( iV, jV, kV )  +  1.0_KDR
              S ( iaV_i ( 1 ), iaV_i ( 2 ), iaV_i ( 3 ) )  &
                =  S ( iaV_i ( 1 ), iaV_i ( 2 ), iaV_i ( 3 ) )  +  1.0_KDR

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
                =  S ( iV, jV, kV )  +  1.0_KDR
              S ( iaV_i ( 1 ), iaV_i ( 2 ), iaV_i ( 3 ) )  &
                =  S ( iaV_i ( 1 ), iaV_i ( 2 ), iaV_i ( 3 ) )  +  1.0_KDR

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
        
  end procedure DetectShocks_CSL_Kernel

  
  module procedure ClearBoundary_CSL_Kernel 

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
        !$OMP& private ( iV, jV, kV )
        !$OMP& firstprivate ( lV, uV )
        do kV = lV ( 3 ), uV ( 3 ) 
          do jV = lV ( 2 ), uV ( 2 )
            do iV = lV ( 1 ), uV ( 1 )

                    S ( iV, jV, kV )  =  0.0_KDR
               S_I_iD ( iV, jV, kV )  =  0.0_KDR
              DF_I_iD ( iV, jV, kV )  =  0.0_KDR
              DF_I_jD ( iV, jV, kV )  =  0.0_KDR
              DF_I_kD ( iV, jV, kV )  =  0.0_KDR

            end do !-- iV
          end do !-- jV
        end do !-- kV
        !$OMP  end OMP_TARGET_DIRECTIVE
      else
        !$OMP  parallel do collapse ( 3 ) &
        !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
        !$OMP& private ( iV, jV, kV ) &
        !$OMP& firstprivate ( lV, uV )
        do kV = lV ( 3 ), uV ( 3 ) 
          do jV = lV ( 2 ), uV ( 2 )
            do iV = lV ( 1 ), uV ( 1 )

                    S ( iV, jV, kV )  =  0.0_KDR
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
        !$OMP& private ( iV, jV, kV )
        !$OMP& firstprivate ( lV, uV )
        do kV = lV ( 3 ), uV ( 3 ) 
          do jV = lV ( 2 ), uV ( 2 )
            do iV = lV ( 1 ), uV ( 1 )

                    S ( iV, jV, kV )  =  0.0_KDR
               S_I_iD ( iV, jV, kV )  =  0.0_KDR
              DF_I_iD ( iV, jV, kV )  =  0.0_KDR
              DF_I_jD ( iV, jV, kV )  =  0.0_KDR
              DF_I_kD ( iV, jV, kV )  =  0.0_KDR

            end do !-- iV
          end do !-- jV
        end do !-- kV
        !$OMP  end OMP_TARGET_DIRECTIVE
      else
        !$OMP  parallel do collapse ( 3 ) &
        !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
        !$OMP& private ( iV, jV, kV ) &
        !$OMP& firstprivate ( lV, uV )
        do kV = lV ( 3 ), uV ( 3 ) 
          do jV = lV ( 2 ), uV ( 2 )
            do iV = lV ( 1 ), uV ( 1 )

                    S ( iV, jV, kV )  =  0.0_KDR
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
       
  end subroutine ClearBoundary_CSL_Kernel


end submodule FluidFeatures_P__Kernel
