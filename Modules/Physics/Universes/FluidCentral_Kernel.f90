#include "Preprocessor"

submodule ( FluidCentral_Template ) FluidCentral_Kernel

  use Basics
  
  implicit none
  
contains


  module procedure ComposePillarsPack_2_Kernel
    
    integer ( KDI ) :: &
      oO, &      !-- oOutgoing, oIncoming
      iC, kC, &  !-- iCell, etc.
      iS, &      !-- iSelected
      nVariables !
    integer ( KDI ), dimension ( 3 ) :: &
      lC, uC

    nVariables  =  size ( iaS )

    lC  =  nGL + 1
    uC  =  nGL + nCB

    !$OMP  parallel do collapse ( 2 ) &
    !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
    !$OMP& private ( iC, kC, iS, oO ) firstprivate ( nVariables )  
    do kC = lC ( 3 ), uC ( 3 )
      do iC = lC ( 1 ), uC ( 1 )
        if ( Crsn_2 ( iC, lC ( 2 ), kC ) < 2.0_KDR ) &
          cycle
        oO  =  oOC_2 ( iC, kC )
        Outgoing ( oO + 1 )  =  Crsn_2 ( iC, lC ( 2 ), kC )
        oO  =  oO + 1
        Outgoing ( oO + 1 : oO + nCB ( 2 ) ) &
          =  Vol ( iC, lC ( 2 ) : uC ( 2 ), kC )
        oO  =  oO + nCB ( 2 )
        do iS = 1, nVariables
          Outgoing ( oO + 1 : oO + nCB ( 2 ) ) &
            =  SV ( iC, lC ( 2 ) : uC ( 2 ), kC, iaS ( iS ) )
          oO  =  oO + nCB ( 2 )
        end do !-- iS
      end do !-- iC
    end do !-- kC
    !$OMP  end parallel do
  
  end procedure ComposePillarsPack_2_Kernel
  
  
  module procedure ComposePillarsPack_3_Kernel
  
    integer ( KDI ) :: &
      oO, &      !-- oOutgoing, oIncoming
      iC, jC, &  !-- iCell, etc.
      iS, &      !-- iSelected
      nVariables !
    integer ( KDI ), dimension ( 3 ) :: &
      lC, uC
      
    nVariables = size ( iaS )
    
    lC  =  nGL + 1
    uC  =  nGL + nCB
    
    !$OMP  parallel do collapse ( 2 ) &
    !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
    !$OMP& private ( iC, jC, iS, oO ) firstprivate ( nVariables )  
    do jC = lC ( 2 ), uC ( 2 )
      do iC = lC ( 1 ), uC ( 1 )
        if ( Crsn_3 ( iC, jC, lC ( 3 ) ) < 2.0_KDR ) &
          cycle
        oO  =  oOC_3 ( iC, jC )
        Outgoing ( oO + 1 ) = Crsn_3 ( iC, jC, lC ( 3 ) )
        oO = oO + 1
        Outgoing ( oO + 1 : oO + nCB ( 3 ) ) &
          =  Vol ( iC, jC, lC ( 3 ) : uC ( 3 ) )
        oO = oO + nCB ( 3 )
        do iS = 1, nVariables
          Outgoing ( oO + 1 : oO + nCB ( 3 ) ) &
            =  SV ( iC, jC, lC ( 3 ) : uC ( 3 ), iaS ( iS ) )
          oO = oO + nCB ( 3 )
        end do !-- iS
      end do !-- iC
    end do !-- jC
    !$OMP  end parallel do
  
  end procedure ComposePillarsPack_3_Kernel
  
  
  module procedure ComposePillarsUnpack_2_Kernel
    
    integer ( KDI ) :: &
      oI, &          !-- oOutgoing, oIncoming
      oP, &          !-- oPillar
      oC, &          !-- oCell
      iC, jC, kC, &  !-- iCell, etc.
      iS, &          !-- iSelected
      iG, &          !-- iGroup (of pillars/processes)
      iB, &          !-- iBrick
      iR, &          !-- iRank
      iPS, &         !-- iPillarSegment
      iP, &          !-- iPillar
      nVariables

    oP = 0
    oI = 0
    nVariables = size ( CoarsenPillar_2, dim = 2 )
    do iG = 1, nGroups 
      oC = 0
      do iB = 1, nBricks ( 2 )
        iR = ( iG - 1 ) * nBricks ( 2 )  +  iB  -  1
        do iPS = 1, nSegmentsFrom_2 ( iR )
          iP = oP + iPS
          nCoarsen_2 ( iP ) &
            = min ( int ( Incoming ( oI + 1 ) + 0.5_KDR ), nCells ( 2 ) )
          oI = oI + 1
          do iS = 1, nVariables
            CoarsenPillar_2 ( oC + 1 : oC + nCB ( 2 ), iS, iP )  &
              =  Incoming ( oI + 1 : oI + nCB ( 2 ) )
            oI = oI + nCB ( 2 )
          end do !-- iS
        end do !-- iPS
        oC = oC + nCB ( 2 )
      end do !-- iB
      oP = oP + nSegmentsFrom_2 ( iR )
    end do !-- iG

  end procedure ComposePillarsUnpack_2_Kernel
  
  
  module procedure ComposePillarsUnpack_3_Kernel
    
    integer ( KDI ) :: &
      oI, &          !-- oOutgoing, oIncoming
      oP, &          !-- oPillar
      oC, &          !-- oCell
      iC, jC, kC, &  !-- iCell, etc.
      iS, &          !-- iSelected
      iG, &          !-- iGroup (of pillars/processes)
      iB, &          !-- iBrick
      iR, &          !-- iRank
      iPS, &         !-- iPillarSegment
      iP, &          !-- iPillar
      nVariables

    oP = 0
    oI = 0
    nVariables = size ( CoarsenPillar_3, dim = 2 )
    do iG = 1, nGroups
      oC = 0
      do iB = 1, nBricks ( 3 )
        iR = ( iG - 1 ) * nBricks ( 3 )  +  iB  -  1
        do iPS = 1, nSegmentsFrom_3 ( iR )
          iP = oP + iPS
          nCoarsen_3 ( iP ) &
            = min ( int ( Incoming ( oI + 1 ) + 0.5_KDR ), nCells ( 3 ) )
          oI = oI + 1
          do iS = 1, nVariables
            CoarsenPillar_3 ( oC + 1 : oC + nCB ( 3 ), iS, iP )  &
              =  Incoming ( oI + 1 : oI + nCB ( 3 ) )
            oI = oI + nCB ( 3 )
          end do !-- iS
        end do !-- iPS
        oC = oC + nCB ( 3 )
      end do !-- iB
      oP = oP + nSegmentsFrom_3 ( iR )
    end do !-- iG

  end procedure ComposePillarsUnpack_3_Kernel
  

  module procedure CoarsenPillarsKernel

    integer ( KDI ) :: &
      oC, &  !-- oCell
      iP, &  !-- iPillar
      iA, &  !-- iAverage
      iV, &  !-- iVariable
      nPillars, &
      nVariables, &
      nAverages
    real ( KDR ) :: &
      Volume_A, &
      Density_A

      nPillars  =  size ( CP, dim = 3 ) 
    nVariables  =  size ( CP, dim = 2 )

    do iP  =  1, nPillars
      oC  =  0
      nAverages  =  size ( CP, dim = 1 )  /  nCoarsen ( iP )
      associate ( Volume => CP ( :, 1, iP ) )
      do iA  =  1, nAverages
        Volume_A  =  sum ( Volume ( oC + 1 : oC + nCoarsen ( iP ) ) )
        do iV = 2, nVariables
          associate ( Density => CP ( :, iV, iP ) )
          Density_A = sum ( Volume ( oC + 1 : oC + nCoarsen ( iP ) ) &
                            *  Density ( oC + 1 : oC + nCoarsen ( iP ) ) )
          Density ( oC + 1 : oC + nCoarsen ( iP ) ) &
            =  Density_A / Volume_A
          end associate !-- Density
        end do !-- iV
        oC  =  oC + nCoarsen ( iP )
      end do !-- iA
      end associate !-- Volume
    end do !-- iP

  end procedure CoarsenPillarsKernel


  module procedure DecomposePillarsPack_2_Kernel
    
    integer ( KDI ) :: &
      oO, &           !-- oOutgoing
      oP, &          !-- oPillar
      oC, &          !-- oCell
      iS, &          !-- iSelected
      iG, &          !-- iGroup (of pillars/processes)
      iB, &          !-- iBrick
      iR, &          !-- iRank
      iPS, &         !-- iPillarSegment
      iP, &          !-- iPillar
      nVariables

    nVariables  =  size ( CoarsenPillar_2, dim = 2 )

    oP  =  0
    oO  =  0
    do iG  =  1, nGroups
      oC  =  0
      do iB  =  1, nBricks ( 2 )
        iR  =  ( iG - 1 ) * nBricks ( 2 )  +  iB  -  1
        do iPS = 1, nSegmentsFrom_2 ( iR )
          iP  =  oP + iPS
          do iS = 2, nVariables
            Outgoing ( oO + 1 : oO + nCB ( 2 ) )  &
              =  CoarsenPillar_2 ( oC + 1 : oC + nCB ( 2 ), iS, iP )
            oO  =  oO + nCB ( 2 )
          end do !-- iS
        end do !-- iPS
        oC  =  oC + nCB ( 2 )
      end do !-- iB
      oP  =  oP + nSegmentsFrom_2 ( iR )
    end do !-- iG
      
  end procedure DecomposePillarsPack_2_Kernel


  module procedure DecomposePillarsPack_3_Kernel
    
    integer ( KDI ) :: &
      oO, &           !-- oOutgoing
      oP, &          !-- oPillar
      oC, &          !-- oCell
      iS, &          !-- iSelected
      iG, &          !-- iGroup (of pillars/processes)
      iB, &          !-- iBrick
      iR, &          !-- iRank
      iPS, &         !-- iPillarSegment
      iP, &          !-- iPillar
      nVariables

    nVariables  =  size ( CoarsenPillar_3, dim = 2 )

    oP  =  0
    oO  =  0
    do iG  =  1, nGroups
      oC  =  0
      do iB  =  1, nBricks ( 3 )
        iR  =  ( iG - 1 ) * nBricks ( 3 )  +  iB  -  1
        do iPS = 1, nSegmentsFrom_3 ( iR )
          iP  =  oP + iPS
          do iS = 2, nVariables
            Outgoing ( oO + 1 : oO + nCB ( 3 ) )  &
              =  CoarsenPillar_3 ( oC + 1 : oC + nCB ( 3 ), iS, iP )
            oO  =  oO + nCB ( 3 )
          end do !-- iS
        end do !-- iPS
        oC  =  oC + nCB ( 3 )
      end do !-- iB
      oP  =  oP + nSegmentsFrom_3 ( iR )
    end do !-- iG
      
  end procedure DecomposePillarsPack_3_Kernel


  module procedure DecomposePillarsUnpack_2_Kernel
    
    integer ( KDI ) :: &
      oI, &      !-- oIncoming
      iC, kC, &  !-- iCell, etc.
      iS, &          !-- iSelected
      nVariables
    integer ( KDI ), dimension ( 3 ) :: &
      lC, uC

    nVariables  =  size ( iaS )

    lC  =  nGL + 1
    uC  =  nGL + nCB

    oI  =  0
    do kC  =  lC ( 3 ), uC ( 3 )
      do iC  =  lC ( 1 ), uC ( 1 )
          if ( Crsn_2 ( iC, lC ( 2 ), kC ) < 2.0_KDR ) &
            cycle
          do iS = 1, nVariables
            SV ( iC, lC ( 2 ) : uC ( 2 ), kC, iaS ( iS ) )  &
              =  Incoming ( oI + 1 : oI + nCB ( 2 ) )
            oI = oI + nCB ( 2 )
            if ( ( iS == iMomentum_2 .or. iS == iMomentum_3 ) &
                 .and. R ( iC, lC ( 2 ), kC )  <  RadiusPolarMomentum ) &
              SV ( iC, lC ( 2 ) : uC ( 2 ), kC, iaS ( iS ) ) = 0.0_KDR
          end do !-- iS
      end do !-- iC
    end do !-- kC

  end procedure DecomposePillarsUnpack_2_Kernel


  module procedure DecomposePillarsUnpack_3_Kernel
    
    integer ( KDI ) :: &
      oI, &      !-- oIncoming
      iC, jC, &  !-- iCell, etc.
      iS, &          !-- iSelected
      nVariables
    integer ( KDI ), dimension ( 3 ) :: &
      lC, uC

    nVariables  =  size ( iaS )

    lC  =  nGL + 1
    uC  =  nGL + nCB

    oI  =  0
    do jC  =  lC ( 2 ), uC ( 2 )
      do iC  =  lC ( 1 ), uC ( 1 )
          if ( Crsn_3 ( iC, jC, lC ( 3 ) ) < 2.0_KDR ) &
            cycle
          do iS = 1, nVariables
            SV ( iC, jC, lC ( 3 ) : uC ( 3 ), iaS ( iS ) )  &
              =  Incoming ( oI + 1 : oI + nCB ( 3 ) )
            oI = oI + nCB ( 3 )
            if ( ( iS == iMomentum_2 .or. iS == iMomentum_3 ) &
                 .and. R ( iC, jC, lC ( 3 ) )  <  RadiusPolarMomentum ) &
              SV ( iC, jC, lC ( 3 ) : uC ( 3 ), iaS ( iS ) ) = 0.0_KDR
          end do !-- iS
      end do !-- iC
    end do !-- jC

  end procedure DecomposePillarsUnpack_3_Kernel


end submodule FluidCentral_Kernel
