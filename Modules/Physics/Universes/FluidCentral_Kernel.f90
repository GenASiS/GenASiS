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
      
    nVariables = size ( iaS )
  
    oO = 0
    
    do kC = 1, nCB ( 3 )
      do iC = 1, nCB ( 1 )
        if ( Crsn_2 ( iC, 1, kC ) < 2.0_KDR ) &
          cycle
        Outgoing ( oO + 1 ) = Crsn_2 ( iC, 1, kC )
        oO = oO + 1
        Outgoing ( oO + 1 : oO + nCB ( 2 ) ) &
          =  Vol ( iC, 1 : nCB ( 2 ), kC )
        oO = oO + nCB ( 2 )
        do iS = 1, nVariables
          Outgoing ( oO + 1 : oO + nCB ( 2 ) ) &
            =  SV ( iC, 1 : nCB ( 2 ), kC, iaS ( iS ) )
          oO = oO + nCB ( 2 )
        end do !-- iS
      end do !-- iC
    end do !-- kC
  
  end procedure ComposePillarsPack_2_Kernel
  
  
  module procedure ComposePillarsPack_3_Kernel
  
    integer ( KDI ) :: &
      oO, &      !-- oOutgoing, oIncoming
      iC, jC, &  !-- iCell, etc.
      iS, &      !-- iSelected
      nVariables !
      
    nVariables = size ( iaS )
    
    
    oO = 0
    do jC = 1, nCB ( 2 )
      do iC = 1, nCB ( 1 )
        if ( Crsn_3 ( iC, jC, 1 ) < 2.0_KDR ) &
          cycle
        Outgoing ( oO + 1 ) = Crsn_3 ( iC, jC, 1 )
        oO = oO + 1
        Outgoing ( oO + 1 : oO + nCB ( 3 ) ) &
          =  Vol ( iC, jC, 1 : nCB ( 3 ) )
        oO = oO + nCB ( 3 )
        do iS = 1, nVariables
          Outgoing ( oO + 1 : oO + nCB ( 3 ) ) &
            =  SV ( iC, jC, 1 : nCB ( 3 ), iaS ( iS ) )
          oO = oO + nCB ( 3 )
        end do !-- iS
      end do !-- iC
    end do !-- jC
  
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
      iP             !-- iPillar
  
    oP = 0
    oI = 0
    do iG = 1, nGroups 
      oC = 0
      do iB = 1, nBricks ( 2 )
        iR = ( iG - 1 ) * nBricks ( 2 )  +  iB  -  1
        do iPS = 1, nSegmentsFrom_2 ( iR )
          iP = oP + iPS
          nCoarsen_2 ( iP ) &
            = min ( int ( Incoming ( oI + 1 ) + 0.5_KDR ), nCells ( 2 ) )
          oI = oI + 1
          do iS = 1, CoarsenPillar_2 ( iP ) % nVariables
            CoarsenPillar_2 ( iP ) % Value ( oC + 1 : oC + nCB ( 2 ), iS )  &
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
      iP             !-- iPillar
    
    oP = 0
    oI = 0
    do iG = 1, nGroups
      oC = 0
      do iB = 1, nBricks ( 3 )
        iR = ( iG - 1 ) * nBricks ( 3 )  +  iB  -  1
        do iPS = 1, nSegmentsFrom_3 ( iR )
          iP = oP + iPS
          nCoarsen_3 ( iP ) &
            = min ( int ( Incoming ( oI + 1 ) + 0.5_KDR ), nCells ( 3 ) )
          oI = oI + 1
          do iS = 1, CoarsenPillar_3 ( iP ) % nVariables
            CoarsenPillar_3 ( iP ) % Value ( oC + 1 : oC + nCB ( 3 ), iS )  &
              =  Incoming ( oI + 1 : oI + nCB ( 3 ) )
            oI = oI + nCB ( 3 )
          end do !-- iS
        end do !-- iPS
        oC = oC + nCB ( 3 )
      end do !-- iB
      oP = oP + nSegmentsFrom_3 ( iR )
    end do !-- iG

  end procedure ComposePillarsUnpack_3_Kernel
  

end submodule FluidCentral_Kernel
