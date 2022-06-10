#include "Preprocessor"

submodule ( Coarsening_C__Form ) Coarsening_C__Kernel

  use Basics
  
  implicit none
  
contains 


  module procedure ComputeKernel

    integer ( KDI ) :: &
      iBC, &  !-- iBlockCoarsen
      iS, &   !-- iSelected
      iF, &   !-- iField
      iTh_1, iTh_2, &
      iPh_1, iPh_2, &
      iRad
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice  =  .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice  =  UseDeviceOption
      
    if ( UseDevice ) then

      !$OMP OMP_TARGET_DIRECTIVE parallel do collapse ( 2 ) &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP private ( iF, iTh_1, iTh_2, iPh_1, iPh_2, iRad )    
      do iBC  =  1, nBC
        do iS  =  1, size ( iaS )

          iF  =  iaS ( iS )

          iTh_1  =  oC ( 2 )  +  iTh ( 1, iBC )
          iTh_2  =  oC ( 2 )  +  iTh ( 2, iBC )
          iPh_1  =  oC ( 3 )  +  iPh ( 1, iBC )
          iPh_2  =  oC ( 3 )  +  iPh ( 2, iBC )
          iRad   =  oC ( 1 )  +  iR  (    iBC )

          FS_4D ( iRad, iTh_1 : iTh_2, iPh_1 : iPh_2, iF )  &
            =  sum ( dV_3D ( iRad, iTh_1 : iTh_2, iPh_1 : iPh_2 )  &
                     *  FS_4D ( iRad, iTh_1 : iTh_2, iPh_1 : iPh_2, iF ) )  &
               /  sum ( dV_3D ( iRad, iTh_1 : iTh_2, iPh_1 : iPh_2 ) )

        end do !-- iS
      end do !-- iBC
      !$OMP end OMP_TARGET_DIRECTIVE parallel do

    else 

      !$OMP parallel do collapse ( 2 ) &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP private ( iF, iTh_1, iTh_2, iPh_1, iPh_2, iRad )    
      do iBC  =  1, nBC
        do iS  =  1, size ( iaS )

          iF  =  iaS ( iS )

          iTh_1  =  oC ( 2 )  +  iTh ( 1, iBC )
          iTh_2  =  oC ( 2 )  +  iTh ( 2, iBC )
          iPh_1  =  oC ( 3 )  +  iPh ( 1, iBC )
          iPh_2  =  oC ( 3 )  +  iPh ( 2, iBC )
          iRad   =  oC ( 1 )  +  iR  (    iBC )

          FS_4D ( iRad, iTh_1 : iTh_2, iPh_1 : iPh_2, iF )  &
            =  sum ( dV_3D ( iRad, iTh_1 : iTh_2, iPh_1 : iPh_2 )  &
                     *  FS_4D ( iRad, iTh_1 : iTh_2, iPh_1 : iPh_2, iF ) )  &
               /  sum ( dV_3D ( iRad, iTh_1 : iTh_2, iPh_1 : iPh_2 ) )

        end do !-- iS
      end do !-- iBC
      !$OMP end parallel do

    end if !-- UseDevice

  end procedure ComputeKernel


end submodule Coarsening_C__Kernel

