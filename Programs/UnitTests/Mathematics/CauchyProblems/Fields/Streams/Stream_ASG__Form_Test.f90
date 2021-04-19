program Stream_ASG__Form_Test

  !-- Stream_AtlasSingleGrid__Form_Test

  use Basics
  use Manifolds
  use FieldSets
  use Streams

  implicit none

  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Atlas_SG_Form ), allocatable :: &
    A
  type ( FieldSet_ASG_Form ), allocatable :: &
    FSA
  type ( Stream_ASG_Form ), allocatable :: &
    SA

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Stream_ASG__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( A )
  call A % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator, &
           PeriodicOption = [ .true., .true., .true. ] )

  allocate ( FSA )
  call FSA % Initialize ( A )

  allocate ( SA )
  call SA % Initialize ( A, GIS )
  call SA % AddFieldSet ( FSA )

  call   A % Show ( )
  call FSA % Show ( )
  call  SA % Show ( )

  call SetFieldSet ( FSA )
  call GIS % Open ( GIS % ACCESS_CREATE )
  call SA % Write ( )
  call GIS % Close ( )

  deallocate ( FSA )
  deallocate ( A )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )

contains


  subroutine SetFieldSet ( FSA )

    class ( FieldSet_ASG_Form ), intent ( inout ) :: &
      FSA

    integer ( KDI ) :: &
      iS, &          !-- iSelected
      iF, &          !-- iField
      iC, jC, kC     !-- iCell, etc.
    integer ( KDI ), dimension ( 3 ) :: &
      oC
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      F_3D  !-- Field

    associate ( FSG  =>  FSA % FieldSet_G )

    call Show ( 'Set FieldSet' )
    call Show ( FSG % Name, 'FieldSet' )
    call Clear ( FSG % FieldSet % Value )

    select type ( G  =>  FSG % Chart )
    class is ( Grid_S_Form )

    associate ( nCB  =>  G % nCellsBrick )

    do iS  =  1, FSG % nFields
      iF  =  FSG % iaSelected ( iS )
      associate ( F  =>  FSG % FieldSet % Value ( :, iF ) )
      call G % SetFieldPointer ( F, F_3D )

      oC  =  ( G % iaBrick  -  1 )  *  nCB
      do kC  =  1,  nCB ( 3 )
        do jC  =  1,  nCB ( 2 )
          do iC  =  1,  nCB ( 1 )
            F_3D ( iC, jC, kC )  &
              =  iS  *  (    1.e0  *  ( oC ( 1 )  +  iC  -  1 )  &
                          +  1.e2  *  ( oC ( 2 )  +  jC  -  1 )  &
                          +  1.e4  *  ( oC ( 3 )  +  kC  -  1 ) )
          end do !-- iC
        end do !-- jC
      end do !-- kC

      end associate !-- F, etc.
    end do !-- iF

    call FSG % ExchangeGhostData ( )

    end associate !-- nCB
    end select !-- G
    end associate !-- FSG

    nullify ( F_3D )

  end subroutine SetFieldSet


end program Stream_ASG__Form_Test
