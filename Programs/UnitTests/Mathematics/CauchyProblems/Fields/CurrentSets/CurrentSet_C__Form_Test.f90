program CurrentSet_C__Form_Test

  !-- CurrentSet_Chart__Form_Test

  use Basics
  use Manifolds
  use Streams
  use Geometries
  use CurrentSets

  implicit none

  type ( MeasuredValueForm ) :: &
    DensityUnit
  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Chart_GS_Form ), allocatable :: &
    C
  type ( Stream_C_Form ), allocatable :: &
    SC
  type ( Geometry_F_C_Form ), allocatable :: &
    GC
  type ( CurrentSet_C_Form ), allocatable :: &
    CSC

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'CurrentSet_C__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( C )
  call C % Initialize &
         ( PeriodicOption = [ .true., .true., .true. ] )

  allocate ( SC )
  call SC % Initialize ( C, GIS )

  allocate ( GC )
  call GC % Initialize ( C )

  DensityUnit  =  UNIT % MASS_DENSITY_MKS

  allocate ( CSC )
  call CSC % Initialize ( GC, DensityUnitOption = DensityUnit )
  call CSC % SetStream ( SC )

  call   C % Show ( )
  call CSC % Show ( )
  call  SC % Show ( )

  call SetWave ( CSC, GC )

  call GIS % Open ( GIS % ACCESS_CREATE )
  call SC % Write ( )
  call GIS % Close ( )

  deallocate ( CSC )
  deallocate ( GC )
  deallocate ( SC )
  deallocate ( C )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )


contains


  subroutine SetWave ( CSC, GC )

    class ( CurrentSet_C_Form ), intent ( inout ) :: &
      CSC
    class ( Geometry_F_C_Form ), intent ( in ), target :: &
      GC

    integer ( KDI ), dimension ( 3 ) :: &
      nWavelengths
    real ( KDR ) :: &
      Offset, &
      Amplitude, &
      Speed
    real ( KDR ), dimension ( 3 ) :: &
      Wavenumber

    select type ( C  =>  CSC % Chart )
    class is ( Chart_GS_Form )

    nWavelengths  =  0
    nWavelengths ( 1 : C % nDimensions )  =  1
    call PROGRAM_HEADER % GetParameter ( nWavelengths, 'nWavelengths' )

    Offset     =  2.0_KDR
    Amplitude  =  1.0_KDR
    Speed      =  1.0_KDR
    call PROGRAM_HEADER % GetParameter ( Offset, 'Offset' )
    call PROGRAM_HEADER % GetParameter ( Amplitude, 'Amplitude' )
    call PROGRAM_HEADER % GetParameter ( Speed, 'Speed' )

    associate ( BoxSize  =>  C % MaxCoordinate  -  C % MinCoordinate )
    where ( BoxSize  >  0.0_KDR )
      Wavenumber  =  nWavelengths / BoxSize
    elsewhere
      Wavenumber  =  0.0_KDR
    end where
    end associate !-- BoxSize

    associate &
      (     X  =>   GC % Storage_FSC % Storage &
                       % Value ( :, GC % CENTER_U_1 ), &
            Y  =>   GC % Storage_FSC % Storage &
                       % Value ( :, GC % CENTER_U_2 ), &
            Z  =>   GC % Storage_FSC % Storage &
                       % Value ( :, GC % CENTER_U_3 ), &
          Rho  =>  CSC % Storage_FSC % Storage &
                       % Value ( :, CSC % DENSITY_DEFAULT ), &
            V  =>  CSC % VelocityDefault_U, &
            K  =>  Wavenumber, &
        Abs_K  =>  sqrt ( dot_product ( Wavenumber, Wavenumber ) ), &
        TwoPi  =>  2.0_KDR  *  CONSTANT % PI )

    Rho  =  Offset  &
            +  Amplitude  &
               *  sin ( TwoPi * (    K ( 1 ) * X  &
                                  +  K ( 2 ) * Y  &
                                  +  K ( 3 ) * Z  ) )

    V ( 1 )  =  Speed  *  K ( 1 )  /  Abs_K
    V ( 2 )  =  Speed  *  K ( 2 )  /  Abs_K
    V ( 3 )  =  Speed  *  K ( 3 )  /  Abs_K
    
    end associate !-- Rho, etc.
    end select !-- C

  end subroutine SetWave


end program CurrentSet_C__Form_Test
