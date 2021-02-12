program TableGeneration_EOS

  use HDF5
  use GenASiS
  implicit none

  type ( EOS_P_HN_OConnorOtt_Form ), allocatable :: &
    EOS_In, &
    EOS_Out

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'TableGeneration_NuLib' )

  call ReadTable ( EOS_In )
  call CoarsenTable ( EOS_In, EOS_Out )
  call WriteTable ( EOS_Out, EOS_In )
  deallocate ( PROGRAM_HEADER )

contains

  
  subroutine ReadTable ( EOS_In )

    type ( EOS_P_HN_OConnorOtt_Form ), intent ( out ), allocatable :: &
      EOS_In

    integer ( KDI ) :: &
      i
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaSelected

    call Show ( 'Reading original table' )

    allocate ( EOS_In )
    call EOS_In % Initialize ( )

    allocate ( iaSelected ( EOS_In % N_VARIABLES ) )
    iaSelected = [ ( i, i = 1, size ( iaSelected ) ) ]
    call EOS_In % SelectVariables ( iaSelected, iaSelected )

  end subroutine ReadTable


  subroutine CoarsenTable ( EOS_In, EOS_Out )

    type ( EOS_P_HN_OConnorOtt_Form ), intent ( inout ) :: &
      EOS_In
    type ( EOS_P_HN_OConnorOtt_Form ), intent ( out ), allocatable :: &
      EOS_Out

    integer ( KDI ) :: &
      iD, iT, iY, iV, &
      iF, &
      DensityFactor
    real ( KDR ) :: &
      dLogD, dLogT, dY
    type ( StorageForm ) :: &
      Fluid

    call Show ( 'Preparing coarsened table' )

    allocate ( EOS_Out )

    associate &
      (      nD  =>  EOS_Out % nDensity, &
             nT  =>  EOS_Out % nTemperature, &
             nY  =>  EOS_Out % nElectronFraction, &
        MinLogD  =>  EOS_Out % MinLogDensity, &
        MaxLogD  =>  EOS_Out % MaxLogDensity, &
        MinLogT  =>  EOS_Out % MinLogTemperature, &
        MaxLogT  =>  EOS_Out % MaxLogTemperature, &
           MinY  =>  EOS_Out % MinElectronFraction, &
           MaxY  =>  EOS_Out % MaxElectronFraction, &
             nV  =>  EOS_Out % N_VARIABLES )

    EOS_Out % EnergyShift  =  EOS_In % EnergyShift

    DensityFactor  =  2
    call PROGRAM_HEADER % GetParameter ( DensityFactor, 'DensityFactor' )

    nD  =  EOS_In % nDensity  /  DensityFactor
    nT  =  EOS_In % nTemperature
    nY  =  EOS_In % nElectronFraction

    MinLogD  =  EOS_In % MinLogDensity
    MaxLogD  =  EOS_In % MaxLogDensity
    MinLogT  =  EOS_In % MinLogTemperature
    MaxLogT  =  EOS_In % MaxLogTemperature
      MinY  =  EOS_In % MinElectronFraction
      MaxY  =  EOS_In % MaxElectronFraction

    call Show ( 'Original table' )
    call Show ( EOS_In % nDensity, 'nDensity' )
    call Show ( EOS_In % nTemperature, 'nTemperature' )
    call Show ( EOS_In % nElectronFraction, 'nElectronFraction' )
    call Show ( EOS_In % MinLogDensity, 'MinLogDensity' )
    call Show ( EOS_In % MaxLogDensity, 'MaxLogDensity' )
    call Show ( EOS_In % MinLogTemperature, 'MinLogTemperature' )
    call Show ( EOS_In % MaxLogTemperature, 'MaxLogTemperature' )
    call Show ( EOS_In % MinElectronFraction, 'MinElectronFraction' )
    call Show ( EOS_In % MaxElectronFraction, 'MaxElectronFraction' )

    call Show ( 'Coarsened table' )
    call Show ( EOS_Out % nDensity, 'nDensity' )
    call Show ( EOS_Out % nTemperature, 'nTemperature' )
    call Show ( EOS_Out % nElectronFraction, 'nElectronFraction' )
    call Show ( EOS_Out % MinLogDensity, 'MinLogDensity' )
    call Show ( EOS_Out % MaxLogDensity, 'MaxLogDensity' )
    call Show ( EOS_Out % MinLogTemperature, 'MinLogTemperature' )
    call Show ( EOS_Out % MaxLogTemperature, 'MaxLogTemperature' )
    call Show ( EOS_Out % MinElectronFraction, 'MinElectronFraction' )
    call Show ( EOS_Out % MaxElectronFraction, 'MaxElectronFraction' )

    allocate ( EOS_Out % LogDensity ( nD ) )
    allocate ( EOS_Out % LogTemperature ( nT ) )
    allocate ( EOS_Out % ElectronFraction ( nY ) )
    associate &
      ( LogD  =>  EOS_Out % LogDensity, &
        LogT  =>  EOS_Out % LogTemperature, &
           Y  =>  EOS_Out % ElectronFraction )

    dLogD = ( MaxLogD - MinLogD ) / dble ( nD - 1 )
    do iD = 1, nD
      LogD ( iD )  =  MinLogD + dble ( iD - 1 ) * dLogD
    end do
    call Show ( LogD, 'Table densities' )

    dLogT = ( MaxLogT - MinLogT ) / dble ( nT - 1 )
    do iT = 1, nT
      LogT ( iT )  =  MinLogT + dble ( iT - 1 ) * dLogT
    end do
    call Show ( LogT, 'Table temperatures' )

    dY = ( MaxY - MinY ) / dble ( nY - 1 )
    do iY = 1, nY
      Y ( iY )  =  MinY + dble ( iY - 1 ) * dY
    end do
    call Show ( Y, 'Table electron fractions' )

    allocate ( EOS_Out % Table ( nD, nT, nY, nV ) )

    call Fluid % Initialize ( [ nD * nT * nY, nV + 3 ] )
    call Show ( shape ( Fluid % Value ), 'Fluid shape' )
    associate &
      (   D  =>  Fluid % Value ( :, nV + 1 ), &
          T  =>  Fluid % Value ( :, nV + 2 ), &
         Ye  =>  Fluid % Value ( :, nV + 3 ) )

    iF  =  0
    do iY  =  1, nY
      do iT  =  1, nT
        do iD  =  1, nD
          iF  =  iF + 1
           D ( iF )  =  10.0_KDR ** LogD ( iD )
           T ( iF )  =  10.0_KDR ** LogT ( iT )
          Ye ( iF )  =  Y ( iY )
        end do !-- iD
      end do !-- iT
    end do !-- iY

    call EOS_In % ComputeFromTemperature ( Fluid, nV + [ 1, 2, 3 ] )

    iF  =  0
    do iY  =  1, nY
      do iT  =  1, nT
        do iD  =  1, nD
          iF  =  iF + 1
          EOS_Out % Table ( iD, iT, iY, 1 : nV )  &
            =  Fluid % Value ( iF, 1 : nV )
        end do !-- iD
      end do !-- iT
    end do !-- iY
    
    end associate !-- D, etc.
    end associate !-- LogD, etc.
    end associate !-- nD, etc.

  end subroutine CoarsenTable


  subroutine WriteTable ( EOS_Out, EOS_In )

    type ( EOS_P_HN_OConnorOtt_Form ), intent ( inout ), allocatable :: &
      EOS_Out, &
      EOS_In

    integer ( KDI ) :: &
      error, &
      accerr
    integer ( HID_T ) :: &
      file_id, &
      dset_id, dspace_id
    integer ( HSIZE_T ), dimension ( 1 ) :: &
      dims1
    integer ( HSIZE_T ), dimension ( 3 ) :: &
      dims3
    character ( LDF ) :: &
      Filename

    Filename = '../Parameters/LS_Coarse.h5'

    call Show ( 'Writing coarsened table' )
    call Show ( Filename, 'Filename' )

    call h5open_f ( error )

    call h5fcreate_f ( trim ( adjustl ( Filename ) ), H5F_ACC_TRUNC_F, &
                       file_id, error )

    associate &
      ( nrho          => EOS_Out % nDensity, &
        ntemp         => EOS_Out % nTemperature, &
        nye           => EOS_Out % nElectronFraction, &
        nvars         => EOS_Out % N_VARIABLES, &
        energy_shift  => EOS_Out % EnergyShift )

    end associate !-- nrho, etc.

    deallocate ( EOS_Out )
    deallocate ( EOS_In )

  end subroutine WriteTable


end program TableGeneration_EOS
