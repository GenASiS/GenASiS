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
    character ( LDF ) :: &
      EOS_Filename

    call Show ( 'Reading original table' )

    allocate ( EOS_In )
    
    EOS_Filename &
      = '../Parameters/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5'
    call PROGRAM_HEADER % GetParameter ( EOS_Filename, 'EOS_Filename' )

    call EOS_In % Initialize ( EOS_Filename )

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
      (      nV  =>  EOS_Out % N_VARIABLES, &  !-- Automatically defined
             nD  =>  EOS_Out % nDensity, &
             nT  =>  EOS_Out % nTemperature, &
             nY  =>  EOS_Out % nElectronFraction, &
        MinLogD  =>  EOS_Out % MinLogDensity, &
        MaxLogD  =>  EOS_Out % MaxLogDensity, &
        MinLogT  =>  EOS_Out % MinLogTemperature, &
        MaxLogT  =>  EOS_Out % MaxLogTemperature, &
           MinY  =>  EOS_Out % MinElectronFraction, &
           MaxY  =>  EOS_Out % MaxElectronFraction, &
        E_Shift  =>  EOS_Out % EnergyShift )

    E_Shift  =  EOS_In % EnergyShift

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
    call Show ( LogD, 'Table densities', &
                IgnorabilityOption = CONSOLE % INFO_3 )

    dLogT = ( MaxLogT - MinLogT ) / dble ( nT - 1 )
    do iT = 1, nT
      LogT ( iT )  =  MinLogT + dble ( iT - 1 ) * dLogT
    end do
    call Show ( LogT, 'Table temperatures', &
                IgnorabilityOption = CONSOLE % INFO_3 )

    dY = ( MaxY - MinY ) / dble ( nY - 1 )
    do iY = 1, nY
      Y ( iY )  =  MinY + dble ( iY - 1 ) * dY
    end do
    call Show ( Y, 'Table electron fractions', &
               IgnorabilityOption = CONSOLE % INFO_3 )

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
          !-- Variable 1: log pressure
          EOS_Out % Table ( iD, iT, iY, 1 )  &
            =  log10 ( Fluid % Value ( iF, 1 ) )
          !-- Variable 2: log ( energy + e_shift )
          EOS_Out % Table ( iD, iT, iY, 2 )  &
            =  log10 ( Fluid % Value ( iF, 2 ) + E_Shift )
          !-- All other variables
          do iV = 3, nV
            EOS_Out % Table ( iD, iT, iY, iV )  &
              =  Fluid % Value ( iF, iV )
          end do !-- iV
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
      iC, &   !-- iCharacter
      iV, &   !-- iVariable
      Error, &
      DataRank
    integer ( HID_T ) :: &
      FileHandle
      
    character ( LDF ) :: &
      Filename
    character ( LDL ) :: &
      nD_String, &
      nT_String, &
      nYe_String
    character ( LDL ), dimension ( 17 ) :: &
      VariableName
      
    iC = index ( EOS_In % Filename, '_' )
    
    write ( nD_String,  fmt = '( i3 )' ) EOS_Out % nDensity
    write ( nT_String,  fmt = '( i3 )' ) EOS_Out % nTemperature
    write ( nYe_String, fmt = '( i3 )' ) EOS_Out % nElectronFraction
    
    Filename = EOS_In % Filename ( 1 : iC ) // 'Coarsen_'  &
               // trim ( adjustl ( nD_String ) ) // 'r_'  &
               // trim ( adjustl ( nT_String ) ) // 't_'  &
               // trim ( adjustl ( nYe_String ) ) // 'y.h5'
    
    EOS_Out % Filename = Filename

    call Show ( 'Writing coarsened table' )
    call Show ( Filename, 'Filename' )

    call H5OPEN_F ( Error )

    call H5FCREATE_F ( trim ( adjustl ( Filename ) ), H5F_ACC_TRUNC_F, &
                      FileHandle, error )
    
    call WriteVariable_0D_HDF5 &
           ( EOS_Out % nDensity, 'pointsrho', FileHandle, Error )
    call WriteVariable_0D_HDF5 &
           ( EOS_Out % nTemperature, 'pointstemp', FileHandle, Error )
    call WriteVariable_0D_HDF5 &
           ( EOS_Out % nElectronFraction, 'pointsye', FileHandle, Error )
    
    call WriteVariable_0D_HDF5 &
           ( EOS_Out % EnergyShift, 'energy_shift', FileHandle, Error )
        
    call WriteVariable_1D_HDF5 &
           ( EOS_Out % LogDensity, 'logrho', FileHandle, Error )
    call WriteVariable_1D_HDF5 &
           ( EOS_Out % LogTemperature, 'logtemp', FileHandle, Error )
    call WriteVariable_1D_HDF5 &
           ( EOS_Out % ElectronFraction, 'ye', FileHandle, Error )
           
    VariableName = [ 'logpress                       ', & 
                     'logenergy                      ', &
                     'entropy                        ', &
                     'munu                           ', &
                     'cs2                            ', &
                     'dedt                           ', &
                     'dpdrhoe                        ', &
                     'dpderho                        ', &
                     'muhat                          ', &
                     'mu_e                           ', &
                     'mu_p                           ', &
                     'mu_n                           ', &
                     'Xa                             ', &
                     'Xh                             ', &
                     'Xn                             ', &
                     'Xp                             ', &
                     'gamma                          ' ]
    
    do iV = 1, size ( VariableName )
      if ( iV == size ( VariableName ) ) then
        !-- Some variable in EOS is skipped / not used
        call WriteVariable_3D_HDF5 &
             ( EOS_Out % Table ( :, :, :,  19 ), &
               trim ( VariableName ( iV ) ), FileHandle, Error )
      else
        call WriteVariable_3D_HDF5 &
             ( EOS_Out % Table ( :, :, :,  iV ), &
               trim ( VariableName ( iV ) ), FileHandle, Error )
      end if
    end do
    
    
    call H5FCLOSE_F ( FileHandle, Error )
    call H5CLOSE_F ( Error )

    deallocate ( EOS_Out )
    deallocate ( EOS_In )

  end subroutine WriteTable
  
  
  subroutine WriteVariable_0D_HDF5 ( Value, Name, Location, Error )
  
    class ( * ), intent ( in ), target :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    integer ( HID_T ), intent ( in ) :: &
      Location
    integer ( KDI ) , intent ( out ) :: &
      Error
      
    integer ( HID_T ) :: &
      Set_ID, &
      Space_ID
    integer ( HSIZE_T ), dimension ( rank ( Value ) ) :: &
      Dimensions
      
    Dimensions = shape ( Value )
    call H5SCREATE_SIMPLE_F ( rank ( Value ), Dimensions, Space_ID, Error )
    call H5DCREATE_F ( Location, trim ( Name ), H5T_NATIVE_DOUBLE, &
                       Space_ID, Set_ID, Error )
    
    select type ( Value )
    type is ( real ( KDR ) )
      call H5DWRITE_F ( Set_ID, H5T_NATIVE_DOUBLE, Value, Dimensions, Error )
    type is ( integer ( KDI ) )
      call H5DWRITE_F ( Set_ID, H5T_NATIVE_INTEGER, Value, Dimensions, Error )
    end select
     
    call H5DCLOSE_F ( Set_ID, Error )
    call H5SCLOSE_F ( Space_ID, Error )
  
  end subroutine WriteVariable_0D_HDF5
  
  
  subroutine WriteVariable_1D_HDF5 ( Value, Name, Location, Error )
  
    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    integer ( HID_T ), intent ( in ) :: &
      Location
    integer ( KDI ) , intent ( out ) :: &
      Error
      
    integer ( HID_T ) :: &
      Set_ID, &
      Space_ID
    integer ( HSIZE_T ), dimension ( rank ( Value ) ) :: &
      Dimensions
      
    Dimensions = shape ( Value )
    call H5SCREATE_SIMPLE_F ( rank ( Value ), Dimensions, Space_ID, Error )
    call H5DCREATE_F ( Location, trim ( Name ), H5T_NATIVE_DOUBLE, &
                       Space_ID, Set_ID, Error )
    call H5DWRITE_F ( Set_ID, H5T_NATIVE_DOUBLE, Value, Dimensions, Error )
    call H5DCLOSE_F ( Set_ID, Error )
    call H5SCLOSE_F ( Space_ID, Error )
  
  end subroutine WriteVariable_1D_HDF5
  
  
  subroutine WriteVariable_3D_HDF5 ( Value, Name, Location, Error )
  
    real ( KDR ), dimension ( :, :, : ), intent ( in ), target :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    integer ( HID_T ), intent ( in ) :: &
      Location
    integer ( KDI ) , intent ( out ) :: &
      Error
      
    integer ( HID_T ) :: &
      Set_ID, &
      Space_ID
    integer ( HSIZE_T ), dimension ( rank ( Value ) ) :: &
      Dimensions
      
    Dimensions = shape ( Value )
    call H5SCREATE_SIMPLE_F ( rank ( Value ), Dimensions, Space_ID, Error )
    call H5DCREATE_F ( Location, trim ( Name ), H5T_NATIVE_DOUBLE, &
                       Space_ID, Set_ID, Error )
    call H5DWRITE_F ( Set_ID, H5T_NATIVE_DOUBLE, Value, Dimensions, Error )
    call H5DCLOSE_F ( Set_ID, Error )
    call H5SCLOSE_F ( Space_ID, Error )
  
  end subroutine WriteVariable_3D_HDF5


end program TableGeneration_EOS
