program Storage_Form_Test

  use OMP_LIB
  use Specifiers
  use ArrayOperations
  use ArrayArrays
  use Storage_Form

  implicit none

  integer ( KDI ) :: &
    i
  real ( KDR ) :: &
    StartTime, &
    TotalTime, &
    DataSize_GB
  character ( LDL ), dimension ( 6 ) :: &
    Variable &
      = [ 'Variable_1', 'Variable_2', 'Variable_3', 'Variable_4', &
          'Variable_5', 'Variable_6'  ], &
    UnitString &
      = [ 'Unit_1', 'Unit_2', 'Unit_3', 'Unit_4', 'Unit_5', 'Unit_6'  ]
  type ( MeasuredValueForm ), dimension ( 6 ) :: &
    VariableUnit
  type ( Integer_1D_Form ), dimension ( 1 ) :: &
    VectorIndices
  type ( StorageForm ), dimension ( 6 ) :: &
    S

  do i = 1, size ( VariableUnit ) 
    call VariableUnit ( i ) % Initialize &
           ( UnitString ( i ), UnitString ( i ), 1.0_KDR )
  end do

  call VectorIndices ( 1 ) % Initialize ( [ 2, 3 ] )

  !-- InitializeAllocate
  call S ( 1 ) % Initialize &
         ( ValueShape = [ 4, 6 ], VectorIndicesOption = VectorIndices, &
           UnitOption = VariableUnit, VectorOption = [ 'Vector_1' ], &
           VariableOption = Variable, NameOption = 'Storage_1')
  call PrintStorage ( S ( 1 ) )

  !-- InitializeClone
  call S ( 2 ) % Initialize &
         ( S ( 1 ), NameOption = 'Storage_2', &
           iaSelectedOption = [ 1, 2, 3, 5 ] )
  call S ( 2 ) % AllocateDevice ( )
  call PrintStorage ( S ( 2 ) )

  !-- InitializeClone, take 2
  call S ( 3 ) % Initialize &
         ( S ( 2 ), NameOption = 'Storage_3', &
           iaSelectedOption = [ 1, 3, 6 ] )
  call S ( 3 ) % AllocateDevice ( )
  call PrintStorage ( S ( 3 ) )
  
  do i = 1, S ( 3 ) % nVariables 
    call random_number ( S ( 3 ) % Value ( :, S ( 3 ) % iaSelected ( i ) ) )
    print*, 'Host Initial', S ( 3 ) % Value ( :, S ( 3 ) % iaSelected ( i ) )
  end do
  call S ( 3 ) % UpdateDevice ( )
  
  call Clear ( S ( 3 ) % Value )
  do i = 1, S ( 3 ) % nVariables 
    print*, 'Host Cleared', S ( 3 ) % Value ( :, S ( 3 ) % iaSelected ( i ) )
  end do
  
  call S ( 3 ) % UpdateHost ( )
  do i = 1, S ( 3 ) % nVariables 
    print*, 'Host From Device', S ( 3 ) % Value ( :, S ( 3 ) % iaSelected ( i ) )
  end do
  
  call S ( 4 ) % Initialize &
         ( ValueShape = [ 256**3, 16 ], PinnedOption = .true. )
  call random_number ( S ( 4 ) % Value )
  call S ( 4 ) % AllocateDevice ( )
         
  
  call S ( 5 ) % Initialize &
         ( S ( 4 ), &
           iaSelectedOption = [ 1, 3, 5, 7, 9, 11, 13, 15 ] )
  call S ( 5 ) % AllocateDevice ( )
  
  DataSize_GB = 1.0_KDR * S ( 4 ) % nValues * S ( 4 ) % nVariables * 8 &
                / 1.0e9_KDR 
  print*, 'Data size (GB)', DataSize_GB
  
  StartTime = OMP_GET_WTIME ( )
  call S ( 5 ) % UpdateDevice ( )
  TotalTime = OMP_GET_WTIME ( ) - StartTime
  print*, 'Overlay Storage Data Transfer'
  print*, 'Device Error                   :', S ( 5 ) % ErrorDevice
  print*, 'Host-to-Device Time (s)        :', TotalTime
  print*, 'Host-to-Device Bandwith (GB/s) :', DataSize_GB / TotalTime / 2
  print*, ''
  
  StartTime = OMP_GET_WTIME ( )
  call S ( 4 ) % UpdateDevice ( )
  TotalTime = OMP_GET_WTIME ( ) - StartTime
  print*, 'Primary Storage Data Transfer'
  print*, 'Device Error                   :', S ( 4 ) % ErrorDevice
  print*, 'Host-to-Device Time (s)        :', TotalTime
  print*, 'Host-to-Device Bandwith (GB/s) :', DataSize_GB / TotalTime
  print*, ''
  
  StartTime = OMP_GET_WTIME ( )
  call S ( 5 ) % UpdateHost ( )
  TotalTime = OMP_GET_WTIME ( ) - StartTime
  print*, 'Overlay Storage Data Transfer'
  print*, 'Device Error                   :', S ( 5 ) % ErrorDevice
  print*, 'Device-to-Host Time (s)        :', TotalTime
  print*, 'Device-to-Host Bandwith (GB/s) :', DataSize_GB / TotalTime / 2
  print*, ''

  StartTime = OMP_GET_WTIME ( )
  call S ( 4 ) % UpdateHost ( )
  TotalTime = OMP_GET_WTIME ( ) - StartTime
  print*, 'Primary Storage Data Transfer'
  print*, 'Device Error                   :', S ( 4 ) % ErrorDevice
  print*, 'Device-to-Host Time (s)        :', TotalTime
  print*, 'Device-to-Host Bandwith (GB/s) :', DataSize_GB / TotalTime
  print*, ''
  
  
  !-- Re-initialize S ( 4 )
  call S ( 6 ) % Initialize &
         ( ValueShape = [ 4, 2 ], PinnedOption = .true., &
           NameOption = 'Storage 6', VariableOption = Variable ( 1 : 2 ) )
  S ( 6 ) % Value = 1.0_KDR
  call S ( 6 ) % AllocateDevice ( )
  call S ( 6 ) % UpdateDevice ( )
  S ( 6 ) % Value = - huge ( 1.0_KDR )
  
  call Compute ( S ( 6 ) % Value ( :, 1 ) )
  
  print*
  print*, 'After Compute:', S ( 6 ) % Value
  
  !call PrintStorage ( S ( 6 ) )
  
  call S ( 6 ) % UpdateHost ( )
  
  print*
  print*, 'After UpdateHost:', S ( 6 ) % Value
  !call PrintStorage ( S ( 6 ) )
  

contains


  subroutine PrintStorage ( S )

    use Specifiers
    use Storage_Form
    
    type ( StorageForm ) :: &
      S

    integer ( KDI ) :: &
      i
    integer ( KBI ) :: &
      Address
    character ( LDN ) :: &
      IndexLabel
    character ( LDB ) :: &
      Buffer

    print *
    print *, 'S % nValues = ', S % nValues
    print *, 'S % nVariables = ', S % nVariables
    print *, 'S % nVectors = ', S % nVectors
    print *, 'S % lName = ', S % lName
    print *, 'S % lVariable = ', S % lVariable
    print *, 'S % lVector = ', S % lVector
    print *, 'S % iaSelected = ', S % iaSelected
    print *, 'S % AllocatedValue = ', S % AllocatedValue
    print *, 'S % Name = ', trim ( S % Name )

    do i = 1, S % nVariables
      print *, &
        'S % Variable (', i, ') = ', &
        trim ( S % Variable ( S % iaSelected ( i ) ) )
    end do

    do i = 1, S % nVectors
      print *, &
        'S % Vector (', i, ') = ', trim ( S % Vector ( i ) )
    end do

    do i = 1, S % nVariables
      print *, &
        'S % Unit (', i, ') % Unit = ', &
        trim ( S % Unit ( S % iaSelected ( i ) ) % Unit )
    end do

    do i = 1, S % nVectors
      print *, &
        'S % VectorIndices (', i, ') = ', S % VectorIndices ( i ) % Value
    end do
    
    if ( allocated ( S % D_Selected ) ) then
    
      print*, 'S % D_Selected: ' 
      do i = 1, S % nVariables
        write ( IndexLabel, fmt = '( i7 )' ) i
        Address = transfer ( S % D_Selected ( i ), 1_KBI )
        write ( Buffer, fmt = ' ( z64 ) ' ) Address
        Buffer = '0x' //  adjustl ( Buffer )
        print &
          '(a38, a32)', &
          '( ' // trim ( adjustl ( IndexLabel ) ) // ' ) = ', Buffer
      end do
    
    end if

  end subroutine PrintStorage
  
  
  subroutine Compute ( A )
    
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      A
      
    integer ( KDI ) :: &
      iV
    
    !$OMP target teams distribute parallel do   
    do iV = 1, size ( A )
      A ( iV ) = A ( iV ) * iV
    end do
    !$OMP end target teams distribute parallel do   
  
  end subroutine
  

end program Storage_Form_Test
