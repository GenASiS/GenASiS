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
  type ( QuantityForm ), dimension ( 6 ) :: &
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
  call PrintStorage &
         ( S ( 2 ), &
           DescriptionOption = 'InitializeClone from ' &
                                // trim ( S ( 1 ) % Name ) )

  !-- Initialize
  call S ( 3 ) % Initialize &
         ( ValueShape = [ 2, 3 ], &
           NameOption = 'Storage_3' )
  call S ( 3 ) % AllocateDevice ( )
  call PrintStorage &
         ( S ( 3 ), &
           DescriptionOption = 'InitializeClone from ' &
                                // trim ( S ( 2 ) % Name ) )
  
  print*, 'Setting values of ' // trim ( S ( 3 ) % Name )
  do i = 1, S ( 3 ) % nVariables 
    call random_number ( S ( 3 ) % Value ( :, S ( 3 ) % iaSelected ( i ) ) )
    print*, 'Host Initial', S ( 3 ) % Value ( :, S ( 3 ) % iaSelected ( i ) )
  end do
  print*
  
  print*, 'Calling UpdateDevice ( ) of ' // trim ( S ( 3 ) % Name )
  call S ( 3 ) % UpdateDevice ( )
  print*
  
  print*, 'Clearing host values of ' // trim ( S ( 3 ) % Name )
  call Clear ( S ( 3 ) % Value )
  do i = 1, S ( 3 ) % nVariables 
    print*, 'Host Cleared', S ( 3 ) % Value ( :, S ( 3 ) % iaSelected ( i ) )
  end do
  print*
  
  print*, 'Calling UpdateHost ( ) of ' // trim ( S ( 3 ) % Name )
  call S ( 3 ) % UpdateHost ( )
  do i = 1, S ( 3 ) % nVariables 
    print*, 'Host From Device', S ( 3 ) % Value ( :, S ( 3 ) % iaSelected ( i ) )
  end do
  print*
  
  call S ( 4 ) % Initialize &
         ( ValueShape = [ 256**3, 16 ], PinnedOption = .true., &
           NameOption = 'Storage_4 Pinned' )
  call random_number ( S ( 4 ) % Value )
  call S ( 4 ) % AllocateDevice ( )
  call PrintStorage ( S ( 4 ) )
         
  call S ( 5 ) % Initialize &
         ( S ( 4 ), &
           iaSelectedOption = [ 1, 3, 5, 7, 9, 11, 13, 15 ], &
           NameOption = 'Storage_5' )
  call S ( 5 ) % AllocateDevice ( )
  call PrintStorage ( S ( 5 ), 'InitializeClone from ' // S ( 4 ) % Name )
  
  
  DataSize_GB = 1.0_KDR * S ( 4 ) % nValues * S ( 4 ) % nVariables * 8 &
                / 1.0e9_KDR 
  print*, 'Data size (GB)', DataSize_GB
  
  StartTime = OMP_GET_WTIME ( )
  call S ( 5 ) % UpdateDevice ( )
  TotalTime = max ( OMP_GET_WTIME ( ) - StartTime, epsilon ( 1.0_KDR ) )
  print*, 'Overlay Storage Data Transfer'
  print*, 'Device Error                   :', S ( 5 ) % ErrorDevice
  print*, 'Host-to-Device Time (s)        :', TotalTime
  print*, 'Host-to-Device Bandwith (GB/s) :', DataSize_GB / TotalTime / 2
  print*, ''
  
  StartTime = OMP_GET_WTIME ( )
  call S ( 4 ) % UpdateDevice ( )
  TotalTime = max ( OMP_GET_WTIME ( ) - StartTime, epsilon ( 1.0_KDR ) )
  print*, 'Primary Storage Data Transfer'
  print*, 'Device Error                   :', S ( 4 ) % ErrorDevice
  print*, 'Host-to-Device Time (s)        :', TotalTime
  print*, 'Host-to-Device Bandwith (GB/s) :', DataSize_GB / TotalTime
  print*, ''
  
  StartTime = OMP_GET_WTIME ( )
  call S ( 5 ) % UpdateHost ( )
  TotalTime = max ( OMP_GET_WTIME ( ) - StartTime, epsilon ( 1.0_KDR ) )
  print*, 'Overlay Storage Data Transfer'
  print*, 'Device Error                   :', S ( 5 ) % ErrorDevice
  print*, 'Device-to-Host Time (s)        :', TotalTime
  print*, 'Device-to-Host Bandwith (GB/s) :', DataSize_GB / TotalTime / 2
  print*, ''

  StartTime = OMP_GET_WTIME ( )
  call S ( 4 ) % UpdateHost ( )
  TotalTime = max ( OMP_GET_WTIME ( ) - StartTime, epsilon ( 1.0_KDR ) )
  print*, 'Primary Storage Data Transfer'
  print*, 'Device Error                   :', S ( 4 ) % ErrorDevice
  print*, 'Device-to-Host Time (s)        :', TotalTime
  print*, 'Device-to-Host Bandwith (GB/s) :', DataSize_GB / TotalTime
  print*, ''
  
  call S ( 4 ) % ReassociateHost ( AssociateVariablesOption = .false. )
  call PrintStorage ( S ( 4 ) )
  call S ( 4 ) % ReassociateHost ( AssociateVariablesOption = .true. )
  call PrintStorage ( S ( 4 ) )
  
!  !-- Re-initialize S ( 4 )
!  call S ( 6 ) % Initialize &
!         ( ValueShape = [ 4, 2 ], PinnedOption = .true., &
!           NameOption = 'Storage 6', VariableOption = Variable ( 1 : 2 ) )
!  S ( 6 ) % Value = 1.0_KDR
!  call S ( 6 ) % AllocateDevice ( )
!  call S ( 6 ) % UpdateDevice ( )
!  S ( 6 ) % Value = - huge ( 1.0_KDR )
!  
!  call Compute ( S ( 6 ) % Value ( :, 1 ) )
!  
!  print*
!  print*, 'After Compute:', S ( 6 ) % Value
!  
!  call PrintStorage ( S ( 6 ) )
!  
!  call S ( 6 ) % UpdateHost ( )
!  
!  print*
!  print*, 'After UpdateHost:', S ( 6 ) % Value
!  call PrintStorage ( S ( 6 ) )
  

contains


  subroutine PrintStorage ( S, DescriptionOption )

    use Specifiers
    use Storage_Form
    
    type ( StorageForm ), intent ( in ) :: &
      S
    character ( * ), intent ( in ), optional :: &
      DescriptionOption

    integer ( KDI ) :: &
      i
    integer ( KBI ) :: &
      Address
    character ( LDN ) :: &
      IndexLabel
    character ( LDB ) :: &
      Buffer
      
    print*
    
    if ( trim ( S % Name ) /= '' ) then
      print*, '===== Printing Storage : ' // trim ( S % Name ) // ' ====='
    else
      print*, '===== Printing Storage ===== '
    end if
    
    if ( present ( DescriptionOption ) ) then
      print*, trim ( DescriptionOption )
    end if
    
    print *
    print *, 'S % nValues         = ', S % nValues
    print *, 'S % nVariables      = ', S % nVariables
    print *, 'S % nVectors        = ', S % nVectors
    print *, 'S % lName           = ', S % lName
    print *, 'S % lVariable       = ', S % lVariable
    print *, 'S % lVector         = ', S % lVector
    print *, 'S % iaSelected      = ', S % iaSelected
    print *, 'S % AllocatedValue  = ', S % AllocatedValue
    print *, 'S % AllocatedDevice = ', S % AllocatedDevice
    print *, 'S % Pinned          = ', S % Pinned
    
    print *, 'S % Name            = ', trim ( S % Name )

    do i = 1, S % nVariables
      if ( trim ( S % Variable ( S % iaSelected ( i ) ) ) /= '' ) &
        print *, &
          'S % Variable (', i, ')  = ', &
          trim ( S % Variable ( S % iaSelected ( i ) ) )
    end do

    do i = 1, S % nVectors
      if ( trim ( S % Vector ( i ) ) /= '' ) &
        print *, &
          'S % Vector (', i, ')    = ', trim ( S % Vector ( i ) )
    end do

    do i = 1, S % nVariables
      if ( trim ( S % Unit ( S % iaSelected ( i ) ) % Unit ) /= '' ) &
        print *, &
          'S % Unit (', i, ') % Unit = ', &
          trim ( S % Unit ( S % iaSelected ( i ) ) % Unit )
    end do

    do i = 1, S % nVectors
      print *, &
        'S % VectorIndices (', i, ') = ', S % VectorIndices ( i ) % Value
    end do
    
    ! if ( allocated ( S % D_Selected ) ) then
    
    !   print*, 'S % D_Selected: ' 
    !   do i = 1, S % nVariables
    !     write ( IndexLabel, fmt = '( i7 )' ) i
    !     Address = transfer ( S % D_Selected ( i ), 1_KBI )
    !     write ( Buffer, fmt = ' ( z64 ) ' ) Address
    !     Buffer = '0x' //  adjustl ( Buffer )
    !     print &
    !       '(a38, a32)', &
    !       '( ' // trim ( adjustl ( IndexLabel ) ) // ' ) = ', Buffer
    !   end do
    
    ! end if
    
    print*

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
