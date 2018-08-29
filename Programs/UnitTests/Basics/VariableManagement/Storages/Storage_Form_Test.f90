program Storage_Form_Test

  use Specifiers
  use ArrayOperations
  use ArrayArrays
  use Storage_Form

  implicit none

  integer ( KDI ) :: &
    i
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
  type ( StorageForm ), dimension ( 5 ) :: &
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


end program Storage_Form_Test
