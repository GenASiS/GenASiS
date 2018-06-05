program PackedStorage_Form_Test

  use Specifiers
  use ArrayArrays
  use Storage_Form
  use PackedStorage_Form

  implicit none

  integer ( KDI ) :: &
    i
  integer ( KDI ), dimension ( 2 ) :: &
    UnpackedIndex
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
  type ( StorageForm ), dimension ( 2 ) :: &
    S
  type ( PackedStorageForm ) :: &
    PS

  do i = 1, size ( VariableUnit ) 
    call VariableUnit ( i ) % Initialize &
           ( UnitString ( i ), UnitString ( i ), 1.0_KDR ) 
  end do

  call VectorIndices ( 1 ) % Initialize ( [ 2, 3 ] )

  call S ( 1 ) % Initialize &
         ( ValueShape = [ 4, 6 ], VectorIndicesOption = VectorIndices, &
           UnitOption = VariableUnit, VectorOption = [ 'Vector_1' ], &
           VariableOption = Variable, NameOption = 'Storage_1')
  call PrintStorage ( S ( 1 ) )

  call S ( 2 ) % Initialize &
         ( S ( 1 ), NameOption = 'Storage_2', &
           iaSelectedOption = [ 2, 3, 6 ] )
  call PrintStorage ( S ( 2 ) )

  UnpackedIndex = [ 1, 3 ]
  call PS % Initialize &
         ( UnpackedIndex, S ( 2 ) % nValues, S ( 2 ) % nVariables )
  call PrintPackedStorage ( PS )

contains


  subroutine PrintStorage ( S )

    use Specifiers
    use Storage_Form

    type ( StorageForm ) :: &
      S

    integer ( KDI ) :: &
      i

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

  end subroutine PrintStorage


  subroutine PrintPackedStorage ( PS )

    use Specifiers
    use PackedStorage_Form

    type ( PackedStorageForm ) :: &
      PS

    integer ( KDI ) :: &
      i

    print *
    print *, 'PackedStorage'
    print *, 'PS % nValues = ', PS % nValues
    print *, 'PS % nVariables = ', PS % nVariables  

    print *, 'PS % iaUnpacked ( 1 ) = ', PS % iaUnpacked ( 1 )
    do i = 2, PS % nValues
      print *, &
        'PS % iaUnpacked (', i, ') = ', PS % iaUnpacked ( i )
    end do

  end subroutine PrintPackedStorage


end program PackedStorage_Form_Test


