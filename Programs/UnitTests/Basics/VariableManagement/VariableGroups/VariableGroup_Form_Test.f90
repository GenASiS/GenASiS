program Storage_Form_Test

  use Specifiers
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
  type ( StorageForm ), dimension ( 3 ) :: &
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
  call S ( 2 ) % Initialize ( S ( 1 ), NameOption = 'Storage_2' )
  call PrintStorage ( S ( 2 ) )

  !-- InitializeClone, take 2
  call S ( 3 ) % Initialize &
         ( S ( 1 ), NameOption = 'Storage_3', &
           iaSelectedOption = [ 2, 3, 6 ] )
  call PrintStorage ( S ( 3 ) )
  
contains


  subroutine PrintStorage ( S )

    use Specifiers
    use Storage_Form

    integer ( KDI ) :: &
      i
    type ( StorageForm ) :: &
      S

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


end program Storage_Form_Test

