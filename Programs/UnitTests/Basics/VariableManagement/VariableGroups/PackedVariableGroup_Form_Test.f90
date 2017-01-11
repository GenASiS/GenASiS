program PackedVariableGroup_Form_Test

  use Specifiers
  use ArrayArrays
  use VariableGroup_Form
  use PackedVariableGroup_Form

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
  type ( VariableGroupForm ), dimension ( 2 ) :: &
    VG
  type ( PackedVariableGroupForm ) :: &
    PVG

  do i = 1, size ( VariableUnit ) 
    call VariableUnit ( i ) % Initialize &
           ( UnitString ( i ), UnitString ( i ), 1.0_KDR ) 
  end do

  call VectorIndices ( 1 ) % Initialize ( [ 2, 3 ] )

  call VG ( 1 ) % Initialize &
         ( ValueShape = [ 4, 6 ], VectorIndicesOption = VectorIndices, &
           UnitOption = VariableUnit, VectorOption = [ 'Vector_1' ], &
           VariableOption = Variable, NameOption = 'VariableGroup_1')
  call PrintVariableGroup ( VG ( 1 ) )

  call VG ( 2 ) % Initialize &
         ( VG ( 1 ), NameOption = 'VariableGroup_2', &
           iaSelectedOption = [ 2, 3, 6 ] )
  call PrintVariableGroup ( VG ( 2 ) )

  UnpackedIndex = [ 1, 3 ]
  call PVG % Initialize &
         ( UnpackedIndex, VG ( 2 ) % nValues, VG ( 2 ) % nVariables )
  call PrintPackedVariableGroup ( PVG )

contains


  subroutine PrintVariableGroup ( VG )

    use Specifiers
    use VariableGroup_Form

    type ( VariableGroupForm ) :: &
      VG

    integer ( KDI ) :: &
      i

    print *
    print *, 'VG % nValues = ', VG % nValues
    print *, 'VG % nVariables = ', VG % nVariables
    print *, 'VG % nVectors = ', VG % nVectors
    print *, 'VG % lName = ', VG % lName
    print *, 'VG % lVariable = ', VG % lVariable
    print *, 'VG % lVector = ', VG % lVector
    print *, 'VG % iaSelected = ', VG % iaSelected
    print *, 'VG % AllocatedValue = ', VG % AllocatedValue
    print *, 'VG % Name = ', trim ( VG % Name )

    do i = 1, VG % nVariables
      print *, &
        'VG % Variable (', i, ') = ', &
        trim ( VG % Variable ( VG % iaSelected ( i ) ) )
    end do

    do i = 1, VG % nVectors
      print *, &
        'VG % Vector (', i, ') = ', trim ( VG % Vector ( i ) )
    end do

    do i = 1, VG % nVariables
      print *, &
        'VG % Unit (', i, ') % Unit = ', &
        trim ( VG % Unit ( VG % iaSelected ( i ) ) % Unit )
    end do

    do i = 1, VG % nVectors
      print *, &
        'VG % VectorIndices (', i, ') = ', VG % VectorIndices ( i ) % Value
    end do

  end subroutine PrintVariableGroup


  subroutine PrintPackedVariableGroup ( PVG )

    use Specifiers
    use PackedVariableGroup_Form

    type ( PackedVariableGroupForm ) :: &
      PVG

    integer ( KDI ) :: &
      i

    print *
    print *, 'PackedVariableGroup'
    print *, 'PVG % nValues = ', PVG % nValues
    print *, 'PVG % nVariables = ', PVG % nVariables  

    print *, 'PVG % iaUnpacked ( 1 ) = ', PVG % iaUnpacked ( 1 )
    do i = 2, PVG % nValues
      print *, &
        'PVG % iaUnpacked (', i, ') = ', PVG % iaUnpacked ( i )
    end do

  end subroutine PrintPackedVariableGroup


end program PackedVariableGroup_Form_Test


