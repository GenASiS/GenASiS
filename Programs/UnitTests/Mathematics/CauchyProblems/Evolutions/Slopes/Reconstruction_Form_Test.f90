program Reconstruction_Form_Test

  use Basics
  use Manifolds
  use Fields
  use Slopes

  implicit none

  integer ( KDI ) :: &
    iD, &
    nCompute
  character ( 1 ), dimension ( 3 ) :: &
    D  =  [ 'X', 'Y', 'Z' ]
  logical ( KDL ) :: &
    DeviceMemory, &
    PinnedMemory, &
    DevicesCommunicate
  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( TimerForm ), pointer :: &
    T
  type ( Atlas_SCG_Form ), allocatable :: &
    A
  type ( FieldSetForm ), allocatable :: &
    FS
  type ( FieldSetForm ), dimension ( : ), allocatable :: &
    FS_I, &
     D_IL,  D_IR
  type ( StreamForm ), allocatable :: &
    S
  type ( Geometry_F_Form ), allocatable :: &
    G
  type ( ReconstructionForm ), allocatable :: &
    R_0, R_1, R_2

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Reconstruction_Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( A )
  call A % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator )

  DeviceMemory  =  OffloadEnabled ( )  .and.  NumberOfDevices ( ) >= 1 
  call PROGRAM_HEADER % GetParameter ( DeviceMemory, 'DeviceMemory' )

  PinnedMemory        =  DeviceMemory
  DevicesCommunicate  =  DeviceMemory
  call PROGRAM_HEADER % GetParameter &
         ( PinnedMemory, 'PinnedMemory' )
  call PROGRAM_HEADER % GetParameter &
         ( DevicesCommunicate, 'DevicesCommunicate' )

  allocate ( FS )
  call FS % Initialize &
         ( A, &
           DeviceMemoryOption = DeviceMemory, &
           PinnedMemoryOption = PinnedMemory, &
           DevicesCommunicateOption = DevicesCommunicate )
  do iD  =  1, 3
    call FS % SetBoundaryConditionsFace &
           ( [ 'PERIODIC', 'PERIODIC' ], iC = 1, iD = iD )
  end do !-- iD

  associate ( nD  =>  A % Chart_GS % nDimensions )

  allocate ( FS_I ( nD ) )
  allocate ( D_IL ( nD ), D_IR ( nD ) )
  do iD = 1, nD
    call FS_I ( iD ) % Initialize &
           ( A, NameOption = 'Fields_I_' // D ( iD ) )
    call D_IL ( iD ) % Initialize &
          ( A, NameOption = 'Difference_IL_' // D ( iD ) )
    call D_IR ( iD ) % Initialize &
          ( A, NameOption = 'Difference_IR_' // D ( iD ) )
  end do !-- iD

  allocate ( S )
  call S % Initialize ( A, GIS )
  call S % AddFieldSet ( FS )
  do iD = 1, nD
    call S % AddFieldSet ( FS_I  ( iD ) )
    call S % AddFieldSet (  D_IL ( iD ) )
    call S % AddFieldSet (  D_IR ( iD ) )
  end do !-- iD

  allocate ( G )
  call G % Initialize &
         ( A, &
           DeviceMemoryOption = DeviceMemory, &
           PinnedMemoryOption = PinnedMemory, &
           DevicesCommunicateOption = DevicesCommunicate )
  call G % SetStream ( S )

  allocate ( R_0 )
  allocate ( R_1 )
  allocate ( R_2 )
  call CONSOLE % SetVerbosity ( 'INFO_2' )
  call R_0 % Initialize &
         ( G, FS, SuffixOption = 'Rcnstrctn_0', OrderOption = 0 )
  call R_1 % Initialize &
         ( G, FS, SuffixOption = 'Rcnstrctn_1', OrderOption = 1 )
  call R_2 % Initialize &
         ( G, FS, SuffixOption = 'Rcnstrctn_2', OrderOption = 2 )
  call CONSOLE % SetVerbosity ( 'INFO_1' )
  !-- Get all Reconstruction timers initialized
  T  =>  R_0 % Timer ( Level = 1 )
  T  =>  R_1 % Timer ( Level = 1 )
  T  =>  R_2 % Timer ( Level = 1 )

  call  A   % Show ( )
  call FS   % Show ( )
  call  G   % Show ( )
  call CONSOLE % SetVerbosity ( 'INFO_2' )
  call  R_0 % Show ( )
  call  R_1 % Show ( )
  call  R_2 % Show ( )
  call CONSOLE % SetVerbosity ( 'INFO_1' )
  call  S   % Show ( )

  call SetWave ( FS, G )
  do iD = 1, nD
    call SetReference ( FS_I ( iD ), G, iD )
  end do !-- iD
  
  call FS % UpdateDevice ( )
  call G  % UpdateDevice ( )

  nCompute  =  1000
  call PROGRAM_HEADER % GetParameter ( nCompute, 'nCompute' )

  call TestReconstruction ( R_0, S, D_IL, D_IR, FS_I )
  call TestReconstruction ( R_1, S, D_IL, D_IR, FS_I )
  call TestReconstruction ( R_2, S, D_IL, D_IR, FS_I )

  end associate !-- nD

  call CONSOLE % SetVerbosity ( 'INFO_2' )
  deallocate ( R_2 )
  deallocate ( R_1 )
  deallocate ( R_0 )
  call CONSOLE % SetVerbosity ( 'INFO_1' )
  deallocate ( G )
  deallocate ( S )
  deallocate ( D_IR, D_IL )
  deallocate ( FS_I )
  deallocate ( FS )
  deallocate ( A )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )


contains


  subroutine SetWave ( FS, G, X_Option, Y_Option, Z_Option )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS
    class ( Geometry_F_Form ), intent ( in ), target :: &
      G
    real ( KDR ), dimension ( : ), intent ( in ), target, optional :: &
      X_Option, Y_Option, Z_Option

    integer ( KDI ), dimension ( 3 ) :: &
      nWavelengths
    real ( KDR ) :: &
      Offset, &
      Amplitude
    real ( KDR ), dimension ( 3 ) :: &
      Wavenumber
    real ( KDR ), dimension ( : ), pointer :: &
      X, Y, Z

    select type ( A  =>  FS % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS )
    associate &
      (  GV  =>   G % Storage_GS % Value, &
        FSV  =>  FS % Storage_GS % Value )

    nWavelengths  =  0
    nWavelengths ( 1 : C % nDimensions )  =  1
    call PROGRAM_HEADER % GetParameter ( nWavelengths, 'nWavelengths' )

    Offset     =  2.0_KDR
    Amplitude  =  1.0_KDR
    call PROGRAM_HEADER % GetParameter ( Offset, 'Offset' )
    call PROGRAM_HEADER % GetParameter ( Amplitude, 'Amplitude' )

    associate ( BoxSize => C % MaxCoordinate - C % MinCoordinate )
    where ( BoxSize > 0.0_KDR )
      Wavenumber = nWavelengths / BoxSize
    elsewhere
      Wavenumber = 0.0_KDR
    end where
    end associate !-- BoxSize

    if ( present ( X_Option ) ) then
      X  =>  X_Option
    else
      X  =>  GV ( :, G % CENTER_U_1 )
    end if

    if ( present ( Y_Option ) ) then
      Y  =>  Y_Option
    else
      Y  =>  GV ( :, G % CENTER_U_2 )
    end if

    if ( present ( Z_Option ) ) then
      Z  =>  Z_Option
    else
      Z  =>  GV ( :, G % CENTER_U_3 )
    end if

    associate &
      (     N  =>  FSV ( :, 1 ), &
            K  =>  Wavenumber, &
        TwoPi  =>  2.0_KDR  *  CONSTANT % PI )

    N  =  Offset  &
          +  Amplitude  &
             *  sin ( TwoPi * ( K ( 1 ) * X  +  K ( 2 ) * Y  +  K ( 3 ) * Z  ) )

    end associate !-- N, etc.
    end associate !-- GV, etc.
    end associate !-- C
    end select !-- A
    
  end subroutine SetWave


  subroutine SetReference ( FS_I, G, iD )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS_I
    class ( Geometry_F_Form ), intent ( in ) :: &
      G
    integer ( KDI ), intent ( in ) :: &
      iD

    associate ( GV  =>  G % Storage_GS % Value )

    associate &
      ( C  =>  GV ( :, G % CENTER_U ( iD ) ), &
        W  =>  GV ( :, G % WIDTH_U  ( iD ) ) )

    select case ( iD )
    case ( 1 )
      call SetWave ( FS_I, G, X_Option  =  C  -  0.5 * W )
    case ( 2 )
      call SetWave ( FS_I, G, Y_Option  =  C  -  0.5 * W )
    case ( 3 )
      call SetWave ( FS_I, G, Z_Option  =  C  -  0.5 * W )
    end select !-- iD

    end associate !-- C, etc.
    end associate !-- GV

  end subroutine SetReference


  subroutine TestReconstruction ( R, S, D_IL, D_IR, FS_I )

    class ( ReconstructionForm ), intent ( inout ) :: &
      R
    type ( StreamForm ), intent ( inout ) :: &
      S
    type ( FieldSetForm ), dimension ( : ), intent ( inout ) :: &
      D_IL, D_IR
    type ( FieldSetForm ), dimension ( : ), intent ( in ) :: &
      FS_I

    integer ( KDI ) :: &
      iD, &  !-- iDimension
      iC     !-- iCompute
    type ( TimerForm ), pointer :: &
      T

    select type ( A  =>  R % FieldSet % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( nD  =>  A % Chart_GS % nDimensions )

    do iD  =  1, nD

      call Show ( 'Computing reconstruction' )
      call Show ( R % Name, 'Reconstruction' )
      call Show ( iD, 'iDimension' )
      call Show ( nCompute, 'nCompute' )

      T  =>  R % Timer ( Level = 1 )
      call T % Start ( )
      do iC  =  1,  nCompute
        call R % Compute ( iC = 1, iD = iD )
      end do 
      call T % Stop ( )
      
      call R % Output_IL % UpdateHost ( )
      call R % Output_IR % UpdateHost ( )

      call CompareFieldSets ( R % Output_IL, FS_I ( iD ), iD )
      call CompareFieldSets ( R % Output_IR, FS_I ( iD ), iD )

      associate & 
        ( FV_I   =>  FS_I ( iD ) % Storage_GS % Value, &
          OV_IL  =>  R % Output_IL % Storage_GS % Value, &
          OV_IR  =>  R % Output_IR % Storage_GS % Value, &
          DV_IL  =>  D_IL ( iD ) % Storage_GS % Value, &
          DV_IR  =>  D_IR ( iD ) % Storage_GS % Value )

      DV_IL  =  OV_IL  -  FV_I
      DV_IR  =  OV_IR  -  FV_I

      end associate !-- OV_IL, etc.

    end do !-- iD

    T  =>  S % TimerWrite ( Level = 1 )
    call T % Start ( )
    call GIS % Open ( GIS % ACCESS_CREATE )
    call S % Write ( )
    call GIS % Close ( )
    call T % Stop ( )

    end associate !-- nD
    end select !-- A

  end subroutine TestReconstruction


  subroutine CompareFieldSets ( FS, FS_R, iD )

    class ( FieldSetForm ), intent ( in ) :: &
      FS, &
      FS_R
    integer ( KDI ), intent ( in ) :: &
      iD

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iF     !-- iField
    type ( CollectiveOperation_R_Form ) :: &
      CO 

    call Show ( 'Comparing FieldSet with reference' )
    call Show ( FS % Name, 'FieldSet' )
    call Show ( iD, 'iDimension' )

    select type ( A  =>  FS % Atlas )
      class is ( Atlas_SCG_Form )
    associate ( C  =>  A % Chart_GS )

    associate ( nF  =>  FS % nFields )

    associate ( Cr  =>  C % Communicator )
    call CO % Initialize &
           ( Cr, nOutgoing = [ 2 * nF ], nIncoming = [ 2 * nF ] )
    end associate !-- Cr

    do iS  =  1, nF
      iF  =  FS % iaSelected ( iS )
      associate &
        ( F_R  =>  FS_R % Storage_GS % Value ( :, iF ), &
          F    =>  FS   % Storage_GS % Value ( :, iF ) )

      !-- proper cells only
      CO % Outgoing % Value ( iS )  &
        =  sum ( pack ( abs ( F  -  F_R ), mask = C % ProperCell ) )
      CO % Outgoing % Value ( nF + iS )  &
        =  sum ( pack ( abs ( F_R ), mask = C % ProperCell ) )

!      !-- with ghost cells
!      CO % Outgoing % Value ( iF )  &
!        =  sum ( abs ( F  -  F_R ) )

      end associate !-- F_R, etc.
    end do !-- iF

    call CO % Reduce ( REDUCTION % SUM )
    call Show (    CO % Incoming % Value (      1 :      nF )  &
                /  CO % Incoming % Value ( nF + 1 : nF + nF ), 'L1 Error' )

    end associate !-- nF
    end associate !-- C
    end select !-- A

  end subroutine CompareFieldSets


end program Reconstruction_Form_Test
