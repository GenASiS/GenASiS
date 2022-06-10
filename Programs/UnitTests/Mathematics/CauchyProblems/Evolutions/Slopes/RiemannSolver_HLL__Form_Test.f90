program RiemannSolver_HLL__Form_Test

  !-- RiemannSolver_HartenLaxVanLeer__Form_Test

  use Basics
  use Manifolds
  use Fields
  use Slopes

  implicit none

  integer ( KDI ) :: &
    nCompute
  type ( GridImageStreamForm ), allocatable :: &
    GIS, &
    GIS_SD  !-- StageDimension
  type ( Atlas_SCG_Form ), allocatable :: &
    A
  type ( StreamForm ), allocatable :: &
    S, &
    S_SD  !-- StageDimension
  type ( Geometry_F_Form ), allocatable :: &
    G
  type ( CurrentSetForm ), allocatable :: &
    CS
  type ( DivergencePart_CS_Form ), allocatable :: &
    DP
  type ( RiemannSolver_HLL_Form ), allocatable :: &
    RS

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'RiemannSolver_HLL__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  allocate ( GIS_SD )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )
  call GIS_SD % Initialize &
         ( 'RiemannSolver_SD_Test', &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( A )
  call A % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( S )
  allocate ( S_SD )
  call S    % Initialize ( A, GIS )
  call S_SD % Initialize ( A, GIS_SD, NameOption = trim ( S % Name ) // '_SD' )

  allocate ( G )
  call G % Initialize ( A )
  call G % SetStream ( S )

  allocate ( CS )
  call CS % Initialize( G )
  call CS % SetStream ( S )

  allocate ( DP )
  call DP % Initialize ( CS )

  allocate ( RS )
  call CONSOLE % SetVerbosity ( 'INFO_2' )
  call RS % Initialize ( CS )
  call CONSOLE % SetVerbosity ( 'INFO_1' )
  call S % AddFieldSet ( RS )

  call  A % Show ( )
  call  G % Show ( )
  call CS % Show ( )
  call DP % Show ( )
  call CONSOLE % SetVerbosity ( 'INFO_2' )
  call RS % Show ( )
  call CONSOLE % SetVerbosity ( 'INFO_1' )
  call  S % Show ( )

  call SetWave ( CS, G )

  nCompute  =  1000
  call PROGRAM_HEADER % GetParameter ( nCompute, 'nCompute' )

  associate ( C  =>  A % Chart_GS )

  call TestRiemannSolver ( iD = 1 )
  if ( C % nDimensions  >  1 ) &
    call TestRiemannSolver ( iD = 2 )
  if ( C % nDimensions  >  2 ) &
    call TestRiemannSolver ( iD = 3 )

  call TestStageDimension ( C % nDimensions )

  end associate !-- C

  call CONSOLE % SetVerbosity ( 'INFO_2' )
  deallocate ( RS )
  call CONSOLE % SetVerbosity ( 'INFO_1' )
  deallocate ( DP )
  deallocate ( CS )
  deallocate ( G )
  deallocate ( S_SD )
  deallocate ( S )
  deallocate ( A )
  deallocate ( GIS_SD )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )


contains


  subroutine SetWave ( CS, G )

    class ( CurrentSetForm ), intent ( inout ) :: &
      CS
    class ( Geometry_F_Form ), intent ( in ), target :: &
      G

    integer ( KDI ), dimension ( 3 ) :: &
      nWavelengths
    real ( KDR ) :: &
      Offset, &
      Amplitude, &
      Speed
    real ( KDR ), dimension ( 3 ) :: &
      Wavenumber

    select type ( A  =>  CS % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS )

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
      (  GV  =>   G % Storage_GS % Value, &
        CSV  =>  CS % Storage_GS % Value )
    associate &
      (     X  =>   GV ( :,  G % CENTER_U_1 ), &
            Y  =>   GV ( :,  G % CENTER_U_2 ), &
            Z  =>   GV ( :,  G % CENTER_U_3 ), &
          Rho  =>  CSV ( :, CS % DENSITY_CS ), &
          V_1  =>  CSV ( :, CS % VELOCITY_CS_U_1 ), &
          V_2  =>  CSV ( :, CS % VELOCITY_CS_U_2 ), &
          V_3  =>  CSV ( :, CS % VELOCITY_CS_U_3 ), &
            K  =>  Wavenumber, &
        Abs_K  =>  sqrt ( dot_product ( Wavenumber, Wavenumber ) ), &
        TwoPi  =>  2.0_KDR  *  CONSTANT % PI )

    Rho  =  Offset  &
            +  Amplitude  &
               *  sin ( TwoPi * (    K ( 1 ) * X  &
                                  +  K ( 2 ) * Y  &
                                  +  K ( 3 ) * Z  ) )

    V_1  =  Speed  *  K ( 1 )  /  Abs_K
    V_2  =  Speed  *  K ( 2 )  /  Abs_K
    V_3  =  Speed  *  K ( 3 )  /  Abs_K
    
    end associate !-- Rho, etc.
    end associate !-- GV, etc.
    end associate !-- C
    end select !-- A

  end subroutine SetWave


  subroutine TestRiemannSolver ( iD )

    integer ( KDI ), intent ( in ) :: &
      iD

    integer ( KDI ) :: &
      iC  !-- iCompute
    type ( TimerForm ), pointer :: &
      T

    call Show ( 'RiemannSolver computation' )
    call Show ( RS % Name, 'RiemannSolver' )
    call Show ( iD, 'iDimension' )
    call Show ( nCompute, 'nCompute' )

    do iC  =  1,  nCompute

      T  =>  RS % Timer_P ( Level = 1 )
      call T % Start ( )
      call RS % Prepare ( iC = 1, iD = iD, T_Option = T )
      call T % Stop ( )

      T  =>  RS % Timer_C ( Level = 1 )
      call T % Start ( )
      call RS % Compute ( DP, iC = 1, iD = iD, T_Option = T )
      call T % Stop ( )

    end do

    T  =>  S % TimerWrite ( Level = 1 )
    call T % Start ( )
    call GIS % Open ( GIS % ACCESS_CREATE )
    call S % Write ( )
    call GIS % Close ( )
    call T % Stop ( )

  end subroutine TestRiemannSolver


  subroutine TestStageDimension ( nD )

    integer ( KDI ), intent ( in ) :: &
      nD

    call Show ( 'RiemannSolver computation with StageDimension' )
    call Show ( RS % Name, 'RiemannSolver' )

    call CONSOLE % SetVerbosity ( 'INFO_2' )
    call RS % SetStream ( S_SD, nS = 1 )  !-- nStages = 1
    call CONSOLE % SetVerbosity ( 'INFO_1' )

    call RS % Prepare ( iC = 1, iD = 1, iS_Option = 1 )
    call RS % Compute ( DP, iC = 1, iD = 1, iS_Option = 1 )
    if ( nD  >  1 ) then
      call RS % Prepare ( iC = 1, iD = 2, iS_Option = 1 )
      call RS % Compute ( DP, iC = 1, iD = 2, iS_Option = 1 )
    end if
    if ( nD  >  2 ) then
      call RS % Prepare ( iC = 1, iD = 3, iS_Option = 1 )
      call RS % Compute ( DP, iC = 1, iD = 3, iS_Option = 1 )
    end if

    call GIS_SD % Open ( GIS_SD % ACCESS_CREATE )
    call S_SD % Write ( )
    call GIS_SD % Close ( )

  end subroutine TestStageDimension


end program RiemannSolver_HLL__Form_Test
