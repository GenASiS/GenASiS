program DivergenceContribution_CS__Form_Test

  use Basics
  use Manifolds
  use FieldSets
  use Geometries
  use CurrentSets

  implicit none

  integer ( KDI ) :: &
    nCompute
  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Atlas_SCG_Form ), allocatable :: &
    A
  type ( StreamForm ), allocatable :: &
    S
  type ( FieldSetForm ), allocatable :: &
    FS_F  !-- FieldSet_Fluxes
  type ( Geometry_F_Form ), allocatable :: &
    G
  type ( CurrentSetForm ), allocatable :: &
    CS

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'DivergenceContribution_CS__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( A )
  call A % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( S )
  call S % Initialize ( A, GIS )

  allocate ( G )
  call G % Initialize ( A )
  call G % SetStream ( S )

  allocate ( CS )
  call CS % Initialize ( G )
  call CS % SetStream ( S )

  CS % nDivergenceContributions  =  1
  allocate ( CS % DivergenceContribution ( 1 ) )
  allocate ( DivergenceContribution_CS_Form &
             :: CS % DivergenceContribution ( 1 ) % Element )
  select type ( DC  =>  CS % DivergenceContribution ( 1 ) % Element )
    class is ( DivergenceContribution_CS_Form )
  call DC % Initialize ( CS )
  end select !-- FS

  allocate ( FS_F )
  call FS_F % Initialize &
         ( A, &
           FieldOption = CS % Balanced, &
           NameOption = 'Fluxes', &
           nFieldsOption = CS % nBalanced )
  call S % AddFieldSet ( FS_F )

  call    A % Show ( )
  call    G % Show ( )
  call   CS % Show ( )
  call FS_F % Show ( )
  call    S % Show ( )

  call SetWave ( CS, G )

  nCompute  =  1000
  call PROGRAM_HEADER % GetParameter ( nCompute, 'nCompute' )

  call TestFluxes ( CS, FS_F, S, iD = 1 )
  call TestFluxes ( CS, FS_F, S, iD = 2 )
  call TestFluxes ( CS, FS_F, S, iD = 3 )

  deallocate ( FS_F )
  deallocate ( CS )
  deallocate ( G )
  deallocate ( S )
  deallocate ( A )
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


  subroutine TestFluxes ( CS, FS_F, S, iD )

    class ( CurrentSetForm ), intent ( inout ) :: &
      CS
    class ( FieldSetForm ), intent ( inout ) :: &
      FS_F
    class ( StreamForm ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iD

    integer ( KDI ) :: &
      iC  !-- iCompute
    type ( TimerForm ), pointer :: &
      T

    associate ( DC  =>  CS % DivergenceContribution ( 1 ) % Element )

    call Show ( 'DivergenceContribution_CS computation' )
    call Show ( DC % Name, 'DivergenceContribution_CS' )
    call Show ( iD, 'iDimension' )
    call Show ( nCompute, 'nCompute' )

    T  =>  DC % Timer ( LevelOption = 1 )
    call T % Start ( )
    do iC  =  1,  nCompute
      call DC % ComputeFluxes ( FS_F, CS, iC = 1, iD = iD )
    end do
    call T % Stop ( )

    T  =>  S % TimerWrite ( LevelOption = 1 )
    call T % Start ( )
    call GIS % Open ( GIS % ACCESS_CREATE )
    call S % Write ( )
    call GIS % Close ( )
    call T % Stop ( )

    end associate !-- DC

  end subroutine TestFluxes


end program DivergenceContribution_CS__Form_Test
