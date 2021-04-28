program Reconstruction_C__Form_Test

  !-- Reconstruction_Chart__Form_Test

  use Basics
  use Manifolds
  use Fields
  use Slopes

  implicit none

  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Chart_GS_Form ), allocatable :: &
    C
  type ( FieldSet_C_Form ), allocatable :: &
    FSC
  type ( Stream_C_Form ), allocatable :: &
    SC
  type ( Geometry_F_C_Form ), allocatable :: &
    GC
  type ( Reconstruction_C_Form ), allocatable :: &
    RC_0, &
    RC_1, &
    RC_2

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Reconstruction_C__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( C )
  call C % Initialize &
         ( PeriodicOption = [ .true., .true., .true. ] )

  allocate ( FSC )
  call FSC % Initialize ( C )

  allocate ( GC )
  call GC % Initialize ( C )

  allocate ( RC_0 )
  call RC_0 % Initialize &
         ( GC, FSC, &
           NameOption = 'Reconstruction_0', &
           StreamedOption = .true., &
           OrderOption = 0 )

  allocate ( RC_1 )
  call RC_1 % Initialize &
         ( GC, FSC, &
           NameOption = 'Reconstruction_1', &
           StreamedOption = .true., &
           OrderOption = 1 )

  allocate ( RC_2 )
  call RC_2 % Initialize &
         ( GC, FSC, &
           NameOption = 'Reconstruction_2', &
           StreamedOption = .true., &
           OrderOption = 2 )

  allocate ( SC )
  call SC % Initialize ( C, GIS )
  call SC % AddFieldSet ( FSC )

  call   C   % Show ( )
  call FSC   % Show ( )
  call  GC   % Show ( )
  call  RC_0 % Show ( )
  call  RC_1 % Show ( )
  call  RC_2 % Show ( )
  call  SC   % Show ( )

  call SetWave ( FSC, GC )

  call GIS % Open ( GIS % ACCESS_CREATE )
  call SC % Write ( )
  call GIS % Close ( )

  deallocate ( RC_2 )
  deallocate ( RC_1 )
  deallocate ( RC_0 )
  deallocate ( GC )
  deallocate ( FSC )
  deallocate ( C )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )


contains


  subroutine SetWave ( FSC, GC, X_Option, Y_Option, Z_Option )

    class ( FieldSet_C_Form ), intent ( inout ) :: &
      FSC
    class ( Geometry_F_C_Form ), intent ( in ), target :: &
      GC
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

    select type ( C  =>  FSC % Chart )
    class is ( Chart_GS_Form )

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
      X  =>  GC % Storage_FSC % Storage % Value ( :, GC % CENTER_U_1 )
    end if

    if ( present ( Y_Option ) ) then
      Y  =>  Y_Option
    else
      Y  =>  GC % Storage_FSC % Storage % Value ( :, GC % CENTER_U_2 )
    end if

    if ( present ( Z_Option ) ) then
      Z  =>  Z_Option
    else
      Z  =>  GC % Storage_FSC % Storage % Value ( :, GC % CENTER_U_3 )
    end if

    associate &
      (     N  =>  FSC % Storage_FSC % Storage % Value ( :, 1 ), &
            K  =>  Wavenumber, &
        TwoPi  =>  2.0_KDR  *  CONSTANT % PI )

    N  =  Offset  &
          +  Amplitude  &
             *  sin ( TwoPi * ( K ( 1 ) * X  +  K ( 2 ) * Y  +  K ( 3 ) * Z  ) )

    end associate !-- N, etc.
    end select !-- C

  end subroutine SetWave


end program Reconstruction_C__Form_Test
