program Reconstruction_C__Form_Test

  !-- Reconstruction_Chart__Form_Test

  use Basics
  use Manifolds
  use Fields
  use Slopes

  implicit none

  integer ( KDI ) :: &
    iD
  character ( 1 ), dimension ( 3 ) :: &
    D  =  [ 'X', 'Y', 'Z' ]
  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Chart_GS_Form ), allocatable :: &
    C
  type ( FieldSet_C_Form ), allocatable :: &
    FSC
  type ( FieldSet_C_Form ), dimension ( : ), allocatable :: &
    FSC_I, &
     DC_IL,  DC_IR
  type ( Stream_C_Form ), allocatable :: &
    SC
  type ( Geometry_F_C_Form ), allocatable :: &
    GC
  type ( Reconstruction_C_Form ), allocatable :: &
    RC_0, RC_1, RC_2

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Reconstruction_C__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( C )
  call C % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator, &
           PeriodicOption = [ .true., .true., .true. ] )

  allocate ( FSC )
  call FSC % Initialize ( C )

  associate ( nD  =>  C % nDimensions )

  allocate ( FSC_I ( nD ) )
  allocate (  DC_IL ( nD ),  DC_IR ( nD ) )
  do iD = 1, nD
    call FSC_I ( iD ) % Initialize &
           ( C, NameOption = 'Fields_I_' // D ( iD ) )
    call DC_IL ( iD ) % Initialize &
          ( C, NameOption = 'Difference_IL_' // D ( iD ) )
    call DC_IR ( iD ) % Initialize &
          ( C, NameOption = 'Difference_IR_' // D ( iD ) )
  end do !-- iD

  allocate ( GC )
  call GC % Initialize ( C )

  allocate ( RC_0 )
  allocate ( RC_1 )
  allocate ( RC_2 )
  call RC_0 % Initialize &
         ( GC, FSC, &
           NameOption = 'Reconstruction_0', &
           StreamedOption = .true., &
           OrderOption = 0 )
  call RC_1 % Initialize &
         ( GC, FSC, &
           NameOption = 'Reconstruction_1', &
           StreamedOption = .true., &
           OrderOption = 1 )
  call RC_2 % Initialize &
         ( GC, FSC, &
           NameOption = 'Reconstruction_2', &
           StreamedOption = .true., &
           OrderOption = 2 )

  allocate ( SC )
  call SC % Initialize ( C, GIS )
  call SC % AddFieldSet ( FSC )
  do iD = 1, nD
    call SC % AddFieldSet ( FSC_I  ( iD ) )
    call SC % AddFieldSet (  DC_IL ( iD ) )
    call SC % AddFieldSet (  DC_IR ( iD ) )
  end do !-- iD

  call   C   % Show ( )
  call FSC   % Show ( )
  call  GC   % Show ( )
  call  RC_0 % Show ( )
  call  RC_1 % Show ( )
  call  RC_2 % Show ( )
  call  SC   % Show ( )

  call SetWave ( FSC, GC )
  do iD = 1, nD
    call SetReference ( FSC_I ( iD ), GC, iD )
  end do !-- iD

  call TestReconstruction ( RC_0, SC, DC_IL, DC_IR, FSC_I )
  call TestReconstruction ( RC_1, SC, DC_IL, DC_IR, FSC_I )
  call TestReconstruction ( RC_2, SC, DC_IL, DC_IR, FSC_I )

  end associate !-- nD

  deallocate ( RC_2 )
  deallocate ( RC_1 )
  deallocate ( RC_0 )
  deallocate ( GC )
  deallocate ( DC_IR, DC_IL )
  deallocate ( FSC_I )
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


  subroutine SetReference ( FSC_I, GC, iD )

    class ( FieldSet_C_Form ), intent ( inout ) :: &
      FSC_I
    class ( Geometry_F_C_Form ), intent ( in ) :: &
      GC
    integer ( KDI ), intent ( in ) :: &
      iD

    associate &
      ( C  =>  GC % Storage_FSC % Storage % Value ( :, GC % CENTER_U ( iD ) ), &
        W  =>  GC % Storage_FSC % Storage % Value ( :, GC % WIDTH_U ( iD ) ) )

    select case ( iD )
    case ( 1 )
      call SetWave ( FSC_I, GC, X_Option  =  C  -  0.5 * W )
    case ( 2 )
      call SetWave ( FSC_I, GC, Y_Option  =  C  -  0.5 * W )
    case ( 3 )
      call SetWave ( FSC_I, GC, Z_Option  =  C  -  0.5 * W )
    end select !-- iD

    end associate !-- C, etc.

  end subroutine SetReference


  subroutine TestReconstruction ( RC, SC, DC_IL, DC_IR, FSC_I )

    class ( Reconstruction_C_Form ), intent ( inout ) :: &
      RC
    type ( Stream_C_Form ), intent ( inout ) :: &
      SC
    type ( FieldSet_C_Form ), dimension ( : ), intent ( inout ) :: &
      DC_IL, DC_IR
    type ( FieldSet_C_Form ), dimension ( : ), intent ( in ) :: &
      FSC_I

    associate ( nD  =>  RC % FieldSet_C % Chart % nDimensions )

    do iD  =  1, nD

      call RC % Compute ( iD )
      call CompareFieldSets ( RC % Output_IL_C, FSC_I ( iD ), iD )
      call CompareFieldSets ( RC % Output_IR_C, FSC_I ( iD ), iD )

      associate & 
        ( FV_I   =>  FSC_I ( iD ) % Storage_FSC % Storage % Value, &
          OV_IL  =>  RC % Output_IL_C % Storage_FSC % Storage % Value, &
          OV_IR  =>  RC % Output_IR_C % Storage_FSC % Storage % Value, &
          DV_IL  =>  DC_IL ( iD ) % Storage_FSC % Storage % Value, &
          DV_IR  =>  DC_IR ( iD ) % Storage_FSC % Storage % Value )

      DV_IL  =  OV_IL  -  FV_I
      DV_IR  =  OV_IR  -  FV_I

      end associate !-- OV_IL, etc.

    end do !-- iD

    call GIS % Open ( GIS % ACCESS_CREATE )
    call SC % Write ( )
    call GIS % Close ( )

    end associate !-- nD

  end subroutine TestReconstruction


  subroutine CompareFieldSets ( FSC, FSC_R, iD )

    class ( FieldSet_C_Form ), intent ( in ) :: &
      FSC, &
      FSC_R
    integer ( KDI ), intent ( in ) :: &
      iD

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iF     !-- iField
    type ( CollectiveOperation_R_Form ) :: &
      CO 

    call Show ( 'Comparing FieldSet with reference' )
    call Show ( FSC % Name, 'FieldSet' )
    call Show ( iD, 'iDimension' )

    associate ( nF  =>  FSC % nFields )

    select type ( C  =>  FSC % Chart )
    class is ( Chart_GS_Form )

    associate ( Cmm  =>  C % Communicator )
    call CO % Initialize &
           ( Cmm, nOutgoing = [ 2 * nF ], nIncoming = [ 2 * nF ] )
    end associate !-- Cmm

    do iS  =  1, nF
      iF  =  FSC % iaSelected ( iS )
      associate &
        ( F_R  =>  FSC_R % Storage_FSC % Storage % Value ( :, iF ), &
          F    =>  FSC   % Storage_FSC % Storage % Value ( :, iF ) )

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
    end select !-- G

    call CO % Reduce ( REDUCTION % SUM )
    call Show (    CO % Incoming % Value (      1 :      nF )  &
                /  CO % Incoming % Value ( nF + 1 : nF + nF ), 'L1 Error' )

    end associate !-- nF

  end subroutine CompareFieldSets


end program Reconstruction_C__Form_Test
