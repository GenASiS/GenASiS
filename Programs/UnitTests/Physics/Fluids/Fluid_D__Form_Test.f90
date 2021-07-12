program Fluid_D__Form_Test

  !-- Fluid_Dust__Form_Test

  use Basics
  use Mathematics
  use Gravitations
  use Fluids

  implicit none

  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Atlas_SCG_Form ), allocatable :: &
    A
  type ( StreamForm ), allocatable :: &
    S
  type ( Gravitation_G_Form ), allocatable :: &
    G
  type ( Units_F_Form ), dimension ( : ), allocatable :: &
    U
  type ( Fluid_D_Form ), allocatable :: &
    F

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Fluid_D__Form_Test', DimensionalityOption = '2D' )

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

  allocate ( U ( 1 ) )
  call U ( 1 ) % Initialize ( TypeOption = 'MKS' )

  allocate ( F )
  call F % Initialize ( G, U )
  call F % SetStream ( S )

  call A       % Show ( )
  call U ( 1 ) % Show ( )
  call F       % Show ( )
  call S       % Show ( )

!   call SetWave ( CSC, G )

  call GIS % Open ( GIS % ACCESS_CREATE )
  call S % Write ( )
  call GIS % Close ( )

  deallocate ( F )
  deallocate ( U )
  deallocate ( G )
  deallocate ( S )
  deallocate ( A )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )


! contains


!   subroutine SetWave ( CSC, GC )

!     class ( Fluid_D_Form ), intent ( inout ) :: &
!       CSC
!     class ( Geometry_F_C_Form ), intent ( in ), target :: &
!       GC

!     integer ( KDI ), dimension ( 3 ) :: &
!       nWavelengths
!     real ( KDR ) :: &
!       Offset, &
!       Amplitude, &
!       Speed
!     real ( KDR ), dimension ( 3 ) :: &
!       Wavenumber

!     select type ( C  =>  CSC % Chart )
!     class is ( Chart_GS_Form )

!     nWavelengths  =  0
!     nWavelengths ( 1 : C % nDimensions )  =  1
!     call PROGRAM_HEADER % GetParameter ( nWavelengths, 'nWavelengths' )

!     Offset     =  2.0_KDR
!     Amplitude  =  1.0_KDR
!     Speed      =  1.0_KDR
!     call PROGRAM_HEADER % GetParameter ( Offset, 'Offset' )
!     call PROGRAM_HEADER % GetParameter ( Amplitude, 'Amplitude' )
!     call PROGRAM_HEADER % GetParameter ( Speed, 'Speed' )

!     associate ( BoxSize  =>  C % MaxCoordinate  -  C % MinCoordinate )
!     where ( BoxSize  >  0.0_KDR )
!       Wavenumber  =  nWavelengths / BoxSize
!     elsewhere
!       Wavenumber  =  0.0_KDR
!     end where
!     end associate !-- BoxSize

!     associate &
!       (     X  =>   GC % Storage_FSC % Storage &
!                        % Value ( :, GC % CENTER_U_1 ), &
!             Y  =>   GC % Storage_FSC % Storage &
!                        % Value ( :, GC % CENTER_U_2 ), &
!             Z  =>   GC % Storage_FSC % Storage &
!                        % Value ( :, GC % CENTER_U_3 ), &
!           Rho  =>  CSC % Storage_FSC % Storage &
!                        % Value ( :, CSC % DENSITY_DEFAULT ), &
!             V  =>  CSC % VelocityDefault_U, &
!             K  =>  Wavenumber, &
!         Abs_K  =>  sqrt ( dot_product ( Wavenumber, Wavenumber ) ), &
!         TwoPi  =>  2.0_KDR  *  CONSTANT % PI )

!     Rho  =  Offset  &
!             +  Amplitude  &
!                *  sin ( TwoPi * (    K ( 1 ) * X  &
!                                   +  K ( 2 ) * Y  &
!                                   +  K ( 3 ) * Z  ) )

!     V ( 1 )  =  Speed  *  K ( 1 )  /  Abs_K
!     V ( 2 )  =  Speed  *  K ( 2 )  /  Abs_K
!     V ( 3 )  =  Speed  *  K ( 3 )  /  Abs_K
    
!     end associate !-- Rho, etc.
!     end select !-- C

!   end subroutine SetWave


end program Fluid_D__Form_Test
