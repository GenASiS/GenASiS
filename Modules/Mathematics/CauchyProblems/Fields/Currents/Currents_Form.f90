module Currents_Form

  use Basics
  use FieldSets
  use Geometries

  implicit none
  private

!     integer ( KDI ), private, parameter :: &
!       N_FIELDS_CS    = 0, &
!       N_VECTORS_CS   = 0, &
!       N_PRIMITIVE_CS = 0, &
!       N_BALANCED_CS  = 0

!   type, public, extends ( FieldSetForm ) :: CurrentsForm
!     !-- Fields and vectors
!     integer ( KDI ) :: &
!       N_FIELDS_CS    = N_FIELDS_CS, &
!       N_VECTORS_CS   = N_VECTORS_CS
!     !-- Defaults not generally used
!     integer ( KDI ) :: &
!       DENSITY_DEFAULT      = 0, &
!       VELOCITY_DEFAULT_U_1 = 0, &
!       VELOCITY_DEFAULT_U_2 = 0, &
!       VELOCITY_DEFAULT_U_3 = 0
!     integer ( KDI ), dimension ( 3 ) :: &
!       VELOCITY_DEFAULT_U = 0
!     !-- Primtive and Balanced
!     integer ( KDI ) :: &
!       nPrimitive = 0, &
!       nBalanced  = 0
!     integer ( KDI ) :: &
!       N_PRIMITIVE_CS = N_PRIMITIVE_CS, &
!       N_BALANCED_CS  = N_BALANCED_CS
!     integer ( KDI ), dimension ( : ), allocatable :: &
!       iaPrimitive, &
!       iaBalanced
!     character ( LDL ), dimension ( : ), allocatable :: &
!       Primitive, &
!       Balanced
!     class ( Geometry_F_orm ), pointer :: &
!       Geometry_C => null ( )
!   contains
!     procedure, private, pass :: &
!       InitializeAllocate_CS
!     generic, public :: &
!       Initialize => InitializeAllocate_CS
!     procedure, public, pass :: &
!       SetVelocityConstant
!     procedure, public, pass :: &
!       SetVelocityLinear
!     procedure, public, pass ( C ) :: &
!       SetStream
!     procedure, private, pass :: &
!       Show_FS
!     procedure, public, pass :: &
!       ComputeFromInitial
!     procedure, public, pass :: &
!       ComputeFromConserved
!     procedure, public, pass ( C ) :: &
!       ComputeFluxes
!     procedure, public, pass ( C ) :: &
!       ComputeEigenspeeds
!     final :: &
!       Finalize
!   end type CurrentsForm

!     private :: &
!       ComputeFluxesKernel, &
!       ComputeEigenspeedsKernel

!     interface

!       module subroutine ComputeFluxesKernel ( D, V_Dim, F_D, UseDeviceOption )
!         use Basics
!         implicit none
!         real ( KDR ), dimension ( : ), intent ( in ) :: &
!           D, &
!           V_Dim
!         real ( KDR ), dimension ( : ), intent ( out ) :: &
!           F_D
!         logical ( KDL ), intent ( in ), optional :: &
!           UseDeviceOption
!       end subroutine ComputeFluxesKernel

!       module subroutine ComputeEigenspeedsKernel &
!                ( V_Dim, EF_P, EF_M, UseDeviceOption )
!         use Basics
!         implicit none
!         real ( KDR ), dimension ( : ), intent ( in ) :: &
!           V_Dim
!         real ( KDR ), dimension ( : ), intent ( out ) :: &
!           EF_P, EF_M
!         logical ( KDL ), intent ( in ), optional :: &
!           UseDeviceOption
!       end subroutine ComputeEigenspeedsKernel

!     end interface

! contains


!   subroutine InitializeAllocate_CS &
!                ( C, GC, FieldOption, VectorOption, NameOption, UnitOption, &
!                  DensityUnitOption, VectorIndicesOption, iaPrimitiveOption, &
!                  iaBalancedOption, nFieldsOption, IgnorabilityOption )

!     class ( CurrentsForm ), intent ( inout ) :: &
!       C
!     class ( Geometry_F_orm ), intent ( in ), target :: &
!       GC
!     character ( * ), dimension ( : ), intent ( in ), optional :: &
!       FieldOption, &
!       VectorOption
!     character ( * ), intent ( in ), optional :: &
!       NameOption
!     type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
!       UnitOption
!     type ( MeasuredValueForm ), intent ( in ), optional :: &
!       DensityUnitOption
!     type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
!       VectorIndicesOption
!     integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
!       iaPrimitiveOption, &
!       iaBalancedOption
!     integer ( KDI ), intent ( in ), optional :: &
!       nFieldsOption, &
!       IgnorabilityOption

!     integer ( KDI ) :: &
!       iP, &  !-- iPrimitive
!       iB, &  !-- iBalanced
!       iF, &  !-- iField
!       iV, &  !-- iVector
!       nFields, &
!       nVectors
!     type ( Integer_1D_Form ), dimension ( : ), allocatable :: &
!       VectorIndices
!     type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
!       Unit
!     character ( LDL ) :: &
!       Name
!     character ( LDL ), dimension ( : ), allocatable :: &
!       Field, &
!       Vector

!     if ( C % Type  ==  '' ) &
!       C % Type  =  'a Currents_C' 
    
!     Name  =  'Currents'
!     if ( present ( NameOption ) ) &
!       Name  =  NameOption

!     C % Geometry_C  =>  GC

!     associate &
!       ( DeviceMemory  =>  GC % Storage_FSC % DeviceMemory, &
!         PinnedMemory  =>  GC % Storage_FSC % PinnedMemory, &
!         DevicesCommunicate  =>  GC % GhostExchange_FSC % DevicesCommunicate ) 

!     !-- Field indices

!     if ( present ( nFieldsOption ) ) then
!       nFields  =  nFieldsOption
!     else

!       C % DENSITY_DEFAULT       =  1
!       C % VELOCITY_DEFAULT_U_1  =  2
!       C % VELOCITY_DEFAULT_U_2  =  3
!       C % VELOCITY_DEFAULT_U_3  =  4

!       nFields  =  C % N_FIELDS_CS  +  4

!       C % VELOCITY_DEFAULT_U  =  [ C % VELOCITY_DEFAULT_U_1, &
!                                      C % VELOCITY_DEFAULT_U_2, &
!                                      C % VELOCITY_DEFAULT_U_3 ]

!     end if

!     !-- Field names

!     if ( present ( FieldOption ) ) then
!       allocate ( Field, source = FieldOption )
!     else
!       allocate ( Field ( nFields ) )
!       Field ( C % N_FIELDS_CS + 1 )  =  'Density'
!       Field ( C % N_FIELDS_CS + 2 )  =  'Velocity_U_1'
!       Field ( C % N_FIELDS_CS + 3 )  =  'Velocity_U_2'
!       Field ( C % N_FIELDS_CS + 4 )  =  'Velocity_U_3'
!     end if !-- FieldOption

!     !-- Units

!     if ( present ( UnitOption ) ) then
!       allocate ( Unit, source = UnitOption )
!     else
!       allocate ( Unit ( nFields ) )
!       if ( present ( DensityUnitOption ) ) &
!         Unit ( C % DENSITY_DEFAULT )  =  DensityUnitOption
!     end if !-- UnitOption

!     !-- Vector indices

!     if ( present ( VectorIndicesOption ) ) then
!       nVectors  =  size ( VectorIndicesOption )
!       allocate ( VectorIndices ( nVectors ) )
!       do iV  =  C % N_VECTORS_CS + 1,  nVectors 
!         call VectorIndices ( iV ) % Initialize ( VectorIndicesOption ( iV ) )
!       end do !-- iV
!     else
!       nVectors  =  C % N_VECTORS_CS + 1
!       allocate ( VectorIndices ( nVectors ) )
!       call VectorIndices ( C % N_VECTORS_CS + 1 ) % Initialize &
!              ( C % VELOCITY_DEFAULT_U )
!     end if

!     !-- Vector names

!     if ( present ( VectorOption ) ) then
!       allocate ( Vector, source = VectorOption )
!     else
!       allocate ( Vector ( nVectors ) )
!       Vector ( C % N_VECTORS_CS + 1 )  =  'Velocity'
!     end if !-- FieldOption

!     !-- Primitive fields

!     if ( present ( iaPrimitiveOption ) ) then
!       C % nPrimitive  =  size ( iaPrimitiveOption )
!       allocate ( C % iaPrimitive, source = iaPrimitiveOption )
!     else
!       C % nPrimitive  =  C % N_PRIMITIVE_CS + 1
!       allocate ( C % iaPrimitive ( C % nPrimitive ) )
!       C % iaPrimitive  =  [ C % DENSITY_DEFAULT ]
!     end if !-- iaPrimitiveOption

!     associate ( nP  =>  C % nPrimitive )
!     allocate ( C % Primitive ( nP ) )
!     do iP  =  1, nP
!       iF  =  C % iaPrimitive ( iP )
!       C % Primitive ( iP )  =  Field ( iF )
!     end do !-- iP
!     end associate !-- nP

!     !-- Balanced fields

!     if ( present ( iaBalancedOption ) ) then
!       C % nBalanced  =  size ( iaBalancedOption )
!       allocate ( C % iaBalanced, source = iaBalancedOption )
!     else
!       C % nBalanced  =  C % N_BALANCED_CS + 1
!       allocate ( C % iaBalanced ( C % nBalanced ) )
!       C % iaBalanced  =  [ C % DENSITY_DEFAULT ]
!     end if !-- iaBalancedOption

!     associate ( nB  =>  C % nBalanced )
!     allocate ( C % Balanced ( nB ) )
!     do iB  =  1, nB
!       iF  =  C % iaBalanced ( iB )
!       C % Balanced ( iB )  =  Field ( iF )
!     end do !-- iB
!     end associate !-- nP

!     !-- FieldSet

!     call C % FieldSetForm % Initialize &
!            ( GC % Chart, &
!              FieldOption = Field, &
!              VectorOption = Vector, &
!              NameOption = Name, &
!              DeviceMemoryOption = DeviceMemory, &
!              PinnedMemoryOption = PinnedMemory, &
!              DevicesCommunicateOption = DevicesCommunicate, &
!              UnitOption = Unit, &
!              VectorIndicesOption = VectorIndices, &
!              nFieldsOption = nFields, &
!              IgnorabilityOption = IgnorabilityOption )

!     end associate !-- DeviceMemory, etc.

!   end subroutine InitializeAllocate_CS


!   subroutine SetVelocityConstant ( C, Direction, Speed )

!     class ( CurrentsForm ), intent ( inout ) :: &
!       C
!     real ( KDR ), dimension ( : ), intent ( in ) :: &
!       Direction
!     real ( KDR ), intent ( in ) :: &
!       Speed

!     associate &
!       ( CSV  =>  C % Storage_FSC % Storage % Value )
!     associate &
!       (     V_1  =>  CSV ( :, C % VELOCITY_DEFAULT_U_1 ), &
!             V_2  =>  CSV ( :, C % VELOCITY_DEFAULT_U_2 ), &
!             V_3  =>  CSV ( :, C % VELOCITY_DEFAULT_U_3 ), &
!             K    =>  Direction, &
!         Abs_K    =>  sqrt ( dot_product ( Direction, Direction ) ) )

!     V_1  =  Speed  *  K ( 1 )  /  Abs_K
!     V_2  =  Speed  *  K ( 2 )  /  Abs_K
!     V_3  =  Speed  *  K ( 3 )  /  Abs_K
    
!     end associate !-- V_1, etc.
!     end associate !-- CSV

!   end subroutine SetVelocityConstant


!   subroutine SetVelocityLinear ( C, Speed, Length )

!     class ( CurrentsForm ), intent ( inout ) :: &
!       C
!     real ( KDR ), intent ( in ) :: &
!       Speed, &
!       Length

!     associate &
!       ( GC  =>  C % Geometry_C )
!     associate &
!       ( CSV  =>  C % Storage_FSC % Storage % Value, &
!          GV  =>   GC % Storage_FSC % Storage % Value )
!     associate &
!       (    V_1  =>  CSV ( :, C % VELOCITY_DEFAULT_U_1 ), &
!            V_2  =>  CSV ( :, C % VELOCITY_DEFAULT_U_2 ), &
!            V_3  =>  CSV ( :, C % VELOCITY_DEFAULT_U_3 ), &
!            X_1  =>   GV ( :,  GC % CENTER_U_1 ) )

!     V_1  =  Speed  *  ( X_1 / Length )
!     V_2  =  0.0_KDR
!     V_3  =  0.0_KDR
    
!     end associate !-- V_1, etc.
!     end associate !-- CSV, etc.
!     end associate !-- GC

!   end subroutine SetVelocityLinear


!   subroutine SetStream ( SC, C )

!     class ( StreamForm ), intent ( inout ) :: &
!       SC
!     class ( CurrentsForm ), intent ( in ) :: &
!       C

!     integer ( KDI ), dimension ( : ), allocatable :: &
!       iaSelected

!     if ( C % DENSITY_DEFAULT  >  0 ) then
!       allocate ( iaSelected ( 4 ) )
!       iaSelected ( 1 )  =  C % DENSITY_DEFAULT
!       iaSelected ( 2 )  =  C % VELOCITY_DEFAULT_U_1
!       iaSelected ( 3 )  =  C % VELOCITY_DEFAULT_U_2
!       iaSelected ( 4 )  =  C % VELOCITY_DEFAULT_U_3
!     else
!       allocate ( iaSelected ( 0 ) )
!     end if

!     call SC % AddFieldSet &
!            ( C, &
!              iaSelectedOption  =  iaSelected )

!   end subroutine SetStream


!   subroutine Show_FS ( FSC )

!     class ( CurrentsForm ), intent ( in ) :: &
!       FSC

!     call FSC % FieldSetForm % Show ( )

!     call Show ( FSC %  nPrimitive,  'nPrimitive', FSC % IGNORABILITY )
!     call Show ( FSC % iaPrimitive, 'iaPrimitive', FSC % IGNORABILITY )
!     call Show ( FSC %   Primitive,   'Primitive', FSC % IGNORABILITY )

!     call Show ( FSC %  nBalanced,  'nBalanced', FSC % IGNORABILITY )
!     call Show ( FSC % iaBalanced, 'iaBalanced', FSC % IGNORABILITY )
!     call Show ( FSC %   Balanced,   'Balanced', FSC % IGNORABILITY )

!   end subroutine Show_FS


!   subroutine ComputeFromInitial ( C )

!     class ( CurrentsForm ), intent ( inout ) :: &
!       C

!   end subroutine ComputeFromInitial


!   subroutine ComputeFromConserved ( C )

!     class ( CurrentsForm ), intent ( inout ) :: &
!       C

!   end subroutine ComputeFromConserved


!   subroutine ComputeFluxes ( FSC, C, iD )

!     class ( FieldSetForm ), intent ( inout ) :: &
!       FSC
!     class ( CurrentsForm ), intent ( in ) :: &
!       C
!     integer ( KDI ), intent ( in ) :: &
!       iD  !-- iDimension
    
!     integer ( KDI ) :: &
!       iDensity

!     if ( C % DENSITY_DEFAULT > 0 ) then

!       call Search ( C % iaBalanced, C % DENSITY_DEFAULT, iDensity )

!       associate &
!         ( FSS  =>  FSC % Storage_FSC % Storage, &
!           CSS  =>  C % Storage_FSC % Storage, &
!           DeviceMemory  =>  FSC % Storage_FSC % DeviceMemory )
!       associate &
!         ( F_D      =>  FSS % Value ( :, iDensity ), &
!             D      =>  CSS % Value ( :, C % DENSITY_DEFAULT ), & 
!             V_Dim  =>  CSS % Value ( :, C % VELOCITY_DEFAULT_U ( iD ) ) ) 
 
!       call ComputeFluxesKernel &
!              ( D, V_Dim, F_D, UseDeviceOption = DeviceMemory )
  
!       end associate !-- F_D, etc.
!       end associate !-- FSS, etc.

!     end if !-- Density default

!   end subroutine ComputeFluxes


!   subroutine ComputeEigenspeeds ( FSC, C, iaEigenspeeds, iD )

!     class ( FieldSetForm ), intent ( inout ) :: &
!       FSC
!     class ( CurrentsForm ), intent ( in ) :: &
!       C
!     integer ( KDI ), dimension ( : ), intent ( in ) :: &
!       iaEigenspeeds
!     integer ( KDI ), intent ( in ) :: &
!       iD  !-- iDimension

!     if ( C % DENSITY_DEFAULT > 0 ) then

!       associate &
!         ( FSS  =>  FSC % Storage_FSC % Storage, &
!           CSS  =>  C % Storage_FSC % Storage, &
!           DeviceMemory  =>  C % Storage_FSC % DeviceMemory )
!       associate &
!         ( EF_P    =>  FSS % Value ( :, iaEigenspeeds ( 1 ) ), &
!           EF_M    =>  FSS % Value ( :, iaEigenspeeds ( 2 ) ), &
!            V_Dim  =>  CSS % Value ( :, C % VELOCITY_DEFAULT_U ( iD ) ) ) 
 
!       call ComputeEigenspeedsKernel &
!              ( V_Dim, EF_P, EF_M, UseDeviceOption = DeviceMemory )
  
!       end associate !-- EF_P, etc.
!       end associate !-- FSS, etc.

!     end if !-- Density default

!   end subroutine ComputeEigenspeeds


!   impure elemental subroutine Finalize ( C )

!     type ( CurrentsForm ), intent ( inout ) :: &
!       C

!     nullify ( C % Geometry_C )

!     if ( allocated ( C % Balanced ) ) &
!       deallocate ( C % Balanced )
!     if ( allocated ( C % Primitive ) ) &
!       deallocate ( C % Primitive )
!     if ( allocated ( C % iaBalanced ) ) &
!       deallocate ( C % iaBalanced )
!     if ( allocated ( C % iaPrimitive ) ) &
!       deallocate ( C % iaPrimitive )

!   end subroutine Finalize

  
end module Currents_Form
