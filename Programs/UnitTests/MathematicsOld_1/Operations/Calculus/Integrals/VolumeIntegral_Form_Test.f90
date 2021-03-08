module VolumeIntegral_Form_Test__Form

  use Basics
  use Manifolds
  use VolumeIntegral_Form

  implicit none
  private

  character ( LDF ), public, parameter :: &
    ProgramName = 'VolumeIntegral_Form_Test'

  type, public :: VolumeIntegral_Form_Test_Form
    type ( GridImageStreamForm ), allocatable :: &
      GridImageStream
    type ( Atlas_SC_CC_Form ), allocatable :: &
      Atlas
    type ( Storage_ASC_Form ), allocatable :: &
      Integrand_ASC
    type ( VolumeIntegralForm ), allocatable :: &
      VolumeIntegral
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type VolumeIntegral_Form_Test_Form

contains


  subroutine Initialize ( VIFT, Name )

    class ( VolumeIntegral_Form_Test_Form ), intent ( inout ), target :: &
      VIFT
    character ( * ), intent ( in ) :: &
      Name

    real ( KDR ), dimension ( 1 ) :: &
      Integral
    type ( Real_1D_Form ), dimension ( 1 ) :: &
      Integrand
    class ( StorageForm ), pointer :: &
      I

    !-- Stream

    allocate ( VIFT % GridImageStream )
    associate ( GIS => VIFT % GridImageStream )
    call GIS % Initialize &
           ( Name, CommunicatorOption = PROGRAM_HEADER % Communicator )
    
    !-- Atlas, etc.

    allocate ( VIFT % Atlas )
    associate ( A => VIFT % Atlas )
    call A % Initialize ( Name, PROGRAM_HEADER % Communicator )
    call A % CreateChart_CC ( )  
    call A % SetGeometry ( )

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
       
    !-- Integrand
    
    allocate ( VIFT % Integrand_ASC )
    associate ( IA => VIFT % Integrand_ASC )
    call IA % Initialize &
           ( A, 'Integrand', nFields = 1, VariableOption = [ 'Integrand' ], &
             WriteOption = .true. )

    I => IA % Storage ( )
    I % Value = 1.0_KDR
    call Integrand ( 1 ) % Initialize ( I % Value ( :, 1 ) )

    !-- VolumeIntegral

    allocate ( VIFT % VolumeIntegral )
    associate ( VI => VIFT % VolumeIntegral )

    call VI % Compute ( C, Integrand, Integral )
    call Show ( Integral ( 1 ), '*** VolumeIntegral', nLeadingLinesOption = 1 )
    call Show ( 4.0_KDR / 3.0_KDR  *  CONSTANT % PI  &
                *  C % MaxCoordinate ( 1 ) ** 3, &
                '*** Expected', nTrailingLinesOption = 1 )

    !-- Write

    call A % OpenStream ( GIS, '1', iStream = 1 )

    call GIS % Open ( GIS % ACCESS_CREATE )
    call A % Write ( iStream = 1 )
    call GIS % Close ( )

    !-- Cleanup

    end associate !-- VI
    end associate !-- IA
    end select !-- C
    end associate !-- A
    end associate !-- GIS
    nullify ( I )
    
  end subroutine Initialize


  subroutine Finalize ( VIFT )

    type ( VolumeIntegral_Form_Test_Form ), intent ( inout ) :: &
      VIFT

    deallocate ( VIFT % VolumeIntegral )
    deallocate ( VIFT % Integrand_ASC )
    deallocate ( VIFT % Atlas )
    deallocate ( VIFT % GridImageStream )

  end subroutine Finalize


end module VolumeIntegral_Form_Test__Form



program VolumeIntegral_Form_Test

  use Basics
  use VolumeIntegral_Form_Test__Form

  implicit none
  
  type ( VolumeIntegral_Form_Test_Form ), allocatable :: &
    VIFT
    
  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( ProgramName )
    
  allocate ( VIFT )
  call VIFT % Initialize ( PROGRAM_HEADER % Name )
  deallocate ( VIFT )

  deallocate ( PROGRAM_HEADER )

end program VolumeIntegral_Form_Test
