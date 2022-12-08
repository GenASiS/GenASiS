module BundleTranspose_Form

  use Basics
  use Manifolds

  implicit none
  private

  type, public :: BundleTransposeForm
    integer ( KDI ) :: &
      IGNORABILITY = 0
    type ( CommunicatorForm ), allocatable :: &
      Communicator
    type ( MessageIncoming_1D_R_Form ), allocatable :: &
      Incoming_F_S, &
      Incoming_S_F
    type ( MessageOutgoing_1D_R_Form ), allocatable :: &
      Outgoing_F_S, &
      Outgoing_S_F
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      StartTranspose_F_S
    ! procedure, public, pass :: &
    !   FinishTranspose_F_S
    ! procedure, public, pass :: &
    !   StartTranspose_S_F
    ! procedure, public, pass :: &
    !   FinishTranspose_S_F
    final :: &
      Finalize
  end type BundleTransposeForm

    private :: &
      Start_F_S_ASCG_ASCG

      private :: &
        LoadMessage_F_S_ASCG_ASCG

    integer ( KDI ), private, parameter :: &
      TAG_RECEIVE_F_S  = 999, &
      TAG_RECEIVE_S_F  = 998, &
      TAG_SEND_F_S     = TAG_RECEIVE_F_S, &
      TAG_SEND_S_F     = TAG_RECEIVE_S_F


contains


  subroutine Initialize ( BT )

    class ( BundleTransposeForm ), intent ( inout ) :: &
      BT

    BT % IGNORABILITY  =  CONSOLE % INFO_4

  end subroutine Initialize


  subroutine StartTranspose_F_S ( BT, B, S_F, S_S, DevicesCommunicate )

    class ( BundleTransposeForm ), intent ( inout ) :: &
      BT
    class ( Bundle_H_Form ), intent ( in ) :: &
      B
    class ( StorageForm ), dimension ( : ), intent ( in ) :: &
      S_F, &  !-- Storage_Fiber
      S_S     !-- Storage_Section
    logical ( KDL ), intent ( in ) :: &
      DevicesCommunicate

    call Show ( 'Starting bundle transpose (fiber to section)', &
                BT % IGNORABILITY )
    call Show ( S_F % Name, 'FieldSet_Fiber', BT % IGNORABILITY )
    call Show ( S_S % Name, 'FieldSet_Section', BT % IGNORABILITY )

    select type ( B )
    class is ( Bundle_ASCG_ASCG_Form )

      call Start_F_S_ASCG_ASCG &
             ( BT, B, B % Portal_F_S, S_F, DevicesCommunicate )

    class default
      call Show ( 'Bundle type not recognized', CONSOLE % ERROR )
      call Show ( 'BundleTransposeForm', 'module', CONSOLE % ERROR )
      call Show ( 'StartTranspose_F_S', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select  !-- B

  end subroutine StartTranspose_F_S


  impure elemental subroutine Finalize ( BT )

    type ( BundleTransposeForm ), intent ( inout ) :: &
      BT

    if ( allocated ( BT % Outgoing_S_F ) ) &
      deallocate ( BT % Outgoing_S_F )
    if ( allocated ( BT % Outgoing_F_S ) ) &
      deallocate ( BT % Outgoing_F_S )

    if ( allocated ( BT % Incoming_S_F ) ) &
      deallocate ( BT % Incoming_S_F )
    if ( allocated ( BT % Incoming_F_S ) ) &
      deallocate ( BT % Incoming_F_S )

    if ( allocated ( BT % Communicator ) ) &
      deallocate ( BT % Communicator )

  end subroutine Finalize


  subroutine Start_F_S_ASCG_ASCG &
               ( BT, B, PH, S_F, DevicesCommunicate )

    class ( BundleTransposeForm ), intent ( inout ) :: &
      BT
    class ( Bundle_ASCG_ASCG_Form ), intent ( in ) :: &
      B
    type ( PortalHeaderForm ), intent ( in ) :: &
      PH
    class ( StorageForm ), dimension ( : ), intent ( in ) :: &
      S_F  !-- Storage_Fiber
    logical ( KDL ), intent ( in ) :: &
      DevicesCommunicate

    integer ( KDI ) :: &
      iT  !-- iTarget

    !-- Allocate on first use

    if ( .not. allocated ( BT % Communicator ) ) then
      allocate ( BT % Communicator )
      call BT % Communicator % Initialize ( B % Communicator )
    end if

    associate &
      ( Communicator        =>  BT % Communicator, &
        nCellsSectionsFrom  =>  PH % nChunksFrom, &
        nFibersBinsTo       =>  PH % nChunksTo, &
        nFields             =>  size ( S_F ) )
 
    if ( .not. allocated ( BT % Incoming_F_S ) &
         .and. .not. allocated ( BT % Outgoing_F_S ) ) then
    
      allocate ( BT % Incoming_F_S )
      allocate ( BT % Outgoing_F_S )

      call BT % Incoming_F_S % Initialize &
             ( Communicator, &
               spread ( TAG_RECEIVE_F_S, dim = 1, ncopies = PH % nSources ), &
               PH % Source, &
               nCellsSectionsFrom * nFields )
      call BT % Outgoing_F_S % Initialize &
             ( Communicator, &
               spread ( TAG_SEND_F_S, dim = 1, ncopies = PH % nTargets ), &
               PH % Target, &
               nFibersBinsTo * nFields )
    
      if ( DevicesCommunicate ) then
        call BT % Incoming_F_S % AllocateDevice ( )
        call BT % Outgoing_F_S % AllocateDevice ( )
      end if 
    
    end if  !-- allocated faces
    
    !-- Post Receives

    call BT % Incoming_F_S % Receive ( )

    !-- Post Sends

    do iT  =  1,  PH % nTargets

      call LoadMessage_F_S_ASCG_ASCG &
             ( B, S_F, BT % Outgoing_F_S % Message ( iT ), DevicesCommunicate, &
               nFibersBinsTo ( iT ), iT )
      
      call BT % Outgoing_F_S % Send ( iT )

    end do !-- iT

    !-- Cleanup

    end associate  !-- Communicator, etc.

  end subroutine Start_F_S_ASCG_ASCG


  subroutine LoadMessage_F_S_ASCG_ASCG &
               ( B, S_F, OutgoingMessage, DevicesCommunicate, nFibersBins, iT )

    class ( Bundle_ASCG_ASCG_Form ), intent ( in ) :: &
      B
    class ( StorageForm ), dimension ( : ), intent ( in ) :: &
      S_F  !-- Storage_Fiber
    type ( MessageOutgoing_R_Form ), intent ( in ) :: &
      OutgoingMessage
    logical ( KDL ), intent ( in ) :: &
      DevicesCommunicate
    integer ( KDI ) :: &
      nFibersBins, &
      iT  !-- iTarget

    integer ( KDI ) :: &
      iFld, &        !-- iField
      iFbr, &        !-- iFiber
      iB, jB, kB, &  !-- iBin, etc.
      nBins, &
      oV, &          !-- oValue
      iV             !-- iValue
    real ( KDR ), dimension ( :, :, :, : ), pointer :: &
      FV  !-- FieldValue

    associate &
      ( nFlds  =>  size ( S_F ), &
        iFbrF  =>  B % iFiberExchangeFirst ( iT ), &
        iFbrL  =>  B % iFiberExchangeLast  ( iT ), &
           nB  =>  B % Chart_GS_Fiber % nCells, &
         iaBF  =>  B % iaBinExchangeFirst ( iT ) % Value, &
         iaBL  =>  B % iaBinExchangeLast ( iT ) % Value, &
          OMV  =>  OutgoingMessage % Value, &
            C  =>  B % Chart_GS_Fiber )

    do iFld  =  1,  nFlds

      call C % SetFieldPointer ( S_F ( iFld ) % Value, FV )
      nBins  =  nFibersBins  /  ( iFbrL - iFbrF + 1 )

      !-- FIXME: parallelize over fibers, private oV and iV
      do iFbr  =  iFbrF,  iFbrL

        oV  =  ( iFbr - 1 ) * nBins 
        iV  =  0

        kBin: do kB  =  1,  nB ( 3 )
          if ( kB  <  iaBF ( 3 ) ) &
            cycle kBin

          jBin: do jB  =  1,  nB ( 2 )
            if ( kB  ==  iaBF ( 3 ) .and. jB  <  iaBF ( 2 ) ) &
              cycle jBin

            iBin: do iB  =  1,  nB ( 1 )
              if ( kB  ==  iaBF ( 3 ) .and. jB  ==  iaBF ( 2 )  &
                   .and. iB  <  iaBF ( 1 ) ) &
                cycle iBin

              iV  =  iV  +  1
              OMV ( oV + iV )  =  FV ( iB, jB, kB, iFbr ) 

            end do iBin !-- iB
          end do jBin !-- jB
        end do kBin !-- kB
      end do !-- iFbr
    end do !-- iFld

    end associate !-- nFlds, etc.

  end subroutine LoadMessage_F_S_ASCG_ASCG


end module BundleTranspose_Form
