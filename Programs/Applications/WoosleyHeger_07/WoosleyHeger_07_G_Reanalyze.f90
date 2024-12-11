program WoosleyHeger_07_G_Reanalyze

  !-- WoosleyHeger_07_Grey

  use GenASiS
  use WoosleyHeger_07_NM__Form

  implicit none

  type ( StorageForm ), allocatable :: &
    Measures_CC
  type ( WoosleyHeger_07_NM_Form ), allocatable :: &
    WH

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'WoosleyHeger_07_G', DimensionalityOption = '1D' )

  allocate ( WH )
  call WH % Initialize_NM ( 'GREY', PROGRAM_HEADER % Name )
  WH % Integrator % ReanalyzeCheckpoint => ReanalyzeCheckpoint
  call WH % Reanalyze ( )
  deallocate ( WH )

  deallocate ( PROGRAM_HEADER )


contains

  subroutine ReanalyzeCheckpoint ( I, CI, nR, iR )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    type ( CurveImageForm ), intent ( inout ), allocatable :: &
      CI
    integer ( KDI ), intent ( in ) :: &
      nR, &
      iR
    
    integer ( KDI ) :: &
      iM  !-- iMeasure
    type ( QuantityForm ), dimension ( : ), allocatable :: &
      SeriesUnit
    character ( LDL ), dimension ( : ), allocatable :: &
      SeriesName
    
    select type ( M => WH % Measures )
    class is ( Measures_R_CC_C_Form )
    
    if ( .not. allocated ( Measures_CC ) ) then
   
      allocate ( SeriesName ( M % nMeasures ) )
      allocate ( SeriesUnit ( M % nMeasures ) )
      do iM  =  1,  M % nMeasures
        SeriesName ( iM )  =  M % Name ( iM )
        SeriesUnit ( iM )  =  M % Unit ( iM )
      end do !-- iS

      allocate ( Measures_CC )
      call Measures_CC % Initialize & 
             ( [ nR, M % nMeasures ], &
               VariableOption = SeriesName, UnitOption = SeriesUnit, &
               NameOption =  'Measures_CC', &
               ClearOption = .true. )
      if ( allocated ( CI ) ) &
        call CI % AddStorage ( Measures_CC )
   
    end if
    
    call I % Analyze ( I % IGNORABILITY )
    
    associate &
      (  MCC_V =>  Measures_CC % Value, &
            MV =>  M % Value, &
            nM =>  M % nMeasures )
    do iM  =  1,  nM
      MCC_V ( iR, iM )  =  MV ( iM )
    end do !-- iM
    end associate !-- SMV, etc.
    
    end select
  
  end subroutine ReanalyzeCheckpoint


end program WoosleyHeger_07_G_Reanalyze
