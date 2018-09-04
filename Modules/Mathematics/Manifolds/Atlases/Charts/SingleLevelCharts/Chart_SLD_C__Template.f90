!-- Chart_SLD_C represents a distributed single-level chart using spherical 
!   coordinates and proportional radial spacing.

module Chart_SLD_C__Template

  !-- Chart_SingleLevelDistributed_Central_Template

  use Basics
  use AtlasBasics
  use ChartBasics
  use Chart_SLD__Form

  implicit none
  private

  type, public, extends ( Chart_SLD_Form ), abstract :: Chart_SLD_C_Template
    integer ( KDI ) :: &
      nCellsPolar, &
      nPillars_2, &
      nPillars_3
  end type Chart_SLD_C_Template

end module Chart_SLD_C__Template
