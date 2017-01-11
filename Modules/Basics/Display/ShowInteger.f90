!-- This wrapper exists for use by CONSOLE_Singleton.

subroutine ShowInteger &
             ( Integer, Description, IgnorabilityOption, DisplayRankOption, &
               nLeadingLinesOption, nTrailingLinesOption )

    !-- Convention on argument order violated because the Integer being
    !   "Show"n is more important than the Description.

    use VariableManagement
    use Show_Command

    implicit none

    integer ( KDI ), intent ( in ) :: &
      Integer
    character ( * ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption

    call Show ( Integer, Description, IgnorabilityOption, DisplayRankOption, &
                nLeadingLinesOption, nTrailingLinesOption )

end subroutine ShowInteger
