!-- This wrapper exists for use by CONSOLE_Singleton.

subroutine ShowMessage &
             ( Character, IgnorabilityOption, DisplayRankOption, &
               nLeadingLinesOption, nTrailingLinesOption )

    use VariableManagement
    use KIND_DEFAULT_Singleton
    use Show_Command

    implicit none

    character ( * ), intent ( in ) :: &
      Character
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption

    call Show ( Character, IgnorabilityOption, DisplayRankOption, &
                nLeadingLinesOption, nTrailingLinesOption )

end subroutine ShowMessage
