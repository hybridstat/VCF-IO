@echo off
rem Please consult "https://www.activestate.com/blog/2010/10/how-install-cpan-modules-activeperl"
rem prior to using this build script if installed outside CPAN or Active State ppm
rem Thanks to the Stack Overflow answer below for checking admin rights
rem https://stackoverflow.com/questions/4051883/batch-script-how-to-check-for-admin-rights#11995662

goto check_Permissions

:check_Permissions
    rem echo Administrative permissions required. Detecting permissions...

    net session >nul 2>&1
    if %errorLevel% == 0 (
        perl -MCPAN -e "install qw(IO::Uncompress::Gunzip List::Util Tie::IxHash::Easy Test::Exception Test::File)"

		perl Makefile.PL
		dmake
		dmake test
		dmake install
    ) else (
        echo Failure: Current permissions inadequate.
    )

    pause >nul
