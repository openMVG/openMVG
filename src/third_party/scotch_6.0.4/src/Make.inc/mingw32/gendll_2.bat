REM   (C) 2008 Yves Secretan (yves.secretan@ete.inrs.ca)
REM   This software is governed by the CeCILL-C license under French law
REM   and abiding by the rules of distribution of free software. You can
REM   use, modify and/or redistribute the software under the terms of the
REM   CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
REM   URL: "http://www.cecill.info".
REM
REM   To be executed in a DOS window.
REM
REM   This file will create the PT-Scotch DLL. It must be adapted to reflect
REM   your environment, in particular library path and library name.
REM
set VC_PATH="C:\Program Files\Microsoft Visual Studio .NET\VC7\BIN"
%VC_PATH%\vcvars32.bat && %VC_PATH%\lib.exe /def:..\..\..\lib\libptscotch.def /out:..\..\..\lib\libptscotch.dll
