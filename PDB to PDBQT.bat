@echo off
REM Convert all .pdb ligands to .pdbqt using prepare_ligand4.py
REM Handles spaces in folder names

REM --- CONFIG ---
set "PREP_LIG=C:\Program Files (x86)\MGLTools-1.5.7\Lib\site-packages\AutoDockTools\Utilities24\prepare_ligand4.py"
set "INPUT_DIR=C:\Users\Jason\Documents\University\Year 4\PHAY60\ivrPruned\ivr\9ety_site2\pubchemligands"
set "OUTPUT_DIR=C:\Users\Jason\Documents\University\Year 4\PHAY60\ivrPruned\ivr\9ety_site2\pubchemligandsPDBQT"

REM --- create output folder if it doesn't exist ---
if not exist "%OUTPUT_DIR%" mkdir "%OUTPUT_DIR%"

REM --- loop over all .pdb files ---
for %%F in ("%INPUT_DIR%\*.pdb") do (
    echo Converting %%~nxF ...
    python "%PREP_LIG%" -l "%%~fF" -o "%OUTPUT_DIR%\%%~nF.pdbqt" -A hydrogens
)

echo.
echo All ligands converted to PDBQT in %OUTPUT_DIR%
pause
