@ECHO OFF
SET QSimovVersion=5.0.5
:: Removing old wheels
del "dist\qsimov_Mowstyl-*.whl"
:: Removing old sources
del "dist\qsimov-Mowstyl-*.tar.gz"

:: Cleaning Python3.9
py -3.9 -m pip uninstall qsimov_Mowstyl -y
rmdir "qsimov_Mowstyl.egg-info" /S /Q
:: Building Python3.9
py -3.9 -m build
:: Installing Python3.9
py -3.9 -m pip install --user dist/qsimov_Mowstyl-%QSimovVersion%-py3-none-any.whl

pause
