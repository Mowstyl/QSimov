@ECHO OFF
SET QSimovVersion=5.2.0
:: Removing old wheels
del "dist\qsimov-*.whl"
:: Removing old sources
del "dist\qsimov-*.tar.gz"

:: Cleaning Python3.12
py -3.12 -m pip uninstall qsimov -y
rmdir "qsimov_Mowstyl.egg-info" /S /Q
:: Building Python3.12
py -3.12 -m build
:: Installing Python3.12
py -3.12 -m pip install --user dist/qsimov-%QSimovVersion%-py3-none-any.whl

pause
