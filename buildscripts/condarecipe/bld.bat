conda config --add channels poliastro
"%PYTHON%" setup.py install
if errorlevel 1 exit 1
