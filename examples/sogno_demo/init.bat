@echo off

set PYTHON_DIR=C:\Python36
set SOGNO_CIGRE_SE=Sogno_cigre_se.py
set SOGNO_CIGRE_SE_CIM=Sogno_cigre_se_cim.py
set DUMMY_SIMULATOR=dummy_simulator.py


start CMD /K "title Dummy_Simulator & python %DUMMY_SIMULATOR%"
TIMEOUT /T 3
start CMD /K "title Sogno_cigre_se_cim & python %SOGNO_CIGRE_SE_CIM%"
start CMD /K "title Sogno_cigre_se & python %SOGNO_CIGRE_SE%"
