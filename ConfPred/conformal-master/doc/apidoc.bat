@ECHO OFF

REM sphinx-apidoc --full --force -o . -H "Orange - Conformal Prediction" -A "Biolab" -V "1.0" ../cp
REM sphinx-apidoc --force --separate -o . ../cp ../cp/tests.py
call make.bat html
