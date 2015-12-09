setlocal
(SET PATH=.;D:\Rtools31\bin;D:\perl64\bin;D:\Rtools31\gcc-4.6.3\bin;D:\miktex29\miktex\bin\x64;D:\Program Files\R\R-3.1.3\bin\x64;C:\windows;C:\windows\system32)
set package="E:\cggreen\Research\dissertation2\chapters\calibration\HardinRockeExtension\HardinRockeExtension"
R CMD INSTALL --library="D:\Program Files\R\R-3.1.3\site-library" %package%
R CMD INSTALL --library="D:\Program Files\R\R-3.0.2\site-library" %package%
R CMD INSTALL --library="D:\Program Files\R\R-2.15.2\site-library" %package% 
d:\cygwin\bin\sleep 5
endlocal
