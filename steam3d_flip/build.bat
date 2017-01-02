@Set Path=C:\Windows\Microsoft.NET\Framework64\v4.0.30319;%PATH%
msbuild /p:Configuration=Release
Release\steam3d_flip.exe %1
set tm=%time: =0% 
set dt=%date%
set dn=mov_%date:~-10,4%%date:~-5,2%%date:~-2,2%%tm:~0,2%%tm:~3,2%%tm:~6,2%
mkdir "%dn%"
move render_*.bmp %dn%
ffmpeg -y -i "%dn%\render_%%02d.bmp" -pix_fmt yuv420p "%dn%\steam3d.mp4"
rem mogrify -format eps "%dn%\*.bmp"