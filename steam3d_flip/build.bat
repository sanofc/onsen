del render*.bmp
del steam3d.mp4
rmdir /s/q output
mkdir output
@Set Path=C:\Windows\Microsoft.NET\Framework64\v4.0.30319;%PATH%
msbuild /p:Configuration=Release
Release\steam3d_flip.exe %1
set tm=%time: =0% 
set dt=%date%
set fn=steam3d_%date:~-10,4%%date:~-5,2%%date:~-2,2%%tm:~0,2%%tm:~3,2%%tm:~6,2%.mp4
ffmpeg -y -i render_%%02d.bmp -pix_fmt yuv420p "%fn%"
"%fn%"	
