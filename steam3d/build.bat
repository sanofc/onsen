del render*.bmp
del steam3d.mp4
rmdir /s/q output
mkdir output
@Set Path=C:\Windows\Microsoft.NET\Framework64\v4.0.30319;%PATH%
msbuild /p:Configuration=Release
Release\steam3d.exe
ffmpeg -y -i render_%%02d.bmp -pix_fmt yuv420p steam3d.mp4
steam3d.mp4
