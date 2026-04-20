ffmpeg -r 50  -start_number 1 -i frame-%%06d.png  -filter_complex "[0:v] fps=20,scale=640:-1,split [a][b];[a] palettegen [p];[b][p] paletteuse" -vcodec libx264 -pix_fmt yuv420p gistemp_frames_with_coast.mp4

ffmpeg -r 50  -start_number 1 -i frame-%%06d.png  -filter_complex "[0:v] fps=20,scale=480:-1,split [a][b];[a] palettegen [p];[b][p] paletteuse"  gistemp_frames_with_coast.gif
