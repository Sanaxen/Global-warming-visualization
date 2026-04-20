ffmpeg -r 10  -start_number 1 -i frame-%%06d.png  -filter_complex "[0:v] fps=20,scale=640:-1,split [a][b];[a] palettegen [p];[b][p] paletteuse" -vcodec libx264 -pix_fmt yuv420p spiral_top_then_tilt_v2_frames.mp4

ffmpeg -r 10  -start_number 1 -i frame-%%06d.png  -filter_complex "[0:v] fps=20,scale=640:-1,split [a][b];[a] palettegen [p];[b][p] paletteuse"  spiral_top_then_tilt_v2_frames.gif
ffmpeg -r 10  -start_number 1 -i frame-%%06d.png  -filter_complex "[0:v] fps=20,scale=320:-1,split [a][b];[a] palettegen [p];[b][p] paletteuse"  spiral_top_then_tilt_v2_frames_mini.gif
