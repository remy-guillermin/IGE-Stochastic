# FILM

After creating the frames with [film_std.py](film_std.py), one can create the movie associated with 
```bash
ffmpeg -framerate 30 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p {name}.mp4
```