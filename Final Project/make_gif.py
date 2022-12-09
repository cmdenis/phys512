import matplotlib.pyplot as plt
import os
import imageio as img


# Function to make the gif given the path to figures to assemble and output path
def makegif(frames_path, output_name, fps = 30):
    # List of images
    img_list = os.listdir(frames_path)
    img_list.sort()
    try:
        img_list.remove(".DS_Store")
    except:
        1==1

    # List of frames
    frames = []

    # Append all frames to list
    for frame in img_list:
        frames.append(img.v2.imread(frames_path+"/"+frame))

    # Save as gif
    img.mimsave(output_name, frames, fps = fps)


# 1) one particle
makegif("figs/single_particle", "gifs/1_particle.gif")

# 2) Two particles
makegif("figs/two_particles", "gifs/2_particles.gif")

# 3) Many particles
makegif("figs/periodic", "gifs/3_periodic.gif")
makegif("figs/non_periodic", "gifs/3_non_periodic.gif")

# 4) RK Many Particles
makegif("figs/rk_periodic", "gifs/rk4_periodic.gif")
makegif("figs/rk_nonperiodic", "gifs/rk4_nonperiodic.gif")