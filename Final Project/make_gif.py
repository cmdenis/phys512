import matplotlib.pyplot as plt
import os
import imageio as img


# 2) Two particles

# List of images
img_list = os.listdir("figs/two_particles")
img_list.sort()
img_list.remove(".DS_Store")

# List of frames
frames = []

for frame in img_list:
    frames.append(img.v2.imread("figs/two_particles/"+frame))

img.mimsave("2_particles.gif", frames, fps = 30)



