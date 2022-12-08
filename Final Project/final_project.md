# Final Project

By Christian Denis for PHYS-512

## Part 1

Here we use an $n$-body simulation to show the behavior of a single free particle. This is more of a test than anything. Nothing too interesting is happening here.

We expect our particle to stay put, and indeed this is what is happening as demonstrated by this uneventful GIF.

![2_two_particles](gifs/1_particle.gif)

Maybe there's some stuff happening on the smaller scale that is not displayed by the GIF however... Let's make a plot of the XY positions over time to see if there's stuff we're missing:

![2_two_particles](figs/single_position.jpg)

Looks pretty static.

## Part 2

Here we use the same $n$-body simulation to show two particles orbiting each other. This result is more interesting than the previous one, but is also sort of a test to make sure our simulation doesn't present any problems. The following GIF shows the simulation:

![2_two_particles](gifs/2_particles.gif)

We see that the two particles are indeed in orbit.

An interesting observation we can make with this is the energy of the system at different points in time. If we plot the kinetic, potential and total energies of the system we get:

![2_two_particles](figs/2particles_energy.jpg)

It appears as though the energy is somewhat constant over the long term, however, it seems to be oscillating through the simulation, this is somewhat surprizing and I have a hard time coming up with an explanation...

It is worth pointing out that for smaller grid sizes there is one somewhat undesirable feature of the simulation is the fast that the particles keep drifting. I think this might be caused by the fact that the potential exists at the center of each cell in a grid which does not exactly correspond to the location of the particles.

Initially I was wondering wether this might be caused by the periodic boundary conditions, but after checking with and without periodic boundary conditions, the same effect was observed in both case, meaning that there it's unlikely due to that.

## Part 3

Here we simulate periodic and non-periodic systems of $n$ particles (here $n$ is about 600, in both cases).

![2_two_particles](gifs/3_periodic.gif)

![2_two_particles](gifs/3_non_periodic.gif)

If we plot the total energy as a function of time, we obtain

![2_two_particles](figs/per_vs_nonper.jpg)

## Part 4

![rk_periodic](gifs/rk4_periodic.gif)

![rk_nonperiodic](gifs/rk4_nonperiodic.gif)

If we plot the total energy as a function of time, we obtain

![2_two_particles](figs/rk_per_vs_nonper.jpg)