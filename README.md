# InsolationExplorer.jl

This is an interactive visualization of the Earth's orbit around the sun, where
we can change the eccentricity (degree to which the orbit is elliptical), the
obliquity or tilt (the angle of the Earth's spin axis with respect to the orbit
normal) and the longitude of perihelion with respect to the moving equinox.

We use the [Julia](https://julialang.org) programming languate with the
[Makie](https://docs.makie.org/stable/) plotting library for speedy interactive
exploration.

## Getting Started

[Install Julia](https://julialang.org/install/), then add `GLMakie`.

Download this package, then navigate to the folder, launch Julia and
execute the following code line-by-line:

```julia
# activate the current directory as the active package
using Pkg
Pkg.activate(".")

# load the packages
using GLMakie
GLMakie.activate!()
using InsolationExplorer

# create the plot
f = explore_insolation()
```

This results in the following image:

![](explore_insolation.png)

## Explore the Insolation!

On the left, we plot the insolation at the top of the atmosphere in watts per square metre (colour) as a function of the Earth's latitude and the true solar longitude.
Below are sliders to change the orbital configuration!
On the right, we see a 3d visualization of the orbit, the Earth, and if you zoom out or reduce the scale factor, the Sun.
Drag the left mouse button to rotate around, right mouse button to shift around. Scroll to zoom in/out.

To get started, I would play around with the obliquity, then exaggerate eccentricity, and then slide the longitude of perihelion around.
