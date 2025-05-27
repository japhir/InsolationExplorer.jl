# InsolationExplorer
<!-- [![Build Status](https://github.com/japhir/InsolationExplorer.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/japhir/InsolationExplorer.jl/actions/workflows/CI.yml?query=branch%3Amain) -->

`InsolationExplorer.jl` is an interactive visualization of the Earth's orbit
around the sun, where we can change the eccentricity (degree to which the orbit
is elliptical), the obliquity or tilt (the angle of the Earth's spin axis with
respect to the orbit normal) and the longitude of perihelion with respect to
the moving equinox.

## Installation

We use the [Julia](https://julialang.org) programming language with the
[Makie](https://docs.makie.org/stable/) plotting library for speedy interactive
exploration.

1. [Install Julia](https://julialang.org/install/).
2. Open the [REPL](https://docs.julialang.org/en/v1/stdlib/REPL/) and
   then copy-paste the below and hit enter.

```julia
using Pkg; Pkg.add("https://github.com/japhir/InsolationExplorer.jl")
```

Running this will take a while, because it will also install `GLMakie`, the
plotting library and all other dependencies. However, after the first install
and launch, everything should be a lot faster!


## Explore the Insolation

Load the package and create the plot
```julia
using InsolationExplorer
i = explore_insolation()
```

This results in the following image:

![](explore_insolation.png)

On the left, we plot the insolation at the top of the atmosphere in watts per
square metre (colour) as a function of the Earth's latitude and the true solar
longitude. Below are sliders to change the orbital configuration!

On the right, we see a 3d visualization of the orbit, the Earth, and if you
zoom out or reduce the scale factor, the Sun. Drag the left mouse button to
rotate around, drag the right mouse button to shift around. Scroll to zoom
in/out.

To get started, I would play around with the obliquity, then exaggerate
eccentricity, and then slide the longitude of perihelion around. Note that the
ranges for the sliders allow for all possible values. If you double-click on
the slider it will reset it to Earth's modern values.

## Explore Solution

In the previous plot it is easy to change eccentricity and obliquity to values
that have never occurred for the Earth. The below plot explores the insolation
as a function of the `ZB18a(1,1)` astronomical solution.

```julia
s = explore_solution()
```

which results in the following image:

![](explore_solution.png)

On the left, we plot the eccentricity, climatic precession, and obliquity. Drag
a rectangle/scroll to zoom in the top panel. Control + left click to reset the
zoom. Click in the bottom panel or slide the time slider to select a time.
Press space to toggle auto play and proceed through time. Use the left and
right arrows for precise seeking.

On the top right, we plot the insolation at the top of the atmosphere in watts
per square metre (colour) as a function of the Earth's latitude and the true
solar longitude. On the bottom right, we see a 3d visualization of the orbit,
the Earth, and if you zoom out or reduce the scale factor, the Sun. Drag the
left mouse button to rotate around, right mouse button to shift around. Scroll
to zoom in/out.

Note that the reference frame remains the same as the previous one, where the
vernal equinox is fixed in the positive x-axis direction and the orbit normal
in the positive y-axis direction.

<!-- Also note that for the 3d visualization we always plot the ellipse with a
fixed semimajor axis length of 1, rather than the time-varying semimajor axis
with extrema of 0.9999722469706677 and 1.000035840979405. -->
