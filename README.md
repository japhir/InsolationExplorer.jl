# InsolationExplorer

[![Build Status](https://github.com/japhir/InsolationExplorer.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/japhir/InsolationExplorer.jl/actions/workflows/CI.yml?query=branch%3Amain)

`InsolationExplorer.jl` is an interactive visualization of the Earth's orbit
around the sun, where we can change the eccentricity (degree to which the orbit
is elliptical), the obliquity or tilt (the angle of the Earth's spin axis with
respect to the orbit normal) and the longitude of perihelion with respect to
the moving equinox.

## Installation

We use the [Julia](https://julialang.org) programming languate with the
[Makie](https://docs.makie.org/stable/) plotting library for speedy interactive
exploration.

[Install Julia](https://julialang.org/install/) via `juliaup`.

Run julia and install this package and the plotting library Makie.

```julia
using Pkg
Pkg.add("GLMakie")
Pkg.add("InsolationExplorer")
```

## Getting Started

load the packages
```julia
using GLMakie
GLMakie.activate!()
using InsolationExplorer
```

Calculate insolation for t₀
```julia
insolation(0.016705, 0.4090928042223287, 1.7962486166737615;
                  longitude = pi/2, latitude = deg2rad(65),
                  S0 = 1360.7, H = nothing)
```

Create the plot
```julia
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
