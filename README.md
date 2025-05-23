# InsolationExplorer

[![Build Status](https://github.com/japhir/InsolationExplorer.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/japhir/InsolationExplorer.jl/actions/workflows/CI.yml?query=branch%3Amain)

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

```
using Pkg; Pkg.add("https://github.com/japhir/InsolationExplorer.jl")
```

Running this will take a while, because it will also install `GLMakie`, the
plotting library and all other dependencies. However, after the first install
and launch, everything should be a lot faster!


## Getting Started

Load the package and create the plot
```julia
using InsolationExplorer
f = explore_insolation()
```

This results in the following image:

![](explore_insolation.png)

## Explore the Insolation!

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
