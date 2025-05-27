"""
    explore_solution()

Plot the interactive astronomical solution explorer.

On the left, we plot the eccentricity, climatic precession, and obliquity. Drag
a rectangle/scroll to zoom in the top panel. Control + left click to reset the
zoom. Click in the bottom panel or slide the time slider to select a time.
Press space to toggle auto play and proceed through time. Use the left and
right arrows for precise seeking.

On the top right, we plot the insolation at the top of the atmosphere in watts per
square metre (colour) as a function of the Earth's latitude and the true solar
longitude.
On the bottom right, we see a 3d visualization of the orbit, the Earth, and if you
zoom out or reduce the scale factor, the Sun. Drag the left mouse button to
rotate around, right mouse button to shift around. Scroll to zoom in/out.

We use textures from https://www.solarsystemscope.com/textures/

See also [`explore_insolation`](@ref).
"""
function explore_solution()
    bg_color = RGB(0.2, 0.2, 0.2)
    fg_color = RGB(0.9, 0.9, 0.9)

    # read reference data
    ref1_1 = CSV.read("out.dat",
                      DataFrame;
                      delim = ' ',
                      stripwhitespace = true,
                      ignorerepeated = true,
                      header = [:time, :eccentricity, :obliquity, :precession, :lpx, :climatic_precession])

    fins = Figure(; backgroundcolor = bg_color,
                  textcolor = fg_color)
    sg = SliderGrid(fins[3,1],
                    (label = "Time\n[←→]",
                     range = -1e5:1:0,
                     format = "{:.0f} kyr",
                     startvalue = -1e5),
                    (label = L"$\lambda_\odot$: True Solar Longitude",
                     range = 0:360,
                     format = "{:.2f}°",
                     startvalue = 90),
                    (label = "Size factor",
                     range = [1.0, 100.0, 200.0, 1000.0, 2000.0, 3000.0, 4000.0],
                     format = "{:.2f}×",
                     startvalue = 3000.0)
                    )
    for s in sg.sliders
        s.color_active[] = RGBA(1.0, 0.9, 0.4, 1.0)
        s.color_active_dimmed[] = RGBA(0.7, 0.6, 0.4, 1.0)
        s.color_inactive[] = RGBA(0.3, 0.3, 0.3, 1.0)
    end

    # make the Observables available
    tm = sg.sliders[1].value
    lon = sg.sliders[2].value
    scale_factor = sg.sliders[3].value

    gl = GridLayout(fins[4,1], tellwidth = false)
    Label(gl[1,1], "Play/Pause ⏯\n[Space]", color = fg_color)
    tg = Toggle(gl[1,2], active = false,
                halign = :left,
                buttoncolor = fg_color,
                framecolor_inactive = bg_color,
                tellwidth = false, tellheight = true)
    tb = Textbox(gl[1,3], placeholder = "Enter Solar Constant S₀\ndefault = 1367.1 Wm¯²",
            validator = Float64,
            tellwidth = false, tellheight = true, )
    # 1361 # Eocene
    # 1368 # same as La04 webtool
    S0 = Observable(1367.1) # Menviel 2019, same as models
    on(tb.stored_string) do s
        S0[] = parse(Float64, s)
    end

    ax_ecc = Axis(fins[1,1],
                  xlabel = "",
                  ylabel = L"$e\sin{\bar\omega}$: Climatic\\Precession (-)")
    ax_obl = Axis(fins[2,1],
                  xlabel = "",
                  ylabel = L"$\epsilon$: Obliquity (°)")
    linkxaxes!(ax_ecc, ax_obl)
    xlims!(ax_obl, (-3000, 0))

    lines!(ax_ecc, ref1_1.time, ref1_1.climatic_precession,
           color = :gray, label = "Climatic Precession")
    lines!(ax_ecc, ref1_1.time, ref1_1.eccentricity,
           color = :black, label = "Eccentricity")
    axislegend(ax_ecc, position = :lt, merge = true, unique = true, labelcolor = :black)
    lines!(ax_obl, ref1_1.time, rad2deg.(ref1_1.obliquity),
           color = :black, label = "ZB18a(1,1)")

    # add current time
    vlines!(ax_ecc, tm, color = :purple)
    vlines!(ax_obl, tm, color = :purple)

    # make time clickable with left-click
    Makie.deactivate_interaction!(ax_obl, :rectanglezoom)
    spoint = select_point(ax_obl.scene)
    on(spoint) do p
        x, y = p
        set_close_to!(sg.sliders[1], x)
    end
    # allow left and right arrow to move back/forward in time
    on(events(fins).keyboardbutton) do event
        if event.action == Keyboard.press || event.action == Keyboard.repeat
            if event.key == Keyboard.left
                set_close_to!(sg.sliders[1], tm[] - 1)
            elseif event.key == Keyboard.right
                set_close_to!(sg.sliders[1], tm[] + 1)
            elseif event.key == Keyboard.space
                tg.active[] = !tg.active[]
                notify(tg.active)
            end
        end
    end
    # when we toggle the play/pause, it will start to scrub through time
    on(fins.scene.events.tick) do tick
        tg.active[] || return
        tm[] > 0 && return
        tm[] += 5 * tick.delta_time
        # set_close_to!(sg.sliders[1], tm[] + 5 * tick.delta_time)
        notify(tm)
    end

    # subset current time slice
    xslice = lift(tm) do t
        # r = ref1_1[ref1_1.time .== t, :]
        @view ref1_1[partialsortperm(abs.(ref1_1.time .- t), 1), :]
    end
    ecc = lift(xslice) do r
        r.eccentricity[1]
    end
    obl = lift(xslice) do r
         r.obliquity[1]
    end
    lpx = lift(xslice) do r
         mod((r.lpx[1] - pi), 2pi)
    end

    # calculate insolation grid
    lats = deg2rad.(-90:90 |> collect)
    lons = deg2rad.(0:359 |> collect)
    rf = zeros(length(lons), length(lats))
    ins = lift(ecc, obl, lpx, S0) do ecc, obl, lpx, S0
        for (j, lat) in enumerate(lats), (i, lon) in enumerate(lons)
            rf[i,j] = insolation(ecc,
                                 obl,
                                 lpx,
                                 longitude = lon,
                                 latitude = lat,
                                 S0 = S0, H = nothing)
        end
        rf
    end

    # plot insolation grid
    ax_ins = Axis(fins[1,2],
                  xlabel = "True Solar Longitude (°)",
                  ylabel = "Latitude (°)",
                  backgroundcolor = bg_color)
    cb1 = image!(ax_ins,
                 rad2deg.(extrema(lons)), rad2deg.(extrema(lats)),
                 ins,
                 interpolate = false,
                 colormap = :solar)
    # scatter!(ax_ins, lon, lat)
    vlines!(ax_ins, lon)
    # TODO: I would like contours, but the labels bug out when updating the observables
    contour!(ax_ins,
             extrema(lons) ./ pi .* 180, extrema(lats) ./ pi .* 180,
             ins,
             labels = true,
             levels = 100:100:600,
             color = :white)
    Colorbar(fins[1,3], cb1, vertical = true, label = L"Insolation (Wm$^{-2}$)")

    # 3d illutration of eccentricity, lpx, and obliquity
    # note that this is in the coordinate system moving with the orbit plane
    # where the x-axis is aligned to the ascending node
    # and at t0, where the precession angle phi is 0.
    a = 1
    ef = 0.3 # the scale factor for vectors drawn from Earth origin
    ep = lift(ecc,lpx) do e,l
        b = a * sqrt(1 - e^2)
        c = e * a # distance from center to each focus
        focus1 = Point3f(-c, 0, 0)
        focus2 = Point3f(c, 0, 0)
        θ = range(0, 2π, length=360)
        x = a * cos.(θ)
        y = b * sin.(θ)
        z = zero(x)  # Initially on the XY plane
        # shift ellipse so that one focus is at the origin
        ellipse_points = [Point3f(x[i] - c, y[i], z[i]) for i in eachindex(x)]
        # rotate ellipse by 90° so that vernal equinox is at 0
        # Rve = RotMatrix3(AngleAxis(-pi, 0,0,1))
        # rotate ellipse by lpx + pi
        # Rz = RotMatrix3(AngleAxis(deg2rad(mod(l + 180, 360)), 0,0,1))
        Rz = RotMatrix3(AngleAxis(mod(l + pi, 2pi), 0,0,1))
        rotated_ellipse = [Rz * p for p in ellipse_points]
        return rotated_ellipse
    end
    perihelion = lift(ep) do e
        [Point3f(0.0,0.0,0.0), Point3f(e[1])]
    end
    aphelion = lift(ep) do e
        [Point3f(0.0,0.0,0.0), Point3f(e[180])]
    end
    # TODO: find out how to draw 90d lines from sun
    semiminor = lift(ep) do e
        [Point3f(e[90]), Point3f(e[270])]
    end
    f2 = lift(ecc,lpx) do e, lpx
        Point3f(2e*cos(lpx), 2e*sin(lpx), 0.0)
    end
    true_anomaly = lift(lpx, lon) do lpx, sl
        # note that 180 must be added to lpx! (or subtracted)
        mod(sl - rad2deg(lpx) + 180, 360)
    end
    # probably incorrect?!? not used yet
    mean_anomaly = lift(true_anomaly, ecc) do nu, e
        M = atan(sqrt(1-e^2)*sind(nu), e + cosd(nu)) - e * (sqrt(1-e^2)*sind(nu))/(1+e*cosd(nu))
        rad2deg(M)
    end
    # now correctly traces nu
    earth_location = lift(lpx,true_anomaly, ecc) do lpx, nu, e
        # distance
        rho = (1 - e^2) / (1 + e * cosd(nu))
        x = rho * cosd(nu)
        y = rho * sind(nu)
        z = zero(x)
        earth_point = Point3f(x, y, z)
        # rotate around by lpx + 180 so that x-axis corresponds to VE
        Rz = RotMatrix3(AngleAxis(mod(lpx + pi, 2pi), 0,0,1))
        rot_earth = Point3f(Rz * earth_point)
        return rot_earth
    end
    # the equinoxes and solstices
    ve = [Point3f(0.0,0.0,0.0), Point3f(2.0,0.0,0.0)]
    ws = [Point3f(0.0,-2.0,0.0), Point3f(0.0,0.0,0.0)]
    ae = [Point3f(-2.0,0.0,0.0), Point3f(0.0,0.0,0.0)]
    ss = [Point3f(0.0,0.0,0.0), Point3f(0.0,2.0,0.0)]
    o = Point3f(0.0, 0.0, 0.0) # the origin
    n = Vec3f(0.,0.,1.) # the orbit normal
    lan = 255.5932487735715 # SNVec.df.lan[1] + SNVec.OMT
    lanv = Vec3f(cosd(lan), sind(lan), 0.0)
    # lanv = lift(lpx) do lpx
    #     Rz = RotMatrix3(AngleAxis(mod(lpx + 180.0, 360)*pi/180, 0,0,1))
    #     Rz * Vec3f(cosd(lan), sind(lan), 0.0)
    # end
    s = lift(obl) do o # the spin vector
        # the spin vector s is rotated around the ascending node
        # by the precession angle ϕ
        # phi is 0 at t0
        # the ascending node is lan (longitude of ascending node)
        # we first transform lan by OMT to get it into ecliptic coords
        # Rz1 = RotMatrix3(AngleAxis(lan / SNVec.R2D, 0,0,1))
        # Rz1 * Point3f(1.0, 0.0, 0.0)
        # Rz = RotMatrix3(AngleAxis(mod(lpx + 180.0, 360)*pi/180, 0,0,1))
        # Rz = RotMatrix3(AngleAxis(0.5pi, 0,0,1))
        # Rz = [[0.0, -1.0, 0.0] [1.0, 0.0, 0.0] [0.0, 0.0, 1.0]]
        # just experimentally
        # or = 90.0 - o
        # opoint = Vec3f(cosd(or), 0.0, sind(or))
        # Vec3f(Rz * opoint)
        Vec3f(0., -sin(o), cos(o))
    end
    # draw obliquity from Earth point
    s2 = lift(s, earth_location) do s, el
        [Point3f(el), Point3f(el + ef * s)]
    end
    # true solar longitude
    lambda = lift(lon) do tl
        Vec3f(cosd(tl),sind(tl),0.)
    end
    rs = lift(earth_location,lambda) do el,lambda
        [Point3f(el), Point3f(el + ef * lambda)] # + or -?
    end
    # draw r from Earth origin
    # r = lift(lat, L) do phi, L
    #     Rx = RotMatrix3(AngleAxis(phi*pi/180, -1.0, 0.0, 0.0))  # Rotation around x-axis
    #     Rz = RotMatrix3(AngleAxis(L*pi/180+pi/2, 0.0, 0.0, 1.0))  # Rotation around z-axis
    #     Vec3f(Rz * (Rx * Point3f(0.0,-ef,0.0)))
    # end
    # r2 = lift(r, earth_location) do r, el # as a line
    #     [el, el + ef * r]
    # end
    n2 = lift(earth_location) do el
        [Point3f(el), Point3f(el + ef * n)]
    end

    # draw orbit
    # ax = Axis3(fins[1,1], aspect = (1,1,1), # for single image
    # ax = Axis3(fins[1,3], aspect = (1,1,1),
    # limits=(-2.0,2.0,-2.0,2.0,-2.0,2.0),
    # # switch to top view
    # perspectiveness = 0.2,
    # azimuth = -pi/2, # front view
    # elevation = pi/2 # top view
    # )
    ax = LScene(
        fins[2:3,2:3],
        show_axis = false,
        # scenekw = (; lights = [PointLight(RGBf(5,5,5), Point3f(0,0,0))])) # more accurate?
        # but it attenuates too fast, so instead point a directional light at the Earth so it seems to be super far
        scenekw = (; lights = [DirectionalLight(RGBf(1,1,1), earth_location)]))
    # lines!(ax, ws, color = :lightblue, label = "Winter Solstice")
    # lines!(ax, ae, color = :brown, label = "Autumnal Equinox")
    # lines!(ax, ss, color = :gold, label = "Summer Solstice")
    # lines!(ax, [o, n], color = :green, label = L"Orbit normal $\vec{n}$")
    # lines!(ax, @lift([o, $s]), color= :red, label = L"Spin $\vec{s}$")


    # coordinate system with 3d arrows
    arrow_tip_length = 0.06
    arrow_length = 1 - arrow_tip_length # arrows of length 1
    opts = (linewidth = 0.01,
            arrowsize = Vec3f(0.05, 0.05, arrow_tip_length),
            # no shading so they don't show up black when backlit by our
            # directional light
            # fxaa = true, faoo = true# ,
            shading = NoShading
            )
    arrows!(ax, [o], [Vec3f(arrow_length,0,0)]; color = :orange, opts...)
    arrows!(ax, [o], [Vec3f(0,arrow_length,0)]; color = :gold, opts...)
    arrows!(ax, [o], [Vec3f(0,0,arrow_length)]; color = :green, opts...)
    # arrows!(ax, [o], @lift([Point3f(arrow_length * $s)]), color= :red, opts...)
    arrows!(ax, [o], @lift([arrow_length * $s]), color = :red,
            linewidth = 0.01, arrowsize = Vec3f(0.05, 0.05, 0.06),
            shading = NoShading)

    lines!(ax, ep, color=:blue, linewidth=2, label = "Earth Orbit")
    lines!(ax, perihelion, color=:blue, label = "Perihelion")
    lines!(ax, aphelion, color=:blue, linestyle = :dot, label = "Aphelion")
    lines!(ax, @lift([o, $lambda]), color= :gold, label = L"True Solar longitude $\lambda_\odot$")
    # lines!(ax, semiminor, color=:blue, linestyle = :dot, label = "Semiminor")
    # lines!(ax, ve, color = :orange, linewidth = 3, label = "Vernal Equinox")
    # this is probably incorrect??
    # lines!(ax, [o, lanv], color = :black, label = L"Longitude of Ascending Node $\Omega$")
    # arc!(ax, Point2f(0.), 0.3, 0., lan/180*pi)

    # add arc for obliquity
    arc!(ax, Point2f(0.), arrow_length - 0.1, 0.,
         # @lift(-$obl/180*pi),
         @lift(-$obl),
         color = :red,
         transformation =
             (;rotation = Makie.rotation_between(Point(1.0, 0.0, 0.0), Point(0.0, cosd(90), sind(90)))),
         label = "Obliquity ϵ")
    # # add obliquity extremes
    # lines!(ax, [o, Vec3f(RotMatrix3(AngleAxis(minimum(x1_1.obliquity), 1.0, 0.0, 0.0)) * n)], color= :red)
    # lines!(ax, [o, Vec3f(RotMatrix3(AngleAxis(maximum(x1_1.obliquity), 1.0, 0.0, 0.0)) * n)], color= :red)
    # plot an arc around the longitude of perihelion
    arc!(ax, Point2f(0), 0.3, 0, @lift(mod(pi + $lpx, 2pi)), color = :purple, L"\bar\omega")
    arc!(ax, Point2f(0), 0.4, 0, @lift(deg2rad($lon)), color = :gold)

    # points
    scatter!(ax, o, color=:gold, markersize=30, label= L"Sun $\odot$")
    scatter!(ax, f2, color=:red, markersize=5, label=L"f_2")
    scatter!(ax, earth_location, color=:blue, markersize=10, label=L"Earth $\oplus$") # might be corrections now?

    # draw them to scale XD haha
    earth_mesh = GeometryBasics.uv_normal_mesh(Sphere(Point3f(0., 0., 0.), 1.0))

    # textures from https://www.solarsystemscope.com/textures/
    earth_map = #load(Makie.assetpath("earth.png"))
        # load(download("http://upload.wikimedia.org/wikipedia/commons/5/56/Blue_Marble_Next_Generation_%2B_topography_%2B_bathymetry.jpg"))
        load(download("https://www.solarsystemscope.com/textures/download/2k_earth_daymap.jpg"))

    earth_clouds = load(download("https://www.solarsystemscope.com/textures/download/2k_earth_clouds.jpg"))

    cloud_alpha = float.(Gray.(earth_clouds))
    cloud_rgba = RGBA.(1.0,1.0,1.0, cloud_alpha)
    blend(fg::RGBA, bg::RGBA) = begin
        α = alpha(fg)
        RGBA(
            red(fg) * α + red(bg) * (1 - α),
            green(fg) * α + green(bg) * (1 - α),
            blue(fg) * α + blue(bg) * (1 - α),
            1.0
        )
    end
    earth_texture = blend.(cloud_rgba, RGBA.(earth_map, 1.0))

    # earth_specular = Float32.(Gray.(load(download("https://www.solarsystemscope.com/textures/download/2k_earth_specular_map.tif"))))
    # earth_normal # maybe someday?
    stars_texture = load(download("https://www.solarsystemscope.com/textures/download/2k_stars_milky_way.jpg"))
    sun_texture = load(download("https://www.solarsystemscope.com/textures/download/2k_sun.jpg"))


    # earth_latitude = lift(lat) do lat
    #     GeometryBasics.Cylinder(Point3f(0.,0.,0.), Point3f(0.,0.,cosd(lat)), sind(lat))
    # end

    INCT = 7.155
    OMT = 75.594

    Sun = meshscatter!(ax, o,
                       marker = earth_mesh,
                       color = sun_texture,
                       rotation = Vec3f(0., sind(INCT), cosd(OMT)),
                       markersize = @lift($scale_factor * 0.00465047),
                       clip_planes = @lift([Plane3f(-$(ax.scene.camera.view_direction), 0.1)]),
                       # markersize = 0.00465047 * 10
                       shading = NoShading)
    Earth = meshscatter!(ax, earth_location,
                         marker = earth_mesh,
                         color = earth_texture,
                         markersize = @lift($scale_factor * 4.26354E-5),
                         # specular = earth_specular,
                         specular = .1,
                         shininess = 0.,
                         # TODO: also rotate around s by L?
                         rotation = s)
    # TODO: plot earth latitude in 3d thing
    # meshscatter!(ax, earth_location, marker = earth_latitude,
    #              rotation = s, markersize = @lift($scale_factor * 4.26354E-5),
    #              shininess = 0., specular = 0.)
    # arc!(ax, earth_location, )

    # lines from Earth
    # lines!(ax, n2, color = :green)
    lines!(ax, s2, color= :red)
    lines!(ax, @lift([Point3f($earth_location), Point3f($earth_location - ef * $s)]), color= :red)
    # lines!(ax, rs, color = :gold)
    # lines!(ax, r2, color= :darkblue, label = "r vector")
    # band!(ax, r3, color= :darkblue, label = "r vector")
    # lines!(ax, @lift($r3.lower), @lift($r3.upper), color = :darkblue, label = "r vector")

    # text!(ax, "Sun", position= Point3f(0,0,0), fontsize=22, color=:gold)
    text!(ax, L"f_2", position= f2, fontsize=22, color=:red)
    # text!(ax, "Earth", position = earth_location, fontsize=22, color=:blue)
    text!(ax, "VE",
          color = :orange, position = (1.1, 0, 0), fontsize=22)
    text!(ax, "SS",
          color = :gold, position = (0, 1.1, 0), fontsize=22)
    text!(ax, "P",
          color = :blue, position= @lift($perihelion[2] .* 1.1), fontsize=22)
    text!(ax, L"\vec{n}",
          color = :green,
          position= @lift($n * arrow_length * 1.05), fontsize=22)
    text!(ax, L"\vec{s}",
          color = :red,
          position= @lift($s * arrow_length * 1.05), fontsize=22)
    text!(ax, L"\epsilon", color = :red, position = @lift(($n + ($s - n)/2) * arrow_length), fontsize = 22)
    text!(ax, L"\nu", color = :gold,
          position = @lift(0.5 .* Vec3f(cosd($lon/2), sind($lon/2), 0.)),
          fontsize = 22)
    text!(ax, L"\bar\omega", color = :purple,
          position = @lift(0.4 .*
              # this is a little annoying because of the whole +180° thing.
              Vec3f(
                  mod(pi + $lpx, 2pi) < pi ? cos($lpx/2-0.5pi) : -cosd($lpx/2-90),
                  mod(pi + $lpx, 2pi) < pi ? sin($lpx/2-0.5pi) : -sin($lpx/2-0.5pi),
                  0.)),
          fontsize = 22)

    cam3d!(ax.scene,
           # projectiontype = :orthographic,
           )
    cameracontrols(ax.scene).settings.center[] = false
    Stars = meshscatter!(ax, o,
                         marker = earth_mesh,
                         color = stars_texture,
                         rotation = Vec3f(0., sind(INCT), cosd(OMT)),
                         markersize = 1000,
                         clip_planes =  @lift([Plane3f($(ax.scene.camera.view_direction), -0.0001)]),
                         shading = NoShading)

    # TODO: figure out how to set initial camera postition
    # below doesn't seem to do anything...
    # update_cam!(finf.content[4].scene, Vec3d(1,1,1), Vec3d(0), Vec3d(0,0,1))
    # update_cam!(finf.content[4].scene,
    #             cameracontrols(finf.content[4].scene),
    #             2pi, pi, pi/20)
    # update_cam!(finf.content[4].scene)
    fins
    # cameracontrols(finf.content[5].scene).eyeposition[] = Vec3f(2)
    # cameracontrols(finf.content[5].scene).lookat[] = Vec3f(0.)
    # cameracontrols(finf.content[5].scene).upvector[] = Vec3f(0., 0., 1.)
    # cameracontrols(finf.content[5].scene).fov[] = 1
    # update_cam!(finf.content[5].scene)
    # projectiontype = Makie.Orthographic,
    #                zoom_shift_lookat = false,
    #                eyeposition = Vec3d(1.0),
    #                center = false
end
