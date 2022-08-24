using GLMakie
using Makie
using LinearAlgebra:norm


function add_source(densities::Matrix{Float32}, forces::Matrix{Float32}, dt::Float32)
    # add external sources to the density field
    densities .+= dt .* forces
end


function diffuse!(densities::Matrix{Float64}, diff::Float64, dt::Float64)
    densities0 = deepcopy(densities)
    a = dt * diff * prod(size(densities))

    # 20 itertion Gauss-Seidl relaxation
    for iter in 1:10
        for row in 2:size(densities)[1]-1
            for col in 2:size(densities)[2]-1
                densities[row, col] = (densities0[row, col] + a * (densities[row-1, col] + densities[row+1, col] + densities[row, col+1] + densities[row, col-1])) / (1 + 4*a)
            end
        end
    end
end


function advect!(x_v::Matrix{Float64}, y_v::Matrix{Float64}, densities::Matrix{Float64}, dt::Float64, SIZE, DIFF)
    dt=1
    d0 = copy(densities)    # old unchanging densities
    for x in 2:SIZE-1
        for y in 2:SIZE-1
            if x==SIZE/2 && y == SIZE/2
                continue
            end
            old_x = round(Int, x - x_v[y, x] * dt)
            
            if old_x < 1
                old_x = x
            end
            if old_x >= SIZE
                old_x = x
            end
            old_y = round(Int, y - y_v[y, x] * dt)
            if old_y < 1
                old_y = y
            end
            if old_y >= SIZE
                old_y = y
            end
            densities[y, x] = d0[old_y, old_x]
        end
    end
end

function step!(x_v::Matrix{Float64}, y_v::Matrix{Float64}, densities::Matrix{Float64}, SIZE::Int, DIFF::Float64, ∇t::Float64, total::Float64)
    diffuse!(densities, DIFF, ∇t)
    t = time()
    advect!(x_v, y_v, densities, ∇t, SIZE, DIFF)    # advect densities
    x_v0, y_v0 = copy.([x_v, y_v])
    advect!(x_v0, y_v0, x_v, 1., SIZE, DIFF)         # advect x velocities
    t = time()
    advect!(x_v0, y_v0, y_v, 1., SIZE, DIFF)         # advect y velocities
    normalize!(densities, total, SIZE)
end

function gen_velocities(SIZE)
    x_velocities = zeros(Float64, (SIZE, SIZE))
    y_velocities = zeros(Float64, (SIZE, SIZE))
    for x in 1:SIZE
        for y in 1:SIZE
            if !(x == SIZE/2 && y == SIZE/2)
                vec = Float64[SIZE/2 - x, SIZE/2 - y] / norm(Float64[SIZE/2 - x, SIZE/2 - y]) * 5
                rot = [cosd(90) -sind(90)
                       sind(90) cosd(90)]
                vec = rot * vec
            else
                vec = [0, 0]
            end
            x_velocities[y, x] = vec[1]
            y_velocities[y, x] = vec[2]
        end
    end
    return x_velocities, y_velocities
end

function gen_line_densities(SIZE, width)
    densities = zeros(Float64, (SIZE, SIZE))
    start_x, end_x = round(Int, SIZE/2) - round(Int, width/2), round(Int, SIZE/2) + round(Int, width/2)
    for i in start_x:end_x
        densities[:, i] .= 1 / ( (round(Int, abs(i-SIZE/2))) + 1 )
    end
    return densities
end

function normalize!(field::Matrix{Float64}, total::Float64, SIZE::Int)
    current_total = sum(field)
    for row in 1:SIZE
        for col in 1:SIZE
            field[row, col] = (field[row, col] / current_total) * total
        end
    end
end


function format_velocities(x_v, y_v, SIZE)
    points = [Point2(row, col) for col in 1:10:SIZE for row in 1:10:SIZE]
    vectors = [[y_v[i], x_v[i]] for i in 1:10:SIZE^2] .* 1
    strength = [sqrt(v[1]^2 + v[2]^2) for v in vectors]
    return points, vectors, strength
end

function add_density(x, y, densities, amount, radius)
    for row in round(Int, x-radius/2):round(Int, x+radius/2)
        for col in round(Int, y-radius/2):round(Int, y+radius/2)
            dist = √( (row-x)^2 + (col-x)^2 ) + 1
            densities[row, col] += (1 / dist) * amount
        end
    end
end


function main(iterations::Int)
    SIZE = 200
    DIFF = 0.00003
    ∇t = 0.01
    # densities = randn(Float64, (SIZE, SIZE)) .|> abs
    densities = gen_line_densities(SIZE, 20)
    TOTAL = sum(densities)
    x_velocities, y_velocities = gen_velocities(SIZE)
    obs = Observable(densities)
    fig = Figure()
    ax = Axis(fig[1, 1], title="Density Distribution")
    ax.aspect = AxisAspect(1)
    vector_ax = Axis(fig[1, 2], backgroundcolor="black", title = "Vector Field")
    vector_ax.aspect = AxisAspect(1)


    deactivate_interaction!(ax, :rectanglezoom)
    pressed = Observable(false)
    on(events(ax.scene).mousebutton, priority=0) do event
        if event.action == Makie.Mouse.press
            pressed[] = true
        end
        if event.action == Makie.Mouse.release
            pressed[] = false
        end
    end
    on(events(ax.scene).mouseposition) do event
        if pressed[]
            x, y = mouseposition(ax.scene)
            println(x, " ", y)
            add_density(x, y, densities, 1, 10)
            notify(obs)
        end
    end
    
    
    
    points, vectors, strength = format_velocities(x_velocities, y_velocities, SIZE)
    obs_points = Observable(points)
    obs_vectors = Observable(vectors)
    obs_strength = Observable(strength)
    screen = display(fig)
    resize!(screen, 1800, 1000)
    heatmap!(ax, obs)
    arrows!(vector_ax, obs_points, obs_vectors, arrowcolor=obs_strength, linecolor="yellow", arrowsize=5, lengthscale=1)
    t = time()
    sleep(1)
    # for iter in 1:iterations
    #     step!(x_velocities, y_velocities, densities, SIZE, DIFF, ∇t, TOTAL)
    #     points, vectors, strength = format_velocities(x_velocities, y_velocities, SIZE)
    #     obs_points[] = points
    #     obs_vectors[] = vectors
    #     obs_strength[] = strength
    #     notify(obs)
    #     sleep(0.05)
    # end

    record(fig, "fluid_sim", 1:iterations, framerate=24) do frame
        step!(x_velocities, y_velocities, densities, SIZE, DIFF, ∇t, TOTAL)
        points, vectors, strength = format_velocities(x_velocities, y_velocities, SIZE)
        obs_points[] = points
        obs_vectors[] = vectors
        obs_strength[] = strength
        notify(obs)
    end
end
