"""
Solving the diffusion equation with a five point finite difference stencil:

|                                    weightx * u[i-1][j]                                    |
|                                                                                           |
| weighty * u[i][j-1]   (diagx * weightx + diagy * weighty) * u[i][j]   weighty * u[i][j+1] |
|                                                                                           |
|                                    weightx * u[i+1][j]                                    |
"""

using Printf

function init_values!(u0, nx, ny, tin, tout)
    u0 .= tout
    for i = 2:nx-1
        for j = 2:ny-1
            u0[i, j] = tin
        end
    end

    return nothing
end
function step!(u0, u, dt, hx, hy, k0)
    diagx = Float64(-2.0 + hx * hx / (2 * k0 * dt))
    diagy = Float64(-2.0 + hy * hy / (2 * k0 * dt))
    weightx = Float64(k0 * dt / (hx * hx))
    weighty = Float64(k0 * dt / (hy * hy))

    for i = 2:size(u, 1) - 1
        for j = 2:size(u, 2) - 1
            u[i, j] =
                weightx * (u0[i-1, j] + u0[i+1, j] + u0[i, j] * diagx) +
                weighty * (u0[i, j-1] + u0[i, j+1] + u0[i, j] * diagy)
        end
    end

    return nothing
end

ps = (
    T0 = -10.0,
    T1 = 10.0,
    D = 1.0,
    nx = 256,
    ny = 256,
    nxd = 4,
    nyd = 2,
    maxiters = 100000000,
    dt = 0.1,
    tol = 0.1,
    path = "output.dat",
)

dx = 1.0 / ps.nx
dy = 1.0 / ps.ny
dt = 0.25 * min(dx, dy)^2 / ps.D
nxt = ps.nx + 2
nyt = ps.ny + 2

u0 = zeros(nxt, nyt)
u = zero(u0)

init_values!(u0, nxt, nyt, ps.T0, ps.T1)
u .= u0

step!(
        u0,
        u,
        dt,
        dx,
        dy,
        ps.D,
    )

global step = 0
global t = 0.0

for iter = 1:4000
    global step += 1
    global t += dt

    step!(
        u0,
        u,
        dt,
        dx,
        dy,
        ps.D,
    )

    u0 .= u
end

using Plots
contourf(u[2:end-1, 2:end-1])
