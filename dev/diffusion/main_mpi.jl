"""
Solving the diffusion equation with a five point finite difference stencil:

|                                    weightx * u[i-1][j]                                    |
|                                                                                           |
| weighty * u[i][j-1]   (diagx * weightx + diagy * weighty) * u[i][j]   weighty * u[i][j+1] |
|                                                                                           |
|                                    weightx * u[i+1][j]                                    |
"""

using Printf
using MPI, KitDistributed

function init_values!(u0, nx, ny, tin, tout)
    u0 .= tout
    for i = 2:nx-1
        for j = 2:ny-1
            u0[i, j] = tin
        end
    end

    return nothing
end
function step!(u0, u, dt, hx, hy, me, xs, ys, xe, ye, k0)
    mep1 = me + 1

    diagx = Float64(-2.0 + hx * hx / (2 * k0 * dt))
    diagy = Float64(-2.0 + hy * hy / (2 * k0 * dt))
    weightx = Float64(k0 * dt / (hx * hx))
    weighty = Float64(k0 * dt / (hy * hy))

    for i = xs[mep1]:xe[mep1]
        for j = ys[mep1]:ye[mep1]
            u[i, j] =
                weightx * (u0[i-1, j] + u0[i+1, j] + u0[i, j] * diagx) +
                weighty * (u0[i, j-1] + u0[i, j+1] + u0[i, j] * diagy)
        end
    end

    res = 0.0
    for j = ys[mep1]:ye[mep1]
        for i = xs[mep1]:xe[mep1]
            diff = u0[i, j] - u[i, j]
            res += diff^2
            u0[i, j] = u[i, j]
        end
    end

    return nothing
end
function write_to_disk(x, x_domains, y_domains, xcell, ycell, temp1, filename)
    f = open(filename, "w")

    # add first x-boundary
    for k = 1:ycell*y_domains+2
        print(f, @sprintf("%15.11f", temp1))
        print(f, "\t")
    end

    # add internal cells + y-boundaries
    c = 0
    for k = 1:x_domains
        for m = 1:xcell
            for i = 1:y_domains
                for j = 1:ycell
                    c += 1
                    if (i == 1 && j == 1)
                        print(f, @sprintf("%15.11f", temp1))
                        print(f, "\t")
                    end
                    print(
                        f,
                        @sprintf(
                            "%15.11f",
                            x[(i-1)*x_domains*xcell*ycell+(k-1)*xcell*ycell+(j-1)*xcell+m]
                        )
                    )
                    print(f, "\t")
                    if (i == y_domains && j == ycell)
                        print(f, @sprintf("%15.11f", temp1))
                        print(f, "\t")
                    end
                end
            end
        end
    end

    # add second x-boundary
    for k = 1:ycell*y_domains+2
        print(f, @sprintf("%15.11f", temp1))
        if (k < ycell * y_domains + 2)
            print(f, "\t")
        end
    end
    close(f)
end

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
nproc = MPI.Comm_size(comm)
const root = 0

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

if rank == root && nproc != ps.nxd * ps.nyd
    throw("number of processes not equal to number of subdomains")
end

dx = 1.0 / ps.nx
dy = 1.0 / ps.ny
dt = 0.25 * min(dx, dy)^2 / ps.D
nxt = ps.nx + 2 + 2 * (ps.nxd - 1)
nyt = ps.ny + 2 + 2 * (ps.nyd - 1)

u0 = zeros(nxt, nyt)
u = zero(u0)

#
nproc = ps.nxd * ps.nyd
#
xindex = zeros(Int, nproc, 2)
yindex = zero(xindex)

nx1 = ps.nx ÷ ps.nxd
ny1 = ps.ny ÷ ps.nyd

ul = zeros(nx1 * ny1)
ug = zeros(ps.nx * ps.ny)

neighbors = KD.domain_neighbor(rank, nproc, ps.nxd, ps.nyd)

rank == root && @show neighbors

KD.domain_coords!(xindex, yindex, nx1, ny1, ps.nxd, ps.nyd)

init_values!(u0, nxt, nyt, ps.T0, ps.T1)

KD.update_ghost!(
    u0,
    neighbors,
    comm,
    rank,
    xindex[:, 1],
    yindex[:, 1],
    xindex[:, 2],
    yindex[:, 2],
    nx1,
    ny1,
)

global step = 0
global t = 0.0
global converged = false

#Starting time
if (rank == root)
    time_init = time()
end

for iter = 1:4000
    global step += 1
    global t += dt

    step!(
        u0,
        u,
        dt,
        dx,
        dy,
        rank,
        xindex[:, 1],
        yindex[:, 1],
        xindex[:, 2],
        yindex[:, 2],
        ps.D,
    )

    KD.update_ghost!(
        u0,
        neighbors,
        comm,
        rank,
        xindex[:, 1],
        yindex[:, 1],
        xindex[:, 2],
        yindex[:, 2],
        nx1,
        ny1,
    )

    MPI.Barrier(comm)
end

if (rank == 0)
    time_final = time()
    #elapsed time
    elapsed_time = time_final - time_init
    println("Elapsed time = ", elapsed_time)
    println("Steps = ", step)
end

xs = xindex[:, 1]
xe = xindex[:, 2]
ys = yindex[:, 1]
ye = yindex[:, 2]

global i = 1
for j = ys[rank+1]:ye[rank+1]
    ul[(i-1)*nx1+1:i*nx1] = u0[xs[rank+1]:xe[rank+1], j]
    global i += 1
end

ug = MPI.Gather(ul, root, comm)

#write to disk
if rank == root
    cd(@__DIR__)
    write_to_disk(ug, ps.nxd, ps.nyd, nx1, ny1, ps.T1, ps.path)
end

MPI.Finalize()
