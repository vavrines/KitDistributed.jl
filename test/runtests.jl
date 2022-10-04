using MPI, KitDistributed

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
nproc = MPI.Comm_size(comm)
const root = 0

ps = (
    T0 = -10.0,
    T1 = 10.0,
    nx = 256,
    ny = 256,
    nxd = 1, #4,
    nyd = 1, #2,
)

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

KD.domain_coords!(xindex, yindex, nx1, ny1, ps.nxd, ps.nyd)

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

MPI.Finalize()
