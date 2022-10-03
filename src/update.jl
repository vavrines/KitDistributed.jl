"""
Update solution in ghost cells

# Arguments
* ``u``: solution
* ``neighbors``: rank of neighboring subdomains
* ``comm``: communicator
* ``id``: rank
* ``xs, xe``: starting & end indices of subdomains in x
* ``ys, ye``: starting & end indices of subdomains in y
* ``xcell, ycell``: number of cells in the subdomain
"""
function update_ghost!(u, neighbors, comm, id, xs, ys, xe, ye, xcell, ycell)
    # Julia's 1-indexed order of rank
    rankj = id + 1

    rreq =
        MPI.Request[MPI.REQUEST_NULL, MPI.REQUEST_NULL, MPI.REQUEST_NULL, MPI.REQUEST_NULL]

    ghost_boundaries = (
        (xe[rankj] + 1, ys[rankj]:ye[rankj]),
        (xs[rankj] - 1, ys[rankj]:ye[rankj]),
        (xs[rankj]:xe[rankj], ye[rankj] + 1),
        (xs[rankj]:xe[rankj], ys[rankj] - 1),
    )

    is_receiving = Bool[false, false, false, false]

    # send
    neighbors[1] >= 0 &&
        MPI.Isend(u[xe[rankj], ys[rankj]:ye[rankj]], neighbors[1], id + 40, comm)
    neighbors[2] >= 0 &&
        MPI.Isend(u[xs[rankj], ys[rankj]:ye[rankj]], neighbors[2], id + 50, comm)
    neighbors[3] >= 0 &&
        MPI.Isend(u[xs[rankj]:xe[rankj], ye[rankj]], neighbors[3], id + 60, comm)
    neighbors[4] >= 0 &&
        MPI.Isend(u[xs[rankj]:xe[rankj], ys[rankj]], neighbors[4], id + 70, comm)

    recv = Vector{Vector{Float64}}(undef, 4)

    # receive
    if neighbors[1] >= 0
        recv[1] = Vector{Float64}(undef, ycell)
        is_receiving[1] = true
        rreq[1] = MPI.Irecv!(recv[1], neighbors[1], neighbors[1] + 50, comm)
    end
    if neighbors[2] >= 0
        recv[2] = Vector{Float64}(undef, ycell)
        is_receiving[2] = true
        rreq[2] = MPI.Irecv!(recv[2], neighbors[2], neighbors[2] + 40, comm)
    end
    if neighbors[3] >= 0
        recv[3] = Vector{Float64}(undef, xcell)
        is_receiving[3] = true
        rreq[3] = MPI.Irecv!(recv[3], neighbors[3], neighbors[3] + 70, comm)
    end
    if neighbors[4] >= 0
        recv[4] = Vector{Float64}(undef, xcell)
        is_receiving[4] = true
        rreq[4] = MPI.Irecv!(recv[4], neighbors[4], neighbors[4] + 60, comm)
    end

    # block until all active requests in the array reqs are complete
    MPI.Waitall!(rreq)

    for k in eachindex(is_receiving)
        if is_receiving[k]
            u[ghost_boundaries[k][1], ghost_boundaries[k][2]] .= recv[k]
        end
    end

    return nothing
end
