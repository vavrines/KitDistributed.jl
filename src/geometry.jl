"""
Calculate the processor ID linked to the given domain

# Arguments
* ``nproc``: number of processors
* ``nxd, nyd``: number of domains in x and y directions
"""
function position_id(nproc, nxd, nyd)
    ids = Array{Int,2}(undef, nxd, nyd)
    for id = 0:nproc-1
        n_row = (id + 1) % nxd > 0 ? (id + 1) % nxd : nxd
        n_col = ceil(Int, (id + 1) / nxd)
        ids[n_row, n_col] = id
    end

    return ids
end


"""
Calculate the domain index linked to the given processor ID

# Arguments
* ``id``: processor ID
* ``nxd, nyd``: number of domains in x and y directions
"""
function id_position(id, nxd, nyd)
    n_row = (id + 1) % nxd > 0 ? (id + 1) % nxd : nxd
    n_col = ceil(Int, (id + 1) / nxd)

    return n_row, n_col
end


"""
Find the ranks of neighboring subdomains
"""
function neighbor_id(idx, idy, nxd, nyd, posid)
    neighbor_N = idx + 1 <= nxd ? idx + 1 : -1
    neighbor_S = idx - 1 > 0 ? idx - 1 : -1
    neighbor_E = idy + 1 <= nyd ? idy + 1 : -1
    neighbor_W = idy - 1 > 0 ? idy - 1 : -1

    neighbors = Array{Int,1}(undef, 4)
    neighbors[1] = neighbor_N >= 0 ? posid[neighbor_N, idy] : -1
    neighbors[2] = neighbor_S >= 0 ? posid[neighbor_S, idy] : -1
    neighbors[3] = neighbor_E >= 0 ? posid[idx, neighbor_E] : -1
    neighbors[4] = neighbor_W >= 0 ? posid[idx, neighbor_W] : -1

    return neighbors
end


"""
Find the ranks of neighboring subdomains
"""
function domain_neighbor(my_id, nproc, nxd, nyd)
    pos = position_id(nproc, nxd, nyd)
    idx, idy = id_position(my_id, nxd, nyd)

    return neighbor_id(idx, idy, nxd, nyd, pos)
end


"""
Calculate starting & ending indices of subdomains
"""
function domain_coords!(
    xindex::T,
    yindex::T,
    xcell,
    ycell,
    nx_domains,
    ny_domains,
) where {T}

    nproc = nx_domains * ny_domains
    xs = Vector{Int}(undef, nproc)
    xe = Vector{Int}(undef, nproc)
    ys = Vector{Int}(undef, nproc)
    ye = Vector{Int}(undef, nproc)

    ys[1:nx_domains] .= 2
    ye[1:nx_domains] .= ys[1:nx_domains] .+ ycell .- 1
    for i = 1:ny_domains-1
        ys[i*nx_domains+1:(i+1)*nx_domains] =
            ys[(i-1)*nx_domains+1:i*nx_domains] .+ ycell .+ 2
        ye[i*nx_domains+1:(i+1)*nx_domains] =
            ys[i*nx_domains+1:(i+1)*nx_domains] .+ ycell .- 1
    end

    for i = 1:ny_domains
        xs[(i-1)*nx_domains+1] = 2
        xe[(i-1)*nx_domains+1] = xs[(i-1)*nx_domains+1] + xcell - 1
    end
    for i = 1:ny_domains
        for j = 2:nx_domains
            xs[(i-1)*nx_domains+j] = xs[(i-1)*nx_domains+(j-1)] + xcell + 2
            xe[(i-1)*nx_domains+j] = xs[(i-1)*nx_domains+j] + xcell - 1
        end
    end

    xindex[:, 1] .= xs
    xindex[:, 2] .= xe
    yindex[:, 1] .= ys
    yindex[:, 2] .= ye

    return nothing

end
