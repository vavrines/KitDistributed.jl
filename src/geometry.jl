function position_id(nproc::T, nx_domains::T, ny_domains::T) where {T<:Integer}
    id_pos = Array{Int,2}(undef, nx_domains, ny_domains)
    for id = 0:nproc-1
        n_row = (id + 1) % nx_domains > 0 ? (id + 1) % nx_domains : nx_domains
        n_col = ceil(Int, (id + 1) / nx_domains)
        id_pos[n_row, n_col] = id
    end

    return id_pos
end
