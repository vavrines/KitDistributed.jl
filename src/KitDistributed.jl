module KitDistributed

export KD

#using KitBase
using MPI

include("geometry.jl")
include("update.jl")

const KD = KitDistributed

end
