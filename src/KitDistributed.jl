module KitDistributed

export KD

using KitBase
using MPI

include("geometry.jl")

const KD = KitDistributed

end
