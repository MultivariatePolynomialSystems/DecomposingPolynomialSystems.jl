using Oscar, LinearAlgebra
include("../../../../src/utils/utils_cpu.jl")

const SUBSCRIPTS = ['₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉']
const SUBSCRIPT_MAP = Dict([first(string(i)) => SUBSCRIPTS[i+1] for i = 0:9])
const SUBSCRIPT_TO_INT_MAP = Dict([SUBSCRIPTS[i+1] => i for i = 0:9])
map_subscripts(index) = join(SUBSCRIPT_MAP[c] for c in string(index))

function buildvar(var::String, indicies::Int...)
    return var*join(map_subscripts.(indicies), "₋")
end

function buildvars(var::String, sizes...)
    reshape(map(i -> buildvar(var, i...), Iterators.product(sizes...)), prod([length(size) for size in sizes]))
end
