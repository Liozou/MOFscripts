using CrystalNets
using Test, Random
using PeriodicGraphs
using StaticArrays
import Base.Threads

const cifs = joinpath(@__DIR__, "cif")
const crystalnetsdir = normpath(@__DIR__, "..")

@testset "Archive" begin
    Threads.@threads for (genome, id) in collect(CrystalNets.CRYSTAL_NETS_ARCHIVE)
        (id == "sxt" || id == "llw-z") && continue # special case for these known unstable nets
        @test reckognize_topology(topological_genome(PeriodicGraph(genome))) == id
    end
end

@testset "Module" begin
    Threads.@threads for target in
        ["pcu", "afy, AFY", "apc, APC", "bam", "bcf", "cdp", "cnd", "ecb", "fiv",
         "ftd", "ftj", "ins", "kgt", "mot", "moz", "muh", "pbz", "qom", "sig",
         "sma", "sod-f", "sod-h", "utj", "utp"]
       @info "Testing $target"
       graph = PeriodicGraph(CrystalNets.REVERSE_CRYSTAL_NETS_ARCHIVE[target])
       n = PeriodicGraphs.nv(graph)
       for k in 1:50
           r = randperm(n)
           offsets = [SVector{3,Int}([rand(-3:3) for _ in 1:3]) for _ in 1:n]
           graph = swap_axes!(offset_representatives!(graph[r], offsets), randperm(3))
           @test reckognize_topology(topological_genome(graph)) == target
       end
   end

   @test_broken reckognize_topology(topological_genome(PeriodicGraph(CrystalNets.REVERSE_CRYSTAL_NETS_ARCHIVE["sxt"]))) == "sxt"
   @test_broken reckognize_topology(topological_genome(PeriodicGraph(CrystalNets.REVERSE_CRYSTAL_NETS_ARCHIVE["llw-z"]))) == "llw-z"

   @info """The following warnings about guessing bonds are expected."""
   @test reckognize_topology(topological_genome(CrystalNet(parse_chemfile(joinpath(cifs, "Moganite.cif"))))) == "mog"
end


@testset "Executable" begin
    safestdout = deepcopy(stdout)
    safeARGS = deepcopy(ARGS)
    piperead, _ = redirect_stdout()
    empty!(ARGS)

    push!(ARGS, "-g", "3   1 2  0 0 0   1 2  0 0 1   1 2  0 1 0   1 2  1 0 0")
    @test CrystalNets.julia_main() == 0
    @test readuntil(piperead, '\n') == "dia"

    empty!(ARGS)
    path = joinpath(cifs, "ABW.cif")
    push!(ARGS, path)
    @test CrystalNets.julia_main() == 0
    @test readuntil(piperead, '\n') == "sra, ABW"

    empty!(ARGS)
    path = joinpath(cifs, "ABW.cif")
    push!(ARGS, "-a", CrystalNets.arc_location*"rcsr.arc", path)
    @test CrystalNets.julia_main() == 0
    @test readuntil(piperead, '\n') == "sra"

    empty!(ARGS)
    path = joinpath(cifs, "RRO.cif")
    push!(ARGS, path)
    @test CrystalNets.julia_main() == 0
    @test readuntil(piperead, '\n') == "RRO"

    empty!(ARGS)
    path = joinpath(cifs, "RRO.cif")
    push!(ARGS, "-a", CrystalNets.arc_location*"rcsr.arc", path)
    @test CrystalNets.julia_main() == 1
    @test readuntil(piperead, '\n') == "UNKNOWN"

    empty!(ARGS)
    path = joinpath(cifs, "HKUST-1.cif")
    push!(ARGS, "-c", "mof", path)
    @test CrystalNets.julia_main() == 0
    @test readuntil(piperead, '\n') == "tbo"

    empty!(ARGS)
    path = joinpath(cifs, "HKUST-1_2.cif")
    push!(ARGS, "-c", "mof", path)
    @test CrystalNets.julia_main() == 0
    @test readuntil(piperead, '\n') == "tbo"

    empty!(ARGS)
    path = joinpath(cifs, "Diamond.cif")
    push!(ARGS, "-c", "atom", path)
    @test CrystalNets.julia_main() == 0
    @test readuntil(piperead, '\n') == "dia"

    empty!(ARGS)
    push!(ARGS, "--help")
    @test CrystalNets.julia_main() == 0
    @test startswith(readuntil(piperead, '\n'), "usage: CrystalNets")
    @test occursin("CRYSTAL_FILE", readuntil(piperead, '\n'))
    @test occursin("Form B", readuntil(piperead, '\n'))
    @test occursin("Form C", readuntil(piperead, '\n'))
    @test isempty(readuntil(piperead, '\n'))
    @test isempty(readuntil(piperead, '\n'))
    @test readuntil(piperead, '\n') == "Automatic reckognition of crystal net topologies."

    redirect_stdout(safestdout)
    empty!(ARGS)
    append!(ARGS, safeARGS)
end
