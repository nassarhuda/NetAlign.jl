@testset "networkalignment" begin
    S,w,li,lj,A,B,L = load_netalign_problem("example-overlap")
    @test sparse(li,lj,w)==L
    @test nnz(S)==30
    @test size(S) == (12,12)

    S2,lw,ei,ej = netalign_setup(A,B,L)
    @test S2 == S
    @test li == ei
    @test lj == ej
    @test lw == w
end
