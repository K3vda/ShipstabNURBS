using ShipstabNURBS
using Test



@testset "FindSpan" begin
    U = [0 0 0 0 2/5 3/5 3/5 1 1 1 1];
    t1,s1 = FindSpan(0,U)
    t2,s2 = FindSpan(2/5,U)
    t3,s3 = FindSpan(1,U)
  @test t1 == 4
  @test s1 == 3
  @test t2 == 7
  @test s2 == 2
  @test t3 == 7
  @test s3 == 0
end

@testset "BasisFuns" begin
    P = [0 1 3 5 6 8 9;0 2 4 2 0 0 3];
    U = [0 0 0 0 2/5 3/5 3/5 1 1 1 1];
    n = size(P,2) - 1;
    m = length(U) - 1;
    p = m-n-1;

    u = range(0,1,length=100);
    funs = zeros(length(u),n+1);

    for i = 1:length(u)
        funs[i,:] = BasisFuns(n,p,u[i],U);
    end
    @test sum(funs) == 100.0
end
