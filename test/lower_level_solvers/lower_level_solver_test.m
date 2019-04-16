
function tests = lower_level_solver_test
    tests = functiontests(localfunctions);
end

function testROFChambollePock(testCase)
    x = [1 2 3;1 0 3;2 4 1];
    param.alpha = 0.1;
    actSolution = solve_rof_cp_single_gaussian(x,param);
    expSolution = [1.08926538199073 1.87288842165649 2.89629345403946;1.02551762730455 0.277061538550202 2.89629353405119;2.02001516425054 3.72857047338053 1.19409440477625];
    verifyEqual(testCase,actSolution,expSolution,'AbsTol',1e-5);
end