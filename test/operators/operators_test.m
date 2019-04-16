
function tests = operators_test
    tests = functiontests(localfunctions);
end

function testFinDiffOperator_val(testCase)
    x = [1 2 3;1 0 3;2 4 1];
    op = FinDiffOperator(size(x),"fn");
    actSolution = op.val(x);
    expSolution = zeros(size(x,1),size(x,2),2);
    expSolution(:,:,1)=[1 1 0;-1 3 0;2 -3 0];
    expSolution(:,:,2)=[0 -2 0;1 4 -2;0 0 0];
    verifyEqual(testCase,actSolution,expSolution);
end

function testFinDiffOperator_conj(testCase)
    op = FinDiffOperator([3,3],"fn");
    y = zeros(3,3,2);
    y(:,:,1)=[1 1 0;-1 3 0;2 -3 0];
    y(:,:,2)=[0 -2 0;1 4 -2;0 0 0];
    actSolution = op.conj(y);
    expSolution = [-1 2 1;0 -10 5;-1 9 -5];
    verifyEqual(testCase,actSolution,expSolution);
end