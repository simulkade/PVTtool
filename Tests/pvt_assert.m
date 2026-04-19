function pvt_assert(condition, test_name)
% pvt_assert  Throw a named error when condition is false.
%
%   pvt_assert(condition, test_name)
%
%   Used inside test files. If condition evaluates to false (or 0),
%   throws an error with identifier 'PVTtest:failed' so the runner can
%   catch it and record a failure without halting the whole suite.
%
% SEE ALSO: run_all_tests

if ~all(condition)
    error('PVTtest:failed', 'FAIL: %s', test_name);
end
