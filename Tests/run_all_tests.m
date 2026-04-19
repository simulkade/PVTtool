% run_all_tests  Execute the full PVTtool test suite and print a summary.
%
%   Run from the project root after calling PVTinitialize:
%     PVTinitialize()
%     cd Tests
%     run_all_tests
%
%   Each test file is a function that runs assertions and returns nothing.
%   A test file passes if it completes without throwing an error.

test_files = { ...
    'test_classes', ...
    'test_eos', ...
    'test_flash_vle', ...
    'test_stability', ...
    'test_stability_lle', ...
    'test_residual_props', ...
    'test_saturation' ...
};

passed = 0;
failed = 0;
failures = {};

fprintf('\n========================================\n');
fprintf('  PVTtool Test Suite\n');
fprintf('========================================\n');

for i = 1:length(test_files)
    name = test_files{i};
    try
        feval(name);
        fprintf('  PASS  %s\n', name);
        passed = passed + 1;
    catch err
        fprintf('  FAIL  %s\n        %s\n', name, err.message);
        failed = failed + 1;
        failures{end+1} = name; %#ok<AGROW>
    end
end

fprintf('----------------------------------------\n');
fprintf('  %d passed, %d failed\n', passed, failed);
fprintf('========================================\n\n');

if failed > 0
    fprintf('Failed tests: %s\n\n', strjoin(failures, ', '));
end
