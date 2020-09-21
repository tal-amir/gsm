% This script tests the accuracy of our code to compute the GSM function.

% Each task_id generates a different instance.
% We used task_id = 1,...,200 in our evaluation.
if ~exist('task_id', 'var')
    task_id = 1;
end

d = 100000;
k = [10, 50, 100, 200];
gamma = [1e-20,1e-10, 1e-5, 1e-2, 0.2:0.2:1, 2:2:10, 1e2, 1e5, 1e10, 1e20];

is_single = true;

vec_type = {'uniform', 'gaussian'};
seed = task_id;

result = gsm_test(d, k, gamma, vec_type, is_single, seed);
