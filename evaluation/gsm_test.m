function result = gsm_test(d, k, gamma, vec_type, is_single, seed)
%function result = gsm_test(d, k, gamma, vec_type, is_single, seed)
%
% Calculates the GSM function tau_{k,gamma}(z) on a randomly generated vector z in R^d
% and compares the result to the ground truth, computed by quad-precision arithmetic.
%
% Input: 
% d: An integer greater than 2 denoting the dimension of z
% k: An integer between 0 and d
% gamma: Real number in [-inf, inf]
% vec_type: (optional) 'uniform' / 'gaussian'. Default: 'uniform' 
%           'uniform' - Each entry of z is drawn randomly i.i.d. uniform in [0,1] 
%           'gaussian' - Each entry is the absolute value of an i.i.d. N(0,1) variable
% is_single: (optional) Default: false
%            Boolean that tells whether to calculate the GSM using single-precision
%            arithmetic for all intermediate values.
% seed: An integer number to initialize the random number generator.
%       Note: If a seed is not provided, it is taken to be constant (hence, the result is
%             constant, not random).
%
% Output: A struct containing measures of precision. The most relevant ones are:
% diff_mu_rel: Relative error of mu: abs(mu-mu_true)/abs(mu_true)
% diff_mu_abs: Absolute error of mu: abs(mu-mu_true)
% diff_mu_ulp: Relative error of mu in units in the last place (ulp)
% diff_theta: Absolute error of theta: l_infinity norm of (theta - theta_true)

if ~exist('seed','var')
    seed = 1234;
    warning('Random seed not provided. Defaulting to %g. Result is not random.', seed);
end

if ~exist('vec_type','var')
    vec_type = 'uniform';
end

if ~exist('is_single','var')
    is_single = false;
end

if ~iscell(vec_type)
    vec_type = {vec_type};
end

result = [];

if (numel(k) > 1) || (numel(gamma) > 1) || (numel(vec_type) > 1)    
    nTests = numel(k) * numel(gamma) * numel(vec_type);   
    t_total = tic;
    counter = 0;
    
    for k_curr = k
        for gamma_curr = gamma
            for i_vectype = 1:numel(vec_type)
                counter = counter + 1;
                
                vec_type_curr = vec_type{i_vectype};
                result_curr = gsm_test(d, k_curr, gamma_curr, vec_type_curr, is_single, seed);
                
                if (counter == 1) || (result_curr.diff_mu_rel > result.diff_mu_rel)
                    k_worst_mu = k_curr;
                    gamma_worst_mu = gamma_curr;
                    vectype_worst_mu = vec_type_curr;                    
                end
                
                if (counter == 1) || (result_curr.diff_theta > result.diff_theta)
                    k_worst_theta = k_curr;
                    gamma_worst_theta = gamma_curr;
                    vectype_worst_theta = vec_type_curr;                    
                end                
                
                result = worst_result(result_curr, result);               
                
                fprintf('Finished %d/%d tests (%2g%%)\n', counter, nTests, counter/nTests*100);                   
            end
        end
    end
    
    result.k_worst_mu = k_worst_mu;
    result.gamma_worst_mu = gamma_worst_mu;
    result.vectype_worst_mu = vectype_worst_mu;

    result.k_worst_theta = k_worst_theta;
    result.gamma_worst_theta = gamma_worst_theta;
    result.vectype_worst_theta = vectype_worst_theta;
    
    t_total = toc(t_total);    
    result.t_total = t_total;
    
    result.seed = seed;
    
    if is_single
        result.class = 'single';
    else
        result.class = 'double';
    end

    fprintf('\n===============================\n\n');
    
    fprintf('\nTotal test time: %s\n', getTimeStr(t_total));

    fprintf('\nerr mu worst:    %g rel.\t(%g ulp)\n', result.diff_mu_rel, result.diff_mu_ulp);
    fprintf('                 at k=%d\tgamma=%g\tvectype=%s\n\n', k_worst_mu, gamma_worst_mu, vectype_worst_mu);

    fprintf('err theta worst: %g abs.\t(%g left, %g right)\n', result.diff_theta, result.diff_theta_l, result.diff_theta_r);
    fprintf('                 at k=%d\tgamma=%g\tvectype=%s\n\n', k_worst_theta, gamma_worst_theta, vectype_worst_theta);
    
    return;
end

% Here vec_type can only contain 1 cell, so we extract it
if iscell(vec_type)
    vec_type = vec_type{1};
end

if strcmp(vec_type, 'uniform')
    z = rand(d,1);
elseif strcmp(vec_type, 'gaussian')
    z = abs(randn([d,1]));
else
    error('Invalid vector type ''%s''', vec_type);
end

if is_single
    z = double(single(z));
    gamma = double(single(gamma));
end

t = tic;
[mu_gt, theta_gt, info_gt] = gsm_v5_1_quad(z,k,gamma,true);
tElapsed_gt = toc(t);

if is_single
    t = tic;
    [mu_test, theta_test, info_test] = gsm_v5_1_single(single(z),k,gamma,true);
    tElapsed_test = toc(t);
else
    t = tic;
    [mu_test, theta_test, info_test] = gsm_v5_1(z,k,gamma,true);
    tElapsed_test = toc(t);
end

% Analyze result
ulp_diff = @(x,x_gt) abs(x-x_gt) / 2^floor(log2(x_gt)) / eps(class(x));

is_k_trivial = (k == 0) || (k == d);
is_gamma_trivial = (gamma == 0) || isinf(gamma);

Il = info_gt.Il;
Ir = info_gt.Ir;

result = struct();
result.version = info_test.version;
result.seed = seed;

result.d = d;
result.k = k;
result.gamma = gamma;

if ~is_k_trivial
    result.dl = info_gt.dl;
    result.dr = info_gt.dr;
else
    result.dl = nan;
    result.dr = nan;
end

result.tElapsed_gt = tElapsed_gt;
result.tElapsed_test = tElapsed_test;

if k > 0
    result.diff_mu_rel = abs(mu_gt-mu_test) / abs(mu_gt);
    result.diff_mu_abs = abs(mu_gt-mu_test);
    result.diff_mu_ulp = ulp_diff(mu_test, mu_gt);
else
    if mu_gt == mu_test
        result.diff_mu_rel = 0;
        result.diff_mu_abs = 0;
        result.diff_mu_ulp = 0;
    else
        result.diff_mu_ulp = inf;
        result.diff_mu_abs = inf;
        result.diff_mu_ulp = inf;
    end
end

result.diff_theta  = max(abs(theta_gt-theta_test));
result.diff_theta_over_k = result.diff_theta / k;

if ~(is_k_trivial || is_gamma_trivial)
    result.diff_theta_l  = max(abs(theta_gt(Il)-theta_test(Il)));
    result.diff_theta_r  = max(abs(theta_gt(Ir)-theta_test(Ir)));
    
    result.diff_bz_rel = max(abs(info_test.bz(2:end)-info_gt.bz(2:end)) ./ abs(info_gt.bz(2:end)));
    result.diff_bzl_rel = max(abs(info_test.bzl(2:end)-info_gt.bzl(2:end)) ./ abs(info_gt.bzl(2:end)));
    result.diff_bzr_rel = max(abs(info_test.bzr(2:end)-info_gt.bzr(2:end)) ./ abs(info_gt.bzr(2:end)));
    
    result.diff_delta = max(abs(info_test.delta-info_gt.delta) ./ abs(info_gt.delta));
    result.diff_log_alpha = max(abs(info_test.log_alpha-info_gt.log_alpha) ./ abs(info_gt.log_alpha));
else
    result.diff_theta_l  = nan;
    result.diff_theta_r  = nan;
    
    result.diff_bz_rel = nan;
    result.diff_bzl_rel = nan;
    result.diff_bzr_rel = nan;
    
    result.diff_delta = nan;
    result.diff_log_alpha = nan;
end


fprintf('\nd=%d\tk=%d\tgamma=%g\n', d, k, gamma);
fprintf('\nTime:  test: %s\tgt: %s\n', getTimeStr(result.tElapsed_test), getTimeStr(result.tElapsed_gt));
fprintf('\nerr mu:    %g rel.\t(%g ulp)\n', result.diff_mu_rel, result.diff_mu_ulp);
fprintf('err theta: %g abs.\t(%g left, %g right)\n\n', result.diff_theta, result.diff_theta_l, result.diff_theta_r);


end


function rout = worst_result(r1, r2)
%function rout = worst_result(r1, r2)
%
% Given two result structs, with different measures of performance, returns a new struct
% with each field corresponding to the worst (max) corresponding fiend among the two 
% input structs.

if isempty(r2)
    rout = r1;
    return;
end

fns = fieldnames(r1);

rout = struct();

rout.version = r1.version;

for i=1:numel(fns)
    fn = fns{i};

    if strcmp(fn, 'version')
        continue
    end
    
    rout.(fn) = max(r1.(fn), r2.(fn));
end

end



function sOut = getTimeStr(t)
% Returns a string of the format HH:MM:SS.FFF describing the time given in
% t. t should be given in seconds.

if isnan(t)
    sOut = 'nan';
    return
end

sOut = sprintf('%s', datestr(t/24/60/60,'HH:MM:SS.FFF'));

if strcmp(sOut(1:3), '00:')
    sOut = sOut(4:end);
end

end

