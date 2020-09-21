function [times, dvals, kvals] = runSpeedTest()
dvals = [1e3, 1e3, 1e3, 1e4, 1e4, 1e4, 1e5, 1e5, 1e5, 1e6, 1e6, 1e6];
kvals = [ 10, 100, 500,  10, 100, 500,  10, 100, 500,  10, 100, 500];

gamma = 10;
nRepeat = 50;

n = numel(dvals);

times = nan(1,n);

for i=1:n
    d = dvals(i);
    k = kvals(i);
    z = randn([1,d]);

    t_curr = 0;
    
    for j=1:nRepeat
        fprintf('d=%g, k=%g, test %d...', d, k, j);
        z = abs(randn([d,1]));
        t_tmp = tic;
        [mu,theta] = gsm_v5_1(z,k,gamma);
        t_tmp = toc(t_tmp);
        t_curr = t_curr+t_tmp;
        
        fprintf(' %g\n', t_tmp);
    end
    
    t_curr = t_curr / nRepeat;
    
    times(i) = t_curr;
    fprintf('d=%g, k=%g average: %g\n\n', d, k, t_curr);
end

end
