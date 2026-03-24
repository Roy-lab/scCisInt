function run_nmf(DataFile, k, outdir, varargin)

A = importdata(DataFile); 
[U, V] = nmf(A, k, varargin{:}); 

writematrix(U, sprintf('%s/U.txt', outdir));
writematrix(V', sprintf('%s/V.txt', outdir));

[~, U_assign] = max(U, [], 2); 
[~, V_assign] = max(V',[], 2);

writematrix(U_assign, sprintf('%s/U_assign.txt', outdir));
writematrix(V_assign, sprintf('%s/V_assign.txt', outdir));

