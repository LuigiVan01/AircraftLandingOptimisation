function A = gen_rate_lim_mat(l)
N = l;
if N == 1 
    A = [1;-1];
else
    A1 = zeros(N, N);
    N_col = size(A1,2);
    N_row = size(A1,1);
    
    A1(1,1) = 1;
    A1(1,2) = -1;
    for ii = 2:N_row
        for jj = 2:N_col
            A1(ii, jj) = A1(ii-1,jj-1);
        end
       
    end
    A = [A1;-A1];
end

end