function D = kron_columns(A, B)
na = size(A,2); nb = size(B,2);
D = zeros(size(A,1), na*nb);
col = 1;
for i = 1:na
    for j = 1:nb
        D(:,col) = A(:,i) .* B(:,j);
        col = col + 1;
    end
end
end
