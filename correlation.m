function cor = correlation(A)
% Computes the correlation function of A 
% Output is centered at size(A)/2 because of the periodic boundary
% conditions

cor = conv2(A,repmat(rot90(A,2),2,2),'same')/numel(A);
cor = cor - (sum(sum(A))/numel(A))^2;
cor = circshift(cor,ceil(size(A)/2));

end
