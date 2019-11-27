function J = jaccard(A,B)

% function to caculate the Jaccard index (intersection over union)
% assumeing two binarised vectors of the same length. 1 = identical, 0 = no
% overlap.

if min(size(A)) ~= 1 || min(size(B)) ~= 1 || numel(A) ~= numel(B)
    error('A and B must be vectors of the same length');
end

J = sum(A & B)./sum(A | B);