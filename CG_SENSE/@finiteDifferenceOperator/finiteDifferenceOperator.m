function  A = finiteDifferenceOperator(direction)

% usage:
%    A = finiteDifferenceOperator(direction)
%
% direction = array dimension in which the operation is applied

s.adjoint = 0;
s.direction = direction;

A = class(s,'finiteDifferenceOperator');