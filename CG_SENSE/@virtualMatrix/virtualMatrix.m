function A = virtualMatrix(fhandle,varargin)

% function A = virtualMatrix(fhandle,fhandle_adjoint, input_size)
% function A = virtualMatrix(fhandle,input_size)
%
% The specified function handle can be used as if it was a matrix,
% i.e. instead of typing "fhandle(z)", "A*z" can be used.
%
% varargin = additional adjoint function handle and size of input
% or only size of input

s.fhandle = fhandle;
if strcmp(class(varargin{1}),'function_handle')
    s.fhandle_adjoint = varargin{1};
    sz = varargin{2};
    s.has_adjoint = 1;
else
    sz = varargin{1};
    s.has_adjoint = 0;    
end
s.sz = sz;
s.adjoint = 0;

A = class(s,'virtualMatrix');