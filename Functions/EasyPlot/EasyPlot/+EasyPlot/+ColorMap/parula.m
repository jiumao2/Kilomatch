function map = parula(N)

if nargin<1 || isempty(N)
	map = parula;
else
	assert(isscalar(N)&&isreal(N),'First argument must be a real numeric scalar.')
end

map = parula(N);
end