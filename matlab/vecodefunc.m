% Vectorize ODE function (apply per each row).

function z = vecodefunc(f,y)
% TODO: work for varargout
z = cell2mat(arrayfun(@(v) f(y(v,:))',1:rows(y),'UniformOutput',false))';
