function [ s ] = scalar( f, g, dt, Norm )
%scalar product of 2 vectors
%   s = 1/Norm * sum(f(i)*g(i)*dt)
s = 0;
for i = 1 : length(f)
    s = s + f(i)*g(i)*dt;
end
s = s/Norm;

end

