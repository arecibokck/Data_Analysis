function y = sinc(x)
% This function calculates the sinc(x) (so sin(x)/x with 1 for x=0) ... yes
% surprisingly this function is not there by default ...

y = sin(x)./x;
y(x==0) = 1;

end