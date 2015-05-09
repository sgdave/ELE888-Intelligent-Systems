function [value,deriv] = sigmoid(x) 
value=tanh(x);
deriv=(1 - (tanh(x))^2);
end
