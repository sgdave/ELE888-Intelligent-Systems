%% 
clc 
clear
close all
%%
eta = 0.1; 
theta = 0.001;


disp('************* Initial weight values ****************');
wih1=rand(1,3) %Wj1 weight vector input to hidden unit no.1
wih2=rand(1,3) %Wj2 weight vector input to hidden unit no.2
who1=rand(1,3) %Wkj weight vector hidden to output unit

x1=[-1 -1  1  1]; % Input1
x2=[-1  1 -1  1]; % Input2
t =[-1  1  1 -1]; % target (tk)

r=0; % initializing epoch
flag=0;

%%%%%%%% back propogation code %%%%%%%%%%%
while(flag==0)
    r=r+1;           % epoch increment
    wc_j1 = [0 0 0]; % initializing weight coefficients (i/p-->h1)
    wc_j2 = [0 0 0]; % initializing weight coefficients (i/p-->h2)
    wc_k  = [0 0 0]; % initializing weight coefficients (h  -->o/p)
    m=0;

    while(m < 4)
        m = m + 1;
        input = [1 x1(m) x2(m)]  % assume bias = 1 always
        net_j1= wih1*input';      % net activation of h1
        net_j2= wih2*input';      % net activation of h2

        % inputs to output layer [1 y1 y2]
        input_2(m,:)=[1 tanh(net_j1) tanh(net_j2)]
        net_k = who1*input_2(m,:)';% net activation of o/p
        op(m) = tanh(net_k);    % Zk calculate nonlinear output value

        % back propagation
        deltak = (t(m)- op(m)) * (1+tanh((net_k))^2); % Sk       error function
        wc_k_temp = eta*deltak*input_2(m,:);     % n*Sk*yi  updating the wc_k (delta_Wkj)

        delta_j1 = (1+tan((net_j1)^2))*(who1(2)*deltak);  % S1 h_unit1 sensitivity
        delta_j2 = (1+tan((net_j2)^2))*(who1(3)*deltak);  % S2 h_unit2 sensitivity
        

        wc_j1_temp = eta * delta_j1 * input;   % delta_Wj1  n*Sj*xi
        wc_j2_temp = eta * delta_j2 * input;   % delta_Wj2

        wc_k  = wc_k  + wc_k_temp;    % updating the wc_k
        wc_j1 = wc_j1 + wc_j1_temp;
        wc_j2 = wc_j2 + wc_j2_temp;
    end

    wih1 = wih1 + wc_j1 % initial wc+updated wc = new wc(input to h1)
    wih2 = wih2 + wc_j2 % initial wc+updated wc = new wc(input to h2)
    who1 = who1 + wc_k  % initial wc+updated wc = new wc(hidden to output)
    
    
    J(r)= 0.5*sum((t-op).^2);   % Eq(9)
    if r==1
        deltaJ = J(r);          % initialize deltaJ
    else
        deltaJ = abs(J(r-1) - J(r)); % update delta
    end
plot(J)
    if(deltaJ<theta);  % condition to determine convergence
        flag=1;
    end
end