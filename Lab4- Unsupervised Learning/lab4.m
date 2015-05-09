clc
clear all
close all
%% Set Up

I= imread('house.tiff');
X = reshape(I, 256*256, 3);
X = double(X);
R = [0,0,1];G = [0,1,0];B = [1,0,0];
figure, plot3(X(:,1), X(:,2), X(:,3),'.','Color',[1, 0, 0])

mean1 = [.1 .1 .1];

%% a- c=2

M = [.1 .1 .1
     .9 .9 .9];

M = M*256;

pre_M = zeros(size(M));

J = [];
M1_history = [M(1,:)];
M2_history = [M(2,:)];

while (pre_M ~= M)

    pre_M = M;
    
    J1 = (X - repmat(M(1,:), size(X,1), 1));
    J1 = sum(J1.^2, 2);

    J2 = (X - repmat(M(2,:), size(X,1), 1));
    J2 = sum(J2.^2, 2);

    cluster1 = J1 < J2;
    cluster2 = ~cluster1;

    M(1,:) = sum(X(cluster1, :)) / sum(cluster1);
    M(2,:) = sum(X(cluster2, :)) / sum(cluster2);
    
    J = [J sum(min(J1, J2))];
    M1_history = [M1_history; M(1,:)];
    M2_history = [M2_history; M(2,:)];
    
    M
    
    
    
end

% i
close all
plot(J)

% ii
figure
plot3(M1_history(:, 1), M1_history(:, 2), M1_history(:, 3), '-*')
hold all
plot3(M2_history(:, 1), M2_history(:, 2), M2_history(:, 3), '-*')

% iii
figure
X1 = X(cluster1, :);
X2 = X(cluster2, :);

plot3(X1(:,1), X1(:,2), X1(:,3),'.','Color', M(1,:)/256)
hold all
plot3(X2(:,1), X2(:,2), X2(:,3),'.','Color', M(2,:)/256)
% iv
figure
XX = repmat(M(1,:), size(X,1), 1) .* repmat(cluster1, 1, size(X,2));
XX = XX + repmat(M(2,:), size(X,1), 1) .* repmat(cluster2, 1, size(X,2));
XX = reshape(XX, size(I, 1), size(I, 2), 3);

subplot(1,2,1);imshow(I)
subplot(1,2,2);imshow(XX/256)


%% b

c = 5;

initial_M_1 = [
  173.8240   60.4672  169.4464
  162.7648   30.5664  197.1968
  241.9712  155.4688   89.6512
   53.4784  115.2256  169.4720
  181.5808  117.4272  106.5472
];
XX=[];
%M = rand(c,3) * 256;
M = initial_M_1;

pre_M = zeros(size(M));

while (pre_M ~= M)

    pre_M = M;
    
    J = zeros(size(X,1), c);
    
    for i = [1:c]
        j = (X - repmat(M(i,:), size(X,1), 1));
        j = sum(j.^2, 2);
        J(:,i) = j;
    end
    [a, cluster] = min(J, [], 2);
    
    for i = [1:c]
        this_cluster = (cluster==i);
        M(i, :) = sum(X(this_cluster, :)) / sum(this_cluster);

    end
XX = zeros(size(X));
for i = [1:c]
    this_cluster = (cluster==i);
    XX = XX + repmat(M(i,:), size(X,1), 1) .* repmat(this_cluster, 1, size(X,2));
end
XX = reshape(XX, size(I, 1), size(I, 2), 3);
subplot(1,2,1);imshow(I)
subplot(1,2,2);imshow(XX/256)



end

M1 = M
cluster1 = cluster;
%
figure
for i = [1:c]
    this_cluster = (cluster==i);
    Xi = X(this_cluster, :);
    plot3(Xi(:,1), Xi(:,2), Xi(:,3),'.','Color', M(i,:)/256)
    hold all
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
initial_M_2 = [
  215.5339  138.4293  240.5963
  213.2267  222.7049  165.2613
   65.6489   67.7834  122.7426
  157.0459   81.4270  163.6651
  149.0558   30.5189  139.4473
];

M = initial_M_2;

pre_M = zeros(size(M));

while (pre_M ~= M)

    pre_M = M;
    
    J = zeros(size(X,1), c);
    
    for i = [1:c]
        j = (X - repmat(M(i,:), size(X,1), 1));
        j = sum(j.^2, 2);
        J(:,i) = j;
    end
    [a, cluster] = min(J, [], 2);
    
    for i = [1:c]
        this_cluster = (cluster==i);
        M(i, :) = sum(X(this_cluster, :)) / sum(this_cluster);
    end
    
end
M2 = M
cluster2 = cluster;

%
figure
for i = [1:c]
    this_cluster = (cluster==i);
    Xi = X(this_cluster, :);
    plot3(Xi(:,1), Xi(:,2), Xi(:,3),'.','Color', M(i,:)/256)
    hold all
end

%% c

N = size(X,1);
XB1 = 0;

for i = [1:c]
    this_cluster = (cluster1==i);
    Xi = X(this_cluster, :);
    mui_j = sort(sum((M1 - repmat(M1(i,:), c, 1)).^2, 2).^.5);
    XB1 = XB1 + sum(sum((Xi - repmat(M1(i,:), size(Xi,1), 1)).^2, 2).^.5) / mui_j(2);
end

XB1 = XB1 / N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
XB2 = 0;

for i = [1:c]
    this_cluster = (cluster2==i);
    Xi = X(this_cluster, :);
    mui_j = sort(sum((M2 - repmat(M2(i,:), c, 1)).^2, 2).^.5);
    XB2 = XB2 + sum(sum((Xi - repmat(M2(i,:), size(Xi,1), 1)).^2, 2).^.5) / mui_j(2);
end

XB2 = XB2 / N