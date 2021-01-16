%Selecting 10 points from which the sin function will be approximated in the interval [-pi.pi].
x = [-pi,-3*pi/4,-pi/2,-pi/4,0,pi/4,1,pi/2,3*pi/4,pi];
y = [0,-0.707,-1,-0.707,0,0.707,0.841,1,0.707,0];

%Plotting selected points.
subplot(1,3,1);
plot(x,y,'-k.','MarkerSize',10);
ylim([-1,1]);

%Plotting results from lagrange method.
targetX = (-pi:0.01:pi);
resY = Lagrange(x,y,targetX);
subplot(1,3,2);
plot(targetX,resY);
ylim([-1,1]);

%Plotting results for linear approximation via minimum squares.
targetX = (-pi:0.025:pi);
resY = MinimumSquares(x,y,2,targetX);
subplot(1,3,3);
plot(targetX,resY);
ylim([-1,1]);

%function that takes a set of points (x,y) and an vector of target x coordinates
%and calculates the vector's corresponding y coordinates via the occurring lagrange polynomial.
function res = Lagrange(x,y,targetX)
    n = length(x);
    L = ones(n,length(targetX));
    for i=1:1:n
        for j=1:1:n
            if i~=j
                L(i,:) = L(i,:).*(targetX - x(j))/(x(i) - x(j));
            end
        end
    end
    res = 0;
    for i=1:1:n
        res = res + y(i) * L(i,:);
    end
end

function resY = MinimumSquares(x,y,rank,targetX)
    %Initializing.
    rows = zeros(rank,rank);
    cols = zeros(rank,1);

    %Calculating row matrix.
    for i = 1:rank
        for j = 0:rank-1
            for k = 1:length(x)
               rows(i,j+1) = rows(i,j+1) + x(k)^(j+i-1);
            end
        end
    end

    %Calculating column vector
    for i = 1: rank
        for j = 1 : length(y)
            cols(i) = cols(i) + x(j)^(i-1) * y(j);
        end
    end

    %Returning approximations after calculating products.
    res = rows \ cols;
    resY = 0;
    for i=1:1:length(res)
        resY = resY + res(i)*targetX.^(i-1);
    end
end