%Defining data for approximation base.
trainX = [7,8,9,10,11,12,13];
trainY_DEI = [5.07,5.06,5.49,5.5,5.595,5.595,5.55];
trainY_INLOT = [0.115,0.116,0.125,0.123,0.123,0.12,0.127];

%Defining actual values for comparison plots.
totalX = [7,8,9,10,11,12,13,14,15,16,17,18];
realY_DEI = [5.07,5.06,5.49,5.5,5.595,5.595,5.55,5.855,5.855,5.98,5.9,5.975];
realY_INLOT = [0.115,0.116,0.125,0.123,0.123,0.12,0.127,0.15,0.144,0.15,0.15,0.153];

targetX = [14,15,16,17,18];

%Plotting.
subplot(2,2,1);
resY = MinimumSquares(trainX,trainY_DEI,2,targetX);
plot (targetX,resY, '- .r', 'MarkerSize', 15);
hold on
resY = MinimumSquares(trainX,trainY_DEI,3,targetX);
plot (targetX,resY, '- .g', 'MarkerSize', 15);
hold on
resY = MinimumSquares(trainX,trainY_DEI,4,targetX);
plot (targetX,resY, '- .b', 'MarkerSize', 15);
title('Predicted Stock Value ΔΕΗ');
hold off

subplot(2,2,2);
resY = MinimumSquares(trainX,trainY_INLOT,2,targetX);
plot (targetX,resY, '- .r', 'MarkerSize', 15);
hold on
resY = MinimumSquares(trainX,trainY_INLOT,3,targetX);
plot (targetX,resY, '- .g', 'MarkerSize', 15);
hold on
resY = MinimumSquares(trainX,trainY_INLOT,4,targetX);
plot (targetX,resY, '- .b', 'MarkerSize', 15);
title('Predicted Stock Value ΙΝΤΡΑΛΟΤ');
hold off

subplot(2,2,3);
plot(totalX,realY_DEI,'-k.','MarkerSize',15);
ylim([2,6]);
xlim([14,18]);
title('Actual Stock Value ΔΕΗ');

subplot(2,2,4);
plot(totalX,realY_INLOT,'-k.','MarkerSize',15);
title('Actual Stock Value ΙΝΤΡΑΛΟΤ');
ylim([0.1,0.25]);
xlim([14,18]);

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

    %Returning matrix with coefficients
    res = rows \ cols;
    resY = 0;
    for i=1:1:length(res)
        resY = resY + res(i)*targetX.^(i-1);
    end
end