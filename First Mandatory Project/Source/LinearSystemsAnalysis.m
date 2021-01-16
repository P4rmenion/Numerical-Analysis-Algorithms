%Author: Parmenion Charistos
%Date: December 2020

%Description: Code for the third exercise of the first mandatory project for the course of Numerical Analysis,
%             given by the Aristotle university of Thessaloniki. Implementation of the PA = LU analysis method
%             for solving a linear equations system Ax = b. Implementation of a Cholesky decomposition function.
%             Implementation of the Gauss Seidel method for solving a linear equations system of a tridiagonal
%             square matrix A.

%E3

%-----PA = LU system analysis-----

%Solves linear systems of the form Ax = b.
%Matrix A should be square.

A = [3 1 -4 1 ; -5 2 1 -2 ; -1 6 -3 -4 ; -2 1 -4 2]; %Matrix of coefficients.
b = [1 ; -3 ; 2 ; 0];                                %Vector of constants.
[P,A,L,U,x] = SolveSystem(A,b);

A = [1 2 3 4 ; 2 3 4 1 ; 3 4 1 2 ; 4 1 2 3];         %Redefining A as a symmetric, positive - definite matrix (in order to be viable for Cholesky decomposition)
L = Cholesky(A);

e = 0.00005;                    %Desired precision.
[x,N] = Gauss_Seidel(10,e); 
x

%Function Definitions

%--function that calculates the solution vector x via PA = LU analysis--

function [P,A,L,U,x] = SolveSystem(A,b)
    n = size(A,1);                      
    [A,P] = SortDrivers(A);     %Sort A and form permutation matrix P.
    [L,U] = LU_Analysis(A);     %Calculate matrices L*U = A.
    r = P*b;                    %Rename P*b product (for ease of use).
    
    %We now have a system of the form LUx = r.
    %Let y = Ux.
    
    y(1) = r(1)/L(1,1);         %Solve Ly = r for y.
    for i=2:1:n
        sum = 0;
        for j=1:1:i-1
           sum = sum+L(i,j)*y(j); 
        end
        y(i) = (r(i)-sum)/L(i,i);
    end
    
    x(n) = y(n)/U(n,n);         %Solve Ux = y for x.
    for i=n-1:-1:1
        sum=0;
        for j=n:-1:i
            sum = sum+U(i,j)*x(j);
        end
        x(i) = (y(i)-sum)/U(i,i);
    end
    x = transpose(x);
end

%--function that sorts driver elements by max absolute value
%--by swapping rows and forms the according permutation matrix P--

function [A,P] = SortDrivers(A)
    n = size(A,1);
    P = zeros(n,n);
    
    for i=1:1:n
        P(i,i) = 1;  %Create n-Identity matrix P.
    end
    
    for j=1:1:n
        max = 0;
        optimal_row = j;
        for i=j:1:n                   %For each column:        
            if abs(A(i,j)) > max      %Find element with max abs.
                max = abs(A(i,j));    %Update max.
                optimal_row = i;      %Save index of row that contains max.
            end
        end
        if j ~= optimal_row           %If max is not a driver element (on the main diagonal), swap rows to make it so.
            A([j optimal_row],:) = A([optimal_row j],:);
            P([j optimal_row],:) = P([optimal_row j],:);
        end
    end
end


%-----LU Analysis-----

function [L,U] = LU_Analysis(A)       
    n = size(A,1);
    U = zeros(n,n);
    L = zeros(n,n);
    for i=1:1:n         %Using the classic LU analysis algorithm (Source: Introduction to Numerical Analysis by Akrivis - Dougalis)
        for j=1:1:i-1
            sum = 0;
            for k=1:1:j-1
                sum = sum + L(i,k)*U(k,j);
            end
            L(i,j) = (A(i,j) - sum) / U(j,j);
        end
        L(i,i) = 1;
        for j=i:1:n
            sum = 0;
            for k=1:1:i-1
                sum = sum + L(i,k)*U(k,j);
            end
            U(i,j) = A(i,j) - sum;
        end
    end
end


%-----Cholesky Decomposition-----

function L = Cholesky(A)
    n = size(A,1);
    for i=1:1:n     %Using the classic Cholesky decomposition algorithm (Source: Introduction to Numerical Analysis by Akrivis - Dougalis)
       for j=1:1:i-1
           sum = 0; 
           for k=1:1:j-1
                sum = sum + L(i,k)*L(j,k);
           end
           L(i,j) = (A(i,j)-sum) / L(j,j);
       end
       
       sum = 0;
       for k=1:1:i-1
           sum = sum + L(i,k)^2;
       end
       L(i,i) = sqrt(A(i,i) - sum);
    end
end


%-----Gauss Seidel Method-----

function [x,N] = Gauss_Seidel(n,e)
    A = zeros(n,n);       %Defining n-n tridiagonal matrix A. 
    
    for i = 1:1:n
       A(i,i) = 5;
       A(i+1,i) = -2;
       A(i,i+1) = -2;
       b(i) = 1;
    end

    b(1) = 3;             %Defining vector of constants b.
    b(end) = 3;
    b = transpose(b);
    
    x = zeros(n,1);       %Defining vector x0.
    norm = 0;             %Initializing infinity norm of x0;
    N = 0;                %Iterations counter.
    
    current_e = intmax;
    x_next = zeros(n,1);  %Initializing next-approximation vector (optimizing calculation cost).
    
    while current_e > e
        N = N + 1;
        for i=1:1:n       %Using classic Gauss-Seidel algorithm
            sum1 = 0;
            for j=1:1:i-1
                sum1 = sum1 + A(i,j)*x_next(j);
            end
            
            sum2 = 0;
            for j=i+1:1:n
                sum2 = sum2 + A(i,j)*x(j);
            end
               
            x_next(i) = (1/A(i,i))*(b(i) - sum1 - sum2);   %Updating next-approximation vector.
        end
        
        norm_next = 0;
        for k=1:1:n                         %Calculating new infinity norm.
            if abs(x_next(k)) > norm_next
                norm_next = abs(x_next(k));
            end
        end
        
        current_e = norm_next - norm;       %Updating current error.
        x = x_next;
        norm = norm_next;                   %Saving current approximation vector and its norm for the following iteration.
    end
end