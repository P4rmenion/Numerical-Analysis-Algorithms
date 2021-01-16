%Author: Parmenion Charistos
%Date: December 2020

%Description: Code for the fourth exercise of the first mandatory project for the course of Numerical Analysis,
%             given by the Aristotle university of Thessaloniki. Implementation of the power method for calculating 
%             the maximum eigenvector of a sample (adjacency) matrix. Examples on a miniature Google matrix for Page ranking.

%E4

A = [0 1 0 0 0 0 0 0 1 0 0 0 0 0 0;     %Defining matrix A.
     0 0 1 0 1 0 1 0 0 0 0 0 0 0 0;
     0 1 0 0 0 1 0 1 0 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0 0 0 0 1 0 0 0;
     1 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 1 1 0 0 0 0;
     0 0 0 0 0 0 0 0 0 1 1 0 0 0 0;
     0 0 0 1 0 0 0 0 0 0 1 0 0 0 0;
     0 0 0 0 1 1 0 0 0 1 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
     0 0 0 0 0 0 1 1 0 0 1 0 0 0 0;
     0 0 0 0 0 0 0 0 1 0 0 0 0 1 0;
     0 0 0 0 0 0 0 0 0 1 1 0 1 0 1;
     0 0 0 0 0 0 0 0 0 0 0 1 0 1 0];
 
 A = transpose(A);  
 
 n = size(A,1);     %Setting n.
 q = 0.15;          %Setting q.
 G = zeros(n,n);    %Initializing Google matrix.
 
 S = zeros(n,1);    %Calculating sums of columns of matrix A.
 for i=1:1:n
    for j=1:1:n
        S(i) = S(i) + A(j,i);
    end
 end
 
 for i=1:1:n        %Defining Google matrix (according to formula).
     for j=1:1:n 
         G(i,j) = q/n + (A(i,j)*(1-q))/S(j);
     end
 end

 v = Power_Method(G); %Getting eigenvector of max eigenvalue of matrix G. 
 
 %-----Power Method-----
 
 %--Implements algorithm on Google matrix, sequentially calculates the
 %--eigenvector of its max eigenvalue and normalizes for its elements to
 %--add up to 1 (since they refer to probability).
 
 function lambda = Power_Method(G)
    n = size(G,1);
    b = rand(n,1); %Generating random column-vector b.
    
    for i=1:1:n    %Implementing algorithm.
        b_next = G*b;
        j = 1;
        while b(j) == 0
            j = j + 1;
        end
        b_next = (1/b_next(j)) * b_next;
        b = b_next;
    end
    
    sum = 0;       %Calculating sum of elements of calculated eigenvector.
    for i=1:1:n
        sum = sum + b(i);
    end
    
    lambda = (1/sum) * b;  %Normalizing eigenvector.
 end
 
 
 