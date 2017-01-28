%Script to calculate the forces on each member of a truss, the cost of the
%truss, the maximum load, the failing member, and the max load to cost
%ratio

load('FinalDesign_Input.mat')
%creates the A matrix by concatenating two copies of C
A = [C;C];

%initializes necessary variables
%index holds the indeces for the two joints each member connects to
index = zeros(1,2);
[r,c] = size(C);
total_length = 0;
%vector of all straw lengths
radii = zeros(1,c);
%finds what the load withstood by the truss is
for i = 1:length(L)
    if L(i)~=0
        trussLoad = L(i);
    end
end

%begins processing the A matrix to find coefficients
%processes the matrix by moving down each column, then moving over to next
%column
for i = 1:c
   for j = 1:r
       %for each column, finds the two joints the member is connected to
       if A(j,i) == 1
           if index(1) == 0
               index(1) = j;
           else
               index(2) = j;
           end
       end
   end
   %finds the x distance and the y distance between the joints at the ends
   %of the member
   xval = abs(X(index(2)) - X(index(1)));
   yval = abs(Y(index(2)) - Y(index(1)));
   %calculates total distance between the two joints
   %equivalent to length of the straw
   radius = sqrt((xval^2) + (yval^2));
   %stores this radius in the radius vector
   radii(i) = radius;
   %updates the variable for the total length of the straws used
   total_length = total_length + radius;
   %if the second joint is to the right of the first joint, the first joint
   %has a positive coefficient and the second joint has a negative
   %coefficient
   %if the second joint is to the left of the first joint, the first joint
   %has a negative coefficient and the second joint has a positive
   %coefficient
   if X(index(2)) > X(index(1))
       A(index(1),i) = A(index(1), i) * (xval/radius);
       A(index(2), i) = A(index(2), i) * (-xval/radius);
   else
       A(index(1),i) = A(index(1), i) * (-xval/radius);
       A(index(2), i) = A(index(2), i) * (xval/radius);
   end
   %updates the coefficients in order to access the coefficients for the
   %forces in the y direction
   index(1) = index(1) + r;
   index(2) = index(2) + r;
   %performs same calculations as above except for y direction
    if Y(index(2)- r) > Y(index(1)- r)
       A(index(1),i) = A(index(1), i) * (yval/radius);
       A(index(2), i) = A(index(2), i) * (-yval/radius);
   else
       A(index(1),i) = A(index(1), i) * (-yval/radius);
       A(index(2), i) = A(index(2), i) * (yval/radius);
    end
   %resets the index vector for processing the next column
   index(1) = 0;
   index(2) = 0;
end

%concatenates Sx and Sy together
S = [Sx;Sy];
%concatenates A and s together to make the final A matrix
A = horzcat(A, S);
%creates the T matrix which holds the values for forces on each member
T = zeros(22,1);
%solves for the values in the T matrix
T = inv(A) * L;

%prints and processes results
fprintf('EK301 Fall 2016\n')
fprintf('Load: %.2f N \n', trussLoad)
fprintf('Member forces in Newtons: \n')

%determines whether to print out C for compression or T for tension
for i = 1:c
    force = T(i);
    if force < 0
        display_char = 'C';
    else
        display_char = 'T';
    end
    %prints out the member number, force withstood, and whether its in compression or
    %tension
    fprintf('m%i: %.3f (%c) \n', i, abs(force), display_char) 
end
%prints out the reaction forces
fprintf('Reaction forces in Newtons: \n')
fprintf('Sx1: %.2f \n', T(c+1))
fprintf('Sy1: %.2f \n', T(c+2))
fprintf('Sy2: %.2f \n', T(c+3))

%Calculations for max load, cost, and load/cost ratio
numJoints = r;
%calculates total cost according to given formula
cost = (numJoints*10) + total_length;
fprintf('Cost of truss: $%.2f \n', cost)
T = abs(T);
%creates vector with the max load each member could withstand according to
%its length
%uses equation Force = 1277.78 / length^2
maxForces = zeros(1,c);
for i = 1:c
    maxForces(i) = 1277.78 / (radii(i)^2);
end
%creates vector of ratios for the max force possibly withstood by each
%member to the force experienced while under the specific load
radiiRatio = zeros(1,c);
radiiRatio = maxForces ./ T(1:c)';
%finds the index and value of the smallest ratio
%this index indicates which member will fail first, as the small ratio
%indicates it is closest to experiencing max load
mymin =radiiRatio(1);
minindex =1;
for i = 2:c
    if (radiiRatio(i) < mymin) && (radiiRatio(i)~=Inf)
        mymin = radiiRatio(i);
        minindex = i;
    end
end
fprintf('The failing member is member %i \n', minindex)
%calculates max load, error of max load, and maxload to cost ratio
max_load = trussLoad * mymin;
load_error_percent = (643.7125/radii(minindex)^3) /(1277.78/(radii(minindex)^2));
load_error = max_load * load_error_percent;
fprintf('Max Load: %.3f N +/- %.3f N\n', max_load, load_error)
maxloadToCost = max_load / cost;
fprintf('Theoretical max load/cost ratio in N/$: %.4f \n', maxloadToCost);

%extra stuff for report
%printing out lengths
fprintf('Additional Details For Each Member \n')
for i = 1:c
    bucklingStrength = 1277.78/(radii(i)^2);
    uncertainty = 643.7125/(radii(i)^3);
    ForceAtMax = T(i) * (max_load/trussLoad);
    fprintf('m%i: \nLength = %.2f cm \nBuckling Strength = %.3f N\n', i,radii(i), bucklingStrength)
    fprintf('Uncertainty: %.3f N \nForce at Max Load: %.3f N\n', uncertainty, ForceAtMax)
end

