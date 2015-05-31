A = importdata('exactsolution_h_256.txt');
X = A(:,1);
Y =A(:,2);
Z = A(:,3);

tri = delaunay(X,Y);
trisurf(tri,X,Y,Z);
xlabel('X');
ylabel('Y');
title 'Exact Solution for Neumann Problem';

A = importdata('solution_h_256.txt');
X = A(:,1);
Y =A(:,2);
Z = A(:,3);
figure
tri = delaunay(X,Y);
trisurf(tri,X,Y,Z);
colormap([1  1  0; 0  1  1]);
xlabel('X');
ylabel('Y'); 
title 'Approximated Solution for Neumann Problem';

