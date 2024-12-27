// Define the radius of the hole
R = 0.3;
// Define half the length of the square domain
L = 1.0;

// Set the points of the square domain
Point(1) = {L, -L, 0};
Point(2) = {L, L, 0};
Point(3) = {-L, L, 0};
Point(4) = {-L, -L, 0};
// Boundary points of the arc
Point(5) = {-L + R, -L, 0};
Point(6) = {-L, -L + R, 0};
// Mid point of the arc
Point(7) = {-L + Cos(Pi/4) * R, -L + Sin(Pi/4) * R, 0};

// Set the arcs by defining {start, center, end} points
Circle(1) = {5, 4, 7};
Circle(2) = {7, 4, 6};

// Set the straight lines by two points
Line(3) = {6, 3};
Line(4) = {3, 2};
Line(5) = {2, 1};
Line(6) = {1, 5};
Line(7) = {2, 7};

// Select the lines of the surface boundary
Curve Loop(1) = {4, 7, 2, 3};
// Create the surface from loop
Plane Surface(1) = {1};

// Select the lines of the surface boundary ('-' stands for reverse)
Curve Loop(2) = {7, -1, -6, -5};
// Create the surface from loop
Plane Surface(2) = {2};

// Set 3 nodes on each lines using transfinite method
Transfinite Line{1, 2, 3, 4, 5, 6, 7} = 3;

// Set transfinite mesh on the surfaces
Transfinite Surface{1};
Transfinite Surface{2};

// Merge the triangular meshes into quadrangular ones
Recombine Surface{1};
Recombine Surface{2};

// Set the elements as first order
Mesh.ElementOrder = 1;
// Generate the mesh using "Frontal-Delaunay for Quads" algorithm
Mesh.Algorithm = 8;

// EOF
