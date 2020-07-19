//+
SetFactory("OpenCASCADE");
Sphere(1) = {-0, -0, 0, 0.7, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(2) = {-0, 0, 0, 1.0, -Pi/2, Pi/2, 2*Pi};
//+
Cylinder(3) = {0.8, 0, 0, 0.4, 0, 0, 0.2, 2*Pi};
//+
BooleanDifference{ Volume{2}; Delete; }{ Volume{1}; Delete; }
//+
BooleanUnion{ Volume{2}; Delete; }{ Volume{3}; Delete; }
//+
Physical Volume(1) = {1};
//+
Physical Surface(2) = {4};
//+
Physical Surface(3) = {3};
