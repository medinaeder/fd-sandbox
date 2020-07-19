//+
SetFactory("OpenCASCADE");
xl = 3.0;
Cylinder(1) = {0, -0, 0, xl, 0, 0, 0.5, 2*Pi};
//+
Cylinder(2) = {-0, 0, 0, xl, 0, 0, 0.4, 2*Pi};
//+
Cylinder(3) = {xl, 0, 0, -0.1, 0, 0, 0.4, 2*Pi};
Box(4) = {-0, -0, -0.5, xl, 0.5, 1};
//+
BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }
//+
//BooleanUnion{ Volume{1}; Delete; }{ Volume{3}; Delete; }
//BooleanDifference{ Volume{1}; Delete; }{ Volume{4}; }
//+
BooleanUnion{ Volume{1}; Delete; }{ Volume{3}; Delete; }
//+
BooleanDifference{ Volume{5}; Delete; }{ Volume{4}; }
//+
Box(6) = {-0, -0, -0.5, xl, 0.3, 1};
BooleanIntersection{ Volume{6}; Delete; }{ Volume{4}; Delete; }
Coherence;

Physical Volume(1) = {6};
//+
Physical Volume(2) = {5};
//+
Physical Surface(3) = {7, 4};
//+
Physical Surface(4) = {5, 9, 6};
