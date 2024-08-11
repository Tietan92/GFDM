// Gmsh project created on Thu Jul  4 14:25:16 2024
//+
SetFactory("OpenCASCADE");
Rectangle(1) = {-0.5, -0.5, 0, 1, 1, 0};
//+
Rectangle(2) = {-5, -5, 0, 25, 10, 0};
//+
Curve Loop(3) = {2, 3, 4, 1};
//+
Curve Loop(4) = {5, 6, 7, 8};
//+
Plane Surface(3) = {3, 4};
//+
Recursive Delete {
  Surface{1}; 
}
//+
Recursive Delete {
  Surface{2}; 
}
//+
Physical Surface("inner_area", 9) = {3};
//+
Physical Curve("inflow", 10) = {8};
//+
Physical Curve("outflow", 11) = {6};
//+
Physical Curve("upper_border", 12) = {7};
//+
Physical Curve("down_border", 13) = {5};
//+
Physical Curve("square", 14) = {1, 2, 3, 4};
