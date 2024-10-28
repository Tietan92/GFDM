// Gmsh project created on Tue Aug 13 17:59:38 2024
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 1, 0.1, 0};
//+
Physical Surface("inner", 5) = {1};
//+
Physical Curve("border_left", 6) = {4};
//+
Physical Curve("border_right", 7) = {2};
//+
Physical Curve("border_up", 8) = {3};
//+
Physical Curve("border_down", 9) = {1};
