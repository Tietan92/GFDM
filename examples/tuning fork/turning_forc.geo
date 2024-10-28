// Gmsh project created on Mon Sep 23 17:23:36 2024
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 0.01, -Pi, 0};
//+
Point(3) = {0.1, 7, 0, 1.0};
//+
Recursive Delete {
  Point{3}; 
}
//+
Point(3) = {-0.01, 0.07, 0, 1.0};
//+
Point(4) = {0.01, 0.07, 0, 1.0};
//+
Line(2) = {3, 1};
//+
Line(3) = {4, 2};
//+
Point(5) = {-0.015, 0, 0, 1.0};
//+
Point(6) = {-0.015, 0.07, 0, 1.0};
//+
Point(7) = {0.015, 0.07, 0, 1.0};
//+
Point(8) = {0.015, 0.0, 0, 1.0};
//+
Point(9) = {0.015, 0.06, 0, 1.0};
//+
Point(10) = {-0.015, 0.06, 0, 1.0};
//+
Line(4) = {3, 6};
//+
Line(5) = {6, 10};
//+
Line(6) = {10, 5};
//+
Line(7) = {4, 7};
//+
Line(8) = {7, 9};
//+
Line(9) = {9, 8};
//+
Point(11) = {-0.005, -0.03, 0, 1.0};
//+
Point(12) = {0.005, -0.03, 0, 1.0};
//+
Line(10) = {5, 11};
//+
Line(11) = {8, 12};
//+
Recursive Delete {
  Curve{10}; 
}
//+
Recursive Delete {
  Curve{11}; 
}
//+
Point(11) = {-0.015, -0.01, 0, 1.0};
//+
Point(12) = {0.015, -0.01, 0, 1.0};
//+
Point(13) = {0.005, -0.03, 0, 1.0};
//+
Point(14) = {-0.005, -0.03, 0, 1.0};
//+
Line(10) = {5, 11};
//+
Line(11) = {8, 12};
//+
Line(12) = {12, 13};
//+
Line(13) = {11, 14};
//+
Point(15) = {0.005, -0.07, 0, 1.0};
//+
Point(16) = {-0.005, -0.07, 0, 1.0};
//+
Line(14) = {14, 16};
//+
Line(15) = {16, 15};
//+
Line(16) = {15, 13};
//+
Curve Loop(1) = {13, 14, 15, 16, -12, -11, -9, -8, -7, 3, -1, -2, 4, 5, 6, 10};
//+
Plane Surface(1) = {1};
//+
Physical Surface("inner", 17) = {1};
//+
Physical Curve("grip_end", 18) = {15};
//+
Physical Curve("tuning_point_left", 19) = {5};
//+
Physical Curve("tuning_point_right", 20) = {8};
//+
Physical Curve("free_surface", 21) = {4, 7, 2, 3, 6, 9, 1, 12, 10, 11, 13, 16, 14};
