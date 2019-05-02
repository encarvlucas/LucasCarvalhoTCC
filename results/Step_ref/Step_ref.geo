// Gmsh project created on Wed May  1 11:24:06 2019
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 2, 0, 1.0};
//+
Point(2) = {2, 2, 0, 1.0};
//+
Point(3) = {3, 2, 0, 1.0};
//+
Point(4) = {3, 1, 0, 1.0};
//+
Point(5) = {5, 1, 0, 1.0};
//+
Point(6) = {5, 0, 0, 1.0};
//+
Point(7) = {3, 0, 0, 1.0};
//+
Point(8) = {2, 0, 0, 1.0};
//+
Point(9) = {2, 1, 0, 1.0};
//+
Point(10) = {0, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 9};
//+
Line(9) = {9, 10};
//+
Line(10) = {10, 1};
//+
Curve Loop(1) = {2, 3, 4, 5, 6, 7, 8, 9, 10, 1};
//+
Plane Surface(1) = {1};
//+
Characteristic Length {3, 8, 9, 4} = 0.05;
//+
Characteristic Length {3, 4, 9, 8} = 0.1;
//+
Characteristic Length {2, 7} = 0.2;
//+
Characteristic Length {1, 10, 5, 6} = 0.4;
//+
Characteristic Length {4, 9, 3, 8} = 0.1;
//+
Characteristic Length {2, 7} = 0.2;
//+
Characteristic Length {1, 10, 5, 6} = 0.3;
