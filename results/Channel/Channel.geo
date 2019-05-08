//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {8, 0, 0, 1.0};
//+
Point(3) = {8, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Point(5) = {4, 0, 0, 1.0};
//+
Point(6) = {4, 1, 0, 1.0};
//+
Line(1) = {4, 6};
//+
Line(2) = {6, 3};
//+
Line(3) = {3, 2};
//+
Line(4) = {2, 5};
//+
Line(5) = {5, 1};
//+
Line(6) = {1, 4};
//+
Curve Loop(1) = {1, 2, 3, 4, 5, 6};
//+
Characteristic Length {4, 1, 3, 2} = 0.5;
//+
Characteristic Length {6, 5} = 0.2;
//+
Plane Surface(1) = {1};
//+
Characteristic Length {6, 5} = 0.1;
//+
Characteristic Length {6, 5} = 0.1;
//+
Characteristic Length {4, 1, 3, 2} = 0.2;
