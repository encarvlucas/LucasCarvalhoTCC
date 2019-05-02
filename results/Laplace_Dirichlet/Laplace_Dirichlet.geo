//+
Point(1) = {0, 0, 0, 0.2};
//+
Point(2) = {1, 0, 0, 0.2};
//+
Point(3) = {1, 1, 0, 0.2};
//+
Point(4) = {0, 1, 0, 0.2};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Characteristic Length {4, 1, 3, 2} = 0.1;
//+
Characteristic Length {4, 3, 2, 1} = 0.05;
