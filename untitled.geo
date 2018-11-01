//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Physical Line(1) = {3};
//+
Physical Line(2) = {4};
//+
Physical Line(3) = {2};
//+
Physical Line(4) = {1};
//+
Physical Surface(5) = {1};
//+
Recombine Surface {1};
//+
Recombine Surface {1};
//+
Recombine Surface {1};
//+
Recombine Surface {1};
//+
Recombine Surface {1};
