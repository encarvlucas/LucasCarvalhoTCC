SetFactory("OpenCASCADE");
//+
Point(6) = {0.5+-0.475696, 0.16425, 0.0, 1.0};
//+
Point(7) = {0.5+-0.335515, 0.365799, 0.0, 1.0};
//+
Point(8) = {0.5+-0.170879, 0.463447, 0.0, 1.0};
//+
Point(9) = {0.5+-0.0306224, 0.487911, 0.0, 1.0};
//+
Point(10) = {0.5+0.0913384, 0.481878, 0.0, 1.0};
//+
Point(11) = {0.5+0.194936, 0.597761, 0.0, 1.0};
//+
Point(12) = {0.5+0.335177, 0.658783, 0.0, 1.0};
//+
Point(13) = {0.5+0.524198, 0.695431, 0.0, 1.0};
//+
Point(14) = {0.5+0.676666, 0.658877, 0.0, 1.0};
//+
Point(15) = {0.5+0.80475, 0.591832, 0.0, 1.0};
//+
Point(16) = {0.5+0.658328, 0.744159, 0.0, 1.0};
//+
Point(17) = {0.5+0.530244, 0.835554, 0.0, 1.0};
//+
Point(18) = {0.5+0.39609, 0.896459, 0.0, 1.0};
//+
Point(19) = {0.5+0.274126, 0.951272, 0.0, 1.0};
//+
Point(20) = {0.5+0.127821, 0.908513, 0.0, 1.0};
//+
Point(21) = {0.5+-0.0245825, 0.84136, 0.0, 1.0};
//+
Point(22) = {0.5+-0.164765, 0.725434, 0.0, 1.0};
//+
Point(23) = {0.5+-0.292753, 0.609514, 0.0, 1.0};
//+
Point(24) = {0.5+-0.428851, 0.429716, 0.0, 1.0};
//+
Point(25) = {0.5+0.409811, 0.930426, 0.0, 1.0};
//+
Point(26) = {0.5+0.512491, 0.886172, 0.0, 1.0};
//+
Point(27) = {0.5+0.590407, 0.845373, 0.0, 1.0};
//+
Point(28) = {0.5+0.674591, 0.782751, 0.0, 1.0};
//+
Point(29) = {0.5+0.746361, 0.698382, 0.0, 1.0};
//+
Bezier(1) = {6, 7, 8, 9, 10};
//+
Bezier(2) = {10, 11, 12, 13, 14, 15};
//+
Bezier(3) = {15, 29, 28, 27, 26, 25, 19};
//+
Bezier(4) = {19, 20, 21, 22, 23, 24, 6};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};

//+
Physical Point(1) -= {24, 7, 8, 9, 23};
//+
Physical Point(1) -= {24};
//+
Recursive Delete {
  Point{7}; 
}
//+
Recursive Delete {
  Point{8}; 
}
//+
Recursive Delete {
  Point{24}; 
}
//+
Recursive Delete {
  Point{23}; 
}
//+
Recursive Delete {
  Point{22}; 
}
//+
Recursive Delete {
  Point{21}; 
}
//+
Recursive Delete {
  Point{20}; 
}
//+
Recursive Delete {
  Point{18}; 
}
//+
Recursive Delete {
  Point{25}; 
}
//+
Recursive Delete {
  Point{26}; 
}
//+
Recursive Delete {
  Point{27}; 
}
//+
Recursive Delete {
  Point{28}; 
}
//+
Recursive Delete {
  Point{29}; 
}
//+
Recursive Delete {
  Point{16}; 
}
//+
Recursive Delete {
  Point{17}; 
}
//+
Recursive Delete {
  Point{13}; 
}
//+
Recursive Delete {
  Point{14}; 
}
//+
Recursive Delete {
  Point{12}; 
}
//+
Recursive Delete {
  Point{11}; 
}
//+
Recursive Delete {
  Point{9}; 
}
//+
Curve Loop(2) = {4, 1, 2, 3};
//+
Plane Surface(2) = {2};
//+
Characteristic Length {6, 10, 19, 15} = 0.1;
//+
Characteristic Length {10, 6, 19, 15} = 0.1;
//+
Characteristic Length {6, 10, 19, 15} = 0.5;
//+
Characteristic Length {6, 10, 19, 15} = 0.05;
//+
Characteristic Length {6, 10, 19, 15} = 0.02;
//+
Characteristic Length {6, 10, 19, 15} = 0.04;
