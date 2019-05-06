

wall = 0.05;

D = 1.0;

/*
 *
 *        AXI-SYMMETRY           AXI-SYMMETRY           AXI-SYMMETRY
 *  -.-----.-----.-----.-----.-----.-----.-----.-----.-----.-----.-----.--- 
 * 0.5*D |      |                                                |   |
 *       -   1  o -------- o - - _ _                             |   | 5*D    
 *                         2         ^ ^ - - -                   |   |
 *                                             ^ ^ ^ - - - _ _ _ o   -
 *                                                                3
 *                L1 = 10*D            L2 = 128.87*D
 *              |----------|--------------------------------------| 
 * */

L1 = 4*D;
//L2 = 128.87*D;
L2 = 10*D;

k=10000;
Point(k+1) = {0,     -D/2.0,                    0,wall};
Point(k+2) = {L1,    -D/2.0,                    0,wall};
Point(k+3) = {L1+L2, -(0.5+ (L2*4.5*D/128.87)), 0,wall};

Line(k+1) = {k+1, k+2};
Line(k+2) = {k+2, k+3};

Point(k+4) = {0,     +D/2.0,                    0,wall};
Point(k+5) = {L1,    +D/2.0,                    0,wall};
Point(k+6) = {L1+L2, +(0.5+ (L2*4.5*D/128.87)), 0,wall};

Line(k+3) = {k+4, k+5};
Line(k+4) = {k+5, k+6};
Line(k+5) = {k+1, k+4};
Line(k+6) = {k+3, k+6};

Physical Line("wallNoSlip") = {k+4, k+3, k+5, -(k+1), -(k+2)};
Physical Line("wallOutflow") = {-(k+6)};

//+
Curve Loop(1) = {10004, -10006, -10002, -10001, 10005, 10003};
//+
Plane Surface(1) = {1};
//+
Characteristic Length {10004, 10001, 10002, 10005, 10003, 10006} = 0.2;
//+
Characteristic Length {10004, 10001, 10002, 10005} = 0.3;
//+
Characteristic Length {10006, 10003} = 0.4;
