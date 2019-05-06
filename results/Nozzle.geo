// Gmsh project created on Thu Jun  4 09:24:57 2009

wall = 0.2;

D = 1.0;

/* sinusoidal wall configuration
 *
 *                      2*pi
 * Eq.: y(x) = A sin( -------- * X - phase)
 *                     lambda
 *
 * A = amplitude
 * lambda = wavelength (channel length / number of corrugations )
 * X = x coordinate
 * phase = phase displacement
 *
 *       2*pi
 * k = -------- --> wave number
 *      lambda
 * */
A = 0.004;
stretch = 10;
phase = 0.0;
xcf = 0.25*stretch;
Xth = 5.0; // X position of throat
nPoints = 40+1; // total number of points in the sinusoidal line
Printf("nPoints: ",nPoints);
Printf("-------------- Simulator2D:setALEBC() --------------");
Printf("  phase: %f",phase);
Printf("  Y: %f",D/2.0 + A*(phase)*(phase));
Printf("----------------------------------------------------");

k = 10000;
j = 1+k;
// top line
For i In {1:nPoints}
 X = stretch*( (i-1)/(nPoints-1) );
 If( X <= Xth )
  Y = 2.0 + 3.0*(((X-phase)/Xth) - 1.5)*((X-phase)/Xth)*((X-phase)/Xth);
 EndIf
 If( X > Xth )
  Y = 3.0 - (X/Xth)*(6.0-4.5*(X/Xth)+(X/Xth)*(X/Xth));
 EndIf
 
 Point(j) = {X, Y, 0, wall};
 j = j + 1;
 Printf("X: %f, Y: %f",X,Y);
EndFor


j = 1+k;
// lines
For i In {1:nPoints-1}
 Line(j) = {j, j+1};
 j = j + 1;
EndFor

k = newp;
/*  
 *   k+1             k+2  6           3  k+3               k+5
 *    o----------------o---o------------o---o-----------------o
 *    
 *    |           ll          | r | r |           lr          |
 *    |<--------------------->|<->|<->|<--------------------->|
 *    |                        stretch                        |
 *    |<----------------------------------------------------->|
 */

Point(k+1) = {0.0, 0.0, 0.0, wall};
Point(k+2) = {stretch,   0.0, 0.0, wall};


//+
Line(10041) = {10001, 10043};
//+
Line(10042) = {10043, 10044};
//+
Line(10043) = {10044, 10041};
//+
Curve Loop(1) = {10042, 10043, -10040, -10039, -10038, -10037, -10036, -10035, -10034, -10033, -10032, -10031, -10030, -10029, -10028, -10027, -10026, -10025, -10024, -10023, -10022, -10021, -10020, -10019, -10018, -10017, -10016, -10015, -10014, -10013, -10012, -10011, -10010, -10009, -10008, -10007, -10006, -10005, -10004, -10003, -10002, -10001, 10041};
//+
Plane Surface(1) = {1};
