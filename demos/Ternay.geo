L = 25;
H = 40;
// densité du maillage
d = 0.5;

// Définition des points
Point(1) = {-4,0,0,d};
Point(2) = {23,0,0,d};
Point(3) = {14,12,0,d};
Point(4) = {4,40,0,10*d};
Point(5) = {0,40,0,d};
Point(6) = {0,20,0,d};
Point(7) = {48.2,40,0,d};  // Le centre du cercle

// Définitions des lignes
Line(1) = {1,2};
Line(2) = {2,3};
Circle(3) = {4,7,3};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};

// Définition du contour (orienté)
Line Loop(7) = {1,2,-3,4,5,6};

// Définition d'une surface de contour "7"
Plane Surface(8) = {7};

// La surface 8 est affectée du numéro physique "1"
Physical Surface(1) = {8};

// Les lignes 5-6 sont affectées du numéro physique "1"
Physical Line(1) = {5,6};
// La ligne 1 est affectée du numéro physique "2"
Physical Line(2) = {1};

//
// Pour convertir le maillage *.msh en maillage Fenics,
// il faut executer "dolfin-convert Ternay.msh Ternay.xml"

