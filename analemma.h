#ifndef ANALEMMA_H
#define ANALEMMA_H

#include "stdafx.h"

typedef struct Rotation { double a11; double a12; double a13;
                          double a21; double a22; double a23;
                          double a31; double a32; double a33;} Rot;
typedef struct Basis { double b11;
                       double b21;
                       double b31;} Bas;
typedef struct coord {double fi; double lam;} coord;

//Transformatiom rad-tograd and grad-to-rad
double toRad(double A);
double toGrad(double A);
coord toRad(coord A);
coord toGrad(coord A);

//Initialization rotatin-matrix
Rot Rx(double a);
Rot Ry(double a);
Rot Rz(double a);

//Initialization coordinate-basis
Bas basis(coord cord);

//Multiplication rotatin-matrix by coordinate-basis
Bas RotBasMult (Rot rot, Bas bas);

//Equation of time
double equationOfTime(double date);

//Siderial time
double S(double time, double date);

//Transformatiom coordinate
coord eklipicToEcvator(coord ecl);
coord ecvatorToHorisont(coord ecv, double time, double date);

//Generate GNUPLOT-script
int generateGnuScript(std::string filename, int number,
                      int sizex_image = 1000, int sizey_image = 800, int size_font = 12,
                      int sizex_diagramm_from = 0, int sizex_diagramm_to = 0,
                      int sizey_diagramm_from = 0, int sizey_diagramm_to = 0,
                      std::string xlabel = "", std::string ylabel = "",
                      std::string title = "");

#endif //ANALEMMA_H
