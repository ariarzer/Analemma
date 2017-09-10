#include <math.h>
#include <stdio.h>

const double PI = 3.141592;
const double epsilon = 23.43;
const double latitude = 60;
const double T = 365.25;
typedef struct Rotation { double a11; double a12; double a13;
                          double a21; double a22; double a23;
                          double a31; double a32; double a33;} Rot;
typedef struct Basis { double b11;
                       double b21;
                       double b31;} Bas;

typedef struct coord {double fi; double lam;} coord;

double rad(double A)
{
    double result = A * PI / 180.;
    return result;
}

double grad(double A)
{
    double result = A * 180. / PI;
    return result;
}

Rot Rx(double a)
{
    Rot Rx = {1, 0,            0,
              0, cos(rad(a)),  sin(rad(a)),
              0, -sin(rad(a)), cos(rad(a))};
    return Rx;
}

Rot Ry(double a)
{
    Rot Ry = {cos(rad(a)), 0, -sin(rad(a)),
              0,           1, 0,
              sin(rad(a)), 0, cos(rad(a))};
    return Ry;
}

Rot Rz(double a)
{
    Rot Rz = {cos(rad(a)),  sin(rad(a)),  0,
             -sin(rad(a)),  cos(rad(a)),  0,
              0,            0,            1};
    return Rz;
}

Bas basis(coord cord)
{
    Bas result;
    result.b11 = cos(rad(cord.fi)) * cos(rad(cord.lam));
    result.b21 = cos(rad(cord.fi)) * sin(rad(cord.lam));
    result.b31 = sin(rad(cord.fi));
    return result;
}

Bas RotBasMult (Rot rot, Bas bas)
{
    Bas result;
    result.b11 = rot.a11 * bas.b11 + rot.a12 * bas.b21 + rot.a13 * bas.b31;
    result.b21 = rot.a21 * bas.b21 + rot.a22 * bas.b21 + rot.a23 * bas.b31;
    result.b31 = rot.a31 * bas.b11 + rot.a32 * bas.b21 + rot.a33 * bas.b31;
    return result;
}

double equationOfTime(double date)
{
    double E, B;
    B = 2. * PI * sin(date / T);
    E = 7.53 * cos(B) + 1.5 * sin(B) - 9.87 * sin(2 * B);
    return E;
}

double S(double time, double date)
{
    double S;
    S = time + 24. * date / T + equationOfTime(date) / 60. + 12;
    if (S > 12) S = S - 24;
    S = S * 360. / 24.;
    return S;
}

coord eklipicToEcvator(coord ecl)
{
    coord ecv;
    ecv.fi =   asin(RotBasMult(Rx(epsilon), basis(ecl)).b31);
    ecv.lam = atan2(RotBasMult(Rx(epsilon), basis(ecl)).b21 / cos(ecv.fi),
                    RotBasMult(Rx(epsilon), basis(ecl)).b11 / cos(ecv.fi));
    ecv.fi =  grad(ecv.fi);
    ecv.lam = grad(ecv.lam);
    return ecv;
}

coord ecvatorToHorisont(coord ecv, double time, double date)
{
    coord horisont;
    coord tmp;

    tmp.fi =   asin(RotBasMult(Rz(S(time, date)), basis(ecv)).b31);
    tmp.lam = atan2(RotBasMult(Rz(S(time, date)), basis(ecv)).b21 / cos(tmp.fi),
                    RotBasMult(Rz(S(time, date)), basis(ecv)).b11 / cos(tmp.fi));
    tmp.fi =  grad(tmp.fi);
    tmp.lam = grad(tmp.lam);

    horisont.fi =   asin(RotBasMult(Ry((90 - latitude)), basis(tmp)).b31);
    horisont.lam = atan2(RotBasMult(Ry((90 - latitude)), basis(tmp)).b21 / cos(horisont.fi),
                         RotBasMult(Ry((90 - latitude)), basis(tmp)).b11 / cos(horisont.fi));
    horisont.fi =  grad(horisont.fi);
    horisont.lam = grad(horisont.lam);

    return tmp;
}

int main ()
{
    double E, M, V, TVR, e;
    coord ecl, ecv, hor;
    int t, j, time;
    time = 12;
    e = 0.016710;
    TVR = 102;
    E = 0;
    V = 0;
    FILE * f = fopen("analemma.dat","w");
    for( t = 1; t < T+1 ; t++ )
    {
        M = (2 * PI * t / T);
        for( j = 0; j < 2000; j++ )
        {
            E = E - ((E - e * sin(E) - M)/(1 - e * cos(E)));
        }
        V = 2 * atan(sqrt ((1 + e)/(1 - e)) * tan (E / 2));
        ecl.fi = 0;
        ecl.lam = grad(V) + TVR;
        ecv = eklipicToEcvator(ecl);
        hor = ecvatorToHorisont(ecv, time, t);
        fprintf(f, "%f %f \n", hor.fi, hor.lam);
//        fprintf(f, "%f %f \n", ecv.fi, ecv.lam);
    }
    fclose (f);
    return 0;
}


