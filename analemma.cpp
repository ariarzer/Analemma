#include <math.h>
#include <stdio.h>
#include <iostream>
#include <sstream>

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
    return A * PI / 180.;
}

double grad(double A)
{
    return A * 180. / PI;
}

coord grad(coord A)
{
    coord result;
    result.fi = grad(A.fi);
    result.lam = grad(A.lam);
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
    double cgf = rad(cord.fi);
    double cgl = rad(cord.lam);
    result.b11 = cos(cgf) * cos(rad(cord.lam));
    result.b21 = cos(cgf) * sin(rad(cord.lam));
    result.b31 = sin(cgf);
    return result;
}

Bas RotBasMult (Rot rot, Bas bas)
{
    Bas result;
    result.b11 = rot.a11 * bas.b11 + rot.a12 * bas.b21 + rot.a13 * bas.b31;
    result.b21 = rot.a21 * bas.b11 + rot.a22 * bas.b21 + rot.a23 * bas.b31;
    result.b31 = rot.a31 * bas.b11 + rot.a32 * bas.b21 + rot.a33 * bas.b31;
    return result;
}

double equationOfTime(double date)
{
    double E, B;
    B = 2. * PI * (date / T);
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

    Bas fi1 = RotBasMult(Rx(epsilon), basis(ecl));

    ecv.fi =   asin(fi1.b31);
    ecv.lam = atan2(fi1.b21 / cos(ecv.fi), fi1.b11 / cos(ecv.fi));
    ecv.fi =  grad(ecv.fi);
    ecv.lam = grad(ecv.lam);
    return ecv;
}

coord ecvatorToHorisont(coord ecv, double time, double date)
{
    coord horisont;
    coord tmp;

    Bas fi1 = RotBasMult(Rz(S(time, date)), basis(ecv));

    tmp.fi =   asin(fi1.b31);
    tmp.lam = atan2(fi1.b21 / cos(tmp.fi), fi1.b11 / cos(tmp.fi));
    tmp.fi =  grad(tmp.fi);
    tmp.lam = grad(tmp.lam);

    fi1 = RotBasMult(Ry((90 - latitude)), basis(tmp));

    horisont.fi =   asin(fi1.b31);
    horisont.lam = atan2(fi1.b21 / cos(horisont.fi), fi1.b11 / cos(horisont.fi));
    horisont.fi =  grad(horisont.fi);
    horisont.lam = grad(horisont.lam);

    return horisont;
}

void GnuOut(FILE *f, int time){
    fprintf(f, "#set size ratio 1 \n");
    fprintf(f, "set terminal png size 2000,900  font "",20"" \n");
    fprintf(f, "#set autoscale fix \n");
    fprintf(f, "set xzeroaxis \nset yzeroaxis \n");
    fprintf(f, "#set terminal png\n");
    fprintf(f, "#set xrange [-1.2:1.2] \n#set yrange [-1.2:1.2]\n");
//    fprintf(f, "set output 'analemma%d.png' \n", time);
    fprintf(f, "set output 'analemma.png' \n");
    fprintf(f, "plot ");
    }

void GnuOut1(FILE *f, int time){
    fprintf(f, "'analemma%d.dat' u ($2/1):($1/1) w l  notitle, ", time);
    }

void GnuOut(FILE *f){
    fprintf(f,"\n");
    fprintf(f,"pause -1");
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
    FILE * fp = fopen("analemma.gnu","w");

    GnuOut(fp, time);
    double r360 = 0;
    double lamprev = 0;
    for(time = 0; time < 24; time++)
    {
        std::stringstream filename;
        filename<<"analemma"<<time<<".dat";
        FILE * f = fopen(filename.str().c_str(),"w");
        GnuOut1(fp,time);
        for( t = 1; t < T+1 ; t++ )
        {
            M = (2 * PI * t / T);
            for( j = 0; j < 200; j++ )
            {
                E = E - ((E - e * sin(E) - M)/(1 - e * cos(E)));
            }
            V = 2 * atan(sqrt ((1 + e)/(1 - e)) * tan (E / 2));
            ecl.fi = 0;
            ecl.lam = grad(V) + TVR;

            ecv = eklipicToEcvator(ecl);
            hor = ecvatorToHorisont(ecv, time, t);
            printf("%f %f \n", hor.lam,lamprev);
            if((t != 1) && (fabs(hor.lam - lamprev) > 180))
                r360 += 360 * (hor.lam - lamprev > 0 ? -1 : 1);
            lamprev = hor.lam;
            fprintf(f, "%f %f \n", hor.fi, hor.lam + r360);
        }
        fclose (f);
    }
    GnuOut(fp);
    fclose (fp);
    return 0;
}


