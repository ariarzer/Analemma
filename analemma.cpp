#include "analemma.h"

//Transformatiom rad-to-grad and grad-to-rad
double toRad(double A)
{
    return A * M_PI / 180.;
}

double toGrad(double A)
{
    return A * 180. / M_PI;
}

coord toRad(coord A)
{
    coord tmp;
    tmp.fi = toRad(A.fi);
    tmp.lam = toRad(A.lam);
    return tmp;
}

coord toGrad(coord A)
{
    coord tmp;
    tmp.fi = toGrad(A.fi);
    tmp.lam = toGrad(A.lam);
    return tmp;
}

//Initialization rotatin-matrix
Rot Rx(double a)
{
    Rot Rx = {1, 0,            0,
              0, cos(toRad(a)),  sin(toRad(a)),
              0, -sin(toRad(a)), cos(toRad(a))};
    return Rx;
}

Rot Ry(double a)
{
    Rot Ry = {cos(toRad(a)), 0, -sin(toRad(a)),
              0,           1, 0,
              sin(toRad(a)), 0, cos(toRad(a))};
    return Ry;
}

Rot Rz(double a)
{
    Rot Rz = {cos(toRad(a)),  sin(toRad(a)),  0,
              -sin(toRad(a)),  cos(toRad(a)),  0,
              0,            0,            1};
    return Rz;
}

//Initialization coordinate-basis
Bas basis(coord cord)
{
    Bas result;
    coord tmp = toRad(cord);
    result.b11 = cos(tmp.fi) * cos(tmp.lam);
    result.b21 = cos(tmp.fi) * sin(tmp.lam);
    result.b31 = sin(tmp.fi);
    return result;
}

//Multiplication rotatin-matrix by coordinate-basis
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
    B = 2. * M_PI * (date / period);
    E = 7.53 * cos(B) + 1.5 * sin(B) - 9.87 * sin(2 * B);
    return E;
}

double S(double time, double date)
{
    double S;
    S = time + 24. * date / period + equationOfTime(date) / 60. + 12;
    if (S > 12) S = S - 24;
    S = S * 360. / 24.;
    return S;
}

coord eklipicToEcvator(coord ecl)
{
    coord ecv;

    Bas fi1 = RotBasMult(Rx(epsilon), basis(ecl));

    ecv.fi =  asin(fi1.b31);
    ecv.lam = atan2(fi1.b21 / cos(ecv.fi), fi1.b11 / cos(ecv.fi));
    ecv.fi =  toGrad(ecv.fi);
    ecv.lam = toGrad(ecv.lam);
    return ecv;
}

coord ecvatorToHorisont(coord ecv, double time, double date)
{
    coord horisont, tmp;

    Bas fi1 = RotBasMult(Rz(S(time, date)), basis(ecv));

    tmp.fi =  asin(fi1.b31);
    tmp.lam = atan2(fi1.b21 / cos(tmp.fi), fi1.b11 / cos(tmp.fi));
    tmp = toGrad(tmp);

    fi1 = RotBasMult(Ry((90 - latitude)), basis(tmp));

    horisont.fi =  asin(fi1.b31);
    horisont.lam = atan2(fi1.b21 / cos(horisont.fi), fi1.b11 / cos(horisont.fi));
    horisont = toGrad(horisont);

    return horisont;
}

int generateGnuScript(std::string filename, int number,
                      int sizex_image = 1000, int sizey_image = 800, int size_font = 12,
                      int sizex_diagramm_from = 0, int sizex_diagramm_to = 0,
                      int sizey_diagramm_from = 0, int sizey_diagramm_to = 0,
                      std::string xlabel = "", std::string ylabel = "",
                      std::string title = "")
{
    std::ofstream file((filename + ".gnu"), std::ios_base::out);
    if (!file.is_open())
    {
        std::cout << "Error open file" << std::endl;
        return 0;
    }
    file << "set size ratio 1" << std::endl
         << "set terminal png size " << sizex_image << "," << sizey_image
         << " font '" << size_font << "'" << std::endl
         << "set xzeroaxis" << std::endl << "set yzeroaxis" << std::endl
         << "set terminal png" << std::endl;
    if ((sizex_diagramm_from != 0) &&  (sizex_diagramm_to != 0))
        file << "set xrange [" << sizex_diagramm_from << ":" << sizex_diagramm_to << "]" << std::endl;
    if ((sizey_diagramm_from != 0) &&  (sizey_diagramm_to != 0))
        file << "set yrange [" << sizey_diagramm_from << ":" << sizey_diagramm_to << "]" << std::endl;
    if (xlabel != "")
        file << "set xlabel '" << xlabel << "'" << std::endl;
    if (ylabel != "")
        file << "set ylabel '" << ylabel << "'" << std::endl;
    if (title != "")
        file << "set title '" << title << "'" << std::endl << std::endl;
    for (int i = 0; i < number ; i++)
        file << "set output '" << (filename + ".png") << std::endl
             << "plot " << (filename + ".dat") << " u ($" << (i + 1) << "/1):($" << (i + 2) << "1) w l  notitle" << std::endl;
    file << "set output '" << (filename + ".png") << std::endl
         << "plot " << (filename + ".dat") << " u ($" << (1) << "/1):($" << (2) << "1) w l  notitle" << std::endl;
    file << "pause -1";
    file.close();
    return 0;
}
