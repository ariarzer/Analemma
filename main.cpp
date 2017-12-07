#include "stdafx.h"
#include "analemma.h"

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
    time = 12;
    double r360 = 0;
    double lamprev = 0;

    generateGnuScript("analemma", 1, 1000, 800, 12, 0, 0, 0, 0,
                      "Right ascension", "Decline", "Analemma in the horizontal coordinate system");

    std::ofstream file(("analemma.dat"), std::ios_base::out);
    if (!file.is_open())
    {
        std::cout << "Error open file" << std::endl;
        return 0;
    }
    for( t = 1; t < period + 1 ; t++ )
    {
        M = (2 * M_PI * t / period);
        for( j = 0; j < 200; j++ )
        {
            E = E - ((E - e * sin(E) - M)/(1 - e * cos(E)));
        }
        V = 2 * atan(sqrt ((1 + e)/(1 - e)) * tan (E / 2));
        ecl.fi = 0;
        ecl.lam = toGrad(V) + TVR;

        ecv = eklipicToEcvator(ecl);
        hor = ecvatorToHorisont(ecv, time, t);

        if((t != 1) && (fabs(hor.lam - lamprev) > 180))
            r360 += 360 * (hor.lam - lamprev > 0 ? -1 : 1);
        lamprev = hor.lam;

        file <<  hor.lam << ' ' <<  hor.fi + r360 << std::endl;
    }
    file.close();
    return 0;
}
