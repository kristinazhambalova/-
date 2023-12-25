#include <stdio.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include <fstream>

int testTochnost(double E1, double E2, double tochnost) {
    if (fabs(E1 - E2) <= tochnost) { return 1; }
    else { return 0; }
}

double EkscentrAnomIter(double Ei, double M, double e, double tochnost) {
    double EI = Ei;
    Ei = M + e * sin(EI);

    if (testTochnost(Ei, EI, tochnost) == 1) { return Ei; }
    EkscentrAnomIter(Ei, M, e, tochnost);
}

int main() {

    double R, rp, ra, E, M, tochnost, An, t;

    std::cout << "Введите радиус перигея, радиус апогея, радиус планеты, Точность (0.0001)" << std::endl;
    std::cin >> rp >> ra >> R >> tochnost;

    double e = (ra - rp) / (ra + rp + R * 2);

    std::ofstream outTime, outM, outE, outAn;
    outM.precision(200);
    outE.precision(200);
    outTime.precision(200);
    outAn.precision(200);

    outE.open("E.txt", std::ios::app);
    outM.open("M.txt", std::ios::app);
    outTime.open("T.txt", std::ios::app);
    outAn.open("An.txt", std::ios::app);

    if (outAn and outE and outM and outTime) {
        for (t = 0; t <= 400; t++) {
            M = 2 * M_PI * (t / 400);
            E = EkscentrAnomIter(M, M, e, tochnost);
            An = M_PI + 2 * atan(pow((1 + e) / (1 - e), 0.5) * tan((M_PI + E) / 2));

            outTime << t / 400 << std::endl;
            outM << M / M_PI << std::endl;
            outE << E / M_PI << std::endl;
            outAn << An / M_PI << std::endl;
        }
        outTime.close();
        outM.close();
        outE.close();
        outAn.close();
    }
}