#include <iostream>
#include <vector>
#include <cmath>
#include <bitset>
#include "help.h"
#include "diag.h"
#include "Eigen/Dense"

using namespace std;

struct PairingSP
{
public:
    PairingSP(int, bool);
    int p; // from 1
    bool spin; // 0 -- spin up, 1 -- spin down
    double spEnergy; // single particle energy
};

PairingSP::PairingSP(int _p, bool _spin) : p(_p), spin(_spin)
{
    spEnergy = p-1;
}

class Pairing
{
public:
    vector<PairingSP> SPStates;
    void generateSP();
    int PMax;
    double g;
};

void Pairing::generateSP()
{
    g = 0.2;
    PMax = 4; //TODO
    for (int p = 1; p <= PMax; p++)
    {
        SPStates.push_back(PairingSP(p,0));
        SPStates.push_back(PairingSP(p,1));
    }
}

class Configuration
{
public:
    Configuration(int);
    int configuration;
};

Configuration::Configuration(int config)
{
    configuration = config;
}

class System
{
public:
    System(int);
    void generateSPOrbitals();
    void generateConfigurations();
    void printConfigurations();
    Pairing spOrbitals;
    double VMatix[8][8][8][8];//TODO
    void generateVMatrix();
    double deltaE();
    void diagonalization();
    void setG(double);
    vector<double> coefficients();
    double getCoefficient(int,int,int,int);
private:
    vector<Configuration> configurations;
    int A; // number of particles, have to be even
};

void System::setG(double _g)
{
    spOrbitals.g = _g;
    generateVMatrix();
}

System::System(int _A) : A(_A)
{
}

void System::diagonalization()
{
    MatrixXd H;
    generate_H(H, spOrbitals.g);
    double E_GS;
    diag(H, E_GS);
    cout << E_GS << endl;
}

void System::generateSPOrbitals()
{
    spOrbitals.generateSP();
}

void System::generateVMatrix()
{
    for (int a = 0; a < A; a++)
    {
        for (int b = 0; b < A; b++)
        {
            for (int i = A; i < spOrbitals.SPStates.size(); i++)
            {
                for (int j = A; j < spOrbitals.SPStates.size(); j++)
                {
                    if ((spOrbitals.SPStates[a].p == spOrbitals.SPStates[b].p) &&
                            (spOrbitals.SPStates[a].spin != spOrbitals.SPStates[b].spin) &&
                            (spOrbitals.SPStates[i].p == spOrbitals.SPStates[j].p) &&
                            (spOrbitals.SPStates[i].spin != spOrbitals.SPStates[j].spin))
                    {
                        VMatix[a][b][i][j] = -0.5*spOrbitals.g;
                        VMatix[i][j][a][b] = -0.5*spOrbitals.g;
                    }
                    else
                    {
                        VMatix[a][b][i][j] = 0.0;
                        VMatix[i][j][a][b] = 0.0;
                    }
                }
            }
        }
    }
}

vector<double> System::coefficients()
{
    vector<double> _coefficients;
    _coefficients.push_back(getCoefficient(2,3,4,5));
    _coefficients.push_back(getCoefficient(2,3,6,7));
    _coefficients.push_back(getCoefficient(0,1,4,5));
    _coefficients.push_back(getCoefficient(0,1,6,7));
    for (int i = 0; i < _coefficients.size(); i++)
        cout << _coefficients.at(i) << "\t";
    cout << endl;
    return _coefficients;
}

double System::getCoefficient(int a, int b, int i, int j)
{
    double Ea = spOrbitals.SPStates.at(a).spEnergy;
    double Eb = spOrbitals.SPStates.at(b).spEnergy;
    double Ei = spOrbitals.SPStates.at(i).spEnergy;
    double Ej = spOrbitals.SPStates.at(j).spEnergy;
    return VMatix[i][j][a][b] / (Ei+Ej-Ea-Eb);
}

double System::deltaE()
{
    double _deltaE = 0.0;
    for (int a = 0; a < A; a++)
    {
        for (int b = 0; b < A; b++)
        {
            for (int i = A; i < spOrbitals.SPStates.size(); i++)
            {
                for (int j = A; j < spOrbitals.SPStates.size(); j++)
                {
                    double Ea = spOrbitals.SPStates.at(a).spEnergy;
                    double Eb = spOrbitals.SPStates.at(b).spEnergy;
                    double Ei = spOrbitals.SPStates.at(i).spEnergy;
                    double Ej = spOrbitals.SPStates.at(j).spEnergy;
                    _deltaE += VMatix[a][b][i][j] * VMatix[i][j][a][b] / (Ei+Ej-Ea-Eb);
                }
            }
        }
    }
    return 0.25*_deltaE;
}

//void System::generateConfigurations()
//{
//    int config = first(A);
//    int configLast = last(spOrbitals.SPStates.size(),A);
//    configurations.push_back(Configuration(config));
//    do
//    {
//        config = next(config);
//        //check
//        configurations.push_back(Configuration(config));
//    }
//    while(config != configLast);
//}

//void System::printConfigurations()
//{
//    for (int i = 0; i < configurations.size(); i++)
//        cout << bitset<8>(configurations.at(i).configuration) << endl;
//}

int main()
{
    System system(4);
//    system.generateConfigurations();
//    system.printConfigurations();
    system.generateSPOrbitals();

//    for (double g = -2.0; g <= 2.0; g += 0.1)
//    {
//        system.setG(g);
//        cout << g << "\t" << 2-system.spOrbitals.g-system.deltaE() << "\t";
//        system.diagonalization();
//    }
    system.setG(1);
    system.coefficients();
    return 0;
}
