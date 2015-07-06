#include <iostream>
#include <vector>
#include <cmath>
#include <bitset>

using namespace std;

struct PairingSP
{
public:
    PairingSP(int, bool);
    int p; // from 1
    bool spin; //0 -- spin up, 1 -- spin down
};

PairingSP::PairingSP(int _p, bool _spin) : p(_p), spin(_spin)
{
}

class Pairing
{
public:
    vector<PairingSP> SPStates;
    void generateSP();
    int PMax;
};

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

template <class T>
class System
{
public:
    System(int, int);
    void generateSPOrbitals();
    void generateConfigurations();
    void printConfigurations();
private:
    vector<T> SPOrbitals;
    vector<Configuration> configurations;
    int numberSPStates;
    int A; // have to be even
};

System::System(int _A, int _numberSPStates) : A(_A), numberSPStates(_numberSPStates)
{
    //SPOrbitals = vector<SPOrbital>(_numberSPStates);
}

void System::generateSPOrbitals()
{
    generatePairingSP(numberSPStates);
}

void System::generatePairingSP(int numberSPStates)
{
    for (int p = 1; p < numberSPStates/2; p++)
    {
        SPOrbitals.push_back(PairingSP(p,0));
        SPOrbitals.push_back(PairingSP(p,1));
    }
}

void System::generateConfigurations()
{
    int config = first(A);
    int configLast = last(numberSPStates,A);
    configurations.push_back(Configuration(config));
    do
    {
        config = next(config);
        //check
        configurations.push_back(Configuration(config));
    }
    while(config != configLast);
}

void System::printConfigurations()
{
    for (int i = 0; i < configurations.size(); i++)
        cout << bitset<8>(configurations.at(i).configuration) << endl;
}

int main()
{
    cout << "Hello world" << endl;
    int g = 10;
    System system(4,8);
    system.generateConfigurations();
    system.printConfigurations();


    return 0;
}
