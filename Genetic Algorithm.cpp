// Genetic Algorithm.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include<vector>
#include<math.h>
#include<map>
#include<time.h>
#include<stdio.h>
#include<stdlib.h>
constexpr auto a1 = 5678.83;
constexpr auto a2 = 0.212;
constexpr auto mu0 = 1.26e-06;
using namespace std;
long double* function(long double* phi, long double* H, long double* tau, long double* M, short int datapoints, long double* DecisionVariableB, long double* DecisionVariableC, long double* DecisionVariableD, long double* DecisionVariableE, long double* DecisionVariableF, short int popsize)
{
    long double error, A, taupredicted;
    long double* Error = new long double[popsize];
    for (int i = 0; i < popsize; i++)
    {
        error = 0;
        //cout << DecisionVariableB[i] << "  " << DecisionVariableC[i] << "  " << DecisionVariableD[i] << "  " << DecisionVariableE[i] << "  " << DecisionVariableF[i] << "\n";
        for (int j = 0; j < datapoints; j++)
        {
            A = (H[j]) / (a1 * pow((phi[j] * M[j]), a2));
            //DecisionVariableB[i] = -1 * DecisionVariableB[i];
            //DecisionVariableD[i] = -1 * DecisionVariableD[i];
            taupredicted = ((pow(A, DecisionVariableB[i])) / ((pow(A, DecisionVariableB[i])) + 1)) * ((DecisionVariableC[i]) * (pow(phi[j], DecisionVariableD[i])) * mu0 * (pow(M[j], DecisionVariableE[i]) * (pow(H[j], DecisionVariableF[i]))));
            //taupredicted = (A / (A + 1)) * 2.5 * phi[j] * mu0 * (pow(M[j], 0.5)) * (pow(H[j], 1.5));
            //cout << "\nA = " << A;
            //cout <<" \n  taupredicted = "<< taupredicted<<"\n";
            error += abs(tau[j] - taupredicted);
        }
        Error[i] = error;
    }
    return Error;
}
long double DV(long double* Initialpop, int i, int j)
{
    int k = 0;
    long double DecodedValue = 0;
    while (j > i)
    {
        DecodedValue += Initialpop[j - 1] * pow(2, k);
        k++;
        j--;
    }
    return(DecodedValue);
}
void VariableCreation(short int popsize, long double* DecisionVariableB, long double* DecisionVariableC, long double* DecisionVariableD, long double* DecisionVariableE, long double* DecisionVariableF, long double* InitialbinarypopB, long double* InitialbinarypopC, long double* InitialbinarypopD, long double* InitialbinarypopE, long double* InitialbinarypopF, short int l, long double lowerboundB, long double lowerboundC, long double lowerboundD, long double lowerboundE, long double lowerboundF, long double upperboundB, long double upperboundC, long double upperboundD, long double upperboundE, long double upperboundF)
{
    int j = 0;
    for (int i = 0; i < popsize; i++)
    {
        DecisionVariableB[i] = lowerboundB + (((upperboundB - lowerboundB) / (pow(2, l) - 1)) * DV(InitialbinarypopB, j, j + l));
        //cout << DecisionVariableB[i] << "\n";
        DecisionVariableC[i] = lowerboundC + (((upperboundC - lowerboundC) / (pow(2, l) - 1)) * DV(InitialbinarypopC, j, j + l));
        //cout << DecisionVariableC[i] << "\n";
        DecisionVariableD[i] = lowerboundD + (((upperboundD - lowerboundD) / (pow(2, l) - 1)) * DV(InitialbinarypopD, j, j + l));
        //cout << DecisionVariableD[i] << "\n";
        DecisionVariableE[i] = lowerboundE + (((upperboundE - lowerboundE) / (pow(2, l) - 1)) * DV(InitialbinarypopE, j, j + l));
        //cout << DecisionVariableE[i] << "\n";
        DecisionVariableF[i] = lowerboundF + (((upperboundF - lowerboundF) / (pow(2, l) - 1)) * DV(InitialbinarypopF, j, j + l));
        //cout << DecisionVariableF[i] << "\n";
        j += l;
    }
}
void TournamentSelectionOperation(short int popsize, long double* error, long double* DecisionVariableB, long double* DecisionVariableC, long double* DecisionVariableD, long double* DecisionVariableE, long double* DecisionVariableF, long double* InitialbinarypopB, long double* InitialbinarypopC, long double* InitialbinarypopD, long double* InitialbinarypopE, long double* InitialbinarypopF, short int l)
{
    vector<pair<int, long double>> vect;
    long double* copyB = new long double[popsize / 2];
    long double* copyC = new long double[popsize / 2];
    long double* copyD = new long double[popsize / 2];
    long double* copyE = new long double[popsize / 2];
    long double* copyF = new long double[popsize / 2];
    long double* copyerror = new long double[popsize / 2];
    long double* binarycopyB = new long double[l * popsize / 2];
    long double* binarycopyC = new long double[l * popsize / 2];
    long double* binarycopyD = new long double[l * popsize / 2];
    long double* binarycopyE = new long double[l * popsize / 2];
    long double* binarycopyF = new long double[l * popsize / 2];
    for (int i = 0; i < popsize; i++)
    {
        vect.push_back(make_pair(i, error[i]));
    }
    vector<pair<int, long double>> result;
    while (!(vect.empty()))
    {
        srand(time(0));
        int i = rand() % vect.size();
        int j = rand() % vect.size();
        if (i != j)
        {
            if (vect[i].second > vect[j].second)
            {
                result.push_back(vect[j]);
            }
            else
            {
                result.push_back(vect[i]);
            }
            vect.erase(vect.begin() + i);
            if (i < j)
            {
                vect.erase(vect.begin() + j - 1);
            }
            else
            {
                vect.erase(vect.begin() + j);
            }
        }
    }
    /*for (int i = 0; i < result.size(); i++)
    {
        cout << result[i].first << "  " << result[i].second << "\n";
    }*/
    for (int i = 0; i < popsize / 2; i++)
    {
        copyB[i] = DecisionVariableB[result[i].first];
        copyC[i] = DecisionVariableC[result[i].first];
        copyD[i] = DecisionVariableD[result[i].first];
        copyE[i] = DecisionVariableE[result[i].first];
        copyF[i] = DecisionVariableF[result[i].first];
        copyerror[i] = result[i].second;
        for (int j = 0; j < l; j++)
        {
            binarycopyB[i * l + j] = InitialbinarypopB[(result[i].first) * l + j];
            binarycopyC[i * l + j] = InitialbinarypopC[(result[i].first) * l + j];
            binarycopyD[i * l + j] = InitialbinarypopD[(result[i].first) * l + j];
            binarycopyE[i * l + j] = InitialbinarypopE[(result[i].first) * l + j];
            binarycopyF[i * l + j] = InitialbinarypopF[(result[i].first) * l + j];
        }
    }
    int j = 0;
    for (int i = 0; i < popsize; i++)
    {
        DecisionVariableB[i] = copyB[j];
        DecisionVariableC[i] = copyC[j];
        DecisionVariableD[i] = copyD[j];
        DecisionVariableE[i] = copyE[j];
        DecisionVariableF[i] = copyF[j];
        error[i] = copyerror[j];
        for (int k = 0; k < l; k++)
        {
            InitialbinarypopB[i * l + k] = binarycopyB[j * l + k];
            InitialbinarypopC[i * l + k] = binarycopyC[j * l + k];
            InitialbinarypopD[i * l + k] = binarycopyD[j * l + k];
            InitialbinarypopE[i * l + k] = binarycopyE[j * l + k];
            InitialbinarypopF[i * l + k] = binarycopyF[j * l + k];
        }
        j++;
        if (j == popsize / 2)
        {
            j = 0;
        }

    }
}
void SinglePointCrossover(long double* InitialbinarypopB, long double* InitialbinarypopC, long double* InitialbinarypopD, long double* InitialbinarypopE, long double* InitialbinarypopF, short int l, short int popsize)
{
    vector<vector<long double>>vectB(popsize);
    vector<vector<long double>>vectC(popsize);
    vector<vector<long double>>vectD(popsize);
    vector<vector<long double>>vectE(popsize);
    vector<vector<long double>>vectF(popsize);
    for (int i = 0; i < popsize; i++)
    {
        for (int j = 0; j < l; j++)
        {
            vectB[i].push_back(InitialbinarypopB[i * l + j]);
            vectC[i].push_back(InitialbinarypopC[i * l + j]);
            vectD[i].push_back(InitialbinarypopD[i * l + j]);
            vectE[i].push_back(InitialbinarypopE[i * l + j]);
            vectF[i].push_back(InitialbinarypopF[i * l + j]);
        }
    }
    srand(time(0));
    while (!vectB.empty())
    {
        int j = rand() % vectB.size();
        int k = rand() % vectB.size();
        int pc = rand() % 11;
        if (pc <= 9 && j != k)
        {
            int p = rand() % vectB[j].size();
            for (int i = p; i < l; i++)
            {
                long double c;
                c = vectB[j][i];
                vectB[j][i] = vectB[k][i];
                vectB[k][i] = c;
            }
            for (int i = 0; i < l; i++)
            {
                InitialbinarypopB[j * l + i] = vectB[j][i];
                InitialbinarypopB[k * l + i] = vectB[k][i];
            }
        }
        if (j > k)
        {
            vectB.erase(vectB.begin() + j);
            vectB.erase(vectB.begin() + k);
        }
        else if (j < k)
        {
            vectB.erase(vectB.begin() + k);
            vectB.erase(vectB.begin() + j);
        }
    }
    srand(time(0));
    while (!vectC.empty())
    {
        int j = rand() % vectC.size();
        int k = rand() % vectC.size();
        int pc = rand() % 11;
        if (pc <= 9 && j != k)
        {
            int p = rand() % vectC[j].size();
            for (int i = p; i < l; i++)
            {
                long double c;
                c = vectC[j][i];
                vectC[j][i] = vectC[k][i];
                vectC[k][i] = c;
            }
            for (int i = 0; i < l; i++)
            {
                InitialbinarypopC[j * l + i] = vectC[j][i];
                InitialbinarypopC[k * l + i] = vectC[k][i];
            }

        }
        if (j > k)
        {
            vectC.erase(vectC.begin() + j);
            vectC.erase(vectC.begin() + k);
        }
        else if (j < k)
        {
            vectC.erase(vectC.begin() + k);
            vectC.erase(vectC.begin() + j);
        }
    }
    srand(time(0));
    while (!vectD.empty())
    {
        int j = rand() % vectD.size();
        int k = rand() % vectD.size();
        int pc = rand() % 11;
        if (pc <= 9 && j != k)
        {
            int p = rand() % vectD[j].size();
            for (int i = p; i < l; i++)
            {
                long double c;
                c = vectD[j][i];
                vectD[j][i] = vectD[k][i];
                vectD[k][i] = c;
            }
            for (int i = 0; i < l; i++)
            {
                InitialbinarypopD[j * l + i] = vectD[j][i];
                InitialbinarypopD[k * l + i] = vectD[k][i];
            }
        }
        if (j > k)
        {
            vectD.erase(vectD.begin() + j);
            vectD.erase(vectD.begin() + k);
        }
        else if (j < k)
        {
            vectD.erase(vectD.begin() + k);
            vectD.erase(vectD.begin() + j);
        }
    }
    srand(time(0));
    while (!vectE.empty())
    {
        int j = rand() % vectE.size();
        int k = rand() % vectE.size();
        int pc = rand() % 11;
        if (pc <= 9 && j != k)
        {
            int p = rand() % vectE[j].size();
            for (int i = p; i < l; i++)
            {
                long double c;
                c = vectE[j][i];
                vectE[j][i] = vectE[k][i];
                vectE[k][i] = c;
            }
            for (int i = 0; i < l; i++)
            {
                InitialbinarypopE[j * l + i] = vectE[j][i];
                InitialbinarypopE[k * l + i] = vectE[k][i];
            }
        }
        if (j > k)
        {
            vectE.erase(vectE.begin() + j);
            vectE.erase(vectE.begin() + k);
        }
        else if (j < k)
        {
            vectE.erase(vectE.begin() + k);
            vectE.erase(vectE.begin() + j);
        }
    }
    srand(time(0));
    while (!vectF.empty())
    {
        int j = rand() % vectF.size();
        int k = rand() % vectF.size();
        int pc = rand() % 11;
        if (pc <= 9 && j != k)
        {
            int p = rand() % vectF[j].size();
            for (int i = p; i < l; i++)
            {
                long double c;
                c = vectF[j][i];
                vectF[j][i] = vectF[k][i];
                vectF[k][i] = c;
            }
            for (int i = 0; i < l; i++)
            {
                InitialbinarypopF[j * l + i] = vectF[j][i];
                InitialbinarypopF[k * l + i] = vectF[k][i];
            }
        }
        if (j > k)
        {
            vectF.erase(vectF.begin() + j);
            vectF.erase(vectF.begin() + k);
        }
        else if (j < k)
        {
            vectF.erase(vectF.begin() + k);
            vectF.erase(vectF.begin() + j);
        }
    }
}
void BitwiseMutation(long double* InitialbinarypopB, long double* InitialbinarypopC, long double* InitialbinarypopD, long double* InitialbinarypopE, long double* InitialbinarypopF, short int popsize, short int l)
{
    srand(time(0));
    for (int i = 0; i < popsize; i++)
    {
        int j = rand() % 11;
        if (j <= 1)
        {
            int k = rand() % l;
            if (InitialbinarypopB[i * l + k] == 0)
            {
                InitialbinarypopB[i * l + k] = 1;
            }
            else
            {
                InitialbinarypopB[i * l + k] = 0;
            }
        }
    }
    for (int i = 0; i < popsize; i++)
    {
        int j = rand() % 11;
        if (j <= 1)
        {
            int k = rand() % l;
            if (InitialbinarypopC[i * l + k] == 0)
            {
                InitialbinarypopC[i * l + k] = 1;
            }
            else
            {
                InitialbinarypopC[i * l + k] = 0;
            }
        }
    }
    for (int i = 0; i < popsize; i++)
    {
        int j = rand() % 11;
        if (j <= 1)
        {
            int k = rand() % l;
            if (InitialbinarypopD[i * l + k] == 0)
            {
                InitialbinarypopD[i * l + k] = 1;
            }
            else
            {
                InitialbinarypopD[i * l + k] = 0;
            }
        }
    }
    for (int i = 0; i < popsize; i++)
    {
        int j = rand() % 11;
        if (j <= 1)
        {
            int k = rand() % l;
            if (InitialbinarypopE[i * l + k] == 0)
            {
                InitialbinarypopE[i * l + k] = 1;
            }
            else
            {
                InitialbinarypopE[i * l + k] = 0;
            }
        }
    }
    for (int i = 0; i < popsize; i++)
    {
        int j = rand() % 11;
        if (j <= 1)
        {
            int k = rand() % l;
            if (InitialbinarypopF[i * l + k] == 0)
            {
                InitialbinarypopF[i * l + k] = 1;
            }
            else
            {
                InitialbinarypopF[i * l + k] = 0;
            }
        }
    }
}
/*void Survival(short int popsize, long double* DecisionVariableB, long double* DecisinoVariableC, long double* DecisionVariableD, long double* DecisionVariableE, long double* DecisionVariableF, long double* error)
{
    map<long double, int> mp;
    for (int i = 0; i < popsize; i++)
    {
        mp.insert({ error[i],i });
    }
}*/
int main()
{
    short int datapoints, popsize;
    cout << "Enter the number of Data points: ";
    cin >> datapoints;
    long double* phi = new long double[datapoints];
    long double* H = new long double[datapoints];
    long double* tau = new long double[datapoints];
    long double* M = new long double[datapoints];
    ifstream indata;
    ofstream Results;
    long double data;
    indata.open("Data.txt");
    indata >> data;
    for (int i = 0; i < datapoints; i++)                 //Takes Input from File Location Entered.
    {
        phi[i] = data;
        indata >> data;
        H[i] = data;
        indata >> data;
        tau[i] = data;
        indata >> data;
        M[i] = data;
        indata >> data;
    }
    indata.close();
    /*for (int i = 0; i < datapoints; i++)
    {
        cout << phi[i] << "  " << H[i] << "  " << tau[i] << "  " << M[i] << "  "<<"\n";
    } */          //This is for checking whether the input is proper or not!!.
    long double lowerboundB, lowerboundC, lowerboundD, lowerboundE, upperboundB, upperboundC, upperboundD, upperboundE, lowerboundF, upperboundF;
    short int l;
    long double precisionB, precisionC, precisionD, precisionE, precisionF;
    cout << "Enter the String Length(even): ";
    cin >> l;
    cout << "Enter the lower bound for B: ";
    cin >> lowerboundB;
    cout << "Enter the upper bound for B: ";
    cin >> upperboundB;
    cout << "Enter the lower bound for C: ";
    cin >> lowerboundC;
    cout << "Enter the upper bound for C: ";
    cin >> upperboundC;
    cout << "Enter the lower bound for D: ";
    cin >> lowerboundD;
    cout << "Enter the upper bound for D: ";
    cin >> upperboundD;
    cout << "Enter the lower bound for E: ";
    cin >> lowerboundE;
    cout << "Enter the upper bound for E: ";
    cin >> upperboundE;
    cout << "Enter the lower bound for F: ";
    cin >> lowerboundF;
    cout << "Enter the upper bound for F: ";
    cin >> upperboundF;
    precisionB = ((upperboundB - lowerboundB) / (pow(2, l) - 1));
    precisionC = ((upperboundC - lowerboundC) / (pow(2, l) - 1));
    precisionD = ((upperboundD - lowerboundD) / (pow(2, l) - 1));
    precisionE = ((upperboundE - lowerboundE) / (pow(2, l) - 1));
    precisionF = ((upperboundF - lowerboundF) / (pow(2, l) - 1));
    cout << "\nEnter the Population Size: ";
    cin >> popsize;
    long double* Error = new long double[popsize];
    long double* DecisionVariableB = new long double[popsize];
    long double* DecisionVariableC = new long double[popsize];
    long double* DecisionVariableD = new long double[popsize];
    long double* DecisionVariableE = new long double[popsize];
    long double* DecisionVariableF = new long double[popsize];
    long double* InitialbinarypopB = new long double[l * popsize];
    long double* InitialbinarypopC = new long double[l * popsize];
    long double* InitialbinarypopD = new long double[l * popsize];
    long double* InitialbinarypopE = new long double[l * popsize];
    long double* InitialbinarypopF = new long double[l * popsize];
    srand(time(0));
    for (int i = 0; i < (l * popsize); i++)                       //Initializing binary Population
    {
        InitialbinarypopB[i] = rand() % 2;
        InitialbinarypopC[i] = rand() % 2;
        InitialbinarypopD[i] = rand() % 2;
        InitialbinarypopE[i] = rand() % 2;
        InitialbinarypopF[i] = rand() % 2;
    }
    cout << "\n";
    int j = 0;
    VariableCreation(popsize, DecisionVariableB, DecisionVariableC, DecisionVariableD, DecisionVariableE, DecisionVariableF, InitialbinarypopB, InitialbinarypopC, InitialbinarypopD, InitialbinarypopE, InitialbinarypopF, l, lowerboundB, lowerboundC, lowerboundD, lowerboundE, lowerboundF, upperboundB, upperboundC, upperboundD, upperboundE, upperboundF);
    /*for (int i = 0; i < popsize; i++)
    {
        DecisionVariableB[i] = lowerboundB + (((upperboundB - lowerboundB) / (pow(2, l) - 1)) * DV(InitialbinarypopB, j, j + l));
        //cout << DecisionVariableB[i] << "\n";
        DecisionVariableC[i] = lowerboundC + (((upperboundC - lowerboundC) / (pow(2, l) - 1)) * DV(InitialbinarypopC, j, j + l));
        //cout << DecisionVariableC[i] << "\n";
        DecisionVariableD[i] = lowerboundD + (((upperboundD - lowerboundD) / (pow(2, l) - 1)) * DV(InitialbinarypopD, j, j + l));
        //cout << DecisionVariableD[i] << "\n";
        DecisionVariableE[i] = lowerboundE + (((upperboundE - lowerboundE) / (pow(2, l) - 1)) * DV(InitialbinarypopE, j, j + l));
        //cout << DecisionVariableE[i] << "\n";
        DecisionVariableF[i] = lowerboundF + (((upperboundF - lowerboundF) / (pow(2, l) - 1)) * DV(InitialbinarypopF, j, j + l));
        //cout << DecisionVariableF[i] << "\n";
        j += l;
    }*/
    cout << "Enter the number of generations: ";
    short int generations;
    cin >> generations;
    cout << "\nPrecision of B = " << precisionB;
    cout << "\nPrecision of C = " << precisionC;
    cout << "\nPrecision of D = " << precisionD;
    cout << "\nPrecision of E = " << precisionE;
    cout << "\nPrecision of F = " << precisionF;
    long double* error = function(phi, H, tau, M, datapoints, DecisionVariableB, DecisionVariableC, DecisionVariableD, DecisionVariableE, DecisionVariableF, popsize);
    while (generations)
    {
        TournamentSelectionOperation(popsize, error, DecisionVariableB, DecisionVariableC, DecisionVariableD, DecisionVariableE, DecisionVariableF, InitialbinarypopB, InitialbinarypopC, InitialbinarypopF, InitialbinarypopE, InitialbinarypopF, l);
        SinglePointCrossover(InitialbinarypopB, InitialbinarypopC, InitialbinarypopD, InitialbinarypopE, InitialbinarypopF, l, popsize);
        BitwiseMutation(InitialbinarypopB, InitialbinarypopC, InitialbinarypopD, InitialbinarypopE, InitialbinarypopF, popsize, l);
        VariableCreation(popsize, DecisionVariableB, DecisionVariableC, DecisionVariableD, DecisionVariableE, DecisionVariableF, InitialbinarypopB, InitialbinarypopC, InitialbinarypopD, InitialbinarypopE, InitialbinarypopF, l, lowerboundB, lowerboundC, lowerboundD, lowerboundE, lowerboundF, upperboundB, upperboundC, upperboundD, upperboundE, upperboundF);
        error = function(phi, H, tau, M, datapoints, DecisionVariableB, DecisionVariableC, DecisionVariableD, DecisionVariableE, DecisionVariableF, popsize);
        generations--;
    }
    /*for (int i = 0; i < popsize; i++)
    {
        cout << DecisionVariableB[i] << "\n";
        cout << DecisionVariableC[i] << "\n";
        cout << DecisionVariableD[i] << "\n";
        cout << DecisionVariableE[i] << "\n";
        cout << DecisionVariableF[i] << "\n";
    }*/
    cout << "\n\n";
    /*for (int i = 0; i < popsize; i++)
    {
        cout << error[i] << "\n";
    }*/
    long double minimum = error[0];
    short int minval = 0;
    for (int i = 0; i < popsize; i++)
    {
        if (error[i] < minimum)
        {
            minimum = error[i];
            minval = i;
        }
    }
    cout << "The minimum error is: " << error[minval];
    cout << "\n The Variables are: \n";
    cout << "\n B = " << DecisionVariableB[minval];
    cout << "\n C = " << DecisionVariableC[minval];
    cout << "\n D = " << DecisionVariableD[minval];
    cout << "\n E = " << DecisionVariableE[minval];
    cout << "\n F = " << DecisionVariableF[minval];
    long double taupred, er = 0, Aval;
    for (int j = 0; j < datapoints; j++)
    {
        Aval = (H[j]) / (a1 * pow((phi[j] * M[j]), a2));
        taupred = ((pow(Aval, DecisionVariableB[minval])) / ((pow(Aval, DecisionVariableB[minval])) + 1)) * ((DecisionVariableC[minval]) * (pow(phi[j], DecisionVariableD[minval])) * mu0 * (pow(M[j], DecisionVariableE[minval]) * (pow(H[j], DecisionVariableF[minval]))));
        //taupredicted = (A / (A + 1)) * 2.5 * phi[j] * mu0 * (pow(M[j], 0.5)) * (pow(H[j], 1.5));
        //cout << "\nA = " << A;
        cout << " \n  %error = " << abs(taupred - tau[j]) / tau[j];
        er += abs(tau[j] - taupred);
    }
    cout << "\n" << er;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
