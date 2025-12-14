/*   ___           _______. _______  __  .__   __. 
    /   \         /       ||   ____||  | |  \ |  | 
   /  ^  \       |   (----`|  |__   |  | |   \|  | 
  /  /_\  \       \   \    |   __|  |  | |  . `  | 
 /  _____  \  .----)   |   |  |     |  | |  |\   | 
/__/     \__\ |_______/    |__|     |__| |__| \__| @ LNS-CT

*/

#include <iostream>
#include <cmath>
#include <ctime>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TROOT.h"

using namespace std;

double ar(double y, double x)
{
    double angolo = 0;
    double rad = 0.0174532952;
    if (x > 0)
        angolo = atan(y / x);
    if (x < 0 && y > 0)
        angolo = atan(y / x) + 180 * rad;
    if (x < 0 && y < 0)
        angolo = atan(y / x) - 180 * rad;
    if (x > 0 && y == 0)
        angolo = 0;
    if (x == 0 && y > 0)
        angolo = 90 * rad;
    if (x < 0 && y == 0)
        angolo = 180 * rad;
    if (x == 0 && y < 0)
        angolo = -90 * rad;
    return (angolo);
}

// TGXCM VIENE ESTRATTO CON UNA DISTRIBUZIONE < DISTX >
double DISTX(double X)
{
    double distx = 1. - X + X;
    return (distx);
}

// PSIG1 VIENE ESTRATTO CON UNA DISTRIBUZIONE < DISTP >
double DISTP(double X)
{
    double distp = 1. - X + X;
    return (distp);
}

// Modulazione della distribuzione di p3
//      < DISTP3 >
double hulten(double X)
{
    // terza particella = neutrone
    double neutron = 1.007825;
    // terza particella = protone
    double proton = 1.0072056;
    double xm3 = proton * 931.5;
    double gpsn = sqrt(2.25 * 59.8) * pow(sqrt(2.25) + sqrt(59.8), 3);
    double gpsd = pow((2.25 + 2. * X * X / (2. * xm3)), 2) * pow((59.8 + 2. * X * X / (2. * xm3)), 2);
    double gps = gpsn / gpsd;
    double distp3 = 2 * gps;
    return distp3;
}

double n14(double X)
{
    double distp3 = exp(-((pow(X,2))/(2*pow(52.4,2))));
    return distp3;
} 

void sim()
{

    string intro;
    string foutname;
    int NMAX;
    vector<double> ER12;
    ER12.reserve(300000);
    int sel = -1;
    int p3func = -1;

    gRandom = new TRandom3();

    bool selvalid = true;
    bool p3fvalid = true;
    printf("\n\n******** SIMULAZIONE 3 CORPI ********\n");

    while (true)
    {
        if (!selvalid)
            printf("\n Scelta errata, seleziona un numero dalla lista. \n\n");

        printf(
            "\n----Opzioni----\n"
            "0) solo spazio delle fasi con tagli dei rivelatori\n"
            "1) spazio delle fasi per distr. Ps\n"
            "\n  Scelta: ");
        cin >> sel; //da migliorare la gestione
        if (sel >= 0 && sel <= 1)
            break;
        else
            selvalid = false;
    }
    if (sel == 1)
    {
        while (true)
        {
            if (!p3fvalid)
                printf("\n Scelta errata, seleziona un numero dalla lista. \n\n");

            printf(
                "\n\n---- Possibili funzioni per la dist. Ps ----\n"
                "0) Hulten\n"
                "1) 14N Gaussian (NITROX)"
                "\n  Scelta: ");
            cin >> p3func; //da migliorare la gestione
            if (p3func <=1 )
                break;
            else
                p3fvalid = false;
        }
    }

    cout << "Nome file input: ";
    cin >> intro;
    ifstream in(intro.c_str());
    while (!in.is_open())
    {
        cout << "Il file " << intro << " non esiste nella directory, riprova." << endl;
        cout << "Nome file input: ";
        cin >> intro;
        ifstream in(intro.c_str());
    }

    cout << "Nome file output: ";
    cin >> foutname;
    cout << "Numero estrazioni: ";
    cin >> NMAX;

    //FILE *fprova = fopen("report.txt","w+");

    TFile *fout = new TFile(foutname.c_str(), "recreate"); //verificare se il file esiste
    TTree *tree = new TTree("sim", " ");

    double thx, e1lab, e2lab, tet1lab, tet2lab, e12, e13, e23, p3, psig1, thcm, p3f, e3lab, tet3lab;
    tree->Branch("thx", &thx, "thx/D");
    tree->Branch("e1lab", &e1lab, "e1lab/D");
    tree->Branch("e2lab", &e2lab, "e2lab/D");
    tree->Branch("tet1lab", &tet1lab, "tet1lab/D");
    tree->Branch("tet2lab", &tet2lab, "tet2lab/D");
    tree->Branch("e12", &e12, "e12/D");
    tree->Branch("e13", &e13, "e13/D");
    tree->Branch("e23", &e23, "e23/D");
    tree->Branch("p3", &p3, "p3/D");
    tree->Branch("psig1", &psig1, "psig1/D");
    tree->Branch("thcm", &thcm, "thcm/D");
    tree->Branch("p3f", &p3f, "p3f/D");
    tree->Branch("e3lab", &e3lab, "e3lab/D");
    tree->Branch("tet3lab", &tet3lab, "tet3lab/D");

    double PIG = 3.14159265;
    double RAD = PIG / 180.;
    double deg = 1 / RAD;
    double MC2 = 931.5;

    string line;
    vector<double> a;
    int ND, counter = 0;
    double temp;

    while (getline(in, line))
    {
        if (line.empty()) continue;
        if (line[0] == '*') continue;
        istringstream iss(line);
        string dummy;
        if (counter == 9)
        {
            iss >> ND;
        }
        else
        {
            iss >> temp;
            a.push_back(temp);
        }

        counter++;
    }
    in.close();

    for (unsigned int jj = 0; jj < a.size(); jj++)
    {
        cerr << "a[" << jj << "]: " << a[jj] << endl;
    }

    double AMP = a[0];
    double XMP = AMP * MC2;
    double AMT = a[1]; //double XMT=AMT*MC2;
    double AM1 = a[2];
    double XM1 = AM1 * MC2;
    double AM2 = a[3];
    double XM2 = AM2 * MC2;
    double AM3 = a[4];
    double XM3 = AM3 * MC2;
    double AMTR = a[5];
    double XMTR = AMTR * MC2;
    double AMX = AM1 + AM2;
    double XMX = AMX * MC2;

    // EP     = ENERGIA DEL BEAM
    // Q0     = Q-VALUE DEL PROCESSO A DUE CORPI (I STADIO)
    // ETHRES = ENERGIA DEL DECADIMENTO  X --> 1 + 2
    double EP = a[6];
    double Q0 = a[7];
    double ETHRES = a[8];

    // ENERGIE DI SOGLIA DEI RIVELATORI
    double E1TH = a[9];
    double E2TH = a[10];
    // LIMITI DEI RIVELATORI
    vector<double> TETD, TETU;
    double td, tu;
    for (int ind = 0; ind < ND; ind++)
    {
        td = a[2 * ind + 11] * RAD;
        tu = a[2 * ind + 12] * RAD;
        if (td > tu)
        {
            cout << "error in angles !" << endl;
            return;
        }
        else
        {
            TETD.push_back(td);
            TETU.push_back(tu);
        }
    }

    // INTERVALLO DI ENERGIA RELATIVA
    double EREL1 = a[2 * ND + 11];
    double EREL2 = a[2 * ND + 12];
    double DEREL = a[2 * ND + 13];

    // LIMITI PER ESTRAZIONE TXCM (TXMIN E TXMAX)
    //			E PSI1 (0 E PSMAX)

    double TXMIN = a[2 * ND + 14];
    double TXMAX = a[2 * ND + 15];
    double PSMAX = a[2 * ND + 16];
    //double WEIGHT=pow(360,2)/(PSMAX*(TXMAX-TXMIN));

    // LIMITI PER ENERGIE RELATIVE E13 , E23

    double E13D = a[2 * ND + 17];
    double E13U = a[2 * ND + 18];
    double E23D = a[2 * ND + 19];
    double E23U = a[2 * ND + 20];
    int nb = 0;

    double VP = sqrt(2. * EP / XMP);
    double VCM = VP * AMP / (AMP + AMT);
    //double WPCM=VP*AMT/(AMP+AMT);
    double EREAZ = EP * AMT / (AMP + AMT);

    //double q3=Q0-ETHRES;

    double argo, xden, xdenb, xdena, xnum, xnumy, xnumx, vty, vtx, p3y, p3x, v2y, v2x, p2y, p2x, v1y, v1x, p1y, p1x, vpy, vpx, ppy, ppx;
    double p2, p1, pp, v2, v1, vp, stc, ctc, sta, cta, P3RAN, r3, XM23, XM13;
    double V23Y, V23X, V13Y, V13X, TETG3, TET3, V3LABY, V3LABX, V2LABY, V2LABX, V1LABY, V1LABX;
    double TET2, TET1, X02, Y02, X01, Y01, DETEST, ETEST, V2LAB, V1LAB, PSI2, PSI1, PSIG2, DISPRAN, PSIRAE, PSIG1, DISTP3;
    double r2, TETGX, TETX, CTETX, STETX, ENX, VNX, VNX2, WXCM, TXCM, DISTRAN, TGXCM, r1, V2CM, V1CM, QX, EREL, r, QQ;
    int IT, IP, j1, j2, jdev, NESTR;
    int I = 0;

    while (1)
    {
        if (!(fout->IsOpen()))
            break;
        I++;
        if (I > 1)
        {
            cout << "EN. REL. = " << ER12[I - 1] << " no. eventi buoni = " << nb << endl;
        }

        ER12[I] = EREL1 + (I - 1) * DEREL;
        QQ = Q0 - ETHRES - ER12[I];
        if (ER12[I] > EREL2)
        {
            tree->Write();
            tree->Print();
            fout->Write();
            fout->Close();
            return;
        }
        if ((EREAZ + QQ) <= 0.)
        {
            tree->Write();
            tree->Print();
            fout->Write();
            fout->Close();
            return;
        }

        for (NESTR = 1; NESTR <= NMAX; NESTR++)
        {

            r = gRandom->Rndm();
            EREL = EREL1 + (I - 1) * DEREL + DEREL * r;
            QX = Q0 - ETHRES - EREL;
            V1CM = sqrt(2. * EREL / (XM1 * (1. + XM1 / XM2)));
            V2CM = V1CM * XM1 / XM2;

            // ESTRAZIONE TGXCM TRA TXMIN E TXMAX
            // FUNCTION DISTX E' LA DISTRIBUZIONE DI TGXCM
            // CHE NE MODULA L'ESTRAZIONE
            r1 = gRandom->Rndm();
            TGXCM = r1 * (TXMAX - TXMIN) + TXMIN;
            DISTRAN = gRandom->Rndm();
            if (DISTX(TGXCM) < DISTRAN)
                continue;
            TXCM = TGXCM * RAD;

            // CALCOLO ENERGIA E VELOCITA' DEL NUCLEO INTERMEDIO ECCITATO
            //          (ENX)      (VNX)   NEL LAB. SYSTEM

            if ((EREAZ + QX) <= 0.)
            { //se l'e. totale è <0 esco dal programma
                tree->Write();
                tree->Print();
                fout->Write();
                fout->Close();
                return;
            }

            WXCM = sqrt(2. * (EREAZ + QX) / (XMX * (1. + XMX / XM3)));
            VNX2 = pow(VCM, 2) + pow(WXCM, 2) + 2. * VCM * WXCM * cos(TXCM);

            if (VNX2 < 0)
            {
                NESTR = NESTR - 1;
                continue;
            }
            VNX = sqrt(VNX2);
            ENX = .5 * XMX * VNX2;

            // CALCOLO ANGOLO DEL NUCLEO INTERM. ECCITATO
            // (TETX)  NEL LAB SYSTEM
            STETX = WXCM * sin(TXCM) / VNX;
            CTETX = (WXCM * cos(TXCM) + VCM) / VNX;
            TETX = ar(STETX, CTETX);
            TETGX = TETX / RAD;

            // ESTRAZIONE PSIG1 TRA 0 E PSMAX
            // FUNCTION DISTP E' LA DISTRIBUZIONE DI PSIG1
            //     CHE NE MODULA L'ESTRAZIONE
            r2 = gRandom->Rndm();
            PSIG1 = r2 * PSMAX;
            PSIRAE = 180. - (PSIG1 + TETGX);
            DISPRAN = gRandom->Rndm();
            if (DISTP(PSIG1) < DISPRAN)
                continue;
            PSIG2 = PSIG1 - 180.;
            PSI1 = PSIG1 * RAD;
            PSI2 = PSIG2 * RAD;
            V1LAB = sqrt(pow(V1CM, 2) + pow(VNX, 2) - 2 * V1CM * VNX * cos(PIG - PSI1));
            V2LAB = sqrt(pow(V2CM, 2) + pow(VNX, 2) - 2 * V2CM * VNX * cos(PIG + PSI2));
            e1lab = .5 * XM1 * pow(V1LAB, 2);
            e2lab = .5 * XM2 * pow(V2LAB, 2);
            ETEST = e1lab + e2lab - EREL;
            DETEST = ETEST * .01;
            if (!(ENX > (ETEST - DETEST) && ENX < (ETEST + DETEST)))
            {
                cout << "ERROR IN ENERGY BALANCE !  NESTR, E12 = " << NESTR << ", " << EREL << endl;
                cout << e1lab << " " << e2lab << " " << ETEST << endl;
            }

            // CALCOLO DI ANGOLI DI 1 E 2 NEL LAB SYSTEM
            // (TET1 E TET2)
            if (e1lab < E1TH || e2lab < E2TH)
                continue; // CONTROLLO SE LE ENERGIE SONO SOPRA SOGLIA DEI RIVELATORI
            Y01 = V1CM * sin(PSI1 + TETX) + VNX * sin(TETX);
            X01 = V1CM * cos(PSI1 + TETX) + VNX * cos(TETX);
            Y02 = V2CM * sin(PSI2 + TETX) + VNX * sin(TETX);
            X02 = V2CM * cos(PSI2 + TETX) + VNX * cos(TETX);
            TET1 = ar(Y01, X01);
            TET2 = ar(Y02, X02);

            // CONTROLLO SE 1 E 2 ENTRANO NEI RIVELATORI
            jdev = 0;
            for (j1 = 0; j1 < ND; j1++)
            {
                if (TET1 > TETD[j1] && TET1 < TETU[j1])
                    jdev = j1 + 1;
            }
            if (jdev == 0)
                continue;
            //fprintf(fprova,"\njdev: %i\n\n",jdev);
            for (j2 = 0; j2 < ND; j2++)
            {
                //fprintf(fprova,"j2: %i\n",j2);
                if ((j2 + 1) == jdev)
                    continue;
                if (TET2 > TETD[j2] && TET2 < TETU[j2])
                {

                    // CALCOLO E13 , E23  ***
                    V1LABX = V1LAB * cos(TET1);
                    V1LABY = V1LAB * sin(TET1);
                    V2LABX = V2LAB * cos(TET2);
                    V2LABY = V2LAB * sin(TET2);
                    V3LABX = (XMP * VP - XM1 * V1LABX - XM2 * V2LABX) / XM3;
                    V3LABY = -(XM1 * V1LABY + XM2 * V2LABY) / XM3;
                    TET3 = ar(V3LABY, V3LABX);
                    TETG3 = TET3 / RAD;
                    p3 = XM3 * sqrt(pow(V3LABX, 2) + pow(V3LABY, 2));
                    e3lab = 0.5 * XM3 * (pow(V3LABX, 2) + pow(V3LABY, 2));
                    p3f = XM3 * sqrt(pow((V3LABX - VP), 2) + pow(V3LABY, 2));

                    if (TET3 >= 0.)
                    {
                        p3 = p3;
                        p3f = p3f;
                    }
                    if (TET3 < 0.)
                    {
                        p3 = -p3;
                        p3f = -p3f;
                    }

                    // Modulazione di ps con varie funzioni, da implementare eckart e la scelta della modulazione 
                    if (sel == 1)
                    {
                        r3 = gRandom->Rndm();
                        P3RAN = r3;
                        if (p3func == 0)
                            DISTP3 = hulten(p3);
                        else if (p3func == 1)
                            DISTP3 = n14(p3);
                        else
                            DISTP3 = 0;
                        if (DISTP3 < P3RAN)
                            break;
                    }

                    V13X = V1LABX - V3LABX;
                    V13Y = V1LABY - V3LABY;
                    V23X = V2LABX - V3LABX;
                    V23Y = V2LABY - V3LABY;
                    XM13 = XM1 * XM3 / (XM1 + XM3);
                    XM23 = XM2 * XM3 / (XM2 + XM3);
                    e13 = .5 * XM13 * (pow(V13X, 2) + pow(V13Y, 2));
                    e23 = .5 * XM23 * (pow(V23X, 2) + pow(V23Y, 2));

                    if (e13 < E13D || e13 > E13U)
                        break;
                    if (e23 < E23D || e23 > E23U)
                        break;

                    if (!(E13U == E13D && E23U == E23D))
                    {
                        cta = cos(TET1);
                        sta = sin(TET1);
                        ctc = cos(TET2);
                        stc = sin(TET2);
                        vp = sqrt(2. * EP / XMP);
                        v1 = V1LAB;
                        v2 = V2LAB;
                        pp = XMP * vp;
                        p1 = XM1 * v1;
                        p2 = XM2 * v2;
                        ppx = pp;
                        ppy = 0.;
                        vpx = vp;
                        vpy = 0.;
                        p1x = p1 * cta;
                        p1y = p1 * sta;
                        v1x = p1x / XM1;
                        v1y = p1y / XM1;
                        p2x = p2 * ctc;
                        p2y = p2 * stc;
                        v2x = p2x / XM2;
                        v2y = p2y / XM2;
                        p3x = ppx - p1x - p2x;
                        p3y = ppy - p1y - p2y;
                        vtx = -p3x / XMTR;
                        vty = -p3y / XMTR;
                        xnumx = (vpx - vtx) * (v1x - v2x);
                        xnumy = (vpy - vty) * (v1y - v2y);
                        xnum = xnumx + xnumy;
                        xdena = sqrt(pow(vpx - vtx, 2) + pow(vpy - vty, 2));
                        xdenb = sqrt(pow(v1x - v2x, 2) + pow(v1y - v2y, 2));
                        xden = xdena * xdenb;
                        argo = xnum / xden;
                        if (argo > -1. || argo < 1.)
                            thcm = acos(xnum / xden) * deg;
                    }

                    IP = (int)(PSIRAE / 2.) + 1;
                    IT = (int)((TGXCM + 180.) / 3.) + 1;
                    if (!(IP <= 0 || IP > 90 || IT <= 0 || IT > 120))
                    {
                        nb = nb + 1; //eventi buoni
                        thx = TGXCM;
                        tet1lab = TET1 / RAD;
                        tet2lab = TET2 / RAD;
                        e12 = ER12[I];
                        psig1 = PSIG1;
                        tet3lab = TETG3;
                        tree->Fill();
                    }
                    break;
                }
                break;
            }
        }
    }
    //fclose(fprova);
}

#ifndef __CLING__
int main()
{
    sim();
    return 0;
}
#endif
