#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <algorithm>

////////////////////
#define n 10000			 	///Ввод размерности решетки
float Polevar = 0;           // Поле
int  prohod_MC = 1000;		/// Ввод числа шагов Монте-Карло
///////////////////////

using namespace std;

///////////////////////////////////		Ввод переменных
int choice_linear_dimension;			// выбор линейного размера системы
double choice_step_temperature;			// выбор шага температуры
int num_podhod=10;			    		// количвество попыток найти систему связей с мин. кол-м ошибок.
double mintemp=0;						// выбор min температуры
double maxtemp=4;						// выбор max температуры


////////////////////////////////////
void Slspin (short *spin);
void SpinFerr(short *spin);
void Vivodspin(int spin);
int Slsvaz( short *svaz, int num_podhod);
void VivodSSS( short *spin,short *svaz);
void Energy(short *spin,short *svaz,float *E);
float EnergyIJ(int i, short *spin, short *svaz);
void VivodEnergy( float *E);
void CopyE(float *E,float *E1);
int Claster( float *E1);
int ClasterFerr( float *E1);
void ColorOut2( float *E, short *spin);
void ColorOut( float *E, short *spin);
float SumEnergy(float *E);
void ColorOutIJ( float *E, short *spin, int i,int j);
void ColorOutSpin( float *E, short *spin);
int MaxClass(int per, float *E1, int tmp, int *Ochered);
int MaxClassFerr(int per, float *E1, int tmp, int *Ochered);
int sosed_sl(int top);
int sosed_sp(int top);
double SROTKL(double OutputMass, double SRKVOTKL);
double MagnMet( short *spin);
void PhazDiagramm(int var04, int var02, int var00, int var20, int var40, float *E);
void MCL(double mintemp, double maxtemp, short *spin, short *svaz, float *E, float *E1, int NUM_step_mcl);
void MCProhod(double temp,short *spin,short *svaz,float *E);


//////////////////////////////////////////////////// Забивает матрицу спинов 1
void SpinFerr(short *spin)
{
    for(int i=0;i<n;i++)
    {
        spin[i]=1;
    }
}

///////////////////////////////////////////////////////Подсчитывает энегрию системы
void Energy(short *spin,short *svaz,float *E)
{
    E[0] = -svaz[0] * spin[0] * spin[1] - svaz[n-1] * spin[0] * spin[n-1] - Polevar;
    E[n-1] = -svaz[n-1] * spin[n-1] * spin[0] - svaz[n-2] * spin[n-1] * spin[n-2] - Polevar;

    for (int i = 1; i < n; i++)
    {
        E[i] = -svaz[i] * spin[i] * spin[i+1] - svaz[i-1] * spin[i] * spin[i-1]- Polevar;
    }

}
///////////////////////////////////////////////////////Копирует массив с энергией Е в Е1
void CopyE(float *E,float *E1)
{
    for (int i = 0; i < n; i++)
    {
        E1[i] = E[i];

    }
}
///////////////////////////////////////////////////////Копирует массив с энергией Е в Е1
void CopyE(short *E,float *E1)
{
    for (int i = 0; i < n; i++)
    {
        E1[i] = E[i];

    }
}
///////////////////////////////////////////////////////Копирует массив с энергией Е в Е1
void CopyE(float *E,short *E1)
{
    for (int i = 0; i < n; i++)
    {
        E1[i] = E[i];

    }
}
/////////////////////////////////////////////////////////Подсчитывает энергию одного спина
float EnergyIJ(int i,short *spin,short *svaz)
{
    int x = n;
    float Ener = 0;
    if (i == 0)
        Ener = -svaz[0] * spin[0] * spin[1] - svaz[n-1] * spin[0] * spin[n-1] - Polevar;
    if (i == n - 1)
        Ener = -svaz[n-1] * spin[n-1] * spin[0] - svaz[n-2] * spin[n-1] * spin[n-2] - Polevar;
    if (i != 0 && i != n - 1)
        Ener = -svaz[i] * spin[i] * spin[i+1] - svaz[i-1] * spin[i] * spin[i-1]- Polevar;
    return Ener;
}
///////////////////////////////////////////////////////////
void MCProhod(double temp,short *spin,short *svaz,float *E)
{

    float En1,En2;
    int step=0;
    double slch;
    double veroyatnost=0;
    while(step<=(n))
    {
        int slspin=rand()%(n);
        int i = slspin;
        int var = spin[i];
        En1 = EnergyIJ(i, spin, svaz);
        spin[i] *= -1;
        En2 = EnergyIJ(i, spin, svaz);
        slch=((double)rand())/RAND_MAX;
        if (En2 <= En1)
            veroyatnost = 1;
        else
            veroyatnost = exp(-(En2 - En1) / temp);
        if (slch < veroyatnost)
            E[i] = En2;
        else
            spin[i] *= -1;
        step++;
    }
    Energy(spin,svaz,E);
}
///////////////////////////////////////////////////////////Подсчет намагниченности
double MagnMet(short *spin)
{
    double Magnvar=0;
    for (int i = 0; i < n; i++)
    {
        Magnvar += spin[i];

    }
    return abs(Magnvar/(n));
}
////////////////////////////////////////////////Для подсчета макс.кластера для спинового стекла
int Claster(float *E1)
{
    int tmp = 0;
    int t;
    int Max = 0;
    int *Ochered;
    Ochered = new int[n];
    for (int i = 0; i < n; i++)
    {
        if (E1[i] == -2 || E1[i] == 0)
            {
                tmp = 1;
                int top = i;
                t = MaxClass(top, E1, tmp,Ochered);
                if (t >= Max)
                    Max = t;
                t = 0;
            }

    }
    delete [] Ochered;
    Ochered = NULL;
    return Max;
}
//////////////////////////////////////////////////////Для подсчета макс.кластера для ферромагеника
int ClasterFerr(float *E1)
{
    int tmp = 0;
    int t;
    int Max = 0;
    int *Ochered;
    Ochered = new int[n];
    for (int i = 0; i < n; i++)
    {
        if (E1[i] == -2 )
            {
                tmp = 1;
                int top = i;
                t = MaxClassFerr(top, E1, tmp, Ochered);
                if (t >= Max)
                    Max = t;
                t = 0;
            }

    }
    delete [] Ochered;
    Ochered = NULL;
    return Max;
}
//////////////////////////////////////////////////////Для подсчета макс.кластера для спинового стекла
int MaxClass(int per,float *E1,int tmp,int *Ochered)
{
    int sl, sp, top;
    int w = 0;
    int r = 1;
    Ochered[0] = per;
    E1[per] = 0;
    while (w < r)
    {
        top = Ochered[w];

        ////
        sl = sosed_sl(top);
        if (E1[sl] == -2 || E1[sl] == 0)
        {
            Ochered[r] = sl;
            E1[sl] = 3;
            r++;
        }
        ////
        sp = sosed_sp(top);
        if (E1[sp] == -2 || E1[sp] == 0)
        {
            Ochered[r] = sp;
            E1[sp] = 3;
            r++;
        }
        ////
        w++;

    }
    return r;


}
/////////////////////////////////////////////////////Для подсчета макс.кластера для ферромагеника
int MaxClassFerr(int per,float *E1,int tmp, int *Ochered)
{
    int sl, sp, top;
    int w = 0;
    int r = 1;
    Ochered[0] = per;
    E1[per] = 0;
    while (w < r)
    {
        top = Ochered[w];

        ////
        sl = sosed_sl(top);
        if (E1[sl] == -2 )
        {
            Ochered[r] = sl;
            E1[sl] = 3;
            r++;
        }
        ////
        sp = sosed_sp(top);
        if (E1[sp] == -2 )
        {
            Ochered[r] = sp;
            E1[sp] = 3;
            r++;
        }
        ////
        w++;

    }
    return r;
}
//////////////////////////////////////////////////// Служебные функции для посчета мах.кластера

int sosed_sl( int top)
{
    int i = top;
    int sl = 0;
    if (i == 0)
        sl = n-1;
    else
        sl = i - 1;
    return sl;
}
int sosed_sp( int top)
{
    int i = top;
    int sp = 0;
    if (i == n - 1)
        sp = 0;
    else
        sp = i + 1;
    return sp;
}
//////////////////////////////////////////////////////////Суммирует энергию системы
float SumEnergy( float *E)
{
    float Sum = 0;
    for (int i = 0; i < n; i++)
    {
        Sum += E[i];
    }
    return Sum;
}
/////////////////////////////////////////////////////////// Получает поле
void POLE(short *spin,short *svaz,float *E1)
{
    int x = n;
    E1[0] = -svaz[0] * spin[1] - svaz[n-1] * spin[n-1] - Polevar;
    E1[n-1] = -svaz[n-2] * spin[n-2] - svaz[n-1] * spin[0] - Polevar;
    for (int i = 1; i < n - 1; i++)
    {
        E1[i] = -svaz[i] * spin[i+1] - svaz[i-1] * spin[i-1] - Polevar;
    }

}
///////////////////////////////////////////////////////////////
void PhazDiagramm(int *var02,int *var00,int *var20, float *E)
{
    *var02=0;
    *var00=0;
    *var20=0;
    for (int i = 0; i < n; i++)
    {
        if(E[i]==-2)
        {
            *var02+=1;
        } else if (E[i]==0)
        {
            *var00+=1;
        } else if (E[i]==2)
        {
            *var02+=1;
        }

    }

}

static void SROTKL(int size,double Outputmass[],double *SRZnach,double *SRKVOtkl)
{
    *SRZnach = 0;
    *SRKVOtkl = 0;
    for (int i = 0; i < size; i++)
    {
        *SRZnach += Outputmass[i];
    }
    *SRZnach = *SRZnach / size;
    for (int i = 0; i < size; i++)
    {
        Outputmass[i] = pow((Outputmass[i] - *SRZnach), 2);
        *SRKVOtkl += Outputmass[i];
    }
    *SRKVOtkl = sqrt(*SRKVOtkl / size);

}

static void SROTKL(int size,int Outputmass[],double *SRZnach,double *SRKVOtkl)
{
    *SRZnach = 0;
    *SRKVOtkl = 0;
    for (int i = 0; i < size; i++)
    {
        *SRZnach += Outputmass[i];
    }
    *SRZnach = *SRZnach / size;
    for (int i = 0; i < size; i++)
    {
        Outputmass[i] = pow((Outputmass[i] - *SRZnach), 2);
        *SRKVOtkl += Outputmass[i];
    }
    *SRKVOtkl = sqrt(*SRKVOtkl / size);

}

int main(int argc, char **argv)
{
    int rank;
    int size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    /////////для сравнения энергий////////////
    struct
    {
        double znach;
        int   rank;
    } SravnenieEnergiy, SravnenieEnergiyOUT;
    SravnenieEnergiy.rank=rank;
    /////////////////////

    ///////// Создание массива spin[n]////////////
    short *spin=new short[n];
    /////////////////////////////////////

    //Создание массива svaz[n]////////////
    short *svaz = new short[n];
    /////////////////////////////////////////

    ////// Создание массива E[n]////////////

    float *E=new float[n];
    /////////////////////////////////////

    ////// Создание массива E1[n]////////////
    float *E1=new float[n];
    ////////////////////////////////////

    ////////// Задание связей ферромагнетика///////////
    for (int i=0;i<n;i++)
        {
            svaz[i]=1;

        }
    /////////////массивы для вывода//////////////////////////
    double OUTEnergy[500];
    double OUTMagn[500];
    double OUTPP[500];
    double OUTPPF[500];

    int OUTPD02[500];
    int OUTPD00[500];
    int OUTPD20[500];

    int OUTPOLE02[500];
    int OUTPOLE00[500];
    int OUTPOLE20[500];

    double SRZnach=0;
    double SRKVOtkl=0;

    ofstream outE("energy.dat",ios::app);
    ofstream outM("Magn.dat",ios::app);
    ofstream outPP("PP.dat",ios::app);
    ofstream outPPF("PPF.dat",ios::app);
    ofstream outC("C.dat",ios::app);
    ofstream outae("ae.dat",ios::app);

    ofstream outPD02("PD02.dat",ios::app);
    ofstream outPD00("PD00.dat",ios::app);
    ofstream outPD20("PD2.dat",ios::app);

    ofstream outPOLE02("POLE02.dat",ios::app);
    ofstream outPOLE00("POLE00.dat",ios::app);
    ofstream outPOLE20("POLE2.dat",ios::app);
    ////////////////////////////
    SpinFerr(spin);
    Energy(spin,svaz,E);
    CopyE(E,E1);
    /////////////////// время///////////////
    srand((unsigned)time(NULL)+rank);
    //////MC///////////////////////////////


    double SumenergyVar,MagnVar,Teploemkost,Vospriimchivost;
    double predsumenergy=0;
    double predMagn=0;
    double maxSF,PPF,maxS,PP;
    int var02=0;
    int var00=0;
    int var20=0;

    int varPOLE02=0;
    int varPOLE00=0;
    int varPOLE20=0;

    double temp = mintemp;
    float Dt = 0.1;
    int Prohod=0;
    int cicle=1;
    while(temp<maxtemp)
    {

        int Prohod=0;
        while(Prohod<prohod_MC)
        {
            MCProhod(temp,spin,svaz,E);
            Prohod++;
        }
        CopyE(E,E1);
        maxSF = ClasterFerr(E1);
        PPF = maxSF / (n);                /////////////
        CopyE(E,E1);
        maxS = Claster(E1);
        PP= maxS / (n);
        SumenergyVar = SumEnergy(E);     /////////////
        MagnVar=MagnMet(spin);			/////////////
        PhazDiagramm(&var02,&var00,&var20,E);
        CopyE(E,E1);
        POLE(spin,svaz,E1);
        PhazDiagramm(&varPOLE02,&varPOLE00,&varPOLE20,E1);

        //////////////////////////////////////////////////

        cout<<"**********************************"<<'\n';
        cout<<"Temperature: "<<temp<<"          "<<'\n';
        cout<<"The number of atoms with the lowest energy (-4) in the largest cluster: "<<maxS<<"          "<<'\n';
        cout<<"The number of atoms with the lowest energy (-4 and -2) in the largest cluster: "<<maxSF<<"          "<<'\n';
        cout<<"The total energy of the system: "<<SumenergyVar<<"          "<<'\n';
        cout<<"The total order parameter: "<<PP<<"          "<<'\n';
        cout<<"The order parameter in Ferromagnet: "<<PPF<<"          "<<'\n';
        cout<<"Magnetization"<<MagnVar<<"           "<<'\n';
        cout<<"**********************************"<<'\n';

        ///////////////////////////////////////////////////
        SravnenieEnergiy.znach=SumenergyVar;

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&SravnenieEnergiy,&SravnenieEnergiyOUT,1,MPI_DOUBLE_INT, MPI_MAXLOC,MPI_COMM_WORLD);

        cout<<SravnenieEnergiyOUT.znach<<"         of"<<SravnenieEnergiyOUT.rank<<"       from"<<rank<<'\n';

        MPI_Bcast(spin,n,MPI_SHORT,SravnenieEnergiyOUT.rank,MPI_COMM_WORLD);
        CopyE(spin,E1);

        /////////////////параметры для циклов////////////
        cicle--;

            if(cicle==0)
            {
                /////////////////// вывод //////////////////////////////
                MPI_Gather(&SumenergyVar,1,MPI_DOUBLE,OUTEnergy,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
                MPI_Gather(&PP,1,MPI_DOUBLE,OUTPP,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
                MPI_Gather(&PPF,1,MPI_DOUBLE,OUTPPF,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
                MPI_Gather(&MagnVar,1,MPI_DOUBLE,OUTMagn,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

                MPI_Gather(&var02,1,MPI_INT,OUTPD02,1,MPI_INT,0,MPI_COMM_WORLD);
                MPI_Gather(&var00,1,MPI_INT,OUTPD00,1,MPI_INT,0,MPI_COMM_WORLD);
                MPI_Gather(&var20,1,MPI_INT,OUTPD20,1,MPI_INT,0,MPI_COMM_WORLD);


                MPI_Gather(&varPOLE02,1,MPI_INT,OUTPOLE02,1,MPI_INT,0,MPI_COMM_WORLD);
                MPI_Gather(&varPOLE00,1,MPI_INT,OUTPOLE00,1,MPI_INT,0,MPI_COMM_WORLD);
                MPI_Gather(&varPOLE20,1,MPI_INT,OUTPOLE20,1,MPI_INT,0,MPI_COMM_WORLD);



                if (rank==0)
                {
                    SRZnach=0;
                    SRKVOtkl=0;
                    SROTKL(size,OUTEnergy,&SRZnach,&SRKVOtkl);
                    outE<<temp<<'\t'<<SRZnach/(n)<<'\t'<<SRKVOtkl/(n)<<endl;

                    Teploemkost=fabs(SRZnach/(n)-predsumenergy/(n))/Dt;
                    predsumenergy=SRZnach;
                    outC<<temp<<'\t'<<Teploemkost<<endl;

                    SRZnach=0;
                    SRKVOtkl=0;
                    SROTKL(size,OUTMagn,&SRZnach,&SRKVOtkl);
                    outM<<temp<<'\t'<<SRZnach<<'\t'<<SRKVOtkl<<endl;
                    predMagn=SRZnach;

                    SRZnach=0;
                    SRKVOtkl=0;
                    SROTKL(size,OUTPP,&SRZnach,&SRKVOtkl);
                    outPP<<temp<<'\t'<<SRZnach<<'\t'<<SRKVOtkl<<endl;

                    SRZnach=0;
                    SRKVOtkl=0;
                    SROTKL(size,OUTPPF,&SRZnach,&SRKVOtkl);
                    outPPF<<temp<<'\t'<<SRZnach<<'\t'<<SRKVOtkl<<endl;

                    SRZnach=0;
                    SRKVOtkl=0;
                    SROTKL(size,OUTPD02,&SRZnach,&SRKVOtkl);
                    outPD02<<temp<<'\t'<<SRZnach/(n)<<'\t'<<SRKVOtkl/(n)<<endl;

                    SRZnach=0;
                    SRKVOtkl=0;
                    SROTKL(size,OUTPD00,&SRZnach,&SRKVOtkl);
                    outPD00<<temp<<'\t'<<SRZnach/(n)<<'\t'<<SRKVOtkl/(n)<<endl;

                    SRZnach=0;
                    SRKVOtkl=0;
                    SROTKL(size,OUTPD20,&SRZnach,&SRKVOtkl);
                    outPD20<<temp<<'\t'<<SRZnach/(n)<<'\t'<<SRKVOtkl/(n)<<endl;


                    SRZnach=0;
                    SRKVOtkl=0;
                    SROTKL(size,OUTPOLE02,&SRZnach,&SRKVOtkl);
                    outPOLE02<<temp<<'\t'<<SRZnach/(n)<<'\t'<<SRKVOtkl/(n)<<endl;

                    SRZnach=0;
                    SRKVOtkl=0;
                    SROTKL(size,OUTPOLE00,&SRZnach,&SRKVOtkl);
                    outPOLE00<<temp<<'\t'<<SRZnach/(n)<<'\t'<<SRKVOtkl/(n)<<endl;

                    SRZnach=0;
                    SRKVOtkl=0;
                    SROTKL(size,OUTPOLE20,&SRZnach,&SRKVOtkl);
                    outPOLE20<<temp<<'\t'<<SRZnach/(n)<<'\t'<<SRKVOtkl/(n)<<endl;



                }

                Polevar=0.05;
                int secondProhod=0;
                while(secondProhod<prohod_MC/2)
                {
                    MCProhod(temp,spin,svaz,E);
                    Prohod++;
                    secondProhod++;
                }
                MagnVar=MagnMet(spin);			/////////////

                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Gather(&MagnVar,1,MPI_DOUBLE,OUTMagn,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

                if(rank==0)
                {
                    SRZnach=0;
                    SRKVOtkl=0;
                    SROTKL(size,OUTMagn,&SRZnach,&SRKVOtkl);
                    Vospriimchivost=fabs((SRZnach-predMagn)/Polevar);

                    outae<<temp<<'\t'<<Vospriimchivost<<endl;
                }


                CopyE(E1,spin);
                Polevar=0;

                /*
                if(temp>=2.2&&temp<2.3)
                    Dt=0.02;
                else
                    Dt=0.1;
                */////////
                temp+=Dt;
                ////////
                temp=(long)(100*(temp+0.005));
                temp=temp/100;
                ///////

                if (temp==0.4||temp==0.5||temp==0.6)
                {
                    cicle=3;
                }

                else
                        cicle=1;

            }
        //////////////////////////////
    }


    ///////////

    MPI_Finalize();

    return 0;
}


