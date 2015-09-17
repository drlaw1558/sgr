// SimToSM_feh.cxx (Version1, based on Version 6 of SimToSM)
// Convert from scf output files to a table of positional data of satellite debris
// Note that Galactic xyz are Galactocentric left-handed, vGSR is GSR velocity seen from solar position
// Code format originally similar to SimToSgr.cxx, but for a general simulation, not just Sgr
//
// Changes made in v2:
//  - Allows automatic calculation for sims in which satellite is completely unbound by end of sim.
//    NB: requires both R00XS002 and R00XS003 otherwise will get lots of errors.
//  - Eliminates custom plane fitting routines that have never been used.
//  - Allows interactive specification which kind of output file to generate
//
// Changes made in v3:
//  - Allows option to output in a format used by POVray for 3D animation
//  - Generalizing this output format non-trivial: current must configure in program and recompile
//  - Fixed problem that never seeded srand48 for distance scatter.
//
// Changes made in v4:
//  - Added swappable code section to chose peri-peri coloring or apo-apo
//  - Changed coding to indicate that final bound particles have Pcol=-1
//
// Changes made in v5:
//  - Added Sgr coordinate ability.
//
// Changes made in v6:
//  - Added computation of angular momentum, which can be used to discriminate
//    leading and trailing tails.
//
// Changes for feh version (10/21/08):
//  - Twiddled dsgr and rsun specification
//  - Now introduces metallicity distribution and work
//    out corresponding K magnitudes
//  - Keep in mind this will NOT work with the SMM multiple-type file outputs
//  - Eliminated a bunch of custom output formats that had accumulated and
//    were not necessary to keep around
//
// 2/6/09: Note that the scf code has been modified to specify leading/trailing tail
//  Format read from S00X files therefore updated
//
// 4/16/09: Modified from SimToSM_feh to ignore feh stuff temporarily
//
// 5/11/09: Modified to incorporate an auto-run feature

// David R. Law <drlaw@astro.caltech.edu>
// Last modified: 5/11/09

// Have currently made this pretty much like the first metal code I wrote.
// Now I need to change it to read in K mags from some file.

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include "/Users/drlaw/GenCode/Library/Coords.h"
#include "/Users/drlaw/GenCode/Library/SgrCoords.h"
#include "/Users/drlaw/GenCode/Library/Random.cxx"
#include "/Users/drlaw/GenCode/Library/Stats.cxx"

using namespace std;

// Global constants
int nbod=100000;
string PartFile,PartFile2,LogFile,OutpFile,BIpopfile;
int nentries,napo,SystNumber;
double mu,ru,pi,radpdeg,msun,cmpkpc,kmpkpc,secpgyr,arcpyr,G,tu,vu;
double m,xframe,yframe,zframe,vxframe,vyframe,vzframe;
double rsun,rparam,DSGR;
int OutpType=1;
bool destroy=false;

// Global vectors, initialize to length 1e5 to save enough space, will resize later
vector<double> x(nbod); vector<double> y(nbod); vector<double> z(nbod);
vector<double> vx(nbod); vector<double> vy(nbod); vector<double> vz(nbod);
vector<double> tub(nbod); vector<double> d(nbod);
vector<double> dGC(nbod); vector<double> b(nbod); vector<double> l(nbod);
vector<double> Xm(nbod); vector<double> Ym(nbod); vector<double> Zm(nbod);
vector<double> vGSR(nbod); vector<double> Pcolor(nbod);
vector<double> lambda(nbod); vector<double> beta(nbod);
vector<double> Kabs(nbod); vector<double> Kapp(nbod); vector<double> Xstat(nbod);
vector<int> Lmflag(nbod);
vector<double> FeH(nbod); vector<double> RA(nbod); vector<double> Dec(nbod);
vector<int> Popn(nbod);
vector<double> mul(nbod); vector<double> mub(nbod);
// Time stuff
vector<double> t(0); vector<double> tapo(0);
// Function prototypes
void GetInput(); void SetParam(); void ReadFrame(); void ReadPart(); void TransformStandard();
void ReadBIfile();
void SetColors(); void Output(); void CalcNentries(); void CalcNbods();
void OutputFehTrend();
void ReadPartDestroyed();int Lcount(string &File);

int main(int argc, char*argv[])
  {
  if (argc>2)// If enough parameters to automate, do so
    {
    string name=argv[1];
    mu=atof(argv[2]);
    ru=0.;
    SystNumber=1;
    DSGR=28.;
    OutpType=1;

    PartFile = "./";
    PartFile.append(name.c_str());
    PartFile2=PartFile;
    LogFile = PartFile;
    OutpFile = PartFile;
    PartFile.append("S002");
    PartFile2.append("S003");
    LogFile.append("LOG");
    cout << LogFile << " " << OutpFile << " " << PartFile << endl;
    
    ifstream fintest(PartFile2.c_str()); // Determine whether file R00XS003 exists
    if (fintest) destroy=true;          // if it does, satellite was destroyed, calculate accordingly
    fintest.close();
    }
  else
    GetInput();// Otherwise do it interactively

  CalcNentries();
  CalcNbods();
  SetParam();
  ReadFrame();
  if (destroy==false) ReadPart();
  else ReadPartDestroyed();
  SetColors();
  //  ReadBIfile();
  TransformStandard();
  Output();
  // OutputFehTrend();
  }

void ReadBIfile()
  {
  ifstream fin(BIpopfile.c_str());
  string Junk;
  double temp,temp3,temp3b,temp4;
  string temp2;
  getline(fin,Junk,'\n');
  for (int i=0; i <nbod; i++)
    {
    fin >> temp3 >> temp3b >> Junk >> Junk >> temp2 >> Junk >> temp4 >> Junk >> temp;
    Kabs[i]=temp;
    FeH[i]=temp4;
    Xstat[i]=temp3b;
    if (temp2 == "A") Popn[i]=0;// Might be wacky
    if (temp2 == "B") Popn[i]=1;
    if (temp2 == "C1") Popn[i]=2;
    if (temp2 == "C2") Popn[i]=3;
    if (temp2 == "C3") Popn[i]=4;
    if (temp2 == "D") Popn[i]=5;
    if (temp2 == "E") Popn[i]=6;
    }
  }

void GetInput() // 
  {
  string name;
  cout << "Data file location (e.g. 'R001'): ";
  cin >> name;
  //  cout << "BI pop'n file: ";
  //  cin >> BIpopfile;
  cout << "Satellite mass (e.g. 1e9): ";
  cin >> mu;
  cout << "Scale parameter (0.0 for default): ";
  cin >> ru;
  cout << "Which Sgr system (1-4, 1 if unsure)?: ";
  cin >> SystNumber;
  cout << "Distance to Sgr = ";
  cin >> DSGR;

  cout << "What kind of output file?" << endl;
  cout << "(1) Revised standard output (lambda, beta, RA,Dec, x, y, z, dist, K, vGSR, Lmflag, Pcol)" << endl;
  cout << "(2) Custom website output (lambda, beta, l, b, ra,dec, xgc, ygc, zgc, x, y, z, u, v, w, dist, K, vgsr, mul, mub, mua, mud, Lmflag, Pcol)" << endl;
  cout << "(3) POVRAY output." << endl;
  cout << "Type: ";
  cin >> OutpType;

  PartFile = "./";
  PartFile.append(name.c_str());
  PartFile2=PartFile;
  LogFile = PartFile;
  OutpFile = PartFile;
  PartFile.append("S002");
  PartFile2.append("S003");
  LogFile.append("LOG");

  ifstream fintest(PartFile2.c_str()); // Determine whether file R00XS003 exists
  if (fintest) destroy=true;          // if it does, satellite was destroyed, calculate accordingly
  fintest.close();
  }

void CalcNentries() // Calculate how many entries in LogFile before sat destroyed so know where frame defined
  {
  if (destroy==false) // Line counting will be messed up if LogFile weird due to destruction
    {
    int lines=Lcount(LogFile);
    nentries=(lines-3)/9-1;
    }
  else
    {
    int lines=0;
    double Junk; double MassRem=1.;
    ifstream fin(LogFile.c_str());
    for(int i=0; i<18; i++) fin >> Junk;
    lines+=9;
    do
      {
      fin >> Junk >> Junk >> Junk;
      fin >> MassRem;
      for(int i=0; i<14; i++) fin >> Junk;
      lines+=9;
      }
    while (MassRem>0);
    nentries=lines/9-2;
    }
  }

void CalcNbods()
  {
  int lines=Lcount(PartFile);

  nbod=lines-1;
  x.resize(nbod); y.resize(nbod); z.resize(nbod);
  vx.resize(nbod); vy.resize(nbod); vz.resize(nbod);
  tub.resize(nbod); d.resize(nbod); dGC.resize(nbod);
  b.resize(nbod); l.resize(nbod);
  Xm.resize(nbod); Ym.resize(nbod); Zm.resize(nbod);
  lambda.resize(nbod); beta.resize(nbod);
  vGSR.resize(nbod); Pcolor.resize(nbod);
  Kabs.resize(nbod); Kapp.resize(nbod); Lmflag.resize(nbod);
  }

void SetParam() // Define constants and vector sizes
  {
  if (ru==0.0)
    ru=0.9*pow(1e9/mu,(-1./3.));

  pi=3.141592653589793;
  radpdeg=pi/180.;
  rsun=8.0;
  msun=1.989e33;
  cmpkpc=3.085678e21;
  kmpkpc=3.085678e16;
  secpgyr=60*60*24*365*1e9;
  // arcpyr converts from radians/sec to arcsec/year
  arcpyr=(3600/radpdeg)*(60*60*24*365);
  G=6.67e-8;
  G=G*msun/cmpkpc/cmpkpc/cmpkpc;

  tu=ru*sqrt(ru/(mu*G));
  vu=(cmpkpc*ru*1e-5)/tu;
  tu=tu/secpgyr;
  cout << "tu = " << tu << endl;

  t.resize(nentries);
  for (int i=0;i<nentries;i++) t[i]=i;
  }

void ReadFrame()
  {
  napo=0;
  double r0=1e5; double r1=1e5; double Junk; double r2;
  double xtemp,ytemp,ztemp,vxtemp,vytemp,vztemp;
  tapo.push_back(0);
  ifstream fin(LogFile.c_str());
  for(int i=0; i<18; i++) fin >> Junk;
  for (int i=1; i < nentries; i++)
    {
    fin >> t[i-1]; fin >> Junk;
    fin >> Junk; fin >> m;
    fin >> xtemp >> ytemp >> ztemp;
    fin >> vxtemp >> vytemp >> vztemp;
    for (int j=0; j<8; j++) fin >> Junk;
    if (i>10)
      {
      r2=sqrt(xtemp*xtemp+ytemp*ytemp+ztemp*ztemp);
      // Uncomment the next line if you want apo-apo color coding
      if((r2<r1) && (r0<r1)) // (Be sure to comment out this line if not using!)
      // Uncomment the next line if you want peri-peri color coding
      //if((r2>r1) && (r0>r1)) // (Be sure to comment out this line if not using!)
        {
        napo=napo+1;
        tapo.push_back(t[i-2]);
        }
      r0=r1;
      r1=r2;
      }
    }
  fin >> t[nentries-1]; fin >> Junk;
  fin >> Junk; fin >> m;
  fin >> xframe >> yframe >> zframe;
  fin >> vxframe >> vyframe >> vzframe;
  fin.close();
  }

void ReadPart() // Read particles file*/
  {
  ifstream fin(PartFile.c_str());
  double Junk,temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9;
  double Lmomx,Lmomy,Lmomz;
  fin >> Junk >> Junk;
  for (int i=0; i <nbod; i++)
    {
    fin >> Junk;
    fin >> temp1 >> temp2 >> temp3;
    fin >> temp4 >> temp5 >> temp6 >> temp7;
    fin >> Junk >> Junk >> temp9;
    temp1=temp1+xframe; temp2=temp2+yframe; temp3=temp3+zframe;
    temp4=temp4+vxframe; temp5=temp5+vyframe; temp6=temp6+vzframe;
    temp8=sqrt(temp1*temp1+temp2*temp2+temp3*temp3);
    temp1=temp1*ru; temp2=temp2*ru; temp3=temp3*ru;
    temp4=temp4*vu; temp5=temp5*vu; temp6=temp6*vu;
    temp8=temp8*ru;
    x[i]=temp1; y[i]=temp2; z[i]=temp3;
    vx[i]=temp4; vy[i]=temp5; vz[i]=temp6;
    tub[i]=temp7;
    Lmflag[i]=temp9;//Leading/trailing flag
    }
  fin.close();
  }

void ReadPartDestroyed() // Read particles file for destroyed sims.  Not checked for new scf with lead/trail flag
  {
  ifstream fin(PartFile.c_str()); // Read R00XS002 file to get unbound data
  double Junk,temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9;
  fin >> Junk >> Junk;

  for (int i=0; i <nbod; i++)
    {
    fin >> Junk;
    fin >> temp1 >> temp2 >> temp3;
    fin >> temp4 >> temp5 >> temp6 >> temp7;
    fin >> Junk >> Junk >> temp9;
    tub[i]=temp7;
    Lmflag[i]=temp9;
    }
  fin.close();

  ifstream fin2(PartFile2.c_str()); // Read R00XS003 file to get position data, frame now taken into account in part file
  fin2 >> Junk >> Junk;
  for (int i=0; i <nbod; i++)
    {
    fin2 >> Junk;
    fin2 >> temp1 >> temp2 >> temp3;
    fin2 >> temp4 >> temp5 >> temp6 >> Junk;
    temp1=temp1*ru; temp2=temp2*ru; temp3=temp3*ru;
    temp4=temp4*vu; temp5=temp5*vu; temp6=temp6*vu;
    x[i]=temp1; y[i]=temp2; z[i]=temp3;
    vx[i]=temp4; vy[i]=temp5; vz[i]=temp6;
    }
  fin2.close();
  }

// Set a color flag for the unbound time of each particle
// 0=yellow, 1=pink, 2=cyan, 3=green, 4=red
void SetColors()
  {
  for (int i=0; i<nentries;i++) t[i]=t[i]*tu;
  int k;
  double Temp;
  for (int i = 0; i < nbod; i++)
    {
    Temp = tub[i];
    k=0;
    for (int j=napo; j>=1; j--)
      {
      k++;
      if ((tapo[j-1]<Temp)&&(tapo[j]>Temp))
	{
	Pcolor[i]=k;
	}
      }
    if (tapo[napo]<Temp) Pcolor[i]=0;
    }
  // Color code all bound particles as Pcolor=-1
  for (int i=0; i<nbod;i++)
    {
    if (tub[i]<0.0001) Pcolor[i]=-1;
    }
  }

void TransformStandard()
  {
  srand48((unsigned)time(NULL));
  double ltemp,btemp,rtemp,vhel;
  OutpFile.append("SGR");
  if (SystNumber==1) OutpFile.append("1");
  else if (SystNumber==2) OutpFile.append("2");
  else if (SystNumber==3) OutpFile.append("3");
  else if (SystNumber==4) OutpFile.append("4");
  srand48((unsigned)time(NULL));

  for (int i=0; i < nbod; i++)
    {
      //    XYZtoLBR(-x[i],y[i],z[i],l[i],b[i],d[i],rsun);
    XYZtoLBRvel(-x[i],y[i],z[i],-vx[i],vy[i],vz[i],l[i],b[i],d[i],mul[i],mub[i],vhel,rsun);
    Kapp[i]=Kabs[i]+5.*log10(d[i]*1000)-5.;
    LBtoRaDec(l[i],b[i],RA[i],Dec[i]);


    if (SystNumber==1) XYZtoSgr1(-x[i],y[i],z[i],Xm[i],Ym[i],Zm[i],d[i],lambda[i],beta[i],rsun);
    else if (SystNumber==2) XYZtoSgr2(-x[i],y[i],z[i],Xm[i],Ym[i],Zm[i],d[i],lambda[i],beta[i],rsun);
    else if (SystNumber==3) XYZtoSgr3(-x[i],y[i],z[i],Xm[i],Ym[i],Zm[i],d[i],lambda[i],beta[i],rsun);
    else if (SystNumber==4) XYZtoSgr4(-x[i],y[i],z[i],Xm[i],Ym[i],Zm[i],d[i],lambda[i],beta[i],rsun);
    vGSR[i]=((x[i]+rsun)*vx[i]+y[i]*vy[i]+z[i]*vz[i])/sqrt((x[i]+rsun)*(x[i]+rsun)+y[i]*y[i]+z[i]*z[i]);
    }
  }

void Output()
  {
  // Kludge the output for Run 601
    /* ifstream finkludge("Runs/Set600/R601SGR1_clean");
  string Junk;
  for (int i=0;i<nbod;i++)
    {
    finkludge >> Junk >> Junk >> Junk >> Junk >> Junk >> Junk >> Junk >> Junk >> Junk >> Junk;
    finkludge >> Lmflag[i] >> Junk;
    }  */

  if (OutpType==1)
    {
    ofstream fout(OutpFile.c_str());
    for (int i=0; i<nbod; i++)
      {
      fout << lambda[i] << " " << beta[i] << " " << RA[i] << " " << Dec[i] << " " << Xm[i] << " " << Ym[i] << " " << Zm[i] << " ";
      fout << d[i] << " " << "-999" << " " << vGSR[i] << " " << Lmflag[i] << " " <<  Pcolor[i] << endl;
      }
    fout.close();
    }

  // Revised web format on 8/31/09
  else if (OutpType==2)
    {
    OutpFile.append(".WEB");
    ofstream fout(OutpFile.c_str());
    double mua,mud,temp1,temp2;
    for (int i=0; i<nbod; i++)
      {
      LBtoRaDecvel(l[i],b[i],mul[i],mub[i],temp1,temp2,mua,mud);
      fout << lambda[i] << " " << beta[i] << " " << l[i] << " " << b[i] << " " << RA[i] << " " << Dec[i] << " ";
      fout << -x[i] << " " << y[i] << " ";
      fout << z[i] << " " << Xm[i] << " " << Ym[i] << " " << Zm[i] << " " << -vx[i] << " " << vy[i] << " ";
      fout << vz[i] << " " << d[i] << " " << "-999" << " " << vGSR[i] << " ";
      fout << mul[i] << " " << mub[i] << " " << mua << " " << mud << " " << Lmflag[i] << " " << Pcolor[i] << endl;
      }
    fout.close();
    }

  else if (OutpType==3)
    {
    srand48((unsigned)time(NULL));
    OutpFile.append(".3D");
    string OutpRed=OutpFile;
    string OutpOrange=OutpFile;
    OutpRed.append("R");
    OutpOrange.append("O");

    ofstream fout1(OutpRed.c_str());
    ofstream fout2(OutpOrange.c_str());
    fout1 << "#declare sgr1 = union {" << endl;
    fout2 << "#declare sgr2 = union {" << endl;
    for (int i=0; i<nbod; i++)
      {
      if (Pcolor[i]<=4)
	{
        if (drand48()<0.5)
	  {
          fout1 << "  object { ster translate <  " << -x[i] << ",\t" << y[i] << ", \t";
          fout1 << z[i] << "> \t}" << endl;
	  }
	else
	  {
	  fout2 << "  object { ster translate <  " << -x[i] << ",\t" << y[i] << ", \t";
          fout2 << z[i] << "> \t}" << endl;
	  }
	}
      }
    fout1 << "}" << endl;
    fout1.close();
    fout2 << "}" << endl;
    fout2.close();
    }

  }

void OutputFehTrend()  
  {
    int Pcolcut=3;// Green debris and younger, ignore Fe/H=-2


  // Do the leading arm
  string Leadfile=OutpFile;
  Leadfile.append(".leadtrend");
  ofstream fout2(Leadfile.c_str());
  vector<double> tempvec(0);
  for (int j=0;j<360;j+=10)
    {
    for (int i=0;i<nbod;i++)
      {
      if ((lambda[i] > j)&&(lambda[i] <= j+10)&&(Pcolor[i]<=Pcolcut)&&(Lmflag[i]>=0)&&(FeH[i] > -1.999))
	{
        tempvec.push_back(FeH[i]);
	}
      }
    fout2 << j+5 << " " << Average(tempvec) << " " << Disp(tempvec) << endl;
    tempvec.resize(0);
    }
  fout2.close();

  // Do the trailing arm
  string Trailfile=OutpFile;
  Trailfile.append(".trailtrend");
  fout2.open(Trailfile.c_str());
  for (int j=0;j<360;j+=10)
    {
    // Special case- ignore green/cyan debris near the core in trailing tail
    if (j < 100) Pcolcut=1;
    else Pcolcut=3;
    for (int i=0;i<nbod;i++)
      {
      if ((lambda[i] > j)&&(lambda[i] <= j+10)&&(Pcolor[i]<=Pcolcut)&&(Lmflag[i]<=0)&&(FeH[i] > -1.999))
	{
        tempvec.push_back(FeH[i]);
	}
      }
    fout2 << j+5 << " " << Average(tempvec) << " " << Disp(tempvec) << endl;
    tempvec.resize(0);
    }
  fout2.close();
  }

int Lcount(string &File)
  {
  ifstream fin(File.c_str());
  int Count=0;
  char InputLine[512];
  while (fin)
    {
    fin.getline(InputLine,512);
    Count++;
    }
  Count--;
  fin.close();
  return Count;
  }
