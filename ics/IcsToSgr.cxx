#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "/Users/drlaw/GenCode/Library/GenTools.h"
#include "/Users/drlaw/GenCode/Library/Coords.h"
#include "/Users/drlaw/GenCode/Library/SgrCoords.h"

using namespace std;

// Global constants
string InpFile,OutpFile;
int lines,SystNumber,OutpType;
double mu,ru,pi,radpdeg,rsun,msun,cmpkpc,kmpkpc,secpgyr,arcpyr,G,tu,vu;
double bmpole,lmpole,langle,ang1,ang2,ang3;
double m,xframe,yframe,zframe,vxframe,vyframe,vzframe;
double rparam,DSGR,PhiCustom,ThetaCustom,PsiCustom,timespec;
char NewPlane;

// Global vectors
vector<double> x(100000); vector<double> y(100000); vector<double> z(100000);
vector<double> vx(100000); vector<double> vy(100000); vector<double> vz(100000);
vector<double> t(100000); vector<double> d(100000);
vector<double> dGC(100000); vector<double> b(100000); vector<double> l(100000);
vector<double> Xm(100000); vector<double> Ym(100000); vector<double> Zm(100000);
vector<double> lambda(100000); vector<double> beta(100000);
vector<double> Xm2(100000); vector<double> Ym2(100000); vector<double> lambda2(100000);
vector<double> XmGC(100000); vector<double> YmGC(100000); vector<double> ZmGC(100000);
vector<double> lambdaGC(100000); vector<double> betaGC(100000);
vector<double> XmGC2(100000); vector<double> YmGC2(100000); vector<double> lambdaGC2(100000);
vector<double> vrhel(100000); vector<double> vGSR(100000); vector<double> bdot(100000);
vector<double> ldot(100000); vector<double> decl(100000); vector<double> rasc(100000);
vector<double> radot(100000); vector<double> decldot(100000); 
vector<double> bdotLSR(100000); vector<double> ldotLSR(100000); vector<double> vLSR(100000);
vector<double> radotLSR(100000); vector<double> decldotLSR(100000); vector<double> Pcolor(100000);
vector<double> tapo(0);
// Function prototypes
void GetInput(); void SetParam(); void ReadFrame(); void ReadPart(); void Transform();
void SetColors(); void Output(); void XYZtoCustom(double X,double Y,double Z,double &Xm1,double &Ym1,double &Zm1,double &d,double &lambda,double &beta,double &rsun);

int main(int argc, char *argv[])
  {
  if (argc==3)
    {
    string name=argv[1];
    timespec=atof(argv[2]);
    InpFile = "./";
    InpFile.append(name.c_str());
    OutpFile = InpFile;
    InpFile.append(".dat");

    DSGR=28.;
    rsun=8.;
    OutpType=2;
    rparam=rsun/DSGR;

    lines=Lcount(InpFile);
    }
  else GetInput(); 
  SetParam();
  ReadFrame();
  Transform();
  Output();
  }

void GetInput() // Get log file and data file names, calculate nentries
  {
  string name;
  cout << "Data file location (e.g. 'orb001'): ";
  cin >> name;
  InpFile = "./";
  InpFile.append(name.c_str());
  OutpFile = InpFile;
  InpFile.append(".dat");

  cout << "DSGR = ";
  cin >> DSGR;
  cout << "rsun = ";
  cin >> rsun;
  cout << "time (1) = ";
  cin >> timespec;
  rparam=rsun/DSGR;

  cout << "What kind of output file?" << endl;
  cout << "(1) Standard DRL output (t, lambda, beta, ra,dec, x, y, z, dist, vGSR, Pcol)" << endl;
  cout << "(2) Tweaked to include Gal coords." << endl;
  cout << "Type: ";
  cin >> OutpType;

  lines=Lcount(InpFile);
  }

void SetParam() // Define constants and vector sizes
  {
  mu=1e9;
  ru=0.9;
  pi=3.141592653589793;
  radpdeg=pi/180.;
  //DSGR=25.;
  //rsun=rparam*DSGR;
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
  }

void ReadFrame()
  {
  x.resize(lines); y.resize(lines); z.resize(lines);
  vx.resize(lines); vy.resize(lines); vz.resize(lines);
  t.resize(lines);
  double Junk,temp1,temp2,temp3,temp4,temp5,temp6;
  ifstream fin(InpFile.c_str());
  for (int i=0; i < lines; i++)
    {
    fin >> t[i];
    fin >> temp1 >> temp2 >> temp3;
    fin >> temp4 >> temp5 >> temp6;
    fin >> Junk >> Junk >> Junk;
    x[i]=temp1; y[i]=temp2; z[i]=temp3;
    vx[i]=temp4; vy[i]=temp5; vz[i]=temp6;
    }
  fin.close();
  }

void Transform()
  {
    if (OutpType==1) OutpFile.append(".SGR1");
    if (OutpType==2) OutpFile.append(".SGR2");
    if (OutpType==3) OutpFile.append(".SGR3");
    for (int i=0; i < lines; i++)
      {
      XYZtoSgr1(-x[i],y[i],z[i],Xm[i],Ym[i],Zm[i],d[i],lambda[i],beta[i],rsun);
      vGSR[i]=((x[i]+rsun)*vx[i]+y[i]*vy[i]+z[i]*vz[i])/sqrt((x[i]+rsun)*(x[i]+rsun)+y[i]*y[i]+z[i]*z[i]);
      }
  }

void Output()
  {
  ofstream fout(OutpFile.c_str());
  double dgc,Junk1,Junk2,Junk3,Junk4,mul,mub,glon,glat,ra,dec;
  for (int i=0; i<lines; i++)
    {
    XYZtoLBRvel(-x[i],y[i],z[i],-vx[i],vy[i],vz[i],glon,glat,Junk3,mul,mub,Junk4,rsun);
    LBtoRaDec(glon,glat,ra,dec);
    if (fabs(t[i])<timespec)
      {
      if (OutpType==1)
        {
	  fout << t[i] << " " << lambda[i] << " " << beta[i] << " " << ra << " " << dec << " " << Xm[i] << " " << Ym[i] << " " << Zm[i] << " ";
        fout << d[i] << " " << vGSR[i] << " " << endl;
	}
      else if (OutpType==2)
	{
	fout << t[i] << " " << lambda[i] << " " << beta[i] << " " << ra << " " << dec << " " << Xm[i] << " " << Ym[i] << " " << Zm[i] << " ";
	fout << -x[i] << " " << y[i] << " " << z[i] << " ";
        fout << d[i] << " " << vGSR[i] << " " << endl;
	}
      else if (OutpType==3)
	{
	if (fabs(t[i])<timespec)
	  {
	  fout << t[i] << " " << lambda[i] << " " << beta[i] << " " << ra << " " << dec << " " << Xm[i] << " " << Ym[i] << " " << Zm[i] << " ";
	  fout << -x[i] << " " << y[i] << " " << z[i] << " ";
          fout << d[i] << " " << vGSR[i] << " " << endl;
	  }
	}
      }
    }
  fout.close();
  }
