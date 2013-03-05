#ifndef BRANT_LUMINOSITY
#define BRANT_LUMINOSITY


 extern float mgu[67626];
 extern float mgg[67626];
 extern float mgr[67626];
 extern float mgi[67626];
 extern float mgz[67626];
 extern float mgj[67626];
 extern float mgh[67626];
 extern float mgk[67626];
 //extern float fms[6,221];
 extern float fms[6*221];
 extern float tk[221];
 extern float zs[51];
 extern int it;
 extern int iz;
 extern int izz;
 extern float mags[8];
 extern float cmags[8];
 extern float s_UBVRIJHK[8];
 extern float s_ugrizJHK[8];


void Initialize_Luminosity(void);
void Load_Magnitudes(void);
double Luminosity_to_Magnitude_UBVRI(double L, int band);
double Luminosity_to_Magnitude_Sloan(double L, int band);
double Magnitude_to_Luminosity_Sloan(double M, int band);
double Magnitude_to_Luminosity_UBVRI(double M, int band);

void Sloan(float mass, float age, float zmet, float zobs, float mags_sloan[]); //returns sloan mags
void UBVRI(float mass, float age, float zmet, float zobs, float mags_ubvri[]); //returns ubvri mags
void L_Sloan(float mass, float age, float zmet, float zobs, float mags_sloan[]); //returns sloan luminosities
void L_UBVRI(float mass, float age, float zmet, float zobs, float mags_ubvri[]); //returns ubvri luminosities
void Get_Mags(float mass, float age, float zmet, float zobs, float mag[]); //returns sloan mags
int glocate(float xx[], int n, float x);
float interpmag(float mgx[], int i, int j, int k, float fi, float fj, float fk, float mass);
//float interpfrc(float fms[6,221], int i, int k, float fi, float fk, float mass);
float interpfrc(float fms[6*221], int i, int k, float fi, float fk, float mass);
float fext(float x,float tv);
#endif
