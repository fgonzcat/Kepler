/////////////////////////////////
// KEPLER.CC: PLANETARY SYSTEM //
/////////////////////////////////
// Creado: 29-Enero-2008.
// Modificacion: 24-Julio-2011.
// Ultima modificacion: 28-Julio-2016.
// Autor: Felipe Gonzalez.

#include "Kepler3D.h"

/************************* WINDOW ******************************/
unsigned char Buttons[3] = {0};
int mousecoord[2];
int win_width=900, win_height=500;
float ancho=150.0, alto=150.0;
float diagonal=sqrt(ancho*ancho+alto*alto);
float dist=diagonal;
float fovy;
float near, far;
float eye[3], center[3], up[3];
float v[3], diag[3], camdir[3];
float tx = 0;
float ty = 0;
bool paused=true;
GLuint planetID, clusterID, bgstarsID, astID;
/***************************************************************/

/*+++++++++++++++++++++++ Units Conversion +++++++++++++++++++++++++++++++*/
long int step=0;
double dt=0.05;// years
const double G=0.00011859645;       // in AU³/(earthmass*year²)
const double AU=1.0;               // Distance in astronomial units (in AU)
const double sunmass=332946.05;    // sunmass=332946.05 earthmass
const double R0=0.5;               // Reference radius (in AU)


const double earthmass =  1.0,	 	Earth_SMA	=1.0,		Earth_R=4.2587571e-05*AU;	// Earth mass, Earth semimajor axis [SMA], Earth radius [UA]
const double sun = 332.830*earthmass,	Sun_SMA		=0.0000, 	Sun_R=109.125*Earth_R;
const double merc = 0.0552*earthmass,	Merc_SMA	=0.3871,	Merc_R=1;
const double venus = 0.814*earthmass,	Venus_SMA	=0.7233,	Venus_R=1;
const double mars = 0.1074*earthmass,	Mars_SMA	=1.5237,	Mars_R=1;
const double jupiter = 318*earthmass,	Jupiter_SMA	=5.2028,	Jupiter_R=1;
const double saturn =  439*earthmass,	Saturn_SMA	=9.5388,	Saturn_R=1;
const double uranus= 14.56*earthmass,	Uranus_SMA	=19.191,	Uranus_R=1;
const double neptune=17.15*earthmass,	Neptune_SMA	=30.061,	Neptune_R=1;
const double pluto = 0.002*earthmass,	Pluto_SMA	=39.529,	Pluto_R=1;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

//------------------------------------------ PARTICLE SETUP ---------------------------------------------------------//
const int num=3;				// Number of planets
double t=0;
const int num_stars=1000;		// Number of stars on the background
double planets_mass=0;			// Sum of all planet's masses
double ast_mass=0;				// Sum of all asteroids' masses
int num_asteroids=0;			// Pieces of dead planets
int dead_ones=0;				// Planets destroyed
bool not_counted[num];			// Identifies if the planet i has already been counted in the list of dead planets
bool wrt=true;
double tau=1e-2;					// Send data to files each tau seconds
int N=0;							// Increasing number of data
Ptcl p[num];
Ptcl ast[4*num];
Ptcl stars[num_stars];
Vector trayec[1000];
Vector trayec2[1000];

char PL[15], AS[15], DEAD[15], TE[15], AM[15], T[15];
double x[num][2];double a_x[4*num][2];
double y[num][2];double a_y[4*num][2];
double z[num][2];double a_z[4*num][2];
double fx[num][2];double a_fx[4*num][2];
double fy[num][2];double a_fy[4*num][2];
double fz[num][2];double a_fz[4*num][2];
double k1x[num][2];double a_k1x[4*num][2];
double k2x[num][2];double a_k2x[4*num][2];
double k3x[num][2];double a_k3x[4*num][2];
double k4x[num][2];double a_k4x[4*num][2];
double k1y[num][2];double a_k1y[4*num][2];
double k2y[num][2];double a_k2y[4*num][2];
double k3y[num][2];double a_k3y[4*num][2];
double k4y[num][2];double a_k4y[4*num][2];
double k1z[num][2];double a_k1z[4*num][2];
double k2z[num][2];double a_k2z[4*num][2];
double k3z[num][2];double a_k3z[4*num][2];
double k4z[num][2];double a_k4z[4*num][2];


ofstream data("datos.dat");


//------------------------------------- MY FUNCTIONS FOR THIS PROGRAM -----------------------------------------------//
/********************************* Gravitational Potential Energy of Two Particles ***********************************/
double Potential(Ptcl m1, Ptcl m2)
{
 double pot_energy=0;
 if (m1.IsAlive() && m2.IsAlive())
 {
  double d=Dist(m1,m2);
  pot_energy=-G*m1.mass()*m2.mass()/d;
  if(m1==m2){pot_energy=0;}
 }
 return pot_energy;
}
/*********************************************************************************************************************/
/******************************************** Potential Energy of The Whole System ***********************************/
double Pot(Ptcl set_of_planets[num],Ptcl set_of_ast[4*num]){
 double pot=0;
 for(int i=0; i<num; i++){if (set_of_planets[i].IsAlive()){
  for(int j=0; j<num; j++){if(i!=j && set_of_planets[j].IsAlive()){
   pot+=Potential(set_of_planets[i],set_of_planets[j]);
  }}
  for(int j=0; j<4*num; j++){if(set_of_ast[j].IsAlive()){
  pot+=Potential(set_of_planets[i],set_of_ast[j]);
  }}
 }}
 for(int i=0; i<4*num; i++){if (set_of_ast[i].IsAlive()){
  for(int j=0; j<num; j++){if(set_of_planets[j].IsAlive()){
   pot+=Potential(set_of_ast[i],set_of_planets[j]);
  }}
  for(int j=0; j<4*num; j++){if(i!=j && set_of_ast[j].IsAlive()){
   pot+=Potential(set_of_ast[i],set_of_ast[j]);
  }}
 }}
 return pot/2.0;
}
/*********************************************************************************************************************/
/******************************************** Kinetic Energy of The Whole System ***********************************/
double Kin(Ptcl set_of_planets[num],Ptcl set_of_ast[4*num]){
 double kin=0;
 for(int i=0; i<num; i++){if (set_of_planets[i].IsAlive()){kin+=set_of_planets[i].Kinetics();}}
 for(int i=0; i<4*num; i++){if (set_of_ast[i].IsAlive()){kin+=set_of_ast[i].Kinetics();}}
 return kin;
}
/*********************************************************************************************************************/
/******************************************** Angular Momentum of The Whole System ***********************************/
Vector TL(Ptcl set_of_planets[num],Ptcl set_of_ast[4*num]){
 Vector l(0,0,0,0);
 for(int i=0; i<num; i++){if (set_of_planets[i].IsAlive()){l+=set_of_planets[i].L();}}
 for(int i=0; i<4*num; i++){if (set_of_ast[i].IsAlive()){l+=set_of_ast[i].L();}}
 return l;
}
/*********************************************************************************************************************/
/********************************** Gravitational Force between particles (Newton) ***********************************/
Vector GForce(Ptcl m1, Ptcl m2){
 double d=Dist(m1,m2);
 double fact=-G*m1.mass()*m2.mass()/(d*d*d);
 Vector f1=fact*(m1.r()-m2.r());
 if (m1==m2){f1=0*e1;}
 return f1;
}
/*********************************************************************************************************************/
/************************************************ Center of Mass *****************************************************/
Vector COM(Ptcl set_of_planets[num], Ptcl set_of_ast[4*num]){
 double Mp=0, Ma=0;
 Vector centre;
 for (int i=0; i<num; i++){if (set_of_planets[i].IsAlive()){	Mp+=set_of_planets[i].mass();}}
 for (int j=0; j<4*num; j++){if (set_of_ast[j].IsAlive()){	Ma+=set_of_ast[j].mass();}}
 for (int i=0; i<num; i++){if (set_of_planets[i].IsAlive()){	centre=centre+( set_of_planets[i].mass()/(Mp+Ma) )*set_of_planets[i].r();}}
 for (int j=0; j<4*num; j++){if (set_of_ast[j].IsAlive()){	centre=centre+( set_of_ast[j].mass()/(Mp+Ma) )*set_of_ast[j].r();}}
 if (Mp==0){centre=0*e1;}
 return centre;
}
/*********************************************************************************************************************/
/************************************************ Integrator *********************************************************/
void Solver(void){
/****************************************Planets' Motion Equations****************************************************/
//-----------------------------------------------Fixing f------------------------------------------------------------//
for (int i=0; i<num; i++){
 fx[i][0]=x[i][1];
 fy[i][0]=y[i][1];
 fz[i][0]=z[i][1];
 fx[i][1]=0;	fy[i][1]=0;	fz[i][1]=0;
 for (int j=0; j<num; j++){if (p[j].IsAlive()){
  fx[i][1]+=GForce(p[i],p[j])[0]/p[i].mass();
  fy[i][1]+=GForce(p[i],p[j])[1]/p[i].mass();
  fz[i][1]+=GForce(p[i],p[j])[2]/p[i].mass();
 }}
 for (int j=0; j<4*num; j++){if (ast[j].IsAlive()){
  fx[i][1]+=GForce(p[i],ast[j])[0]/p[i].mass();
  fy[i][1]+=GForce(p[i],ast[j])[1]/p[i].mass();
  fz[i][1]+=GForce(p[i],ast[j])[2]/p[i].mass();
 }}
}	
//-------------------------------------------------------------------------------------------------------------------//

//-----------------------------------------------Fixing k1-----------------------------------------------------------//
// k1xij=dt*fxij(x11,x12,...,x21,...)
	for (int i=0; i<num; i++){
		k1x[i][0]=dt*fx[i][0];
		k1y[i][0]=dt*fy[i][0];
		k1z[i][0]=dt*fz[i][0];

		k1x[i][1]=dt*fx[i][1];
		k1y[i][1]=dt*fy[i][1];
		k1z[i][1]=dt*fz[i][1];
	}
//-------------------------------------------------------------------------------------------------------------------//

//-----------------------------------------------Fixing k2-----------------------------------------------------------//
// k2xij=dt*fxij(x11+k1x11/2,x12+k1x12/2,...,x21+k1x21/2,...)
	//**************Next is used to evaluate fxij in xij+k1xij/2************//
	for (int i=0; i<num; i++){
		x[i][0]+=k1x[i][0]/2;	y[i][0]+=k1y[i][0]/2;	z[i][0]+=k1z[i][0]/2;
		x[i][1]+=k1x[i][1]/2;	y[i][1]+=k1y[i][1]/2;	z[i][1]+=k1z[i][1]/2;
		p[i].set_r(x[i][0]*e1+y[i][0]*e2+z[i][0]*e3);
		p[i].set_v(x[i][1]*e1+y[i][1]*e2+z[i][1]*e3);
	}
	//**********************************************************************//
	//*******************************Fixing*********************************//
	for (int i=0; i<num; i++){
		k2x[i][0]=dt*x[i][1];
		k2y[i][0]=dt*y[i][1];
		k2z[i][0]=dt*z[i][1];
		
		k2x[i][1]=0; k2y[i][1]=0; k2z[i][1]=0;
		for (int j=0; j<num; j++){if (p[j].IsAlive()){
			k2x[i][1]+=dt*GForce(p[i],p[j])[0]/p[i].mass();
			k2y[i][1]+=dt*GForce(p[i],p[j])[1]/p[i].mass();
			k2z[i][1]+=dt*GForce(p[i],p[j])[2]/p[i].mass();
		}}
		for (int j=0; j<4*num; j++){if (ast[j].IsAlive()){
			k2x[i][1]+=dt*GForce(p[i],ast[j])[0]/p[i].mass();
			k2y[i][1]+=dt*GForce(p[i],ast[j])[1]/p[i].mass();
			k2z[i][1]+=dt*GForce(p[i],ast[j])[2]/p[i].mass();
		}}
	}
	//***********************************************************************//
	//***********Next is used to undo the evaluation of fxij above***********//
	for (int i=0; i<num; i++){
		x[i][0]-=k1x[i][0]/2;	y[i][0]-=k1y[i][0]/2;	z[i][0]-=k1z[i][0]/2;
		x[i][1]-=k1x[i][1]/2;	y[i][1]-=k1y[i][1]/2;	z[i][1]-=k1z[i][1]/2;
		p[i].set_r(x[i][0]*e1+y[i][0]*e2+z[i][0]*e3);
		p[i].set_v(x[i][1]*e1+y[i][1]*e2+z[i][1]*e3);
	}
	//**********************************************************************//
//-------------------------------------------------------------------------------------------------------------------//

//-----------------------------------------------Fixing k3-----------------------------------------------------------//
// k3xij=dt*fxij(x11+k2x11/2,x12+k2x12/2,...,x21+k2x21/2,...)
	//**************Next is used to evaluate fxij in xij+k2xij/2************//
	for (int i=0; i<num; i++){
		x[i][0]+=k2x[i][0]/2;	y[i][0]+=k2y[i][0]/2;	z[i][0]+=k2z[i][0]/2;
		x[i][1]+=k2x[i][1]/2;	y[i][1]+=k2y[i][1]/2;	z[i][1]+=k2z[i][1]/2;
		p[i].set_r(x[i][0]*e1+y[i][0]*e2+z[i][0]*e3);
		p[i].set_v(x[i][1]*e1+y[i][1]*e2+z[i][1]*e3);
	}
	//**********************************************************************//
	//*******************************Fixing*********************************//
	for (int i=0; i<num; i++){
		k3x[i][0]=dt*x[i][1];
		k3y[i][0]=dt*y[i][1];
		k3z[i][0]=dt*z[i][1];
		
		k3x[i][1]=0; k3y[i][1]=0; k3z[i][1]=0;
		for (int j=0; j<num; j++){if (p[j].IsAlive()){
			k3x[i][1]+=dt*GForce(p[i],p[j])[0]/p[i].mass();
			k3y[i][1]+=dt*GForce(p[i],p[j])[1]/p[i].mass();
			k3z[i][1]+=dt*GForce(p[i],p[j])[2]/p[i].mass();
		}}
		for (int j=0; j<4*num; j++){if (ast[j].IsAlive()){
			k3x[i][1]+=dt*GForce(p[i],ast[j])[0]/p[i].mass();
			k3y[i][1]+=dt*GForce(p[i],ast[j])[1]/p[i].mass();
			k3z[i][1]+=dt*GForce(p[i],ast[j])[2]/p[i].mass();
		}}
	}
	//**********************************************************************//
	//**********Next is used to undo the evaluation of fxij above***********//
	for (int i=0; i<num; i++){
		x[i][0]-=k2x[i][0]/2;	y[i][0]-=k2y[i][0]/2;	z[i][0]-=k2z[i][0]/2;
		x[i][1]-=k2x[i][1]/2;	y[i][1]-=k2y[i][1]/2;	z[i][1]-=k2z[i][1]/2;
		p[i].set_r(x[i][0]*e1+y[i][0]*e2+z[i][0]*e3);
		p[i].set_v(x[i][1]*e1+y[i][1]*e2+z[i][1]*e3);
	}
	//**********************************************************************//
//-------------------------------------------------------------------------------------------------------------------//

//-----------------------------------------------Fixing k4-----------------------------------------------------------//
// k4xij=fxij(x11+k3x11,x12+k3x12,...,x21+k3x21,...)
	//**************Next is used to evaluate fxij in xij+k3xij**************//
	for (int i=0; i<num; i++){
		x[i][0]+=k3x[i][0];	y[i][0]+=k3y[i][0];	z[i][0]+=k3z[i][0];
		x[i][1]+=k3x[i][1];	y[i][1]+=k3y[i][1];	z[i][1]+=k3z[i][1];
		p[i].set_r(x[i][0]*e1+y[i][0]*e2+z[i][0]*e3);
		p[i].set_v(x[i][1]*e1+y[i][1]*e2+z[i][1]*e3);
	}
	//**********************************************************************//
	//*******************************Fixing*********************************//
	for (int i=0; i<num; i++){
		k4x[i][0]=dt*x[i][1];
		k4y[i][0]=dt*y[i][1];
		k4z[i][0]=dt*z[i][1];

		k4x[i][1]=0; k4y[i][1]=0; k4z[i][1]=0;
		for (int j=0; j<num; j++){if (p[j].IsAlive()){
			k4x[i][1]+=dt*GForce(p[i],p[j])[0]/p[i].mass();
			k4y[i][1]+=dt*GForce(p[i],p[j])[1]/p[i].mass();
			k4z[i][1]+=dt*GForce(p[i],p[j])[2]/p[i].mass();
		}}
		for (int j=0; j<4*num; j++){if (ast[j].IsAlive()){
			k4x[i][1]+=dt*GForce(p[i],ast[j])[0]/p[i].mass();
			k4y[i][1]+=dt*GForce(p[i],ast[j])[1]/p[i].mass();
			k4z[i][1]+=dt*GForce(p[i],ast[j])[2]/p[i].mass();
		}}
	}
	//***********************************************************************//
	//***********Next is used to undo the evaluation of fxij above***********//
	for (int i=0; i<num; i++){
		x[i][0]-=k3x[i][0];	y[i][0]-=k3y[i][0];	z[i][0]-=k3z[i][0];
		x[i][1]-=k3x[i][1];	y[i][1]-=k3y[i][1];	z[i][1]-=k3z[i][1];
		p[i].set_r(x[i][0]*e1+y[i][0]*e2+z[i][0]*e3);
		p[i].set_v(x[i][1]*e1+y[i][1]*e2+z[i][1]*e3);
	}
	//**********************************************************************//
//-------------------------------------------------------------------------------------------------------------------//

//------------------------------------------Computing the Solutions--------------------------------------------------//
	for (int i=0; i<num; i++){
		x[i][0]+=(k1x[i][0]+2.0*k2x[i][0]+2.0*k3x[i][0]+k4x[i][0])/6.0;
		x[i][1]+=(k1x[i][1]+2.0*k2x[i][1]+2.0*k3x[i][1]+k4x[i][1])/6.0;
		y[i][0]+=(k1y[i][0]+2.0*k2y[i][0]+2.0*k3y[i][0]+k4y[i][0])/6.0;
		y[i][1]+=(k1y[i][1]+2.0*k2y[i][1]+2.0*k3y[i][1]+k4y[i][1])/6.0;
		z[i][0]+=(k1z[i][0]+2.0*k2z[i][0]+2.0*k3z[i][0]+k4z[i][0])/6.0;
		z[i][1]+=(k1z[i][1]+2.0*k2z[i][1]+2.0*k3z[i][1]+k4z[i][1])/6.0;
	}	
//-------------------------------------------------------------------------------------------------------------------//

//----------------------------------------New Positions and Velocities-----------------------------------------------//
	for (int i=0; i<num; i++){
		p[i].set_r(x[i][0]*e1+y[i][0]*e2+z[i][0]*e3);
		p[i].set_v(x[i][1]*e1+y[i][1]*e2+z[i][1]*e3);
	}
//-------------------------------------------------------------------------------------------------------------------//

//-----------------------------------------------Crushings-----------------------------------------------------------//
/*	for (int i=0; i<num; i++){
		//------------------- Killing Planets --------------------------------//
		for (int j=0; j<num;j++){
			if (i!=j && Crushing(p[i],p[j]) && p[i].IsAlive() && p[j].IsAlive()){
				if 	 (fabs(p[i].rho()-p[j].rho())<1e-10){	p[i].kill(); p[j].kill();}
				else if (p[i].rho()>p[j].rho()){	p[j].kill();}
				else if (p[i].rho()<p[j].rho()){	p[i].kill();}
			}
		}
		for (int j=0; j<4*num;j++){
			if (Crushing(p[i],ast[j]) && p[i].IsAlive() && ast[j].IsAlive()){
				if 	 (fabs(p[i].rho()-ast[j].rho())<1e-10){	p[i].kill(); ast[j].kill(); num_asteroids-=1;}
				else if (p[i].rho()>ast[j].rho()){	ast[j].kill(); num_asteroids-=1;}
				else if (p[i].rho()<ast[j].rho()){	p[i].kill();}
			}
		}
		//--------------------------------------------------------------//
		//--------------------- Creating Asteroids --------------------//
		if (p[i].IsDead() && not_counted[i]){
			num_asteroids+=4;
			p[i].set_tod(t);	p[i].set_pod(p[i].r()); p[i].set_vod(p[i].v());
			for (int j=4*i; j<4*i+4; j++){
				//----------- Asteroids Boundary Conditions -------------//
				ast[j].resurrect();
				ast[j].set_mass(0.25*p[i].mass());
				ast[j].set_size(0.5*p[i].size());
				ast[j].set_red(p[i].R()); ast[j].set_green(p[i].G()); ast[j].set_blue(p[i].B());
				ast[j].set_r(p[i].pod()+.75*p[i].size()*(cos(0.25*M_PI+(j-4*i)*M_PI/2)*e1+sin(0.25*M_PI+(j-4*i)*M_PI/2)*e2));
				ast[j].set_v(p[i].vod().Rotation(dazar(0,2*M_PI),dazar(0,2*M_PI),dazar(0,2*M_PI)));
				//------------------------------------------------------//

				//------------------------------ Initialize a_x[], a_y[], a_z[] -----------------------------//
				a_x[j][0]=ast[j].r()[0];		a_y[j][0]=ast[j].r()[1];		a_z[j][0]=ast[j].r()[2];
				a_x[j][1]=ast[j].v()[0];		a_y[j][1]=ast[j].v()[1];		a_z[j][1]=ast[j].v()[2];
				//-------------------------------------------------------------------------------------------//
			}
		}
		//-------------------------------------------------------------//
	}*/
//-------------------------------------------------------------------------------------------------------------------//
/*********************************************************************************************************************/

/************************************************Counters*************************************************************/
//-------------------------------------------Counting Dead Planets---------------------------------------------------//
	for (int i=0; i<num; i++){
		if (p[i].IsDead() && not_counted[i]){
			not_counted[i]=false;
			dead_ones++;
		}
	}
//-------------------------------------------------------------------------------------------------------------------//

//-------------------------------------------Counting Masses---------------------------------------------------------//
	planets_mass=0; ast_mass=0;
	for (int i=0; i<num; i++){if (p[i].IsAlive())	planets_mass+=p[i].mass();	}
	for (int i=0; i<4*num; i++){if (ast[i].IsAlive())	ast_mass+=ast[i].mass();	}
//-------------------------------------------------------------------------------------------------------------------//
/*********************************************************************************************************************/

/****************************************Asteroids' Motion Ecuations**************************************************/
//---------------------------------------------- Fixing f -----------------------------------------------------------//
	for (int i=0; i<num; i++){if (p[i].IsDead()){

		for (int j=4*i; j<4*i+4; j++){
			a_fx[j][0]=a_x[j][1];
			a_fy[j][0]=a_y[j][1];
			a_fz[j][0]=a_z[j][1];
			a_fx[j][1]=0;	a_fy[j][1]=0;	a_fz[j][1]=0;
			for (int l=0; l<4*num; l++){if (ast[l].IsAlive()){
				a_fx[j][1]+=GForce(ast[j],ast[l])[0]/ast[j].mass();
				a_fy[j][1]+=GForce(ast[j],ast[l])[1]/ast[j].mass();
				a_fz[j][1]+=GForce(ast[j],ast[l])[2]/ast[j].mass();
			}}
			for (int l=0; l<num; l++){if (p[l].IsAlive()){
				a_fx[j][1]+=GForce(ast[j],p[l])[0]/ast[j].mass();
				a_fy[j][1]+=GForce(ast[j],p[l])[1]/ast[j].mass();
				a_fz[j][1]+=GForce(ast[j],p[l])[2]/ast[j].mass();
			}}
		}
//-------------------------------------------------------------------------------------------------------------------//

//---------------------------------------------- Fixing a_k1 --------------------------------------------------------//
// a_k1xij=dt*a_fxij(a_x11,a_x12,...,a_x21,...)
		for (int j=4*i; j<4*i+4; j++){
			a_k1x[j][0]=dt*a_fx[j][0];
			a_k1y[j][0]=dt*a_fy[j][0];
			a_k1z[j][0]=dt*a_fz[j][0];

			a_k1x[j][1]=dt*a_fx[j][1];
			a_k1y[j][1]=dt*a_fy[j][1];
			a_k1z[j][1]=dt*a_fz[j][1];
		}
//-------------------------------------------------------------------------------------------------------------------//

//---------------------------------------------- Fixing a_k2 --------------------------------------------------------//
// a_k2xij=dt*a_fxij(a_x11+a_k1x11/2,a_x12+a_k1x12/2,...,a_x21+a_k1x21/2,...)
		//**************Next is used to evaluate a_fxij in a_xij+a_k1xij/2************//
		for (int j=4*i; j<4*i+4; j++){
			a_x[j][0]+=a_k1x[j][0]/2;	a_y[j][0]+=a_k1y[j][0]/2;	a_z[j][0]+=a_k1z[j][0]/2;
			a_x[j][1]+=a_k1x[j][1]/2;	a_y[j][1]+=a_k1y[j][1]/2;	a_z[j][1]+=a_k1z[j][1]/2;
			ast[j].set_r(a_x[j][0]*e1+a_y[j][0]*e2+a_z[j][0]*e3);
			ast[j].set_v(a_x[j][1]*e1+a_y[j][1]*e2+a_z[j][1]*e3);
		}
		//**********************************************************************//
		//*******************************Fixing*********************************//
		for (int j=4*i; j<4*i+4; j++){
			a_k2x[j][0]=dt*a_x[j][1];
			a_k2y[j][0]=dt*a_y[j][1];
			a_k2z[j][0]=dt*a_z[j][1];
		
			a_k2x[j][1]=0; a_k2y[j][1]=0; a_k2z[j][1]=0;
			for (int l=0; l<4*num; l++){if (ast[l].IsAlive()){
				a_k2x[j][1]+=dt*GForce(ast[j],ast[l])[0]/ast[j].mass();
				a_k2y[j][1]+=dt*GForce(ast[j],ast[l])[1]/ast[j].mass();
				a_k2z[j][1]+=dt*GForce(ast[j],ast[l])[2]/ast[j].mass();
			}}
			for (int l=0; l<num; l++){if (p[l].IsAlive()){
				a_k2x[j][1]+=dt*GForce(ast[j],p[l])[0]/ast[j].mass();
				a_k2y[j][1]+=dt*GForce(ast[j],p[l])[1]/ast[j].mass();
				a_k2z[j][1]+=dt*GForce(ast[j],p[l])[2]/ast[j].mass();
			}}
		}
		//***********************************************************************//
		//***********Next is used to undo the evaluation of fxij above***********//
		for (int j=4*i; j<4*i+4; j++){
			a_x[j][0]-=a_k1x[j][0]/2;	a_y[j][0]-=a_k1y[j][0]/2;	a_z[j][0]-=a_k1z[j][0]/2;
			a_x[j][1]-=a_k1x[j][1]/2;	a_y[j][1]-=a_k1y[j][1]/2;	a_z[j][1]-=a_k1z[j][1]/2;
			ast[j].set_r(a_x[j][0]*e1+a_y[j][0]*e2+a_z[j][0]*e3);
			ast[j].set_v(a_x[j][1]*e1+a_y[j][1]*e2+a_z[j][1]*e3);
		}
	//**********************************************************************//
//-------------------------------------------------------------------------------------------------------------------//

//---------------------------------------------- Fixing a_k3 --------------------------------------------------------//
// a_k3xij=dt*a_fxij(a_x11+a_k2x11/2,a_x12+k2x12/2,...,a_x21+a_k2x21/2,...)
		//**************Next is used to evaluate fxij in xij+k2xij/2************//
		for (int j=4*i; j<4*i+4; j++){
			a_x[j][0]+=a_k2x[j][0]/2;	a_y[j][0]+=a_k2y[j][0]/2;	a_z[j][0]+=a_k2z[j][0]/2;
			a_x[j][1]+=a_k2x[j][1]/2;	a_y[j][1]+=a_k2y[j][1]/2;	a_z[j][1]+=a_k2z[j][1]/2;
			ast[j].set_r(a_x[j][0]*e1+a_y[j][0]*e2+a_z[j][0]*e3);
			ast[j].set_v(a_x[j][1]*e1+a_y[j][1]*e2+a_z[j][1]*e3);
		}
		//**********************************************************************//
		//*******************************Fixing*********************************//
		for (int j=4*i; j<4*i+4; j++){
			a_k3x[j][0]=dt*a_x[j][1];
			a_k3y[j][0]=dt*a_y[j][1];
			a_k3z[j][0]=dt*a_z[j][1];
		
			a_k3x[j][1]=0; a_k3y[j][1]=0; a_k3z[j][1]=0;
			for (int l=0; l<4*num; l++){if (ast[l].IsAlive()){
				a_k3x[j][1]+=dt*GForce(ast[j],ast[l])[0]/ast[j].mass();
				a_k3y[j][1]+=dt*GForce(ast[j],ast[l])[1]/ast[j].mass();
				a_k3z[j][1]+=dt*GForce(ast[j],ast[l])[2]/ast[j].mass();
			}}
			for (int l=0; l<num; l++){if (p[l].IsAlive()){
				a_k3x[j][1]+=dt*GForce(ast[j],p[l])[0]/ast[j].mass();
				a_k3y[j][1]+=dt*GForce(ast[j],p[l])[1]/ast[j].mass();
				a_k3z[j][1]+=dt*GForce(ast[j],p[l])[2]/ast[j].mass();
			}}
		}
		//**********************************************************************//
		//**********Next is used to undo the evaluation of fxij above***********//
		for (int j=4*i; j<4*i+4; j++){
			a_x[j][0]-=a_k2x[j][0]/2;	a_y[j][0]-=a_k2y[j][0]/2;	a_z[j][0]-=a_k2z[j][0]/2;
			a_x[j][1]-=a_k2x[j][1]/2;	a_y[j][1]-=a_k2y[j][1]/2;	a_z[j][1]-=a_k2z[j][1]/2;
			ast[j].set_r(a_x[j][0]*e1+a_y[j][0]*e2+a_z[j][0]*e3);
			ast[j].set_v(a_x[j][1]*e1+a_y[j][1]*e2+a_z[j][1]*e3);
		}
		//**********************************************************************//
//-------------------------------------------------------------------------------------------------------------------//

//---------------------------------------------- Fixing a_k4 --------------------------------------------------------//
// a_k4xij=a_fxij(a_x11+a_k3x11,a_x12+a_k3x12,...,a_x21+a_k3x21,...)
		//**************Next is used to evaluate fxij in xij+k2xij/2************//
		for (int j=4*i; j<4*i+4; j++){
			a_x[j][0]+=a_k3x[j][0];	a_y[j][0]+=a_k3y[j][0];	a_z[j][0]+=a_k3z[j][0];
			a_x[j][1]+=a_k3x[j][1];	a_y[j][1]+=a_k3y[j][1];	a_z[j][1]+=a_k3z[j][1];
			ast[j].set_r(a_x[j][0]*e1+a_y[j][0]*e2+a_z[j][0]*e3);
			ast[j].set_v(a_x[j][1]*e1+a_y[j][1]*e2+a_z[j][1]*e3);
		}
		//**********************************************************************//
		//*******************************Fixing*********************************//
		for (int j=4*i; j<4*i+4; j++){
			a_k4x[j][0]=dt*a_x[j][1];
			a_k4y[j][0]=dt*a_y[j][1];
			a_k4z[j][0]=dt*a_z[j][1];

			a_k4x[j][1]=0; a_k4y[j][1]=0; a_k4z[j][1]=0;
			for (int l=0; l<4*num; l++){if (ast[l].IsAlive()){
				a_k4x[j][1]+=dt*GForce(ast[j],ast[l])[0]/ast[j].mass();
				a_k4y[j][1]+=dt*GForce(ast[j],ast[l])[1]/ast[j].mass();
				a_k4z[j][1]+=dt*GForce(ast[j],ast[l])[2]/ast[j].mass();
			}}
			for (int l=0; l<num; l++){if (p[l].IsAlive()){
				a_k4x[j][1]+=dt*GForce(ast[j],p[l])[0]/ast[j].mass();
				a_k4y[j][1]+=dt*GForce(ast[j],p[l])[1]/ast[j].mass();
				a_k4z[j][1]+=dt*GForce(ast[j],p[l])[2]/ast[j].mass();
			}}
		}
		//***********************************************************************//
		//***********Next is used to undo the evaluation of fxij above***********//
		for (int j=4*i; j<4*i+4; j++){
			a_x[j][0]-=a_k3x[j][0];	a_y[j][0]-=a_k3y[j][0];	a_z[j][0]-=a_k3z[j][0];
			a_x[j][1]-=a_k3x[j][1];	a_y[j][1]-=a_k3y[j][1];	a_z[j][1]-=a_k3z[j][1];
			ast[j].set_r(a_x[j][0]*e1+a_y[j][0]*e2+a_z[j][0]*e3);
			ast[j].set_v(a_x[j][1]*e1+a_y[j][1]*e2+a_z[j][1]*e3);
		}
		//**********************************************************************//
//-------------------------------------------------------------------------------------------------------------------//

//------------------------------------------Computing the Solutions--------------------------------------------------//
		for (int j=4*i; j<4*i+4; j++){
			a_x[j][0]+=(a_k1x[j][0]+2.0*a_k2x[j][0]+2.0*a_k3x[j][0]+a_k4x[j][0])/6.0;
			a_x[j][1]+=(a_k1x[j][1]+2.0*a_k2x[j][1]+2.0*a_k3x[j][1]+a_k4x[j][1])/6.0;
			a_y[j][0]+=(a_k1y[j][0]+2.0*a_k2y[j][0]+2.0*a_k3y[j][0]+a_k4y[j][0])/6.0;
			a_y[j][1]+=(a_k1y[j][1]+2.0*a_k2y[j][1]+2.0*a_k3y[j][1]+a_k4y[j][1])/6.0;
			a_z[j][0]+=(a_k1z[j][0]+2.0*a_k2z[j][0]+2.0*a_k3z[j][0]+a_k4z[j][0])/6.0;
			a_z[j][1]+=(a_k1z[j][1]+2.0*a_k2z[j][1]+2.0*a_k3z[j][1]+a_k4z[j][1])/6.0;
		}
//-------------------------------------------------------------------------------------------------------------------//

//----------------------------------------New Positions and Velocities-----------------------------------------------//
		for (int j=4*i; j<4*i+4; j++){
			ast[j].set_r(a_x[j][0]*e1+a_y[j][0]*e2+a_z[j][0]*e3);
			ast[j].set_v(a_x[j][1]*e1+a_y[j][1]*e2+a_z[j][1]*e3);
		}
//-------------------------------------------------------------------------------------------------------------------//

//-----------------------------------------------Crushings-----------------------------------------------------------//
/*		for (int j=0; j<4*num; j++){
			for (int l=0; l<4*num;l++){
				if (j!=l && Crushing(ast[j],ast[l]) && ast[j].IsAlive() && ast[l].IsAlive()){
					if 	 (fabs(ast[j].rho()-ast[l].rho())<1e-10){	ast[j].kill(); ast[l].kill(); num_asteroids-=2;}
					else if (ast[j].rho()>ast[l].rho()){	ast[l].kill(); num_asteroids-=1;}
					else if (ast[j].rho()<ast[l].rho()){	ast[j].kill(); num_asteroids-=1;}
				}
			}
		}*/
//-------------------------------------------------------------------------------------------------------------------//
	}}
/*********************************************************************************************************************/
}
/*********************************************************************************************************************/
/************************************************** Send Data ********************************************************/
void Send_Data(){
	if (fabs(t-N*tau)<=1e-3){
		data.precision(5);
		data << t << "\t" << Kin(p,ast) << "\t\t" << Pot(p,ast) << "\t\t" << Kin(p,ast)+Pot(p,ast) << "\t" << TL(p,ast) << endl;
		N++;
	}
}
/*********************************************************************************************************************/


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// BOUNDARY CONDITIONS ///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void 	Boundary_Conditions(void){
 randomize();
 // data << "Time" << "\t" << "Kinetics" << "\t" << "Potential"<< "\t" << "Energy" << "\t\t" << "Angular Moment" << endl;
 
/********************************************PLANETS******************************************************************/
// PLANET'S MASSES
/* for (int i=0; i<num; i++){p[i].set_mass(1);}
 p[0].set_mass(1e7);
 p[5].set_mass(300);
 p[3].set_mass(1e5);
 p[6].set_mass(1e-3);

 // PLANET'S RATIO
 for (int i=0; i<num; i++){	p[i].set_size(1);}
 p[0].set_size(3);
 p[5].set_size(0.8);
 p[6].set_size(0.4);

 // INITIAL POSITIONS
 p[0].set_r(5*e1+0*e2);
 p[1].set_r(0*e1+7*e2);
 p[2].set_r(-7*e1+0*e2);
 p[3].set_r(12*e1-0*e2);
 p[4].set_r(0*e1-15*e2);
 p[5].set_r(0*e1+17*e2);
 p[6].set_r(20*e1+20*e2);
 p[7].set_r(-15*e1-15*e2);

 // INICIAL VELOCITIES
 //	for (int i=0; i<num; i++){ p[i].set_v(pow(-1,double(i))*0.5*e2);}
 p[0].set_v(0*e1+0*e2);
 p[1].set_v(-13*e1-3*e2+5*e3);
 p[2].set_v(0*e1-10*e2);
 p[3].set_v(0*e1-15*e2);
 p[4].set_v(10*e1+0*e2);
 p[5].set_v(-5*e1-0*e2);
 p[6].set_v(2*e1-3*e2);
 p[7].set_v(-0*e1-50*e2);*/
	// PLANET'S MASSES
	for (int i=0; i<num; i++){p[i].set_mass(earthmass);}
	p[0].set_mass(100*sunmass);

	// PLANET'S RATIO
	for (int i=0; i<num; i++){	p[i].set_size(5);}
	p[0].set_size(10);

	// INITIAL POSITIONS
	p[0].set_r(10*e1+0*e2);
	p[1].set_r(-50*e1+0*e2);
	p[2].set_r(80*e1+0*e2);

	// INICIAL VELOCITIES    (in AU/year)
//	for (int i=0; i<num; i++){ p[i].set_v(pow(-1,double(i))*0.5*e2);}
	p[0].set_v(0*e1+0*e2);
	p[1].set_v(1*e1+4*e2);
	p[2].set_v(-.5*e1-5*e2);

 // COLORS
 for (int i=0; i<num; i++){ p[i].set_red(0); p[i].set_green(0.8); p[i].set_blue(1);} // Default=blue
 p[0].set_red(1); p[0].set_green(1); p[0].set_blue(0);
 p[1].set_red(1); p[1].set_green(0); p[1].set_blue(0);
 p[2].set_red(1); p[2].set_green(0.5); p[2].set_blue(0);
 p[3].set_red(1); p[3].set_green(0.9); p[3].set_blue(0.5);
 p[4].set_red(.3); p[4].set_green(0.3); p[4].set_blue(.8);
 p[6].set_red(.3); p[6].set_green(1); p[6].set_blue(0);


 /********************************************STARS********************************************************************/
 for (int i=0; i<num_stars; i++)
 {
 stars[i].set_size(dazar(0,2));
 stars[i].set_r(dazar(-2*ancho,2*ancho)*e1+dazar(-2*ancho,2*ancho)*e2+dazar(-2*ancho,2*ancho)*e3);
 }
 /*********************************************************************************************************************/

 /********************************************ASTEROIDS****************************************************************/
 for (int i=0; i<4*num; i++)
 {
 ast[i].set_mass(0);
 ast[i].set_size(0);
 ast[i].set_r(0*e1);
 ast[i].set_v(0*e1);
 ast[i].kill();
 }
 /*********************************************************************************************************************/

 /**************************Initialize x[][], y[][] & z[][]**************************/
 for (int i=0; i<num; i++){
 x[i][0]=p[i].r()[0];		y[i][0]=p[i].r()[1];		z[i][0]=p[i].r()[2];
 x[i][1]=p[i].v()[0];		y[i][1]=p[i].v()[1];		z[i][1]=p[i].v()[2];
 not_counted[i]=true;
 }
 /***********************************************************************************/
}

void drawString(char *s){
     unsigned int i;
     for (i = 0; i < strlen(s); i++){	glutBitmapCharacter (GLUT_BITMAP_HELVETICA_12, s[i]);}
}


void drawStringBig(char *s){
     unsigned int i;
     for (i = 0; i < strlen(s); i++){	glutBitmapCharacter (GLUT_BITMAP_HELVETICA_18, s[i]);}
}

//------------------------------------------GLUT FUNCTIONS-----------------------------------------------------------//

//------------------------------ MakePlanet ------------------------------//
void MakePlanet(GLuint id)
{
 int nlon=16,nlat=10, radius=1;
 int i,j;
 float lon,lat;
 float loninc=2*M_PI/nlon,latinc=M_PI/nlat;
 float x,y,z;
 
 glNewList(id,GL_COMPILE);
 /* South-pole triangular fan */
 glBegin(GL_TRIANGLE_FAN);
 glNormal3f(0,-1,0);
 glVertex3f(0,-radius,0);
 lon = 0;
 lat = -M_PI/2 + latinc;
 y = sin(lat);
 for (i=0; i<=nlon; i++)
 {
  x = cos(lon)*cos(lat);
  z = -sin(lon)*cos(lat);
  glNormal3f(x,y,z);
  glVertex3f(x*radius,y*radius,z*radius);
  lon += loninc;
 }
 glEnd();

 /* Quadrilateral stripes to cover the sphere */
 for (j=1; j<nlat-1; j++)
 {
  lon = 0;
  glBegin(GL_QUAD_STRIP);
  for (i=0; i<=nlon; i++)
  {
   x = cos(lon)*cos(lat);
   y = sin(lat);
   z = -sin(lon)*cos(lat);
   glNormal3f(x,y,z);
   glVertex3f(x*radius,y*radius,z*radius);
   x = cos(lon)*cos(lat+latinc);
   y = sin(lat+latinc);
   z = -sin(lon)*cos(lat+latinc);
   glNormal3f(x,y,z);
   glVertex3f(x*radius,y*radius,z*radius);
   lon += loninc;
  }
  glEnd();
  lat += latinc;
 }

 /* North-pole triangular fan */
 glBegin(GL_TRIANGLE_FAN);
 glNormal3f(0,1,0);
 glVertex3f(0,radius,0);
 y = sin(lat);
 lon = 0;
 for (i=0; i<=nlon; i++)
 {
  x = cos(lon)*cos(lat);
  z = -sin(lon)*cos(lat);
  glNormal3f(x,y,z);
  glVertex3f(x*radius,y*radius,z*radius);
  lon += loninc;
 }
 glEnd();

 glEndList();
}

//------------------------------ MakeStarsBG ------------------------------//
void MakeStarsBG(GLuint id)
{
 glNewList(id,GL_COMPILE);
  glColor3f(1,1,1);
  glPointSize(1.0);
  glBegin(GL_POINTS);
  for (int i=0; i<num_stars; i++)
  {
   glVertex3f(stars[i].r()[0],stars[i].r()[1],stars[i].r()[2]);
  }
  glEnd();
 glEndList();
}

//------------------------------ MakeAsteriods ------------------------------//
void MakeAsteroids(GLuint id)
{
 glNewList(id,GL_COMPILE);
	for (int i=0; i<4*num; i++){if (ast[i].IsAlive()){
		glPushMatrix();
		 glTranslatef(ast[i].r()[0], ast[i].r()[1],ast[i].r()[2]);
		 glRotatef(120*ast[i].R(), 0.0, 1.0, 0.0);
//		 glRotatef(100*t, 0.0, 0.0, 1.0);
 		glColor3d(0.5*ast[i].R(),0.5*ast[i].G(),0.5*ast[i].B()); 
   glCallList(planetID);
		glPopMatrix();
	}}
 glEndList();
}


//****************************************************************************//
//-------------------------------- Cluster -----------------------------------//
//****************************************************************************//
void MakeCluster(GLuint id)
{
 glNewList(id,GL_COMPILE);
 for (int i=0; i<num; i++){if (p[i].IsAlive())
 {
  glPushMatrix();
  glTranslatef(p[i].r()[0], p[i].r()[1],p[i].r()[2]);
  glRotatef(120*p[i].R(), 0.0, 1.0, 0.0);
//	 	glRotatef(100*t, 0.0, 0.0, 1.0);
  glScalef(p[i].size(),p[i].size(),p[i].size());
  glColor3d(p[i].R(),p[i].G(),p[i].B());
  glCallList(planetID);
  glPopMatrix();
 }}
 glEndList();
}
//----------------------------------------------------------------------------//

//--------------------------- RotateVector -----------------------------------//
// It is used to rotate the camera in the "MoveCamera" function
// Method taken from http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/
void RotateVector(float * vect, float * axis, float ang)
{
 float u=axis[0],v=axis[1],w=axis[2];
 float x=vect[0],y=vect[1],z=vect[2];
 float norm2=u*u+v*v+w*w;

 float rv[3];
 rv[0] = u*(u*x+v*y+w*z)+(x*(v*v+w*w)-u*(v*y+w*z))*cos(ang)+sqrt(u*u+v*v+w*w)*(-w*y+v*z)*sin(ang);
 rv[1] = v*(u*x+v*y+w*z)+(y*(u*u+w*w)-v*(u*x+w*z))*cos(ang)+sqrt(u*u+v*v+w*w)*(w*x-u*z)*sin(ang);
 rv[2] = w*(u*x+v*y+w*z)+(z*(u*u+v*v)-w*(u*x+v*y))*cos(ang)+sqrt(u*u+v*v+w*w)*(-v*x+u*y)*sin(ang);

 for (int q=0; q<3; ++q) vect[q]=rv[q]/norm2;

}
//------------------------------------- SetCamera ----------------------------//
// This function will be used in any function that involves a camera movement
void SetCamera(float dist)
{
 for (int q=0;q<3;++q) center[q] = 0+tx*v[q]+ty*up[q];//0.5*diag[q] + tx*v[q] + ty*up[q];
 for (int q=0;q<3;++q) eye[q] = center[q] + dist*camdir[q];
 far=4*fabs(dist);
 near=0.5*fabs(dist);
}

void SetCamera(float theta, float phi)
{
 // 'camdir' is the unitary vector that points from the center of the cell towards the camera
 camdir[0] = sin(theta)*cos(phi);
 camdir[1] = sin(theta)*sin(phi);
 camdir[2] = cos(theta);
 float phz = 0.5*M_PI;
 // 'v' is a 90° rotation of the vector 'camdir' (adding 0.5*M_PI to the "theta" angle in spherical coordinates) along the plane phi=const.
 v[0] = sin(theta+phz)*cos(phi);
 v[1] = sin(theta+phz)*sin(phi);
 v[2] = cos(theta+phz);
 // 'up'='camdir'x'v' vector indicates which direction is up (the direction from the bottom to the top of the viewing volume)
 up[0] = camdir[1]*v[2] - camdir[2]*v[1];
 up[1] = camdir[2]*v[0] - camdir[0]*v[2];
 up[2] = camdir[0]*v[1] - camdir[1]*v[0];

  // Set the far-clip (perspetive mode) always beyond the objects
 far=4*fabs(dist);
 near=0.5*fabs(dist);
 SetCamera(dist);
}

//------------------------------- MoveCamera ---------------------------------//
void MoveCamera(double angx, double angy)
{
 float rotx = M_PI*angx/180.0;
 float roty = M_PI*angy/180.0;

 float horiz[3];
 // 'horiz'='up'x'camdir'='v'
 horiz[0] = up[1]*camdir[2] - up[2]*camdir[1];
 horiz[1] = up[2]*camdir[0] - up[0]*camdir[2];
 horiz[2] = up[0]*camdir[1] - up[1]*camdir[0];

 RotateVector(camdir, horiz, roty);
 RotateVector(camdir, up, rotx);

 RotateVector(v, horiz, roty);
 RotateVector(v, up, rotx);
 
 // 'up'='camdir'x'v'
 up[0] = camdir[1]*v[2] - camdir[2]*v[1];
 up[1] = camdir[2]*v[0] - camdir[0]*v[2];
 up[2] = camdir[0]*v[1] - camdir[1]*v[0];

 SetCamera(dist);
}

//----------------------------- Init -----------------------------------------//
void Inicializa(void)
{
 GLfloat light_diffuse[]   = {1.0, 1.0, 1.0, 1.0};
 GLfloat light_position[] = { 30, 30, 30, 0};
 GLfloat light_position1[] = {0, 0, -1, 0};
 GLfloat light_position2[] = {1, 0, 0, 0};
 GLfloat light_position3[] = {-1, 0, 0, 0};
  
 GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
 GLfloat mat_shininess[] = { 50.0 };
 glShadeModel (GL_SMOOTH);
 
 /* Define normal light */
 glLightfv(GL_LIGHT0,GL_DIFFUSE,light_diffuse);
 glLightfv(GL_LIGHT1,GL_DIFFUSE,light_diffuse);
 glLightfv(GL_LIGHT2,GL_DIFFUSE,light_diffuse);
 glLightfv(GL_LIGHT3,GL_DIFFUSE,light_diffuse);
 glLightfv(GL_LIGHT0,GL_POSITION,light_position);
 glLightfv(GL_LIGHT1,GL_POSITION,light_position1);
 glLightfv(GL_LIGHT2,GL_POSITION,light_position2);
 glLightfv(GL_LIGHT3,GL_POSITION,light_position3);

 /* Define material of objects */
 glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
 glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);

 /* Enable a single OpenGL light */
 glEnable(GL_LIGHTING);
 glEnable(GL_LIGHT0);
 glEnable(GL_LIGHT1);
 glEnable(GL_LIGHT2);
 glEnable(GL_LIGHT3);

 /* Use depth buffering for hidden surface elimination */
 glEnable(GL_DEPTH_TEST);
 
 /* Enable the color material mode */
 glEnable(GL_COLOR_MATERIAL);

 /* Don't mess up the lighting equations */
 glEnable(GL_NORMALIZE);

 near = (GLdouble)( 0.5*(diagonal-0.5*diag[2]) );
 far  = (GLdouble)( 2.0*(diagonal+0.5*diag[2]) );

 double min=diagonal/sqrt(3);
 fovy = (GLfloat)( 0.5*alto/(2*near) );
 fovy = (GLfloat)( 2*atan((float)fovy)/M_PI*180.0 ); // Field of view touches exactly the borders of the cell
 fovy = (GLfloat)( 1.2*fovy); // We open the field of view a little bit more.
 
 for (int q=0; q<3; q++) diag[q]=0.5*(alto+ancho);
 double theta=0, phi=0;
 camdir[0] = sin(theta)*cos(phi);
 camdir[1] = sin(theta)*sin(phi);
 camdir[2] = cos(theta);
 float phz = 0.5*M_PI;
 /* Set 'v', a 90° rotation of the vector 'camdir' (adding 0.5*M_PI to the "theta" angle in spherical coordinates) along the plane phi=const. */
 v[0] = sin(theta+phz)*cos(phi); 
 v[1] = sin(theta+phz)*sin(phi);
 v[2] = cos(theta+phz);
 /* Set 'up'='camdir'x'v' vector, that indicates which direction is up (the direction from the bottom to the top of the viewing volume) */
 up[0] = camdir[1]*v[2] - camdir[2]*v[1];
 up[1] = camdir[2]*v[0] - camdir[0]*v[2];
 up[2] = camdir[0]*v[1] - camdir[1]*v[0];
 float cm = 0.0;
 for (int q=0;q<3;++q) cm += up[q]*up[q];
 for (int q=0;q<3;++q) up[q] /= sqrt(cm);
 

 /* SetCamera(dist) */
 MoveCamera(0,0);
 
 /* Make a planet */
 planetID = glGenLists(1);
 MakePlanet(planetID);

 /* Make the stars background */
 bgstarsID = glGenLists(1);
 MakeStarsBG(bgstarsID);
 
  /* Make the asteriords */
 astID = glGenLists(1);
 MakeAsteroids(astID);
  
 /* Make the whole cluster of planets */
 clusterID = glGenLists(1);
 MakeCluster(clusterID);
 
}


//----------------------------- Reshape ---------------------------------------//
void Reshape(int w, int h)
{
  win_width=w;
  win_height=h;
  // set the GL viewport to match the full size of the window
  glViewport(0, 0, (GLsizei)w, (GLsizei)h);
  float aspect = w/(float)h;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  // we math the aspect ratio (w/h) with the viewport ratio in both perspectives
  gluPerspective(fovy,aspect,near,far);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

//----------------------------- Mouse ---------------------------------------//
void Mouse(int button,int state,int x,int y){

 mousecoord[0] = x;
 mousecoord[1] = y;

 switch(button)
 {
  case GLUT_LEFT_BUTTON:
   Buttons[0] = ((GLUT_DOWN==state)?1:0);
   break;
  case GLUT_MIDDLE_BUTTON:
   Buttons[1] = ((GLUT_DOWN==state)?1:0);
   break;
  case GLUT_RIGHT_BUTTON:
   Buttons[2] = ((GLUT_DOWN==state)?1:0);
   break;
  default:
   break;		
 }
 mousecoord[0] = x;
 mousecoord[1] = y;
 
 glutPostRedisplay();
}


//---------------------------- Mouse Motion ---------------------------------//
void MouseMotion(int x,int y)
{

 // drag and drop
 int dx = x-mousecoord[0];
 int dy = y-mousecoord[1];
 mousecoord[0] = x;
 mousecoord[1] = y;
 
 if( Buttons[1] )            // Boton del centro (scroll)
 {
  dist -= (float) 1.1f * dx;
  MoveCamera(0,0);
 }
 else if( Buttons[0] )      // Boton izquierdo
 {
  if (sqrt(dx*dx+dy*dy) < 0.2) return;
  MoveCamera(-0.5*dx,-0.5*dy);
 }
 else if( Buttons[2] )     // Boton derecho
 {
  tx -= (float) 0.1f * dx;
  ty += (float) 0.1f * dy;
  MoveCamera(0,0);
 }
 Reshape(win_width,win_height);
 glutPostRedisplay();
}





//----------------------------- Keyboard -------------------------------------//
void Keyboard(unsigned char key, int x, int y)
{
 switch (key){
 case 27:
 case 'q': {exit(0); data.close();}
  break;
 case 'p':
 case ' ': paused = !paused;
  break;
 case 's':{wrt=!wrt;}
  break;
 case 'r':{dt=-dt;}
  break;
 case 'a':{ tx=0; ty=0; SetCamera(0,0);}
  break;
 }
}
void Up_Down(int key, int x, int y){
	switch(key) {
		case GLUT_KEY_UP:{
				p[0].set_size(p[0].size()+.1);
				if (p[0].size()>=30){p[0].set_red(1),p[0].set_green(0),p[0].set_blue(0);}
		}
		break;
		case GLUT_KEY_DOWN:{	p[0].set_size(p[0].size()-.1);}
		break;
		case GLUT_KEY_RIGHT: dt+=1e-3;   // dt in years
		break;
		case GLUT_KEY_LEFT: dt-=1e-3;    // dt in years
		break;
		default:break;
	}
	glutPostRedisplay();
}


void DrawStatistics(void){
	/**************** Writing in the window ****************/
	glPushMatrix();
	glDisable(GL_LIGHTING);
	glTranslatef(0,0,-dist);
	glColor3f (0.0, 1.0, 0.0);
// glScaled(dist,dist,dist);
	sprintf(PL,"%s %u","Planets: ",num-dead_ones);					glRasterPos3f (-ancho,0,-dist); 	drawString(PL);
	sprintf(AS,"%s %u","Asteroids: ",num_asteroids);				glRasterPos3f (-ancho/2, 2.8,0); 	drawString(AS);
	sprintf(DEAD,"%s %u","Dead Planets: ",dead_ones);				glRasterPos3f (-5.5, 2.6,0); 	drawString(DEAD);
	sprintf(TE,"%s %4.4f","Total Energy: ", Kin(p,ast)+Pot(p,ast));	glRasterPos3f (-5.5, 2.4,0); 	drawString(TE);
	sprintf(AM,"%s %4.4f","Angular Momentum: ",TL(p,ast).Mod());	glRasterPos3f (-5.5, 2.2,0); 	drawString(AM);
	sprintf(T,"%s %4.2f %s","Time: ",t," [years]");					glRasterPos3f (-5.5, 2.0,0); 	drawString(T);
 
	glEnable(GL_LIGHTING);
	glPopMatrix();
	/********************************************************/
}

void Dibuja(void)
{
// glClearColor(1,1,1,1);
 glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
 glLoadIdentity();

 glPushMatrix();
 gluLookAt(
 (GLfloat)eye[0],(GLfloat)eye[1],(GLfloat)eye[2],
 (GLfloat)center[0],(GLfloat)center[1],(GLfloat)center[2],
 (GLfloat)up[0],(GLfloat)up[1],(GLfloat)up[2]);
 
 
 //|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//
 //-------------------------------------------------- OBJECTS ------------------------------------------------------//
 //|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//
 //***********************************PLANE***************************************//
 glPushMatrix();
  glColor3f(0,1,0);
  glBegin(GL_LINES);
   for(float i=-ancho/2;i<=ancho/2;i+=10){	glVertex3f(i,-alto/2,0);	glVertex3f(i,alto/2,0);}
   for(float i=-alto/2; i<=alto/2; i+=10){	glVertex3f(ancho/2,i,0);	glVertex3f(-ancho/2,i,0);}
  glEnd();
 glPopMatrix();

	//*******************************CENTER OF MASS*********************************//
	glPushMatrix();
	glPointSize(5);
	glColor3f(1,0,0);
	glBegin(GL_POINTS);
		glVertex3f(COM(p,ast)[0], COM(p,ast)[1], COM(p,ast)[2]);
	glEnd();
	glPopMatrix();

 // STARS
 glCallList(bgstarsID);
 // PLANETS
 glCallList(clusterID);
 // ASTEROIDS
 glCallList(astID);
 // TRAJECTORIES
 glColor3f(p[1].R(),p[1].G(),p[1].B());
 glBegin(GL_LINES);
 for (int q=0;q<N;q++)
 {
  glVertex3f(trayec[q][0],trayec[q][1],trayec[q][2]);
 }
 glEnd();
 glColor3f(p[3].R(),p[3].G(),p[3].B());
 glBegin(GL_LINES);
 for (int q=0;q<N;q++)
 {
  glVertex3f(trayec2[q][0],trayec2[q][1],trayec2[q][2]);
 }
 glEnd();


 glPopMatrix();
	//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//
	//---------------------------------------------- END OF OBJECTS ---------------------------------------------------//
	//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//
//	if (wrt) DrawStatistics();
 glutPostRedisplay();
 glutSwapBuffers();
}
//------------------------------------------------------------------------------------------------//

void Idle(void)
{
//	Send_Data();
 if(!paused)
 {
  Solver();
  t+=dt;
  MakeCluster(clusterID);
  if (step%10==0) {trayec[N]=p[1].r(); trayec2[N]=p[3].r();N++;}
  step++;
 }
}

//************************************MAIN****************************************//
int main(int argc, char** argv)
{
 // Principal Window
 glutInit(&argc, argv);
 glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
 glutInitWindowSize(win_width,win_height);
 glutInitWindowPosition(0,0);
 glutCreateWindow("Kepler 3D");
 // Program's Functions:
 Boundary_Conditions();
 Inicializa();

 glutIdleFunc(Idle);
 glutDisplayFunc(Dibuja);
 glutKeyboardFunc(Keyboard);
 glutSpecialFunc(Up_Down);
 glutMouseFunc(Mouse);
 glutReshapeFunc(Reshape);
 glutMotionFunc(MouseMotion);

 glutMainLoop();
 return 0;
}
//********************************************************************************//
