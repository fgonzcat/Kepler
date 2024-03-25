/////////////////////////////////
// KEPLER.CC: PLANETARY SYSTEM //
/////////////////////////////////
// Creado: 7-Enero-2008.
// Ultima modificacion: 27-Julio-2016.
// Autor: Felipe Gonzalez.

#include <GL/glut.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "class-ptcl.h"
#include "class-vect.h"
#include "class-azar.h"

using namespace std;

const double G=0.00011859645;       // in AU³/(earthmass*year²)
const double earthmass=1.0;         // sunmass=332946.05 earthmass
const double sunmass=332946.05;    // sunmass=332946.05 earthmass
const double R0=0.5;               // Reference radius (in AU)
double dt=0.01;
//---------------------------WINDOW-----------------------------//
double ancho=1300.0;
double alto=700.0;
double xmax=ancho/2, ymax=alto/2, xmin=-ancho/2, ymin=-alto/2;
double fac=0.02;
//--------------------------------------------------------------//

//-----------------------PARTICLE SETUP-------------------------//
const int num=8;               // Number of particles
const int num_stars=5000;      // Number of stars on the background
GLint star_display_list;
GLint all_stars_display_list;
int dead_ones=0;				// Planets destroyed
bool not_counted[num];
double t=0;
//double UA=1.495E11;
//double mt=5.98E24, ratio_t=6378000;						// Earth: mass[kg], radio[m]
//double sun=331.1*mt, radio_sun=109.125*ratio_t;
//double merc=0.06*mt, radio_merc=;
//double venus=0.82*mt;
//double mars=0.11*mt;
//double jupiter=318*mt;
//double saturn=95.1*mt;
//double saturn=14.6*mt;
//double uranus=17.2*mt;
//double pluto=0.002*mt;
Ptcl p[num];
Ptcl stars[num_stars];
double x[num][2];
double y[num][2];
double z[num][2];
double fx[num][2];
double fy[num][2];
double fz[num][2];
double k1x[num][2];
double k2x[num][2];
double k3x[num][2];
double k4x[num][2];
double k1y[num][2];
double k2y[num][2];
double k3y[num][2];
double k4y[num][2];
double k1z[num][2];
double k2z[num][2];
double k3z[num][2];
double k4z[num][2];

//--------------------------------------------------------------//

//--------------------------------------INTERNAL FUNCTIONS---------------------------------------//
void circulo(double radio, double x, double y)
{
 int puntos=20;
 glBegin(GL_POLYGON);
 for(int i=0;i<puntos;i++)
  glVertex3f(x+radio*cos((2*M_PI*i)/puntos),y+radio*sin((2*M_PI*i)/puntos),0);
 glEnd();
}

void star(double radio)
{
  glPushMatrix();
 for(int i=0; i<6; i++)
 {
  glRotatef(60.0,  0.0, 0.0, 1.0);
  glBegin(GL_POLYGON);
   glVertex3f(-0.5*radio, -0.5*radio , 0.0);
   glVertex3f( 0.5*radio, -0.5*radio , 0.0);
   glVertex3f(0.0*radio, (0.5*sqrt(3)*radio+0.5*radio) , 0.0);
  glEnd();
 }
  glPopMatrix();
}


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


// Potential Energy of the Whole System
double Pot(Ptcl set_of_ptcls[num]){
	double pot=0;
	for(int i=0; i<num; i++){
		for(int j=0; j<num; j++){
			if(i!=j){	pot+=Potential(set_of_ptcls[i],set_of_ptcls[j]);}
		}
	}
	return pot/2.0;
}

// Kinetic Energy of the Whole System
double Kin(Ptcl set_of_ptcls[num]){
	double kin=0;
	for(int i=0; i<num; i++){
		kin+=set_of_ptcls[i].Kinetics();
	}	
	return kin;
}

Vector COM(Ptcl set_of_ptcls[num]){
	double M=0;
	Vector centre;
	for (int i=0; i<num; i++){	M+=set_of_ptcls[i].mass();}
	for (int i=0; i<num; i++){	centre=centre+(set_of_ptcls[i].mass()*set_of_ptcls[i].r())/M;}
	return centre;
}

// Angular momentum of the Whole System
Vector TL(Ptcl set_of_ptcls[num]){
	Vector l(0,0,0,0);
	for(int i=0; i<num; i++){
		l=l+set_of_ptcls[i].L();
	}	
	return l;
}
/********************************** Gravitational GForce between particles (Newton) ***********************************/
Vector GForce(Ptcl m1, Ptcl m2)
{
 double d=Dist(m1,m2);
 double fact=-G*m1.mass()*m2.mass()/(d*d*d);
 Vector f1=fact*(m1.r()-m2.r());
 if (m1==m2){f1=Vector(0,0,0);}
 return f1;
}
/*********************************************************************************************************************/

//-------------------------------------------------------------------------------------------//

//--------------------------------------GLUT FUNCTIONS---------------------------------------//
void Inicializa(void)
{
glClear(GL_COLOR_BUFFER_BIT);	// Limpia la pantalla (2D)
// glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);	// (3D)
glMatrixMode(GL_PROJECTION);
glLoadIdentity();
glOrtho(-fac*ancho, fac*ancho, -fac*alto, fac*alto, -10,10);
//glEnable(GL_DEPTH_TEST);									// (3D)
}

GLint MakeStar(double radius)
{
 GLint starDL;

 // Create the id for the list
 starDL = glGenLists(1);

 // start list
 glNewList(starDL,GL_COMPILE);
  star(radius);
 glEndList();

 return(starDL);
}

GLint MakeStars()
{
 GLint all_starsDL;

 // Create the id for the list
 all_starsDL = glGenLists(1);

 // start list
 glNewList(all_starsDL,GL_COMPILE);
 for (int i=0; i<num_stars; i++)
 {
  glPushMatrix();
  glTranslatef( stars[i].r()[0], stars[i].r()[1], 0.0);
  glCallList(star_display_list);
  glPopMatrix();
 }
 glEndList();

 return(all_starsDL);
}




void Boundary_Conditions(void)
{
 randomize();
 // PLANET'S MASSES
 for (int i=0; i<num; i++){p[i].set_mass(earthmass);}
 p[0].set_mass(sunmass);
 p[3].set_mass(10*earthmass);
 p[5].set_mass(5*earthmass);
 p[6].set_mass(300*earthmass);

 // PLANET'S RADIUS
 for (int i=0; i<num; i++){	p[i].set_size(0.4*R0);}
 p[0].set_size(R0);
 p[4].set_size(0.6*R0);
 p[5].set_size(0.2*R0);
 p[6].set_size(0.8*R0);

 // INITIAL POSITIONS AND VELOCITIES
 p[0].set_r(0,0,0);      p[0].set_v(0,0,0);
 p[1].set_r(1,0,0);      p[1].set_v(0,2*M_PI,0);
 p[2].set_r(-1,-3,0);    p[2].set_v(3,0,0);
 p[3].set_r(3,0,0);      p[3].set_v(1,3,0);
 p[4].set_r(5,0,0);      p[4].set_v(0,3,0);
 p[5].set_r(-5,2,0);      p[5].set_v(1,1,0);
 p[6].set_r(-7,0,0);     p[6].set_v(1,-2,0);
 p[7].set_r(7,4,0);    p[7].set_v(0,1,0);


 // COLORS
 for (int i=0; i<num; i++){ p[i].set_red(0); p[i].set_green(0.8); p[i].set_blue(1);} // Default=blue
 p[0].set_red(1); p[0].set_green(1); p[0].set_blue(0);
 p[1].set_red(1); p[1].set_green(0); p[1].set_blue(0);
 p[2].set_red(1); p[2].set_green(0.5); p[2].set_blue(0);
 p[3].set_red(1); p[3].set_green(0.9); p[3].set_blue(0.5);
 p[4].set_red(.3); p[4].set_green(0.3); p[4].set_blue(.8);
 p[6].set_red(.3); p[6].set_green(1); p[6].set_blue(0);


 /**************************** STARS POSITIONS **************************************/
 for (int i=0; i<num_stars; i++)
 {
  stars[i].set_r(dazar(-200*R0,200*R0), dazar(-200*R0,200*R0),0);
 }
 /***********************************************************************************/
 /**************************Initialize x[][], y[][] & z[][]**************************/
 for (int i=0; i<num; i++){
 x[i][0]=p[i].r()[0];		y[i][0]=p[i].r()[1];		z[i][0]=p[i].r()[2];
 x[i][1]=p[i].v()[0];		y[i][1]=p[i].v()[1];		z[i][1]=p[i].v()[2];
 not_counted[i]=true;
 }

 star_display_list = MakeStar(0.1*R0);
 all_stars_display_list = MakeStars();
}

void Solver(void)
{
 //----------------------Fixing f-------------------//
 for (int i=0; i<num; i++){

 fx[i][0]=x[i][1];
 fy[i][0]=y[i][1];
 fz[i][0]=z[i][1];
 fx[i][1]=0;	fy[i][1]=0;	fz[i][1]=0;
 for (int j=0; j<num; j++){
 fx[i][1]+=GForce(p[i],p[j])[0]/p[i].mass();
 fy[i][1]+=GForce(p[i],p[j])[1]/p[i].mass();
 fz[i][1]+=GForce(p[i],p[j])[2]/p[i].mass();
 }
 }	
 //------------------------------------------------//

 //--------------------Fixing k1-------------------//
 // k1xij=dt*fxij(x11,x12,...,x21,...)
 for (int i=0; i<num; i++){
 k1x[i][0]=dt*fx[i][0];
 k1y[i][0]=dt*fy[i][0];
 k1z[i][0]=dt*fz[i][0];

 k1x[i][1]=dt*fx[i][1];
 k1y[i][1]=dt*fy[i][1];
 k1z[i][1]=dt*fz[i][1];
 }
 //------------------------------------------------//

 //--------------------Fixing k2-------------------//
 // k2xij=dt*fxij(x11+k1x11/2,x12+k1x12/2,...,x21+k1x21/2,...)
 //**************Next is used to evaluate fxij in xij+k1xij/2************//
 for (int i=0; i<num; i++){
 x[i][0]+=k1x[i][0]/2;	y[i][0]+=k1y[i][0]/2;	z[i][0]+=k1z[i][0]/2;
 x[i][1]+=k1x[i][1]/2;	y[i][1]+=k1y[i][1]/2;	z[i][1]+=k1z[i][1]/2;
 p[i].set_r(x[i][0], y[i][0], z[i][0]);
 p[i].set_v(x[i][1], y[i][1], z[i][1]);
 }
 //**********************************************************************//
 //*********Fixing***********//
 for (int i=0; i<num; i++){
 k2x[i][0]=dt*x[i][1];
 k2y[i][0]=dt*y[i][1];
 k2z[i][0]=dt*z[i][1];

 k2x[i][1]=0; k2y[i][1]=0; k2z[i][1]=0;
 for (int j=0; j<num; j++){
 k2x[i][1]+=dt*GForce(p[i],p[j])[0]/p[i].mass();
 k2y[i][1]+=dt*GForce(p[i],p[j])[1]/p[i].mass();
 k2z[i][1]+=dt*GForce(p[i],p[j])[2]/p[i].mass();
 }
 }
 //**************************//
 //***********Next is used to undo the evaluation of fxij above***********//
 for (int i=0; i<num; i++){
 x[i][0]-=k1x[i][0]/2;	y[i][0]-=k1y[i][0]/2;	z[i][0]-=k1z[i][0]/2;
 x[i][1]-=k1x[i][1]/2;	y[i][1]-=k1y[i][1]/2;	z[i][1]-=k1z[i][1]/2;
 p[i].set_r(x[i][0], y[i][0], z[i][0]);
 p[i].set_v(x[i][1], y[i][1], z[i][1]);
 }
 //**********************************************************************//
 //------------------------------------------------//

 //--------------------Fixing k3-------------------//
 // k3xij=dt*fxij(x11+k2x11/2,x12+k2x12/2,...,x21+k2x21/2,...)
 //**************Next is used to evaluate fxij in xij+k2xij/2************//
 for (int i=0; i<num; i++){
 x[i][0]+=k2x[i][0]/2;	y[i][0]+=k2y[i][0]/2;	z[i][0]+=k2z[i][0]/2;
 x[i][1]+=k2x[i][1]/2;	y[i][1]+=k2y[i][1]/2;	z[i][1]+=k2z[i][1]/2;
 p[i].set_r(x[i][0], y[i][0], z[i][0]);
 p[i].set_v(x[i][1], y[i][1], z[i][1]);
 }
 //**********************************************************************//
 //*********Fixing***********//
 for (int i=0; i<num; i++){
 k3x[i][0]=dt*x[i][1];
 k3y[i][0]=dt*y[i][1];
 k3z[i][0]=dt*z[i][1];

 k3x[i][1]=0; k3y[i][1]=0; k3z[i][1]=0;
 for (int j=0; j<num; j++){
 k3x[i][1]+=dt*GForce(p[i],p[j])[0]/p[i].mass();
 k3y[i][1]+=dt*GForce(p[i],p[j])[1]/p[i].mass();
 k3z[i][1]+=dt*GForce(p[i],p[j])[2]/p[i].mass();
 }
 }
 //**************************//
 //***********Next is used to undo the evaluation of fxij above***********//
 for (int i=0; i<num; i++){
 x[i][0]-=k2x[i][0]/2;	y[i][0]-=k2y[i][0]/2;	z[i][0]-=k2z[i][0]/2;
 x[i][1]-=k2x[i][1]/2;	y[i][1]-=k2y[i][1]/2;	z[i][1]-=k2z[i][1]/2;
 p[i].set_r(x[i][0], y[i][0], z[i][0]);
 p[i].set_v(x[i][1], y[i][1], z[i][1]);
 }
 //**********************************************************************//
 //------------------------------------------------//

 //--------------------Fixing k4-------------------//
 // k4xij=fxij(x11+k3x11,x12+k3x12,...,x21+k3x21,...)
 //**************Next is used to evaluate fxij in xij+k2xij/2************//
 for (int i=0; i<num; i++){
 x[i][0]+=k3x[i][0];	y[i][0]+=k3y[i][0];	z[i][0]+=k3z[i][0];
 x[i][1]+=k3x[i][1];	y[i][1]+=k3y[i][1];	z[i][1]+=k3z[i][1];
 p[i].set_r(x[i][0], y[i][0], z[i][0]);
 p[i].set_v(x[i][1], y[i][1], z[i][1]);
 }
 //**********************************************************************//
 //*********Fixing***********//
 for (int i=0; i<num; i++){
 k4x[i][0]=dt*x[i][1];
 k4y[i][0]=dt*y[i][1];
 k4z[i][0]=dt*z[i][1];

 k4x[i][1]=0; k4y[i][1]=0; k4z[i][1]=0;
 for (int j=0; j<num; j++){
 k4x[i][1]+=dt*GForce(p[i],p[j])[0]/p[i].mass();
 k4y[i][1]+=dt*GForce(p[i],p[j])[1]/p[i].mass();
 k4z[i][1]+=dt*GForce(p[i],p[j])[2]/p[i].mass();
 }
 }
 //**************************//
 //***********Next is used to undo the evaluation of fxij above***********//
 for (int i=0; i<num; i++){
 x[i][0]-=k3x[i][0];	y[i][0]-=k3y[i][0];	z[i][0]-=k3z[i][0];
 x[i][1]-=k3x[i][1];	y[i][1]-=k3y[i][1];	z[i][1]-=k3z[i][1];
 p[i].set_r(x[i][0], y[i][0], z[i][0]);
 p[i].set_v(x[i][1], y[i][1], z[i][1]);
 }
 //**********************************************************************//

 //------------------------------------------------//

 //------------Computing the Solutions------------//
 for (int i=0; i<num; i++){
 x[i][0]+=(k1x[i][0]+2.0*k2x[i][0]+2.0*k3x[i][0]+k4x[i][0])/6.0;
 x[i][1]+=(k1x[i][1]+2.0*k2x[i][1]+2.0*k3x[i][1]+k4x[i][1])/6.0;
 y[i][0]+=(k1y[i][0]+2.0*k2y[i][0]+2.0*k3y[i][0]+k4y[i][0])/6.0;
 y[i][1]+=(k1y[i][1]+2.0*k2y[i][1]+2.0*k3y[i][1]+k4y[i][1])/6.0;
 z[i][0]+=(k1z[i][0]+2.0*k2z[i][0]+2.0*k3z[i][0]+k4z[i][0])/6.0;
 z[i][1]+=(k1z[i][1]+2.0*k2z[i][1]+2.0*k3z[i][1]+k4z[i][1])/6.0;
 }	
 //------------------------------------------------//
 //----------New Positions and Velocities----------//
 for (int i=0; i<num; i++){
 p[i].set_r(x[i][0], y[i][0], z[i][0]);
 p[i].set_v(x[i][1], y[i][1], z[i][1]);
 }
 //------------------------------------------------//
 }

 void Tiempo(void){ t+= dt;	glutPostRedisplay();}

 void reshape (int w, int h){
 if (!h)
 return;
 glViewport(0, 0, w, h);
 glMatrixMode(GL_PROJECTION);
 glLoadIdentity();
 ancho=w; alto=h;
 glOrtho(-fac*ancho, fac*ancho, -fac*alto, fac*alto, -10,10);
 glMatrixMode(GL_MODELVIEW);
 glLoadIdentity();
 }

 void Up_Down(int key, int x, int y){
 switch(key) {
 case GLUT_KEY_UP:{
 p[0].set_size(p[0].size()+.1);
 if (p[0].size()>=10*R0){p[0].set_red(1),p[0].set_green(0),p[0].set_blue(0);}
 }
 break;
 case GLUT_KEY_DOWN:{	p[0].set_size(p[0].size()-.1);}
 break;
 case GLUT_KEY_RIGHT: dt+=1e-3/2;
 break;
 case GLUT_KEY_LEFT: dt-=1e-3/2;
 break;
 default:break;
 }
 glutPostRedisplay();
 }

void Keyboard(unsigned char key, int x, int y){
 switch (key)
 {
  case 27: exit(0);
  break;
  case 'q': exit(0);
  break;
  case 'p':{dt=0;}
  break;
  case 'm':{dt=1e-2;}
  break;
  case 'z':{fac-=0.0005;reshape(ancho,alto);}
  break;
  case 'Z':{fac+=0.0005;reshape(ancho,alto);}
  break;
  case 'r':{dt=-dt;}
 }
}

/************************ Renderizar texto *************************/
void drawString(const char *s)
{
 unsigned int i;
 for (i = 0; i < strlen(s); i++){	glutBitmapCharacter (GLUT_BITMAP_HELVETICA_12, s[i]);}
}


void drawStringBig(char *s)
{
 unsigned int i;
 for (i = 0; i < strlen(s); i++){	glutBitmapCharacter (GLUT_BITMAP_HELVETICA_18, s[i]);}
}
/*******************************************************************/


void Dibuja(void){
	glClear(GL_COLOR_BUFFER_BIT);
	Solver();
	/***********************STARS*************************/
	glColor3f(1,1,1);
	glCallList(all_stars_display_list);

/*	for (int i=0; i<num_stars; i++)
	{
	   glPushMatrix();
	   glTranslatef( stars[i].r()[0], stars[i].r()[1], 0.0);
	   glCallList(star_display_list);
	   glPopMatrix();
 	}*/

	/*****************************************************/
	for (int i=0; i<num; i++){
/*		for (int j=0; j<num;j++){
			if (i!=j && destroyed(p[i],p[j])){
				if (p[i].mass()==p[j].mass()){p[i].erase(); p[j].erase(); }
				else if (p[i].mass()>p[j].mass()){p[j].erase();}
				else {p[i].erase();}
			}
		}*/
		
//		glColor3d(cos(i+1),cos(i+.1),sin(i));
		glColor3d(p[i].R(),p[i].G(),p[i].B());
		//**********************Planets************************//
		circulo(p[i].size(),p[i].r()[0],p[i].r()[1]);
		//*****************************************************//
	}
	//**********************************Counting Dead Planets***********************************************//
//	for (int i=0; i<num; i++){if (p[i].confirm_dead() && not_counted[i]){dead_ones++; not_counted[i]=false;}}
	//******************************************************************************************************//
	/****************Center of Mass***************/
	glPointSize(5);
	glBegin(GL_POINTS);
		glColor3f(1,0,0);
		glVertex3f(COM(p)[0],COM(p)[1],COM(p)[2]);
	glEnd();
	/********************************************/
	/************* Writing in the window **************/
	glColor3f (0.0, 1.0, 0.0);
	const char *label1="Dead Planets";
	const char *label2="Total Energy";
	const char *label3="Angular Momentum";
	const char *label4="Time";
	char cad1[15], cad2[15], cad3[15], cad4[15];
	 glPushMatrix();
	sprintf(cad1,"%s %u",": ",dead_ones);	
	glRasterPos3f (-0.9*fac*ancho,       0.9*fac*alto,0); 	drawString(label1);
	glRasterPos3f (-0.9*fac*(ancho-300), 0.9*fac*alto,0); 	drawString(cad1);

	sprintf(cad2,"%s %4.4f (ME AU/year^2)",": ",Kin(p)+Pot(p));
	glRasterPos3f (-0.9*fac*ancho,       0.85*fac*alto,0); 	drawString(label2);
	glRasterPos3f (-0.9*fac*(ancho-300), 0.85*fac*alto,0); 	drawString(cad2);

	sprintf(cad3,"%s %4.4f",": ",TL(p).Mod());
	glRasterPos3f (-0.9*fac*ancho,       0.80*fac*alto,0); 	drawString(label3);
	glRasterPos3f (-0.9*fac*(ancho-300), 0.80*fac*alto,0); 	drawString(cad3);

	sprintf(cad4,"%s %4.2f years",": ",t);
	glRasterPos3f (-0.9*fac*ancho,       0.75*fac*alto,0); 	drawString(label4);
	glRasterPos3f (-0.9*fac*(ancho-300), 0.75*fac*alto,0); 	drawString(cad4);
	 glPopMatrix();
	/************************************/

	glutSwapBuffers();
}
//------------------------------------------------------------------------------------------------//

//************************************MAIN****************************************//
int main(int argc, char** argv){
	glutInit(&argc, argv);
//	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);	//(3D)
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);						//(2D)
	glutInitWindowSize(int(ancho),int(alto));
	glutInitWindowPosition(0,0);
	glutCreateWindow("Kepler 2D");
	// Funciones del Programa:
	Inicializa();
	Boundary_Conditions();
	glutIdleFunc(Tiempo);
	glutDisplayFunc(Dibuja);
	glutKeyboardFunc(Keyboard);
	glutSpecialFunc(Up_Down);
	glutReshapeFunc(reshape);

	glutMainLoop();
	return 0;
}
//********************************************************************************//
