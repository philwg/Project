#include <iostream>
#include <cmath>
#include <vector>
#include <GL/glut.h>

#define ROUNDZERO 0.0000000000001f
#define PT_SIZE 12.0f

using namespace std;

struct Point
{
	private:
		float x, y, z;
		float w;
	
	public:
		Point():Point(0.0f, 0.0f) {}
		Point(float init_x, float init_y):Point(init_x, init_y, 0.0f) {}
		Point(float init_x, float init_y, float init_z):Point(init_x, init_y, init_z, 0.0) {}
		Point(float init_x, float init_y, float init_z, float init_w) {
			setX(init_x);
			setY(init_y);
			setZ(init_z);
			setW(init_w);
		}
		
		~Point() {}
		
		void setX(float new_x) { x = new_x; }
		void setY(float new_y) { y = new_y; }
		void setZ(float new_z) { z = new_z; }
		void setW(float new_w) { w = new_w; }
		
		float getX() { return x; }
		float getY() { return y; }
		float getZ() { return z; }
		float getW() { return w; }
		
		void addX(float delta_X) { x+= delta_X; }
		void addY(float delta_Y) { y+= delta_Y; }
		void addZ(float delta_Z) { z+= delta_Z; }
		void addW(float delta_W) { w+= delta_W; }
		void add(Point P2) {
			addX(P2.getX());
			addY(P2.getY());
			addZ(P2.getZ());
			addW(P2.getW());
		}
		
		static Point milieu(Point P1, Point P2) {
			return Point((P1.getX()+P2.getX())/2.0f,
					  	 (P1.getY()+P2.getY())/2.0f,
					     (P1.getZ()+P2.getZ())/2.0f,
					     (P1.getW()+P2.getW()));
		}
		
		Point times(float k) { 
			Point P;
			P.setX(x*k);
			P.setY(y*k);
			P.setZ(z*k);
			P.setW(w);
			return P;
		}
		
		Point Weighted() { return this->times(this->getW()); }
		
		static Point weightedSum(vector<Point> Pts) {
			Point wS;
			for(vector<Point>::iterator Pt = Pts.begin(); Pt != Pts.end(); ++Pt) {
				wS.add(Pt->Weighted());
			}
			return wS;
		}
		
		static Point barycentre(vector<Point> Pts) {
			Point bary = Point::weightedSum(Pts);
			return bary.times(1.0/bary.getW());
		}
		
		float* to3DVertex() { return new float[3] {getX(), getY(), getZ()}; }
		float* to2DVertex() { return new float[2] {getX(), getY()}; }
};

struct BezierMatrix
{
	private:
	
		vector<int> Values;
		int Dim;
	
		int Pascal(int i, int n) {
			if ((i<0)||(i>n)) return 0.0f;
			else if ((i==0)||(i==n)) return 1.0f;
			else return Pascal(i-1, n-1) + Pascal(i, n-1);
		}
		
		void flush() {
			while (Values.size()>0) Values.pop_back();
			Dim = 0;
		}

	public:

		BezierMatrix():BezierMatrix(3) {}
		BezierMatrix(int n) { updateToDim(n); }
		
		~BezierMatrix() { flush(); }
		
		void updateToDim(int n) {
			flush();
			for(int j=0; j<n; j++) {
				for(int i=0; i<=j; i++) {
					Values.push_back(pow(-1,n-1-j)*Pascal(j,n-1)*Pascal(i,j));
				}
			}
			Dim = n;	
		}
		
		int getDim() { return Dim; }
		
		int getValue(int i, int j) {
			if ((i+j)>=Dim) return 0;
			else return Values[(((i+j)*(i+j+1))/2)+j];
		}
};


/* ------------------------------------------------------------------- */
/*				Les variables globales de la classe			   		   */
/* ------------------------------------------------------------------- */
const float CatmullRomMatrix[4][4] { 0.0f,  1.0f,  0.0f,  0.0f,
					  				-0.5f,  0.0f,  0.5f,  0.0f,
									 1.0f, -2.5f,  2.0f, -0.5f,
					   				-0.5f,  1.5f,  -1.5f, 0.5f};
BezierMatrix BM;
bool BMatrixOn {false};

int selected {-1};
int action 	  {0};

float pickingPrecision {10.0f};
float displayScale      {1.0f};

bool PointsDisplayed {false};
bool BezierDisplayed {false};
bool CatRomDisplayed {false};
bool SplineDisplayed {false};

vector<Point> ctrlPoints;

/*	Définition des variables concernant la fenêtre d'affichage	*/
int dispWdwHeight {1000};	
int dispWdwWidth  {1000};
int dispWdwXpos		{50};
int dispWdwYpos		{50};

/*	Definition des variables de glOrtho (projection orthogonale)	*/
int orthoHeight {dispWdwHeight};
int orthoWidth  {dispWdwWidth};
int zOrthoMin 	{-1};
int zOrthoMax 	{ 1};

/*	Couleurs de certains objets	*/
float bgColor[] { 0.12f, 0.12f, 0.12f, 1.0f };
float lnColor[] { 0.96f, 0.48f, 0.0f };
float ptColor[] { 0.96f, 0.0f,  0.0f };

float BcvColor[]  { 0.96f, 0.96f, 0.0f };
float ScvColor[]  { 0.96f, 0.0f,  0.96f};
float CRcvColor[] { 0.0f,  0.96f, 0.96f};

/*---------------------------------------------------------------------*/
/* 				Prototypes des fonctions & main()					   */
/*---------------------------------------------------------------------*/

int Pascal(int i, int n);
float Bernstein(int i, int n, float t);

Point Bezier(float t);
Point CatmullRomPW(float t, vector<Point> CR_4Pts4);

void renderBitmapString(float x, float y, void *font, const char *string);
void displayPoints();
void displayBezierCurve();
void displayCatRomCurve();
void displaySplineCurve();
void display();

void mouse(int button, int state, int x, int y);
void mousemotion(int x, int y);
void initMenu();
void menu(int item);

void zoomIO(float d);
void keyboard(unsigned char touche, int x, int y);

void idle();
void reshape(int x, int y);


int main(int argc, char** argv) //........................... Fonction de lancement de l'application
{
	/* Initialisation de glut et creation de la fenêtre */
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowPosition(dispWdwXpos, dispWdwYpos);
	glutInitWindowSize(dispWdwWidth, dispWdwHeight);
	glutCreateWindow("Courbes Parametriques et Splines");

	/* Initialisation des styles */
	glClearColor(bgColor[0], bgColor[1], bgColor[2], bgColor[3]);
	glPointSize(6.0f);
	glLineWidth(2.0f);
	glEnable(GL_DEPTH_TEST);
	initMenu();

	/* Enregistrement des fonctions de rappel */
	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutReshapeFunc(reshape);
	glutMouseFunc(mouse);
	glutMotionFunc(mousemotion);
	glutIdleFunc(idle);

	/* Entree dans la boucle principale glut */
	glutMainLoop();
	return 0;
}


/* ------------------------------------------------------------------- */
/*					Fonctions de calcul de courbes				       */
/* ------------------------------------------------------------------- */

int Pascal(int i, int n)
{
	if ((i<0)||(i>n)) return 0;
	else if ((i==0)||(i==n)) return 1;
	else return Pascal(i-1, n-1) + Pascal(i, n-1);
}

float Bernstein(int i, int n, float t)
{
	float value {0.0f};
	if ((t>=0.0f) && (t<=1.0f)) value = (float)Pascal(i, n)*pow(1.0f-t,i)*pow(t, n-i);
	if (value < ROUNDZERO) value = 0.0f;
	return value;
}

Point Bezier(float t) {
	int n = ctrlPoints.size();	
	float Bx {0.0f};
	float By {0.0f};
	float Bz {0.0f};
	for(int i=0; i<n; i++) {
		float coef_Bernstein = Bernstein(i, n-1, t);
		Bx += ctrlPoints[i].getX()*coef_Bernstein;
		By += ctrlPoints[i].getY()*coef_Bernstein;
		Bz += ctrlPoints[i].getZ()*coef_Bernstein;
	}
	Point BezierPoint {Bx, By, Bz};
	return BezierPoint;
}

Point BezierPW(float t, vector<Point> B_ctrls) {
	int d = B_ctrls.size();
	Point B_Pt;
	vector<Point> InterPt;
	for(int i=0; i<d; i++) {
		for(int j=0; j<d; j++) {
			B_ctrls[j].setW(BM.getValue(i,j));
		}
		B_Pt = Point::weightedSum(B_ctrls);
		B_Pt.setW(pow(t,d-1-i));
		InterPt.push_back(B_Pt);
	}
	return Point::weightedSum(InterPt);
}

Point CatmullRomPW(float t, vector<Point> CR_4Pts)
{
	Point CR_Pt;
	vector<Point> InterPt;
	for(int i=0; i<4; i++) {
		for(int j=0; j<4; j++) {
			CR_4Pts[j].setW(CatmullRomMatrix[i][j]);
		}
		CR_Pt = Point::weightedSum(CR_4Pts);
		CR_Pt.setW((double)pow((double)t,i));
		InterPt.push_back(CR_Pt);
	}
	return Point::weightedSum(InterPt);
}


/* ------------------------------------------------------------------- */
/*					Fonctions display d'OpenGL					       */
/* ------------------------------------------------------------------- */

void renderBitmapString(float x, float y, void *font, const char *string)
{
	const char *c;
	glRasterPos2f(x, y);
	for (c = string; *c != '\0'; c++) {
		glutBitmapCharacter(font, *c);
	}
}

void displayPoints()
{
	glPointSize(PT_SIZE);
	char writeBuffer[15] = { '\0' };
	int idPoint {0};
	Point prevPoint;
	for(vector<Point>::iterator ctrlPoint = ctrlPoints.begin(); ctrlPoint != ctrlPoints.end(); ++ctrlPoint) {
		sprintf(writeBuffer, "%d", idPoint++);
		glColor3fv(ptColor);
		renderBitmapString(	ctrlPoint->getX() - 3.0f,
							ctrlPoint->getY() + 10.0f,
							GLUT_BITMAP_HELVETICA_18,
							writeBuffer);
		glBegin(GL_POINTS);
			glVertex2fv(ctrlPoint->to2DVertex());
		glEnd();
		if (idPoint>1) {
			glColor3fv(lnColor);
			glPushAttrib(GL_ENABLE_BIT);
				glLineStipple(4, 0x8888);
				glEnable(GL_LINE_STIPPLE);
				glBegin(GL_LINES);
					glVertex2fv(prevPoint.to2DVertex());
					glVertex2fv(ctrlPoint->to2DVertex());
				glEnd();
			glPopAttrib();
		}
		prevPoint = Point {ctrlPoint->getX(), ctrlPoint->getY()};
	}
}

void displayBezierCurve()
{
	if (!BMatrixOn) {
		Point prevPoint = Bezier(0.0f);
		for (float t=0.01f; t<=1.0f; t+=0.01f) {
			Point crntPoint = Bezier(t);
			glColor3fv(BcvColor);
			glBegin(GL_LINES);
				glVertex2fv(prevPoint.to2DVertex());
				glVertex2fv(crntPoint.to2DVertex());
			glEnd();
			prevPoint = crntPoint;
		}
	}
	else {
		vector<Point> crntCtrls;
		int d = BM.getDim();
		int n = ctrlPoints.size();
		int pwNbr = (n-1)/(d-1);
		for(int i=0; i<pwNbr; i++) {
			for(int j=0; j<d; j++) {
				crntCtrls.push_back(ctrlPoints[(d-1)*i+j]);		
			}
			Point prevPoint = BezierPW(0.0f, crntCtrls);
			for (float t=0.01f; t<=1.0f; t+=0.01f) {
				Point crntPoint = BezierPW(t, crntCtrls);
				glColor3fv(BcvColor);
				glBegin(GL_LINES);
					glVertex2fv(prevPoint.to2DVertex());
					glVertex2fv(crntPoint.to2DVertex());
				glEnd();
				prevPoint = crntPoint;	
			}
			while (crntCtrls.size()>0) crntCtrls.pop_back();			
		}		
	}
}

void displayCatmullRomCurve()
{
	int n = ctrlPoints.size();
	vector<Point> crnt4Points;
	Point prevPoint;
	for(int i=0; i<n-3; i++) {
		for(int j=0; j<4; j++) {
			crnt4Points.push_back(ctrlPoints[i+j]);		
		}
		prevPoint = CatmullRomPW(0.0f, crnt4Points);
		for (float t=0.01f; t<=1.0f; t+=0.01f) {
			Point crntPoint = CatmullRomPW(t, crnt4Points);
			glColor3fv(CRcvColor);
			glBegin(GL_LINES);
				glVertex2fv(prevPoint.to2DVertex());
				glVertex2fv(crntPoint.to2DVertex());
			glEnd();
			prevPoint = crntPoint;	
		}
		while (crnt4Points.size()>0) crnt4Points.pop_back();
	}
}

void displaySplineCurve() {

}

void display()	//------------------------------------------- Gestion de l'affichage
{
	/* Effacement de l'image avec la couleur de fond */
	glClearColor(bgColor[0], bgColor[1], bgColor[2], bgColor[3]);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	/* Objets à afficher */
	glLoadIdentity();
	glPushMatrix();
		glScalef(displayScale,displayScale,1.0f);
		if (PointsDisplayed) displayPoints();
		if (BezierDisplayed) displayBezierCurve();
		if (CatRomDisplayed) displayCatmullRomCurve();
	glPopMatrix();
	
	/* Nettoyage */
	glFlush();
	glutSwapBuffers();
}


/* ------------------------------------------------------------------- */
/*					Fonctions de Gestion de la Souris			       */
/* ------------------------------------------------------------------- */

void mouse(int mButton, int mState, int x, int y)	//------- Gestion du clic de la souris
{
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport); 
	if (mButton == GLUT_LEFT_BUTTON) {
		switch (action) {

			case 11: 	//----------------------------------- On ajoute un Point
				if (mState == GLUT_DOWN) {
					float new_y = (float)(dispWdwHeight - y) * orthoHeight / dispWdwHeight;
					float new_x = (float)x * orthoWidth / dispWdwWidth;
					if (	(new_y >= 10.0f)
						 && (new_y <= (float)orthoHeight - 10.0f)
						 && (new_x >= 10.0f)
						 && (new_x <= (float)orthoWidth - 10.0f))
					{
						Point newPoint { new_x, new_y };
						ctrlPoints.push_back(newPoint);
						PointsDisplayed = true;
					}
				}
				break;

			case 12: 	//----------------------------------- On supprime un point
				if ((mState == GLUT_DOWN) && (ctrlPoints.size()>0)) {
					GLuint selectBuf[2 * ctrlPoints.size()]; 					
					GLuint *pointeurBufferselected, selectedVertexNum; 	
					GLint selectedObjectNb; 							
					glSelectBuffer(2 * ctrlPoints.size(), selectBuf);			
					glRenderMode(GL_SELECT); 							
					glPushMatrix();
						glMatrixMode(GL_PROJECTION);
						glLoadIdentity();
						gluPickMatrix(x, viewport[3] - y, pickingPrecision, pickingPrecision, viewport); 
						glOrtho(0.0f, orthoWidth, 0.0f, orthoHeight, zOrthoMin, zOrthoMax);
						glInitNames();
						glPushName(2);
						glColor3fv(ptColor);
						for (int i=0; i < ctrlPoints.size(); i++) {
							glLoadName(i);							
							glBegin(GL_POINTS);
								glVertex2fv(ctrlPoints.at(i).to2DVertex());
							glEnd();
						}
					glPopMatrix();
					glFlush();
					selectedObjectNb = glRenderMode(GL_RENDER); 
					if (selectedObjectNb) {
						pointeurBufferselected = (GLuint *)selectBuf; 
						pointeurBufferselected += 3; 						// On se place la où le nom est écrit
						selectedVertexNum = *pointeurBufferselected; 		// On récupère le numero du sommet touché
						ctrlPoints.erase(ctrlPoints.begin() + selectedVertexNum);
					}
				}
				break;
				
			case 13: 	//----------------------------------- On insère un point de contrôle
				if ((mState == GLUT_DOWN) && (ctrlPoints.size()>0)) {
					GLuint selectBuf[2 * ctrlPoints.size()]; 					
					GLuint *pointeurBufferselected, selectedVertexNum; 	
					GLint selectedObjectNb; 							
					glSelectBuffer(2 * ctrlPoints.size(), selectBuf);			
					glRenderMode(GL_SELECT); 							
					glPushMatrix();
						glMatrixMode(GL_PROJECTION);
						glLoadIdentity();
						gluPickMatrix(x, viewport[3] - y, pickingPrecision, pickingPrecision, viewport); 
						glOrtho(0.0f, orthoWidth, 0.0f, orthoHeight, zOrthoMin, zOrthoMax);
						glInitNames();
						glPushName(2);
						glColor3fv(ptColor);
						for (int i=0; i < ctrlPoints.size(); i++) {
							glLoadName(i);							
							glBegin(GL_POINTS);
								glVertex2fv(ctrlPoints.at(i).to2DVertex());
							glEnd();
						}
					glPopMatrix();
					glFlush();
					selectedObjectNb = glRenderMode(GL_RENDER); 
					if (selectedObjectNb) {
						pointeurBufferselected = (GLuint *)selectBuf; 
						pointeurBufferselected += 3;
						selectedVertexNum = *pointeurBufferselected;
						vector<Point>::iterator Pts_it = ctrlPoints.begin();
						Point newPoint = Point::milieu(ctrlPoints.at(selectedVertexNum),
													   ctrlPoints.at(selectedVertexNum+1));
						ctrlPoints.insert(Pts_it+selectedVertexNum+1, 1, newPoint);
					}
				}
				break;
				
			case 14:	//----------------------------------- On déplace un point
				if (mState==GLUT_DOWN) {
					GLuint selectBuf[2 * ctrlPoints.size()];
					GLuint *pointeurBufferselected;
					GLint selectedObjectNb;
					glSelectBuffer(2 * ctrlPoints.size(), selectBuf);
					glRenderMode(GL_SELECT);
					glPushMatrix();
						glMatrixMode(GL_PROJECTION);
						glLoadIdentity();
						gluPickMatrix(x, viewport[3] - y, pickingPrecision, pickingPrecision, viewport);
						glOrtho(0.0, orthoWidth, 0.0, orthoHeight, zOrthoMin, zOrthoMax);
						glInitNames();
						glPushName(2);
						glColor3fv(ptColor);
						for (int i=0; i < ctrlPoints.size(); i++) {
							glLoadName(i); 
							glBegin(GL_POINTS);
								glVertex2fv(ctrlPoints.at(i).to2DVertex());
							glEnd();
						}
					glPopMatrix();
					glFlush();
					selectedObjectNb = glRenderMode(GL_RENDER);
					if (selectedObjectNb) {
						pointeurBufferselected = (GLuint *)selectBuf;
						pointeurBufferselected += 3;
						selected = *pointeurBufferselected; 
					}
				}
				break;
							
			default :
				break;
		}
		if (mState==GLUT_UP) selected = -1;
		reshape(viewport[2], viewport[3]);
		glutPostRedisplay();
	}
}

void mousemotion(int x, int y)	//--------------------------- Gestion du mouvement de la souris
{
	if (selected != -1) {
	
		float new_x = (float)x * orthoWidth / dispWdwWidth;
		float new_y = (float)(dispWdwHeight - y) * orthoHeight / dispWdwHeight;	
		
		if (action==14) { 
			
			if (new_y < 10.0f) ctrlPoints.at(selected).setY(10.0f);
			else if (new_y > orthoHeight - 10.0f) ctrlPoints.at(selected).setY(orthoHeight - 10.0f);
			else ctrlPoints.at(selected).setY(new_y);
				
			if (new_x < 10.0f) ctrlPoints.at(selected).setX(10.0f);
			else if (new_x > orthoWidth - 10.0f) ctrlPoints.at(selected).setX(orthoWidth - 10.0f);
			else ctrlPoints.at(selected).setX(new_x);
		}
	}
	glutPostRedisplay();
}

void initMenu()	//------------------------------------------- Initialisation du menu flottant (bouton droit)
{
	int sm_points = glutCreateMenu(menu);
	glutAddMenuEntry("   Add       ", 11);
	glutAddMenuEntry("   Delete    ", 12);
	glutAddMenuEntry("   Insert    ", 13);
	glutAddMenuEntry("   Move      ", 14);
	
	int sm_B_matrix = glutCreateMenu(menu);
	glutAddMenuEntry("   Increase   (D)", 23);
	glutAddMenuEntry("   Decrease  (d)", 24);

	int main_menu = glutCreateMenu(menu);
	glutAddSubMenu	("   CONTROLS      ", sm_points);
	glutAddMenuEntry("   Display/Hide  ", 15);
	glutAddMenuEntry("  ____________   ",  0);
	glutAddMenuEntry("   BEZIER        ",  0);
	glutAddMenuEntry("   Display/Hide  ", 21);
	glutAddMenuEntry("   Rec./Matrix  (m)", 22);
	glutAddSubMenu  ("   PW degree     ", sm_B_matrix);
	glutAddMenuEntry("  ____________   ",  0);
	glutAddMenuEntry("   CATMULL-ROM   ",  0);
	glutAddMenuEntry("   Display/Hide  ", 31);
	glutAddMenuEntry("  ____________   ",  0);
	glutAddMenuEntry("   Cancel        ",  4);
	glutAddMenuEntry("   Exit          ",  5);
	glutAttachMenu(GLUT_RIGHT_BUTTON);
}

void menu(int item)	//--------------------------------------- Gestion des actions du menu flottant
{
	int B_dim = BM.getDim();
	switch (item) {

		case 11:			//------------------------------- Ajouter un point
			action = item;
			glutSetCursor(GLUT_CURSOR_CROSSHAIR);
			break;
		
		case 12:			//------------------------------- Supprimer un point
			action = item;
			glutSetCursor(GLUT_CURSOR_DESTROY);	
			break;

		case 13:			//------------------------------- Insérer un point
			action = item;
			glutSetCursor(GLUT_CURSOR_HELP);	
			break;
			
		case 14:			//------------------------------- Déplacer un point
			action = item;
			glutSetCursor(GLUT_CURSOR_INFO);
			break;
							
		case 15:			//------------------------------- Afficher/Cacher les points de contrôle
			if (PointsDisplayed) PointsDisplayed = false;
			else PointsDisplayed = true;
			break;
			
		case 21:			//------------------------------- Afficher/Cacher la courbe de Bezier
			if (BezierDisplayed) BezierDisplayed = false;
			else BezierDisplayed = true;
			break;

		case 22:			//------------------------------- mode recursif ou matrice par morceaux pour Bezier
			if (BMatrixOn) BMatrixOn = false;
			else BMatrixOn = true;
			break;

		case 23:			//------------------------------- Augmentation de la taille des morceaux
			if ((B_dim<ctrlPoints.size())&&(B_dim<14)) BM.updateToDim(B_dim+1);
			break;
			
		case 24:			//------------------------------- Diminution de la taille des morceaux
			if (B_dim>2) BM.updateToDim(B_dim-1);
			break;
		
		case 31:			//------------------------------- Afficher la courbe de Catmull-Rom
			if (CatRomDisplayed) CatRomDisplayed = false;
			else CatRomDisplayed = true;
			break;
			
		case 4:				//------------------------------- Désactiver l'action en cours
			action = 0;
			selected = -1;
			glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
		break;
		
		case 5:				//------------------------------- Quitter
			exit(0);
			break;
			
		default:
			break;
	}
	glutPostRedisplay();
}


/* ------------------------------------------------------------------- */
/*				Fonction de Gestion des saisies clavier		       	   */
/* ------------------------------------------------------------------- */

void zoomIO(float d)	//----------------------------------- Mise à jour des paramètres de zoom
{
	displayScale *= d;
}

void keyboard(unsigned char touche, int x, int y)	//------- Gestion des saisies clavier
{
	int B_dim = BM.getDim();
	switch(touche) {

		/* Quitter */
		case 'q': 
			exit(0);
			
		/* Zoomer */
		case 'Z': 
			zoomIO(1.1f);
			break;
			
		/* Dézoomer */
		case 'z': 
			zoomIO(1/1.1f);
			break;
		
		/* Activer et désactiver la version matricielle des Beziers */	
		case 'm':
			if (BMatrixOn) BMatrixOn = false;
			else BMatrixOn = true;
			break;
		
		case 'd':
			if (B_dim>2) BM.updateToDim(B_dim-1);
			break;
		
		case 'D':
			if ((B_dim<ctrlPoints.size())&&(B_dim<14)) BM.updateToDim(B_dim+1);
			break;
		
		default :
			break;
	}
	glutPostRedisplay();
}


/* ------------------------------------------------------------------- */
/*							Fonction IDLE						       */
/* ------------------------------------------------------------------- */

void idle()
{
	/*if (true) {	
		glutPostRedisplay();
	}*/
}


/* ------------------------------------------------------------------- */
/*					Fonction de Redimensionnement				       */
/* ------------------------------------------------------------------- */

void reshape(int width, int height)
{
	glViewport(0, 0, width, height); 			
	dispWdwWidth = width;
	dispWdwHeight = height;
	glMatrixMode(GL_PROJECTION); 
	glLoadIdentity(); 
	glOrtho(0.0f, orthoWidth, 0.0f, orthoHeight, zOrthoMin, zOrthoMax); 
	glMatrixMode(GL_MODELVIEW); 
	glLoadIdentity();
}
