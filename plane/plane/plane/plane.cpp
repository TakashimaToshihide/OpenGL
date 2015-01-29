#pragma comment(linker, "/SUBSYSTEM:WINDOWS /ENTRY:mainCRTStartup")
#include<stdlib.h>
#include<gl/glut.h>
#include<math.h>
#include<iostream>
#include<fstream>
#include<string.h>
using namespace std;

#define pi 3.1415926535

char str[256];
double Amp[10000][3] = {};

double f = 1; //�ŗL�U�����̐ݒ�
double tscale = 200; //���ԃX�P�[���̐ݒ� ���t���[���łP��U������悤�ɂ��邩
double vmax = 0.2; //�U���̃X�P�[���𒲐��i����܂�傫���Ɣj�]����悤�ȁj


//�@���x�N�g���̎Z�o
int CalculateNormal(double p1[3], double p2[3],	double p3[3], float n[3])
{
	float v1[3];
	float v2[3];
	float cross[3];
	float length;
	
	int i;
		
	/* v1 = p1 - p2�����߂� */
	for (i = 0; i < 3; i++) {
		v1[i] = p1[i] - p2[i];
	}

	/* v2 = p3 - p2�����߂� */
	for (i = 0; i < 3; i++) {
		v2[i] = p3[i] - p2[i];
	}

	/* �O��v2�~v1�i= cross�j�����߂� */

	for (i = 0; i < 3; i++) {
		cross[i] = v2[(i+1)%3] * v1[(i+2)%3] - v2[(i+2)%3] * v1[(i+1)%3];
	}

	/* �O��v2�~v1�̒���|v2�~v1|�i= length�j�����߂� */
	length = sqrtf(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);

	/* ����|v2�~v1|��0�̂Ƃ��͖@���x�N�g���͋��߂��Ȃ� */
	if (length == 0.0f) {
		return 0;
	}

	/* �O��v2�~v1�𒷂�|v2�~v1|�Ŋ����Ė@���x�N�g��n�����߂� */
	for (i = 0; i < 3; i++) {
		n[i] = cross[i] / length;
	}
	return 1;
}


/* �e�_�̐U���̍ő�l��ǂݍ��ނ̂Ń��[�h�֐��Ōv�Z����K�v�͂Ȃ��C������
double mode(double x){
	double lamda = 0;
	double phai = 0;
	
	lamda = pi/2;

	phai = vmax*x;

	return phai;
}
*/

double frequency(double amp, double t){
	
	double v;
	//v = mode(x)*(a*cos(2*pi*t/tscale)+b*cos(2*pi*t/tscale));
	v = vmax*amp*sin(2*pi*t/tscale);

	return v;
}


//���C�g�̈ʒu

GLfloat light0pos[] = { 0.0, 1.0, -3.0, 1.0 };
GLfloat light1pos[] = { 3.0, 3.0, 3.0, 1.0 };
//GLfloat green[] = { 0.8, 0.8, 0.8, 1.0 };


//���[�h�`��̃f�[�^��ǂݍ���
void Load(double Amp[][3]){

	int i=0;
	
	ifstream fin("shellcs002.dat");
	
	while((fin >> Amp[i][0] >> Amp[i][1] >> Amp[i][2]) != 0){
		i++;
	}

	return;
}

//�`��֐�
static void namisen(double t)
{
	int m,n;
	int w,h;
	int j,i;
	double posx=0;
	double posy=0;
	int R[10000]={};
	int face[10000][3] = {};
	float A[3] = {};
	float B[3] = {};
	double Check[100] ={};
	double color_value;
	double color_ave;
	GLdouble normal[1000][3] = {}; //�@���x�N�g��

	//���b�V����؂���
	//���̕�����
	w = 21;
	//�c�̕�����
	h = 21;

	GLdouble vertex[1000][3] = {};

	for(i=0;i<w*h;i++){
		vertex[i][0] = Amp[i][0];
		vertex[i][1] = Amp[i][1];
		vertex[i][2] = frequency(Amp[i][2],t);
	}

	
	//�O�p�`�|���S���̏ꍇ

	for(m=0; m<h*(w-1); m++){

		//�@���x�N�g�����i�[����ϐ���������
		for(n=0;n<3;n++){
			A[n]=0;
			B[n]=0;
		}

		//���܂肩��[�����̕s�K�v�ȃ|���S�������O
		R[m] = (m+1)%w;

		//�|���S���̒��_���蓖�āA�@���x�N�g���̊��蓖��
		if(R[m] !=0){
			face[2*m][0] = m;
			face[2*m][1] = m+1;
			face[2*m][2] = m+w;

			CalculateNormal(vertex[m], vertex[m+1], vertex[m+w], A);
			normal[2*m][0] = A[0];
			normal[2*m][1] = A[1];
			normal[2*m][2] = A[2];

			face[2*m+1][0] = m+1;
			face[2*m+1][1] = m+w+1;
			face[2*m+1][2] = m+w;

			CalculateNormal(vertex[m+1], vertex[m+w+1], vertex[m+w], B);
			normal[2*m+1][0] = B[0];
			normal[2*m+1][1] = B[1];
			normal[2*m+1][2] = B[2];

		}

	}

	//�|���S���̕`��

	glBegin(GL_TRIANGLES);
	for(j=0;j<h*(w-1);j++){

		/*

		color_ave = ((vertex[face[2*j][0]][2]+vertex[face[2*j][1]][2]+vertex[face[2*j][2]][2]+vertex[face[2*j+1][1]][2])/4/0.2+1)/2;
		color_value = (-cos(4*pi*color_ave)+1)/2;

		if(0<=color_ave && color_ave<0.25){
			glColor3d(0.0,color_value,1.0);
		}

		else if(0.25<=color_ave && color_ave<0.5){
			glColor3d(0.0,1.0,color_value);
		}
		else if(0.5<=color_ave && color_ave<0.75){
			glColor3d(color_value,1.0,0.0);
		}
		else{
			glColor3d(1.0,color_value,0.0);
		}

		*/

		//���܂肩��s�K�v�ȃ|���S�����X�L�b�v
		R[j] = (j+1)%w;
		if(R[j]==0)
			continue;

		else
		glNormal3dv(normal[2*j]);
		for(i=3;--i>=0;){
			color_ave = ((vertex[face[2*j][i]][2])/0.2+1)/2;
			color_value = (-cos(4*pi*color_ave)+1)/2;

				if(0<=color_ave && color_ave<0.25){
					glColor3d(0.0,color_value,1.0);
				}

				else if(0.25<=color_ave && color_ave<0.5){
					glColor3d(0.0,1.0,color_value);
				}
				else if(0.5<=color_ave && color_ave<0.75){
					glColor3d(color_value,1.0,0.0);
				}
				else{
					glColor3d(1.0,color_value,0.0);
				}


			glVertex3dv(vertex[face[2*j][i]]);
		}


		glNormal3dv(normal[2*j+1]);
		for(i=3;--i>=0;){
			color_ave = ((vertex[face[2*j+1][i]][2])/0.2+1)/2;
			color_value = (-cos(4*pi*color_ave)+1)/2;

				if(0<=color_ave && color_ave<0.25){
					glColor3d(0.0,color_value,1.0);
				}

				else if(0.25<=color_ave && color_ave<0.5){
					glColor3d(0.0,1.0,color_value);
				}
				else if(0.5<=color_ave && color_ave<0.75){
					glColor3d(color_value,1.0,0.0);
				}
				else{
					glColor3d(1.0,color_value,0.0);
				}

			glVertex3dv(vertex[face[2*j+1][i]]);
		}

	}


	glEnd();
}

void idle(void)
{
	glutPostRedisplay();
}

void display(void)
{
	static int t = 0;
	double T = 0;
	int i = 0;

	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClear(GL_COLOR_BUFFER_BIT);
	glLoadIdentity();//�������H
	gluLookAt(3.0, 4.0, 5.0, 0.5, 0.5, 0.0, 0.0, 1.0, 0.0);//�J�����̈ʒu��ς���
	 /* �����̈ʒu�ݒ� */
	
	//glLightfv(GL_LIGHT0, GL_POSITION, light0pos);
	//glLightfv(GL_LIGHT1, GL_POSITION, light1pos);//���C�g�̈ʒu
	//glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, green);//�}�`�̐�
	


	glRotated(90.0, 0.0, 0.0, 1.0);
	glTranslated(0.0, -1.0, 0.0);

	glPushMatrix();
	glRotated(180.0, 1.0, 0.0, 0.0);
	glTranslated(0.0, -1.0, 1.0);
	namisen(double(t));
	glPopMatrix();

	glPushMatrix();
	glRotated(-90, 1.0, 0.0, 0.0);
	glTranslated(0.0, 0.0, 1.0);
	namisen(double(t));
	glPopMatrix();

	glPushMatrix();
	glRotated(90.0, 1.0, 0.0, 0.0);
	glTranslated(0.0, -1.0, 0.0);
	namisen(double(t));
	glPopMatrix();
	
	glPushMatrix();
	namisen(double(t));
	glPopMatrix();

	//�^�C���J�E���^�̕\��

	glColor3d(1.0,1.0,1.0);

	T = t/(f*tscale);

	sprintf(str,"t=%f(s)",T);
	glRasterPos2f(-1.2,1.0);
	glColor3d(1.0,1.0,1.0);
	for(i=0; i<strlen(str);i++){
		glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24,str[i]);
	}

	glFlush();

	if(++t >= 100000) t = 0;
}


void resize(int w, int h)
{
	glViewport(0, 0, w, h);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(30.0, (double)w/(double)h, 1.0, 100.0);

	glMatrixMode(GL_MODELVIEW);
}


void mouse(int button, int state, int x, int y)
{
  switch (button) {
  case GLUT_LEFT_BUTTON:
    if (state == GLUT_DOWN) {
      /* �A�j���[�V�����J�n */
      glutIdleFunc(idle);
    }
    else {
      /* �A�j���[�V������~ */
      glutIdleFunc(0);
    }
    break;
  case GLUT_RIGHT_BUTTON:
    if (state == GLUT_DOWN) {
      /* �R�}���� (1�X�e�b�v�����i�߂�) */
        glutPostRedisplay();
    }
    break;
  default:
    break;
  }
}
  
void keyboard(unsigned char key, int x, int y)
{
  switch (key) {
  case 'q':
  case 'Q':
  case '\033':  /* '\033' �� ESC �� ASCII �R�[�h */
    exit(0);
  default:
    break;
  }
}

void init(void)
{
	glClearColor(0.0, 0.0, 0.0, 0.0);
	//glEnable(GL_DEPTH_TEST);
	
	//���ĂȂ��n�̂��߃J�����O�͖��������Ă���
	/*
	glEnable(GL_CULL_FACE);
	glCullFace(GL_FRONT);
	*/

	
	//glEnable(GL_LIGHTING);
	//glEnable(GL_LIGHT0);
	//glEnable(GL_LIGHT1);
	//glLightfv(GL_LIGHT1, GL_DIFFUSE, green);
	//glLightfv(GL_LIGHT1, GL_SPECULAR, green);
	
	
	
}

int main(int argc, char *argv[])
{

	Load(Amp);	//dat�t�@�C���̓ǂݍ���
	glutInit(&argc, argv);
	//glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitDisplayMode(GLUT_RGBA);
	glutCreateWindow(argv[0]);
	glutDisplayFunc(display);
	glutReshapeFunc(resize);
	glutMouseFunc(mouse);
	glutKeyboardFunc(keyboard);
	init();
	glutMainLoop();
	return 0;
}