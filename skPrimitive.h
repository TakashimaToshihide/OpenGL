//skPrimitive.h

//---------------------------------------------------------------------
void calcNormal(float *p1,float *p2,float *p3,float *nn)
{	float len;
	nn[0] = (p2[1]-p1[1])*(p3[2]-p2[2])-(p2[2]-p1[2])*(p3[1]-p2[1]);
	nn[1] = (p2[2]-p1[2])*(p3[0]-p2[0])-(p2[0]-p1[0])*(p3[2]-p2[2]);
	nn[2] = (p2[0]-p1[0])*(p3[1]-p2[1])-(p2[1]-p1[1])*(p3[0]-p2[0]);
	len = (float)sqrt((double)(nn[0]*nn[0]+nn[1]*nn[1]+nn[2]*nn[2]));

	if(len != 0.0)
	{	nn[0] /= len; nn[1] /= len; nn[2] /= len; }
}
//-----------------------------------------------------------------------------
void skWireCube(void)
{
	//���A�����P�̗�����(�����̂̒��S�����_)
	float p[8][3]={{0.5f,0.5f,0.5f},{-0.5f,0.5f,0.5f},{-0.5f,-0.5f,0.5f},
	{0.5f,-0.5f,0.5f},{0.5f,0.5f,-0.5f},{-0.5f,0.5f,-0.5f},{-0.5f,-0.5f,0-0.5f},
	{0.5f,-0.5f,-0.5f}};

	glBegin(GL_QUADS);
		//z����
		glVertex3fv(p[0]); glVertex3fv(p[1]);
		glVertex3fv(p[2]); glVertex3fv(p[3]);
		//x����
		glVertex3fv(p[0]); glVertex3fv(p[3]);
		glVertex3fv(p[7]); glVertex3fv(p[4]);
		//y����
		glVertex3fv(p[0]); glVertex3fv(p[4]);
		glVertex3fv(p[5]); glVertex3fv(p[1]);
		//-x����
		glVertex3fv(p[1]); glVertex3fv(p[5]);
		glVertex3fv(p[6]); glVertex3fv(p[2]);
		//-y����
		glVertex3fv(p[2]); glVertex3fv(p[6]);
		glVertex3fv(p[7]); glVertex3fv(p[3]);
		//-z����
		glVertex3fv(p[4]); glVertex3fv(p[7]);
		glVertex3fv(p[6]); glVertex3fv(p[5]);
	glEnd();
}
//---------------------------------------------------------------------------
void skSolidCube(void)
{   //���A�����P�̗�����
	//�e�_�̍��W(�����̂̒��S�����_)
	float p[8][3]={{0.5f,0.5f,0.5f},{-0.5f,0.5f,0.5f},{-0.5f,-0.5f,0.5f},
	{0.5f,-0.5f,0.5f},{0.5f,0.5f,-0.5f},{-0.5f,0.5f,-0.5f},{-0.5f,-0.5f,0-0.5f},{0.5f,-0.5f,-0.5f}};

	glBegin(GL_QUADS);
		glNormal3f(0.0f,0.0f,1.0f); //z����
		glVertex3fv(p[0]); glVertex3fv(p[1]);
		glVertex3fv(p[2]); glVertex3fv(p[3]);

		glNormal3f(1.0f,0.0f,0.0f); //x����
		glVertex3fv(p[0]); glVertex3fv(p[3]);
		glVertex3fv(p[7]); glVertex3fv(p[4]);

		glNormal3f(0.0f,1.0f,0.0f); //y����
		glVertex3fv(p[0]); glVertex3fv(p[4]);
		glVertex3fv(p[5]); glVertex3fv(p[1]);

		glNormal3f(-1.0f,0.0f,0.0f); //-x����
		glVertex3fv(p[1]); glVertex3fv(p[5]);
		glVertex3fv(p[6]); glVertex3fv(p[2]);

		glNormal3f(0.0f,-1.0f,0.0f); //-y����
		glVertex3fv(p[2]); glVertex3fv(p[6]);
		glVertex3fv(p[7]); glVertex3fv(p[3]);

		glNormal3f(0.0f,0.0f,-1.0f); //-z����
		glVertex3fv(p[4]); glVertex3fv(p[7]);
		glVertex3fv(p[6]); glVertex3fv(p[5]);
	glEnd();
}

//--------------------------------------------------------------------------
void skWireSphere(int Nxy)//�P�ʒ��a�̋�
{
	int i,j,Nz;
//	int np[21][21];//np[17][17];//i,j�_�̒��_�ԍ� (i:x-y�f�ʕ����_  j:�������������_)
	double phai; //x-y���ʂɑ΂���Ίp�i�ܓx�j
	double theta; //x���ɑ΂���Ίp�i�o�x)
	float p[21][21][3];//p[289][3]; //�e�_�̍��W
	float p1[3],p2[3],p3[3],p4[3];

	//if(Nxy > 16) Nxy = 16;
	if(Nxy > 20) Nxy = 20;
	Nz = Nxy;
	//���_�ԍ��A���W
	for(i=0;i<=Nxy;i++)
	{	theta = 2.0*M_PI*(double)i/(double)Nxy;
		for(j=0;j<=Nz;j++)
		{	
			phai = M_PI/2.0-M_PI*(double)j/(double)Nz;
			p[i][j][0] = (float)(0.5*cos(theta)*cos(phai)); //x���W
			p[i][j][1] = (float)(0.5*sin(theta)*cos(phai)); //y
			p[i][j][2] = (float)(0.5*sin(phai));            //z
		}
	}

	//���_����`
	for(i=0;i<Nxy;i++)
		for(j=0;j<Nz-1;j++)
		{ //�ʂ̒��_
			p1[0] = p[i][j][0]; //x���W
			p1[1] = p[i][j][1]; //y���W
			p1[2] = p[i][j][2]; //z���W
			p2[0] = p[i][j+1][0];
			p2[1] = p[i][j+1][1];
			p2[2] = p[i][j+1][2];
			p3[0] = p[i+1][j+1][0];
			p3[1] = p[i+1][j+1][1];
			p3[2] = p[i+1][j+1][2];
			p4[0] = p[i+1][j][0];
			p4[1] = p[i+1][j][1];
			p4[2] = p[i+1][j][2];

			glBegin(GL_QUADS);
				glVertex3fv(p1); glVertex3fv(p2);
				glVertex3fv(p3); glVertex3fv(p4);
			glEnd();
		}
	//��
	j = Nz-1;
	for(i=0;i<Nxy;i++)
		{ //�ʂ̒��_
			p1[0] = p[i][j][0];
			p1[1] = p[i][j][1];
			p1[2] = p[i][j][2];
			p2[0] = p[i][j+1][0];
			p2[1] = p[i][j+1][1];
			p2[2] = p[i][j+1][2];
			p3[0] = p[i+1][j][0];
			p3[1] = p[i+1][j][1];
			p3[2] = p[i+1][j][2];

			glBegin(GL_TRIANGLES);
				glVertex3fv(p1); glVertex3fv(p2);
				glVertex3fv(p3);
			glEnd();
		}
}

//--------------------------------------------------------------------------
void skSolidSphere(int Nxy)//�P�ʒ��a�̋�
{  
	int i,j,Nz;
//	int np[21][21];//i,j�_�̒��_�ԍ� (i:x-y�f�ʕ����_  j:�������������_)
	double phai; //x-y���ʂɑ΂���Ίp�i�ܓx�j
	double theta; //x���ɑ΂���Ίp�i�o�x)
	float p[21][21][3]; //�e�_�̍��W
	float p1[3],p2[3],p3[3],p4[3];

	if(Nxy > 21) Nxy = 20;
	Nz = Nxy;
	//���_�ԍ�,���W
	for(i=0;i<=Nxy;i++)
	{	theta = 2.0*M_PI*(double)i/(double)Nxy;
		for(j=0;j<=Nz;j++)
		{
			phai = M_PI/2.0-M_PI*(double)j/(double)Nz;
			p[i][j][0] = (float)(0.5*cos(theta)*cos(phai)); //x���W
			p[i][j][1] = (float)(0.5*sin(theta)*cos(phai)); //y
			p[i][j][2] = (float)(0.5*sin(phai));            //z
		}
	}
    
	//�e�p�b�`��`��
	for(i=0;i<Nxy;i++)
		for(j=0;j<Nz-1;j++)
		{ //�ʂ̒��_   x���W                       y���W                       z���W
			p1[0] = p[i][j][0]    ; p1[1] = p[i][j][1]    ; p1[2] = p[i][j][2];
			p2[0] = p[i][j+1][0]  ; p2[1] = p[i][j+1][1]  ; p2[2] = p[i][j+1][2];
			p3[0] = p[i+1][j+1][0]; p3[1] = p[i+1][j+1][1]; p3[2] = p[i+1][j+1][2];
			p4[0] = p[i+1][j][0]  ; p4[1] = p[i+1][j][1]  ; p4[2] = p[i+1][j][2];

			glBegin(GL_QUADS);
				glNormal3fv(p1);glVertex3fv(p1); 
				glNormal3fv(p2);glVertex3fv(p2);
				glNormal3fv(p3);glVertex3fv(p3); 
				glNormal3fv(p4);glVertex3fv(p4);
			glEnd();
		}	
	//��
	j = Nz-1;
	for(i=0;i<Nxy;i++)
		{ //�ʂ̒��_   x���W                       y���W                       z���W
			p1[0] = p[i][j][0]   ; p1[1] = p[i][j][1]    ; p1[2] = p[i][j][2];
			p2[0] = p[i][j+1][0] ; p2[1] = p[i][j+1][1]  ; p2[2] = p[i][j+1][2];
			p3[0] = p[i+1][j][0] ; p3[1] = p[i+1][j][1]  ; p3[2] = p[i+1][j][2];

			
			glBegin(GL_TRIANGLES);
				glNormal3fv(p1);glVertex3fv(p1); 
				glNormal3fv(p2);glVertex3fv(p2);
				glNormal3fv(p3);glVertex3fv(p3); 
			glEnd();
		}	
}
//--------------------------------------------------------------------------
void skSolidSemiSphere(int Nxy)//�P�ʒ��a�̔���
{  
	int i,j,Nz;
	double phai; //x-y���ʂɑ΂���Ίp�i�ܓx�j
	double theta; //x���ɑ΂���Ίp�i�o�x)
	float p[21][11][3]; //�e�_�̍��W
	float p1[3],p2[3],p3[3],p4[3];

	if(Nxy > 21) Nxy = 20;
	Nz = Nxy/2;
	//���_�ԍ��A���W
	for(i=0;i<=Nxy;i++)
	{	theta = 2.0*M_PI*(double)i/(double)Nxy;
		for(j=0;j<=Nz;j++)
		{	
			phai = M_PI*(1.0-(double)j/(double)Nz)/2.0;
			p[i][j][0] = (float)(0.5*cos(theta)*cos(phai)); //x���W
			p[i][j][1] = (float)(0.5*sin(theta)*cos(phai)); //y
			p[i][j][2] = (float)(0.5*sin(phai));            //z
		}
	}
    
	//�e�p�b�`�̕`��
	for(i=0;i<Nxy;i++)
		for(j=0;j<Nz;j++)
		{ //�ʂ̒��_   x���W                       y���W                       z���W
			p1[0] = p[i][j][0]    ; p1[1] = p[i][j][1]    ; p1[2] = p[i][j][2];
			p2[0] = p[i][j+1][0]  ; p2[1] = p[i][j+1][1]  ; p2[2] = p[i][j+1][2];
			p3[0] = p[i+1][j+1][0]; p3[1] = p[i+1][j+1][1]; p3[2] = p[i+1][j+1][2];
			p4[0] = p[i+1][j][0]  ; p4[1] = p[i+1][j][1]  ; p4[2] = p[i+1][j][2];

			glBegin(GL_QUADS);
				glNormal3fv(p1);glVertex3fv(p1); 
				glNormal3fv(p2);glVertex3fv(p2);
				glNormal3fv(p3);glVertex3fv(p3); 
				glNormal3fv(p4);glVertex3fv(p4);
			glEnd();
		}	
	//��
	j = Nz;
	glBegin(GL_POLYGON);
		glNormal3f(0.0f,0.0f,-1.0f);//-z����
		for(i=Nxy;i>0;i--) glVertex3fv(p[i][j]);
	glEnd();
}
//--------------------------------------------------------------------------
void skWireSemiSphere(int Nxy)//�P�ʒ��a�̔���
{  
	int i,j,Nz;
	double phai; //x-y���ʂɑ΂���Ίp�i�ܓx�j
	double theta; //x���ɑ΂���Ίp�i�o�x)
	float p[21][11][3]; //�e�_�̍��W
	float p1[3],p2[3],p3[3],p4[3];

	if(Nxy > 21) Nxy = 20;
	Nz = Nxy/2;
	//���_�ԍ��A���W
	for(i=0;i<=Nxy;i++)
	{	theta = 2.0*M_PI*(double)i/(double)Nxy;
		for(j=0;j<=Nz;j++)
		{	
			phai = M_PI*(1.0-(double)j/(double)Nz)/2.0;
			p[i][j][0] = (float)(0.5*cos(theta)*cos(phai)); //x���W
			p[i][j][1] = (float)(0.5*sin(theta)*cos(phai)); //y
			p[i][j][2] = (float)(0.5*sin(phai));            //z
		}
	}
    
	//�e�p�b�`�̕`��
	for(i=0;i<Nxy;i++)
		for(j=0;j<Nz;j++)
		{ //�ʂ̒��_   x���W                       y���W                       z���W
			p1[0] = p[i][j][0]    ; p1[1] = p[i][j][1]    ; p1[2] = p[i][j][2];
			p2[0] = p[i][j+1][0]  ; p2[1] = p[i][j+1][1]  ; p2[2] = p[i][j+1][2];
			p3[0] = p[i+1][j+1][0]; p3[1] = p[i+1][j+1][1]; p3[2] = p[i+1][j+1][2];
			p4[0] = p[i+1][j][0]  ; p4[1] = p[i+1][j][1]  ; p4[2] = p[i+1][j][2];

			glBegin(GL_QUADS);
				glVertex3fv(p1); 
				glVertex3fv(p2);
				glVertex3fv(p3); 
				glVertex3fv(p4);
			glEnd();
		}	
	//��
	j = Nz;
	glBegin(GL_POLYGON);
		for(i=Nxy;i>0;i--) glVertex3fv(p[i][j]);
	glEnd();
}

//----------------------------------------------------------------------
void skSolidPrism(int Nxy)//�P�ʒ��a�A�P�ʍ����̑��p��
{	
	//Nxy--������
	float p[41][3]; //���_���W
	float nn[3];
	double theta0,theta;
	int i,ii;

	if(Nxy>21) { Nxy = 20;}
	theta0 = 2*M_PI/(double)Nxy;
	for(i=0;i<Nxy;i++)
	{   theta = theta0*(double)i;
		p[i][0] = (float)(0.5*cos(theta)); //����x����
		p[i][1] = (float)(0.5*sin(theta)); //������
		p[i][2] = 0.5f;                //������(����)
		p[i+Nxy][0] = (float)(0.5*cos(theta)); //�����x����
		p[i+Nxy][1] = (float)(0.5*sin(theta)); //������
		p[i+Nxy][2] = -0.5f;               //������
	}

	//���
	glBegin(GL_POLYGON);
		glNormal3f(0.0f,0.0f,1.0f);
		for(i=0;i<Nxy;i++) glVertex3fv(p[i]); 
	glEnd();
	//����
	glBegin(GL_POLYGON);
		glNormal3f(0.0f,0.0f,-1.0f); 
		for(i=2*Nxy-1;i>=Nxy;i--) glVertex3fv(p[i]); 
	glEnd();

	//����
	for(i=0;i<Nxy;i++)
	{
		ii = i+1;
		if(ii == Nxy) ii = 0;
		glBegin(GL_QUADS);
			calcNormal(p[i], p[i+Nxy], p[ii], nn);
			glNormal3f(nn[0],nn[1],nn[2]);
			glVertex3fv(p[i]);      glVertex3fv(p[i+Nxy]);
			glVertex3fv(p[ii+Nxy]); glVertex3fv(p[ii]);
		glEnd();
	}
}
//----------------------------------------------------------------------
void skWirePrism(int Nxy)//�P�ʒ��a�A�P�ʍ����̑��p��
{	
	//Nxy--������
	float p[41][3]; //���_���W
	double theta0,theta;
	int i,ii;

	if(Nxy>21) { Nxy = 20;}
	theta0 = 2*M_PI/(double)Nxy;
	for(i=0;i<Nxy;i++)
	{   theta = theta0*(double)i;
		p[i][0] = (float)(0.5*cos(theta)); //����x����
		p[i][1] = (float)(0.5*sin(theta)); //������
		p[i][2] = 0.5f;                //������(����)
		p[i+Nxy][0] = (float)(0.5*cos(theta)); //�����x����
		p[i+Nxy][1] = (float)(0.5*sin(theta)); //������
		p[i+Nxy][2] = -0.5f;               //������
	}

	//���
	glBegin(GL_POLYGON);
		for(i=0;i<Nxy;i++) glVertex3fv(p[i]);
	glEnd();
	//����
	glBegin(GL_POLYGON);
		for(i=2*Nxy-1;i>=Nxy;i--) glVertex3fv(p[i]);
	glEnd();

	//����
	for(i=0;i<Nxy;i++)
	{
		ii = i+1;
		if(ii == Nxy) ii = 0;
		glBegin(GL_QUADS);
			glVertex3fv(p[i]);      glVertex3fv(p[i+Nxy]);
			glVertex3fv(p[ii+Nxy]); glVertex3fv(p[ii]);
		glEnd();
	}
}
//----------------------------------------------------------------------
void skSolidCylinder(int Nxy)//�P�ʒ��a�A�P�ʍ���
{
	//Nxy--������
	float p[41][3]; //���_���W
	double theta0,theta;
	int i,ii;

	if(Nxy>21) { Nxy = 20;}
	theta0 = 2*M_PI/(double)Nxy;
	for(i=0;i<Nxy;i++)
	{   theta = theta0*(double)i;
		p[i][0] = (float)(0.5*cos(theta)); //����x����
		p[i][1] = (float)(0.5*sin(theta)); //������
		p[i][2] = 0.5f;                //������(����)
		p[i+Nxy][0] = (float)(0.5*cos(theta)); //�����x����
		p[i+Nxy][1] = (float)(0.5*sin(theta)); //������
		p[i+Nxy][2] = -0.5f;               //������
	}

	//���
	glBegin(GL_POLYGON);
		glNormal3f(0.0f,0.0f,1.0f);
		for(i=0;i<Nxy;i++) glVertex3fv(p[i]); 
	glEnd();
	//����
	glBegin(GL_POLYGON);
		glNormal3f(0.0f,0.0f,-1.0f); 
		for(i = 2*Nxy-1;i >= Nxy;i--) glVertex3fv(p[i]); 
	glEnd();

	//����
	for(i=0;i<Nxy;i++)
	{
		ii = i+1;
		if(ii == Nxy) ii = 0;
		glBegin(GL_QUADS);
			glNormal3f(p[i][0],p[i][1],0.0f);glVertex3fv(p[i]);
			glNormal3f(p[i][0],p[i][1],0.0f);glVertex3fv(p[i+Nxy]);
			glNormal3f(p[ii][0],p[ii][1],0.0f);glVertex3fv(p[ii+Nxy]); 
			glNormal3f(p[ii][0],p[ii][1],0.0f);glVertex3fv(p[ii]);
		glEnd();
	}
}

//----------------------------------------------------------------------
void skWireCylinder(int Nxy)//�P�ʒ��a�A�P�ʍ���
{	
	//Nxy--������
	float p[41][3]; //���_���W
	double theta0,theta;
	int i,ii;

	if(Nxy>20) { Nxy = 20;}
	theta0 = 2*M_PI/(double)Nxy;
	for(i=0;i<Nxy;i++)
	{   theta = theta0*(double)i;
		p[i][0] = (float)(0.5*cos(theta)); //����x����
		p[i][1] = (float)(0.5*sin(theta)); //������
		p[i][2] = 0.5f;                //������(����)
		p[i+Nxy][0] = (float)(0.5*cos(theta)); //�����x����
		p[i+Nxy][1] = (float)(0.5*sin(theta)); //������
		p[i+Nxy][2] = -0.5f;               //������
	}

	//���
	glBegin(GL_POLYGON);
		for(i=0;i<Nxy;i++) glVertex3fv(p[i]); 
	glEnd();
	//����
	glBegin(GL_POLYGON);
		for(i = 2*Nxy-1;i >= Nxy;i--) glVertex3fv(p[i]); 
	glEnd();

	//����
	for(i=0;i<Nxy;i++)
	{
		ii = i+1;
		if(ii == Nxy) ii = 0;
		glBegin(GL_QUADS);
			glVertex3fv(p[i]);
			glVertex3fv(p[i+Nxy]);
			glVertex3fv(p[ii+Nxy]); 
			glVertex3fv(p[ii]);
		glEnd();
	}
}

//-----------------------------------------------------------------------
void skSolidPyramid(int Nxy)//���p��
{   //����̒��S�𕨑̂̒��S
	//���a�A�����P�̑�`
	float nn[3];//�@������
	float p[21][3]; //���_���W
	double r; //���a
	double theta0,theta;
	int i,ii;

	if(Nxy >21) { Nxy = 20; }
	theta0 = 2.0*M_PI/(double)Nxy;
	r=0.5;
	//���_
	p[0][0] =p[0][1] = 0.0;
	p[0][2] = 1.0;
	for(i=1;i<=Nxy;i++)
	{   theta = theta0*(double)i;
		p[i][0] = (float)(r*cos(theta)); //�����x����
		p[i][1] = (float)(r*sin(theta)); //y����
		p[i][2] = 0.0f;                  //������
	}

	//����
	glBegin(GL_POLYGON);
		glNormal3f(0.0f,0.0f,-1.0f); 
		for(i=Nxy;i>=1;i--) glVertex3fv(p[i]); 
	glEnd();

	for(i=1;i<=Nxy;i++)
	{
		ii = i+1;
		if(ii == Nxy+1) ii = 1;
		glBegin(GL_TRIANGLES);
			calcNormal(p[0],p[i],p[ii],nn);
			glNormal3d(nn[0],nn[1],nn[2]); 
			glVertex3fv(p[0]);
			glVertex3fv(p[i]);
			glVertex3fv(p[ii]);
		glEnd();
	}
}
//-----------------------------------------------------------------------
void skWirePyramid(int Nxy)//���p��
{   //����̒��S�𕨑̂̒��S
	//���a�A�����P�̑�`
	float p[21][3]; //���_���W
	double r; //���a
	double theta0,theta;
	int i,ii;

	if(Nxy >21) { Nxy = 20; }
	theta0 = 2.0*M_PI/(double)Nxy;
	r=0.5;
	//���_
	p[0][0] =p[0][1] = 0.0;
	p[0][2] = 1.0;
	for(i=1;i<=Nxy;i++)
	{   theta = theta0*(double)i;
		p[i][0] = (float)(r*cos(theta)); //�����x����
		p[i][1] = (float)(r*sin(theta)); //y����
		p[i][2] = 0.0f;                  //������
	}

	//����
	glBegin(GL_POLYGON);
		for(i=Nxy;i>=1;i--) glVertex3fv(p[i]); 
	glEnd();

	for(i=1;i<=Nxy;i++)
	{
		ii = i+1;
		if(ii == Nxy+1) ii = 1;
		glBegin(GL_TRIANGLES);
			glVertex3fv(p[0]);
			glVertex3fv(p[i]);
			glVertex3fv(p[ii]);
		glEnd();
	}
}

//-----------------------------------------------------------------------
void skSolidCone(int Nxy)//�~��
{   //����̒��S�𕨑̂̒��S
	//���a�A�����P�̑�`
	float nx1,nx2,ny1,ny2,nxy,nz;//�@������
	float p[21][3]; //���_���W
	double r; //���a
	double theta0,theta,alpha;
	int i,ii;

	if(Nxy >21) { Nxy = 21; }
	theta0 = 2.0*M_PI/(double)Nxy;
	r=0.5;
	//���_���W
	p[0][0] =p[0][1] = 0.0;
	p[0][2] = 1.0;
	//��̍��W
	for(i=1;i<=Nxy;i++)
	{   theta = theta0*(double)i;
		p[i][0] = (float)(r*cos(theta)); //���x����
		p[i][1] = (float)(r*sin(theta)); //y����
		p[i][2] = 0.0f;                  //������
	}

	//��
	glBegin(GL_POLYGON);
		glNormal3f(0.0f,0.0f,-1.0f); 
		for(i=Nxy;i>=1;i--) glVertex3fv(p[i]); 
	glEnd();

	alpha = (float)atan(0.5);
	nxy = (float)cos(alpha);
	nz = (float)sin(alpha);
	for(i=1;i<=Nxy;i++)
	{
		if(i==Nxy) ii=1; else ii=i+1;
		glBegin(GL_POLYGON);
			//��̒��_1
			theta = theta0*(double)i;
			nx1 = nxy*(float)cos(theta);
			ny1 = nxy*(float)sin(theta);
			glNormal3d(nx1,ny1,nz); glVertex3fv(p[i]);
			//��̒��_2
			theta = theta0*(double)ii;
			nx2 = nxy*(float)cos(theta); 
			ny2 = nxy*(float)sin(theta);
			glNormal3d(nx2,ny2,nz); glVertex3fv(p[ii]);
			//���_
			nx1 = (nx1+nx2)/2.0f; 
			ny1 = (ny1+ny2)/2.0f;
			glNormal3d(nx1,ny1,nz); glVertex3fv(p[0]);			
		glEnd();
	}
}
//-----------------------------------------------------------------------
void skWireCone(int Nxy)//�~��
{   //����̒��S�𕨑̂̒��S
	//���a�A�����P�̑�`
	float p[21][3]; //���_���W
	double r; //���a
	double theta0,theta;
	int i,ii;

	if(Nxy >21) { Nxy = 21; }
	theta0 = 2.0*M_PI/(double)Nxy;
	r=0.5;
	//���_���W
	p[0][0] =p[0][1] = 0.0;
	p[0][2] = 1.0;
	//��̍��W
	for(i=1;i<=Nxy;i++)
	{   theta = theta0*(double)i;
		p[i][0] = (float)(r*cos(theta)); //���x����
		p[i][1] = (float)(r*sin(theta)); //y����
		p[i][2] = 0.0f;                  //������
	}

	//��
	glBegin(GL_POLYGON);
		for(i=Nxy;i>=1;i--) glVertex3fv(p[i]); 
	glEnd();

	for(i=1;i<=Nxy;i++)
	{
		if(i==Nxy) ii=1; else ii=i+1;
		glBegin(GL_TRIANGLES);
			//��̒��_1
			glVertex3fv(p[i]);
			//��̒��_2
			glVertex3fv(p[ii]);
			//���_
			glVertex3fv(p[0]);
		glEnd();
	}
}
//-----------------------------------------------------------------------
void skSolidFrustum(int Nxy,float ratio)//�~����
{   //����̒��S�𕨑̂̒��S
	//r1:����̔��a�Cr2:���̔��a
	//���a�A�����P�̑�`,���̕����ނ�ratio�{,
	double nn[3];//�@������
	float p[61][3]; //���_���W
	float r1, r2;
	double theta0, theta, alpha, nxy;
	int i,ii;

	if(Nxy >31) { Nxy = 31; }
	theta0 = 2.0*M_PI/(double)Nxy;
	r1 = 0.5;
	r2 = ratio * r1;
	for(i = 0;i < Nxy;i++)
	{   theta = theta0*(double)i;
		p[i][0] = r2*(float)cos(theta); //����x����
		p[i][1] = r2*(float)sin(theta); //������
		p[i][2] = 1.0f;                 //������
		p[i+Nxy][0] = r1*(float)cos(theta); //�����x����
		p[i+Nxy][1] = r1*(float)sin(theta); //y����
		p[i+Nxy][2] = 0.0f;                  //������
	}

	//���
	glBegin(GL_POLYGON);
		glNormal3f(0.0f,0.0f,1.0f);
		for(i=0;i<Nxy;i++) glVertex3fv(p[i]);
	glEnd();
	//����
	glBegin(GL_POLYGON);
		glNormal3f(0.0f,0.0f,-1.0f);
		for(i=2*Nxy-1;i>=Nxy;i--) glVertex3fv(p[i]);
	glEnd();

	alpha = atan(0.5*(1.0-ratio));
	nxy = cos(alpha);
	nn[2] = sin(alpha);
	for(i=0;i<Nxy;i++)
	{
		if(i == Nxy-1) ii = 0; else ii = i+1;
		glBegin(GL_POLYGON);
			theta = theta0*(double)i;
			nn[0] = nxy*cos(theta);
			nn[1] = nxy*sin(theta);
			glNormal3d(nn[0],nn[1],nn[2]);
			glVertex3fv(p[i]); glVertex3fv(p[i+Nxy]);
			theta = theta0*(double)ii;
			nn[0] = nxy*cos(theta);
			nn[1] = nxy*sin(theta);
			glNormal3d(nn[0],nn[1],nn[2]);
			glVertex3fv(p[ii+Nxy]); glVertex3fv(p[ii]);
		glEnd();
	}
}
//-----------------------------------------------------------------------
void skWireFrustum(int Nxy, float ratio)//�~����
{   //����̒��S�𕨑̂̒��S
	//r1:����̔��a�Cr2:���̔��a
	//���a�A�����P�̑�`,���̕����ނ�ratio�{,
	float p[61][3]; //���_���W
	float r1, r2;
	double theta0, theta;//, alpha;
	int i,ii;

	if(Nxy >31) { Nxy = 31; }
	theta0 = 2.0*M_PI/(double)Nxy;
	r1 = 0.5;
	r2 = ratio * r1;
	for(i = 0;i < Nxy;i++)
	{   theta = theta0*(double)i;
		p[i][0] = r2*(float)cos(theta); //����x����
		p[i][1] = r2*(float)sin(theta); //������
		p[i][2] = 1.0f;                 //������
		p[i+Nxy][0] = r1*(float)cos(theta); //�����x����
		p[i+Nxy][1] = r1*(float)sin(theta); //y����
		p[i+Nxy][2] = 0.0f;                  //������
	}

	//���
	glBegin(GL_POLYGON);
		for(i=0;i<Nxy;i++) glVertex3fv(p[i]); 
	glEnd();
	//����
	glBegin(GL_POLYGON);
		for(i=2*Nxy-1;i>=Nxy;i--) glVertex3fv(p[i]); 
	glEnd();

	for(i=0;i<Nxy;i++)
	{
		if(i == Nxy-1) ii = 0; else ii = i+1;
		glBegin(GL_POLYGON);
			glVertex3fv(p[i]);glVertex3fv(p[i+Nxy]);
			glVertex3fv(p[ii+Nxy]);glVertex3fv(p[ii]);
		glEnd();
	}
}

//--------------------------------------------------------
//���݂̖��������`
void skSolidSquare(void)
{
	float p[4][3];
	p[0][0] = 0.5f;  p[0][1] = 0.5f;  p[0][2] = 0.0f;
	p[1][0] = -0.5f; p[1][1] = 0.5f;  p[1][2] = 0.0f;
	p[2][0] = -0.5f; p[2][1] = -0.5f; p[2][2] = 0.0f;
	p[3][0] = 0.5f;  p[3][1] = -0.5f; p[3][2] = 0.0f;

	glBegin(GL_QUADS);
		glNormal3f(0.0f,0.0f,1.0f);//������
		glVertex3fv(p[0]); glVertex3fv(p[1]);
		glVertex3fv(p[2]); glVertex3fv(p[3]);
	glEnd();
}
//--------------------------------------------------------
//���݂̖��������`
void skWireSquare(void)
{
	float p[4][3];
	p[0][0] = 0.5f;  p[0][1] = 0.5f;  p[0][2] = 0.0f;
	p[1][0] = -0.5f; p[1][1] = 0.5f;  p[1][2] = 0.0f;
	p[2][0] = -0.5f; p[2][1] = -0.5f; p[2][2] = 0.0f;
	p[3][0] = 0.5f;  p[3][1] = -0.5f; p[3][2] = 0.0f;

	glBegin(GL_QUADS);
		glVertex3fv(p[0]); glVertex3fv(p[1]);
		glVertex3fv(p[2]); glVertex3fv(p[3]);
	glEnd();
}

//-------------------------------------------------------------------------
//���݂̖������`�g(�O�������̒����͂P�A�����P�A�@������������,���z=0�j
void skSolidFrame0(double wTop, double wFrame)
{   //wFrame:frame�̕�(0.5�ȉ�)
	//wTop:�㕔�̕��i����1�ɑ΂���j
	float p[8][3];

	//���W
	p[0][0] = 0.0; p[0][1] = -0.5; p[0][2] = 0.0;
	p[1][0] = 0.0; p[1][1] =  0.5; p[1][2] = 0.0;
	p[2][0] = 0.0; p[2][1] =  0.5f*(float)wTop; p[2][2] =  1.0;
	p[3][0] = 0.0; p[3][1] = -0.5f*(float)wTop; p[3][2] =  1.0;
	p[4][0] = 0.0; p[4][1] = -0.5f+(float)wFrame; p[4][2] = (float)wFrame;
	p[5][0] = 0.0; p[5][1] =  0.5f-(float)wFrame; p[5][2] = (float)wFrame;
	p[6][0] = 0.0; p[6][1] =  (0.5f-(float)wFrame)*(float)wTop; p[6][2] =  1.0f-(float)wFrame;
	p[7][0] = 0.0; p[7][1] = (-0.5f+(float)wFrame)*(float)wTop; p[7][2] =  1.0f-(float)wFrame;

	glBegin(GL_QUADS);
		glNormal3f(1.0f,0.0f,0.0f);
		glVertex3fv(p[0]); glVertex3fv(p[4]);
		glVertex3fv(p[7]); glVertex3fv(p[3]);
	glEnd();
	glBegin(GL_QUADS);
		glNormal3f(1.0f,0.0f,0.0f);
		glVertex3fv(p[0]); glVertex3fv(p[1]);
		glVertex3fv(p[5]); glVertex3fv(p[4]);
	glEnd();
	glBegin(GL_QUADS);
		glNormal3f(1.0f,0.0f,0.0f);
		glVertex3fv(p[1]); glVertex3fv(p[2]);
		glVertex3fv(p[6]); glVertex3fv(p[5]);
	glEnd();
	glBegin(GL_QUADS);
		glNormal3f(1.0f,0.0f,0.0f);
		glVertex3fv(p[2]); glVertex3fv(p[3]);
		glVertex3fv(p[7]); glVertex3fv(p[6]);
	glEnd();
}
//-------------------------------------------------------------------------
//���݂̖������`�g(�O�������̒����͂P�A�����P�A�@������������,���z=0�j
void skWireFrame0(double wTop, double wFrame)
{   //wFrame:frame�̕�(0.5�ȉ�)
	//wTop:�㕔�̕��i����1�ɑ΂���j
	float p[8][3];

	//���W
	p[0][0] = 0.0; p[0][1] = -0.5; p[0][2] = 0.0;
	p[1][0] = 0.0; p[1][1] =  0.5; p[1][2] = 0.0;
	p[2][0] = 0.0; p[2][1] =  0.5f*(float)wTop; p[2][2] =  1.0;
	p[3][0] = 0.0; p[3][1] = -0.5f*(float)wTop; p[3][2] =  1.0;
	p[4][0] = 0.0; p[4][1] = -0.5f+(float)wFrame; p[4][2] = (float)wFrame;
	p[5][0] = 0.0; p[5][1] =  0.5f-(float)wFrame; p[5][2] = (float)wFrame;
	p[6][0] = 0.0; p[6][1] =  (0.5f-(float)wFrame)*(float)wTop; p[6][2] =  1.0f-(float)wFrame;
	p[7][0] = 0.0; p[7][1] = (-0.5f+(float)wFrame)*(float)wTop; p[7][2] =  1.0f-(float)wFrame;

	glBegin(GL_QUADS);
		glVertex3fv(p[0]); glVertex3fv(p[4]);
		glVertex3fv(p[7]); glVertex3fv(p[3]);
	glEnd();
	glBegin(GL_QUADS);
		glVertex3fv(p[0]); glVertex3fv(p[1]);
		glVertex3fv(p[5]); glVertex3fv(p[4]);
	glEnd();
	glBegin(GL_QUADS);
		glVertex3fv(p[1]); glVertex3fv(p[2]);
		glVertex3fv(p[6]); glVertex3fv(p[5]);
	glEnd();
	glBegin(GL_QUADS);
		glVertex3fv(p[2]); glVertex3fv(p[3]);
		glVertex3fv(p[7]); glVertex3fv(p[6]);
	glEnd();
}
//----------------------------------------------------------------
void skSolidTorus(int Nm, int Ns, float ratio)//2003.3.6�ύX
{	//Nm:��~�������_��,Ns:�f�ʉ~�������_��
	int i,ii,j,jj;
	double phai,ph1,ph2; //��~�������_��x���ɑ΂���Ίp
	double theta,th1,th2; //�f�ʉ~�������_��x-z���ʂɑ΂���Ίp
	float rr[21], zz[21];
	float p[21][21][3]; //�e�_�̍��W
	float p1[3],p2[3],p3[3],p4[3];
	float r1, r2;

	r1 = 0.5;//��~�������a
	r2 = ratio * r1;//�f�ʔ��a
	//if(r1 < r2) { MessageBox(NULL,"skSolidTorus�ɂ����� r1 < r2 �ƂȂ�܂��� ","���a�̒l",MB_OK);return;}
	if( Nm > 21 ) Nm = 21;// { MessageBox(NULL,"skSolidTorus�ɂ�����Nm>20 �ƂȂ�܂��� ","Nm�̒l",MB_OK);return;}
	if( Ns > 21 ) Ns = 21;//{ MessageBox(NULL,"skSolidTorus�ɂ�����Nc>20 �ƂȂ�܂��� ","Nc�̒l",MB_OK);return;}

	//��{�f��(x-z)�̍��W
	for(j = 0; j < Ns; j++)
	{	theta = M_PI-2.0*M_PI*(double)j/(double)Ns;
		rr[j] = r1 + r2 * (float)cos(theta); //���_����̋���
		zz[j] = r2 * (float)sin(theta);//��
	}
	//���̒f�ʂ̍��W
	for(i = 0; i < Nm; i++)
	{	phai = 2.0*M_PI*(double)i/(double)Nm;
		for(j = 0; j < Ns; j++)
		{	
			p[i][j][0] = rr[j] * (float)cos(phai); //x���W
			p[i][j][1] = rr[j] * (float)sin(phai); //y
			p[i][j][2] = zz[j];   //z
		}
	}

	//���_����`���`��
	for(i = 0; i < Nm; i++){
		ii = i+1;
		if(ii == Nm) ii = 0;
		ph1 = 2.0 * M_PI * (double)i/(double)Nm;
		ph2 = 2.0 * M_PI * (double)ii/(double)Nm;
		for(j = 0;j < Ns;j++)
		{
			jj = j+1;
			if(jj == Ns) jj = 0;
			th1 = M_PI-2.0 * M_PI * (double)j / (double)Ns;
			th2 = M_PI-2.0*M_PI*(double)jj/(double)Ns;
			//�ʂ̒��_   x���W                       y���W                       z���W
			p1[0] = p[i][j][0]  ; p1[1] = p[i][j][1]  ; p1[2] = p[i][j][2];
			p2[0] = p[i][jj][0] ; p2[1] = p[i][jj][1] ; p2[2] = p[i][jj][2];
			p3[0] = p[ii][jj][0]; p3[1] = p[ii][jj][1]; p3[2] = p[ii][jj][2];
			p4[0] = p[ii][j][0] ; p4[1] = p[ii][j][1] ; p4[2] = p[ii][j][2];

			glBegin(GL_QUADS);
				glNormal3d(cos(th1)*cos(ph1),cos(th1)*sin(ph1),sin(th1));glVertex3fv(p1); 
				glNormal3d(cos(th2)*cos(ph1),cos(th2)*sin(ph1),sin(th2));glVertex3fv(p2);
				glNormal3d(cos(th2)*cos(ph2),cos(th2)*sin(ph2),sin(th2));glVertex3fv(p3); 
				glNormal3d(cos(th1)*cos(ph2),cos(th1)*sin(ph2),sin(th1));glVertex3fv(p4);
			glEnd();
		}
	}
}
//--------------------------------------------------------------------
//���݊�(��{�p���ł������ǂ̎������j
//���a��data[0],data[1]
//�v�f�̒����͉�,data[2]����
void skSolidTube(int Nxy,int Nz, float *data)
{
	//3����]�idata�́C�p�x�P�ʂ�deg�j
	//Nxy--xy�f�ʕ�����
	//Nz---������������
	//data---��]�p�x(�ŏ���2��x,y�����̒��a�C2�Ԗڂ���Nz�͗v�f�̒����j
	//xyza---���̒��S�����W�ƕ����i�K�w�\���Ƃ���Ƃ��̖߂�l)
	float p[101][21][3]; //���_���W
	double x0, y0, z0, x, y, z, xx, yy, zz;
	float xyz[2][3];
	float pitch[101];
	float nn[3];
	float a[101][20], b[101][20], c[101][20];//���ʂ̖@��
	double sX, sY;//x,y�������a
	double theta0,theta;
	double alpha0, beta0, gamma0;//1�O�̒��S���̉�]�p�x
	double alpha, beta, gamma;//���݂̐ߖʂ̉�]�p�x
	double pp = M_PI/180.0f;
	int i, ii, k;

	if( Nxy > 20 ) Nxy = 20;
	if( Nz > 100 ) Nz = 100;

	sX = data[0];
	sY = data[1];
//char buf[30];
//wsprintf(buf , "%d, %d",(int)(sX * 1000.0), (int)(sY * 1000.0));
//MessageBox(NULL,buf ,"ccc",MB_OK);

    for(k = 0; k < Nz; k++) pitch[k] = data[k + 2];
	alpha0 = beta0 = gamma0 = 0.0;
	x0 = 0.0; y0 = 0.0; z0 = 0.0;//����̒��S
	//k = 0 �̂Ƃ�
	theta0 = 2*M_PI/(double)Nxy;
	for(i = 0;i < Nxy; i++)
	{   theta = theta0 * (double)i;
		//�e���_�̍��W
		p[0][i][0] =  (float)(0.5*cos(theta) * sX); //x����
		p[0][i][1] =  (float)(0.5*sin(theta) * sY); //������
		p[0][i][2] = 0.0f;                          //������
		//�@������
		a[0][i] = p[0][i][0];
		b[0][i] = p[0][i][1];
		c[0][i] = p[0][i][2];
	}

	for( k = 0; k < Nz; k++){
		alpha0 += data[3*k+Nz+3] * pp;//���S����data[3*k]����x����]�i�ݐς����)
		if(k == Nz) alpha = alpha0 ;//+ data[3*k+3] * pp;
		else alpha = alpha0 + data[3*k+Nz+6] * pp / 2.0;
		beta0 += data[3*k+Nz+4] * pp;
		if(k == Nz) beta = beta0 ;//+ data[3*k+4] * pp;
		else beta = beta0 + data[3*k+Nz+7] * pp / 2.0;
		gamma0 += data[3*k+Nz+5] * pp;
		if(k == Nz) gamma = gamma0 ;//+ data[3*k+5] * pp;
		else gamma = gamma0 + data[3*k+Nz+8] * pp / 2.0;

		//���S�̕��s�ړ���(����x=0,y=0)
		x = 0.0; y = 0.0; z = pitch[k]; //pitch[0]=0�Ƃ��Ă�������
		//x����]
		xx = x;
		yy = -z * (float)sin(alpha0);
		zz = z * (float)cos(alpha0) ;
		//y����]
		x = xx; y = yy; z = zz;
		xx = z * (float)sin(beta0);
		yy = y;
		zz = z * (float)cos(beta0);
		//z����]
		x = xx; y = yy; z = zz;
		xx = x * (float)cos(gamma0) - y * (float)sin(gamma0);
		yy = x * (float)sin(gamma0) + y * (float)cos(gamma0);
		zz = z;
		//���s�ړ�
		x0 += xx; y0 += yy; z0 += zz;

		if(k == Nz){//���̒��S���W
			xyz[1][0] = (float)x0;
			xyz[1][1] = (float)y0;
			xyz[1][2] = (float)z0;
		}
		for(i = 0;i < Nxy; i++)
		{
			theta = theta0 * (double)i;
			//���s�ړ��O�̉�]�O�̊e���_���W�ix-y����)
			x = (float)(0.5 * cos(theta)) * sX;
			y = (float)(0.5 * sin(theta)) * sY;
			z = 0.0f;
			//alpha����x ����]
			xx = x ;
			yy = y * (float)cos(alpha) - z * (float)sin(alpha);
			zz = y * (float)sin(alpha) + z * (float)cos(alpha);
			//beta����������]
			x = xx; y = yy; z = zz;
			xx = x * (float)cos(beta) + z * (float)sin(beta);
			yy = y;
			zz = -x * (float)sin(beta) + z * (float)cos(beta);
			//gamma����z ����]
			x = xx; y = yy; z = zz;
			xx = x * (float)cos(gamma) - y * (float)sin(gamma);
			yy = x * (float)sin(gamma) + y * (float)cos(gamma);
			zz = z;
			//�@�������͉�]����
			a[k][i] = (float)xx;
			b[k][i] = (float)yy;
			c[k][i] = (float)zz;
			//���S���W�������s�ړ�
			p[k][i][0] = (float)(xx + x0);
			p[k][i][1] = (float)(yy + y0);
			p[k][i][2] = (float)(zz + z0);
		}//i
	}//k
/*	//��[�̈ʒu�C��]�p�i�K�w�\���ɂ��邽�߂̖߂�l�j
	xyza[0] = (float)x0;
	xyza[1] = (float)y0;
	xyza[2] = (float)z0;
	xyza[3] = (float)alpha;
	xyza[4] = (float)beta;
	xyza[5] = (float)gamma; */
	//���
	glBegin(GL_POLYGON);
		calcNormal(xyz[1], p[Nz][0], p[Nz][1],nn);
		glNormal3f(nn[0],nn[1],nn[2]);
		for(i = 0;i < Nxy; i++) glVertex3fv(p[Nz][i]);
	glEnd();
	//����
	xyz[0][0] = 0.0; xyz[0][1] = 0.0; xyz[0][2] = 0.0;
	glBegin(GL_POLYGON);
		calcNormal(xyz[0], p[0][1], p[0][0],nn);
		glNormal3f(nn[0],nn[1],nn[2]);
		for(i = Nxy-1;i >= 0; i--) glVertex3fv(p[0][i]);
	glEnd();

	//����
	for(k = 0 ;k < Nz; k++){
		for(i = 0;i < Nxy;i++)
		{
			ii = i + 1;
			if(ii == Nxy) ii = 0;
			glBegin(GL_POLYGON);
				glNormal3f(a[k][i],b[k][i],c[k][i]);glVertex3fv(p[k][i]);
				glNormal3f(a[k][ii],b[k][ii],c[k][ii]);glVertex3fv(p[k][ii]);
				glNormal3f(a[k+1][ii],b[k+1][ii],c[k+1][ii]);glVertex3fv(p[k+1][ii]);
				glNormal3f(a[k+1][i],b[k+1][i],c[k+1][i]);glVertex3fv(p[k+1][i]);
			glEnd();
		}
	}
}
//----------------------------------------------------------------
void skSolidSpring(int Nm, int Ns, int Np, double radius, double ratio, double length)
{
	//Nm:��~�������_��,Ns:�f�ʉ~�������_��,Np:�s�b�`��
	//radius:��~�������a,ratio:�f�ʔ��a/��~�������a,length:�S��
	int i, ii, j, jj, k, k1, k2;
	double phai, ph1, ph2; //��~�������_��x���ɑ΂���Ίp
	double theta, th1, th2; //�f�ʉ~�������_��x-z���ʂɑ΂���Ίp
	float p[21][21][101][3];//p[44100][3]; //�e�_�̍��W
	float p1[3], p2[3], p3[3], p4[3];
	double rr[21], zz[21];
	double pitch, dp, hh;
	double r1, r2;

	r1 = radius; //��~�������a
	r2 = ratio * r1;//�f�ʔ��a
	if( Nm > 20 ) Nm = 20;
	if( Ns > 20 ) Ns = 20;
	if( Np > 100 ) Np = 100;
	pitch = length / (double)Np;
    //�k�񂾂Ƃ��̐���
	if(pitch < 2 * r2) pitch = 2.0f * r2 ;

	dp = pitch / (double)Nm;

	//��{�f��(x-z)�̍��W
	for(j = 0; j < Ns; j++)
	{	theta = M_PI - 2.0 * M_PI*(double)j/(double)Ns;
		rr[j] = r1 + r2 * cos(theta); //���_����̋���
		zz[j] = r2 * sin(theta);//��
	}

	//���̒f�ʂ̍��W
	hh = 0;
	for(k = 0; k < Np; k++)
		for(i = 0; i < Nm; i++)
		{	phai = 2.0 * M_PI * (double)i/(double)Nm;
			for(j = 0; j < Ns; j++)
			{
				p[i][j][k][0] = (float)(rr[j] * cos(phai)); //x���W
				p[i][j][k][1] = (float)(rr[j] * sin(phai)); //y
				p[i][j][k][2] = (float)(zz[j] + hh) ;              //z
			}
			hh += dp;//���S���̍�����dp������
		}

	//�ŏI�[(k=Np-1,i=Nm)
	k = Np - 1; i = Nm;
	for(j = 0; j < Ns; j++){
		phai = 0.0;
		p[i][j][k][0] = (float)(rr[j] * cos(phai)); //x���W
		p[i][j][k][1] = (float)(rr[j] * sin(phai)); //y
		p[i][j][k][2] = (float)(zz[j] + hh) ;            //z
	}

	//���_����`���`��
	for(k = 0; k < Np; k++)
	for(i = 0; i < Nm; i++){
		ii = i+1;
		k1 = k; k2 = k;
		if(ii == Nm) {
			if(k < Np-1) { ii = 0; k2 = k + 1; }
		}
		ph1 = 2.0*M_PI*(double)i / (double)Nm;
		ph2 = 2.0*M_PI*(double)ii / (double)Nm;
		for(j = 0;j < Ns; j++)
		{
			jj = j+1;
			if(jj == Ns) jj = 0;
			th1 = M_PI - 2.0 * M_PI * (double)j / (double)Ns;
			th2 = M_PI - 2.0 * M_PI * (double)jj / (double)Ns;
			//�ʂ̒��_   x���W                       y���W                       z���W
			p1[0] = p[i][j][k1][0]   ; p1[1] = p[i][j][k1][1]   ; p1[2] = p[i][j][k1][2];
			p2[0] = p[i][jj][k1][0]  ; p2[1] = p[i][jj][k1][1]  ; p2[2] = p[i][jj][k1][2];
			p3[0] = p[ii][jj][k2][0] ; p3[1] = p[ii][jj][k2][1] ; p3[2] = p[ii][jj][k2][2];
			p4[0] = p[ii][j][k2][0]  ; p4[1] = p[ii][j][k2][1]  ; p4[2] = p[ii][j][k2][2];

			glBegin(GL_QUADS);
				glNormal3d(cos(th1)*cos(ph1),cos(th1)*sin(ph1),sin(th1));glVertex3fv(p1);
				glNormal3d(cos(th2)*cos(ph1),cos(th2)*sin(ph1),sin(th2));glVertex3fv(p2);
				glNormal3d(cos(th2)*cos(ph2),cos(th2)*sin(ph2),sin(th2));glVertex3fv(p3);
				glNormal3d(cos(th1)*cos(ph2),cos(th1)*sin(ph2),sin(th1));glVertex3fv(p4);
			glEnd();
		}
	}
	//�n�[�ik=0)
	glBegin(GL_POLYGON);
		glNormal3d(0.0,-1.0,0.0);
		for(j = Ns-1; j >= 0; j--) glVertex3fv(p[0][j][0]);
	glEnd();
	//�I�[�ik=Np)
	glBegin(GL_POLYGON);
		glNormal3d(0.0,1.0,0.0);
		for(j = 0; j < Ns; j++) glVertex3fv(p[Nm][j][Np-1]);
	glEnd();

}
//-----------------------------------------------------------------
//��2���֐�
void skSolidSuper1(int Nxy, int Nz, double eps1, double eps2, double p1, double p2, double p3)
{
	//�㉺�̒��S�����_
	int i,j,ip,im,np,npL,npR,npU,npD,k1,k2;
	double ct,phai,theta,z,fz;
	float a[31][31],b[31][31],c[31][31];
	float n1[3],n2[3],n3[3],n4[3];
	double cc;
	float pd[961][3];

	if( Nxy > 30 ) Nxy = 30;
	if( Nz > 30 ) Nz = 30;

	//���a
        float r = 0.5f;
	//�㔼��
        for(j = 0 ;j <= Nz;j++)
	{
                phai = (M_PI/(double)Nz) * ((double)Nz / 2.0 - (double)j);
                if(phai > 0.0) z = (float)(pow(sin(phai),eps1));//z
                else z = - (float)(pow(sin(fabs(phai)), eps1));
                //�`��֐�
                if(z < 0.0) fz = (p2-p3)*z + p2;
                else fz = (p1-p2)*z + p2;
                for (i = 0 ;i <= Nxy / 2;i++)
                {
			k1 = Nxy * j + i;//�������猩�č���
			k2 = Nxy * j + Nxy - i;//�E��
			theta = 2.0 * M_PI * (double)i/(double)Nxy;
			ct = cos(theta);
			if (ct >= 0) { cc = pow(ct, eps2);}
                        else         { cc = -pow(fabs(ct),eps2); }
			if(j == 0 || j == Nz) {
				pd[k1][0] = 0.0f;
				pd[k1][1] = 0.0f;
			}

			else {
				pd[k1][0] = r * (float)(pow(cos(phai),eps1)*cc*fz);
				pd[k1][1] = r * (float)(pow(cos(phai),eps1)*pow(fabs(sin(theta)),eps2)*fz);
			}
                        pd[k2][0] = pd[k1][0];
                        pd[k2][1] = -pd[k1][1];
			pd[k1][2] = r * (float)z;
			pd[k2][2] = r * (float)z;
                }
        }


	//���ʂ̖@������
	for(i = 0;i < Nxy;i++){
 		ip = i+1;
		if(ip == Nxy) ip = 0;
		im = i-1;
		if(i == 0) im = Nxy-1;

		//�^��(Top)
		a[i][0] = 0.0; b[i][0] = 0.0; c[i][0] = 1.0;
		//�^���iBottom)
		a[i][Nz] = 0.0; b[i][Nz] = 0.0; c[i][Nz] = -1.0;

		for(j=1;j<Nz;j++)//�ׂ荇��4�̎O�p�`�̖@���x�N�g���𕽋ω�
		{
			np = j*Nxy+i;//���ړ_
			npL = j*Nxy+im;//����
			npR = j*Nxy+ip;//�E��
			npU = np-Nxy;//��
			npD = np+Nxy;//��
			if(j == 1) {
				n1[0]=0.0; n1[1]=0.0; n1[2]=1.0;//Top
				n2[0]=0.0; n2[1]=0.0; n2[2]=1.0;//Top
				calcNormal(pd[np],pd[npL],pd[npD],n3);//�O���猩�č���
				calcNormal(pd[np],pd[npD],pd[npR],n4);//�E��
			}
			if(j == Nz-1){
				calcNormal(pd[np],pd[npU],pd[npL],n1);//�O���猩�č���
				calcNormal(pd[np],pd[npR],pd[npU],n2);//�E��
				n3[0]=0.0; n3[1]=0.0; n3[2]=-1.0;//Bottom
				n4[0]=0.0; n4[1]=0.0; n4[2]=-1.0;//Bottom
			}
			else {
				calcNormal(pd[np],pd[npU],pd[npL],n1);//�O���猩�č���
				calcNormal(pd[np],pd[npR],pd[npU],n2);//�E��
				calcNormal(pd[np],pd[npL],pd[npD],n3);//�O���猩�č���
				calcNormal(pd[np],pd[npD],pd[npR],n4);//�E��
			}
			a[i][j] = (float)((n1[0]+n2[0]+n3[0]+n4[0])/4.0f);//������
			b[i][j] = (float)((n1[1]+n2[1]+n3[1]+n4[1])/4.0f);//��
			c[i][j] = (float)((n1[2]+n2[2]+n3[2]+n4[2])/4.0f);//��
		}
	}
	//�\��
	//���
	for(i = 0;i < Nxy;i++)
	{	ip = i+1;
		if(ip == Nxy) ip = 0;

		glBegin(GL_TRIANGLES);
			glNormal3d(a[i][0],b[i][0],c[i][0]);glVertex3fv(pd[i]);
			glNormal3d(a[i][1],b[i][1],c[i][1]);glVertex3fv(pd[Nxy+i]);
			glNormal3d(a[ip][1],b[ip][1],c[ip][1]);glVertex3fv(pd[Nxy+ip]);
		glEnd();
	}
	//����
	j = Nz - 1;
	for(i = 0;i < Nxy;i++)
	{	ip = i+1;
		if(ip == Nxy) ip = 0;

		glBegin(GL_TRIANGLES);
			glNormal3d(a[i][j],b[i][j],c[i][j]);glVertex3fv(pd[j*Nxy+i]);
			glNormal3d(a[i][j+1],b[i][j+1],c[i][j+1]);glVertex3fv(pd[(j+1)*Nxy+i]);
			glNormal3d(a[ip][j],b[ip][j],c[ip][j]);glVertex3fv(pd[j*Nxy+ip]);
		glEnd();
	}
	//����(4�p�`�p�b�`�j
	//for(j = 1;j < Nz-1;j++)
	for(j = 0;j < Nz;j++)
		for(i = 0;i < Nxy;i++)
		{	ip = i+1;
			if(ip == Nxy) ip=0;

			glBegin(GL_QUADS);
				glNormal3d(a[i][j],b[i][j],c[i][j]);glVertex3fv(pd[j*Nxy+i]);
				glNormal3d(a[i][j+1],b[i][j+1],c[i][j+1]);glVertex3fv(pd[(j+1)*Nxy+i]);
				glNormal3d(a[ip][j+1],b[ip][j+1],c[ip][j+1]);glVertex3fv(pd[(j+1)*Nxy+ip]);
				glNormal3d(a[ip][j],b[ip][j],c[ip][j]);glVertex3fv(pd[j*Nxy+ip]);
			glEnd();
		}
}
//-----------------------------------------------------------------
//�㔼���̒�2���֐�
void skSolidSuper2(int Nxy, int Nz, double eps1, double eps2, double p1, double p2)
{
	//��̒��S�����_
	int i,j,ip,im,np,npL,npR,npU,npD,k1,k2;
	double ct,phai,theta,z,fz;
	float a[31][31],b[31][31],c[31][31];
	float n1[3],n2[3],n3[3],n4[3];
	double cc;
	float pd[961][3];

	if( Nxy > 30 ) Nxy = 30;
	if( Nz > 30 ) Nz = 30;

	//���a
    float r = 0.5f;
	//�㔼������

    for(j = 0 ;j <= Nz ;j++)
	{
        phai = (M_PI/(double)Nz) * ((double)Nz  - (double)j)/2.0;
		z = (float)(pow(sin(phai),eps1));//z
//		fz = (taper - 1.0) * z + 1.0;
//		fz = 0.5f * (taper + 1.0f + (taper -1.0f) * z);
        //�`��֐�
        fz = (p1-p2)*z + p2;
        for (i = 0 ;i<= Nxy / 2;i++)
        {
			k1 = Nxy * j + i;//�O���猩�ĉE��
			k2 = Nxy * j + Nxy - i;//����
			theta = 2.0*M_PI * (double)i / (double)Nxy;
			ct = cos(theta);
			if (ct >= 0) { cc = pow(ct, eps2); }
            else         { cc = -pow(fabs(ct),eps2); }
			if(j == 0) {
				pd[k1][0] = 0.0f;
				pd[k1][1] = 0.0f;
			}
			else{
 				pd[k1][0] = r * (float)(pow(cos(phai),eps1) * cc * fz);
				pd[k1][1] = r * (float)(pow(cos(phai),eps1) * pow(fabs(sin(theta)),eps2) * fz);
			}
            pd[k2][0] = pd[k1][0];
            pd[k2][1] = -pd[k1][1];
			pd[k1][2] = r * (float)z;
			pd[k2][2] = r * (float)z;
        }
    }

	//���ʂ̖@������
	for(i = 0;i < Nxy; i++){
 	  ip = i + 1;
	  if(ip==Nxy) ip = 0;
	  im=i-1;
	  if(i==0) im = Nxy - 1;

		//�^��(Top)
		a[i][0] = 0.0; b[i][0] = 0.0; c[i][0] = 1.0;
	  for(j=1;j<Nz;j++)
	  {
		np = j*Nxy+i;//���ړ_
		npL = j*Nxy+im;//����
		npR = j*Nxy+ip;//�E��
		npU = np-Nxy;//��
		npD = np+Nxy;//��
		if(j==1) {
			n1[0]=0.0; n1[1]=0.0; n1[2]=1.0;
			n2[0]=0.0; n2[1]=0.0; n2[2]=1.0;
			calcNormal(pd[np],pd[npL],pd[npD],n3);//�O���猩�č���
			calcNormal(pd[np],pd[npD],pd[npR],n4);//�E��
		}
		else {
			calcNormal(pd[np],pd[npU],pd[npL],n1);//�O���猩�č���
			calcNormal(pd[np],pd[npR],pd[npU],n2);//�E��
			calcNormal(pd[np],pd[npL],pd[npD],n3);//�O���猩�č���
			calcNormal(pd[np],pd[npD],pd[npR],n4);//�E��
		}
		a[i][j] = (float)((n1[0]+n2[0]+n3[0]+n4[0])/4.0f);//������
		b[i][j] = (float)((n1[1]+n2[1]+n3[1]+n4[1])/4.0f);//��
		c[i][j] = (float)((n1[2]+n2[2]+n3[2]+n4[2])/4.0f);//��
	  }
	  j = Nz;//��ԉ��̑���(���2�̎O�p�`�̕��ρj
		calcNormal(pd[np],pd[npU],pd[npL],n1);//�O���猩�č���
		calcNormal(pd[np],pd[npR],pd[npU],n2);//�E��
		a[i][j] = (float)((n1[0]+n2[0])/2.0f);
		b[i][j] = (float)((n1[1]+n2[1])/2.0f);
		c[i][j] = (float)((n1[2]+n2[2])/2.0f);
	}

	//�\��
	//���
	for(i = 0;i < Nxy;i++)
	{	ip = i+1;
		if(ip == Nxy) ip = 0;

		glBegin(GL_TRIANGLES);
			glNormal3d(a[i][0],b[i][0],c[i][0]);glVertex3fv(pd[i]);
			glNormal3d(a[i][1],b[i][1],c[i][1]);glVertex3fv(pd[Nxy+i]);
			glNormal3d(a[ip][1],b[ip][1],c[ip][1]);glVertex3fv(pd[Nxy+ip]);
		glEnd();
	}
	//����
	glBegin(GL_POLYGON);
		glNormal3f(0.0f,0.0f,-1.0f);
		for(i=(Nz+1)*Nxy-1;i>=Nz*Nxy;i--) glVertex3fv(pd[i]);
	glEnd();
	//����(4�p�`�p�b�`�j
	for(j = 1;j < Nz;j++)
		for(i = 0;i < Nxy;i++)
		{	ip=i+1;
			if(ip==Nxy) ip=0;

			glBegin(GL_QUADS);
				glNormal3d(a[i][j],b[i][j],c[i][j]);glVertex3fv(pd[j*Nxy+i]);
				glNormal3d(a[i][j+1],b[i][j+1],c[i][j+1]);glVertex3fv(pd[(j+1)*Nxy+i]);
				glNormal3d(a[ip][j+1],b[ip][j+1],c[ip][j+1]);glVertex3fv(pd[(j+1)*Nxy+ip]);
				glNormal3d(a[ip][j],b[ip][j],c[ip][j]);glVertex3fv(pd[j*Nxy+ip]);
			glEnd();
		}
}
//--------------------------------------------------------------------
//�P�����̃e�[�v
void skSolidTape1(int Np, float *data)//�P�ʕ��A
{	//�Z���̃T�C�Y���P,�S��1
	//�������̒��������ƂȂ�悤�ɃZ�����g�̒����͉ρj
	//��{�p���̎�������x��
	//Np---x�������Z����
	float pitch ;
	float p[501][2][3]; //���_���W
	float a[501][2], b[501][2], c[501][2];//�@������
	float nn[3];
	int k;

	if( Np > 500 ) Np = 500;
	pitch = 1.0f / (float)Np;

	//�e���_�̍��W
	for(k = 0; k <= Np; k++)
	{
		//y > 0���̍��W
		p[k][0][0] = (float)k * pitch;
		p[k][0][1] = 0.5f;
		p[k][0][2] = data[k];
		//y < 0 ���̍��W
		p[k][1][0] = (float)k * pitch;
		p[k][1][1] = -0.5f;
		p[k][1][2] = data[k];
	}
	//�e�ʂ̖@������
	for(k = 0; k < Np; k++)
	{
		calcNormal(p[k][0],p[k][1],p[k+1][1],nn);
		a[k][0] = a[k][1] = nn[0];//x����
		b[k][0] = b[k][1] = nn[1];//y����
		c[k][0] = c[k][1] = nn[2];//z����
	}
	a[Np][0] = a[Np][1] = a[Np-1][0];//�I�[��1�O�Ɠ����Ƃ���
	b[Np][0] = b[Np][1] = b[Np-1][0];
	c[Np][0] = c[Np][1] = c[Np-1][0];

	//�\�ʂ̕`��
	for(k = 0 ;k < Np; k++)
	{
		glBegin(GL_POLYGON);
			glNormal3f(a[k][0],b[k][0],c[k][0]);glVertex3fv(p[k][0]);
			glNormal3f(a[k][1],b[k][1],c[k][1]);glVertex3fv(p[k][1]);
			glNormal3f(a[k+1][1],b[k+1][1],c[k+1][1]);glVertex3fv(p[k+1][1]);
			glNormal3f(a[k+1][0],b[k+1][0],c[k+1][0]);glVertex3fv(p[k+1][0]);
		glEnd();
	}
	//���ʂ̕`��
	for(k = 0 ;k < Np; k++)
	{
		glBegin(GL_POLYGON);
			glNormal3f(-a[k][0],-b[k][0],-c[k][0]);glVertex3fv(p[k][0]);
			glNormal3f(-a[k+1][0],-b[k+1][0],-c[k+1][0]);glVertex3fv(p[k+1][0]);
			glNormal3f(-a[k+1][1],-b[k+1][1],-c[k+1][1]);glVertex3fv(p[k+1][1]);
			glNormal3f(-a[k][1],-b[k][1],-c[k][1]);glVertex3fv(p[k][1]);
		glEnd();
	}
}
//--------------------------------------------------------------------
void skSolidTape2(int Np, float *data)
{
	//�Z��(���P,����1/Np,�e�[�v�S��1�j
	//�Z�����g�̒����͌Œ�i�������̒����͕ω�����j
	//��{�p���̎�������x��
	//data[k]�͐����g�̌��z�f�[�^�i�P�ʂ�rad�j
	//Np---�������s�b�`��
	float pitch ;
	float p[501][2][3]; //���_���W
	float nn[3];
	float a[501][2], b[501][2], c[501][2];//�@������
	int k;

	if( Np > 500 ) Np = 500;
	pitch = 1.0f / (float)Np;

	//k = 0(���_�j
	p[0][0][0] = p[0][1][0] = 0.0;      //x����
	p[0][0][1] = 0.5; p[0][1][1] = -0.5;//y����
	p[0][0][2] = p[0][1][2] = 0.0f;//z����

	for(k = 1; k <= Np; k++)
	{
		p[k][0][0] = p[k][1][0] = p[k-1][0][0] + (float)cos(data[k-1]) * pitch;//x����
		p[k][0][1] = 0.5; p[k][1][1] = -0.5;                                   //y����
		p[k][0][2] = p[k][1][2] = p[k-1][0][2] + (float)sin(data[k-1]) * pitch;//z����
	}

	//�e�ʂ̖@������
	for(k = 0; k < Np; k++)
	{
		calcNormal(p[k][0],p[k][1],p[k+1][1],nn);
		a[k][0] = a[k][1] = nn[0];//x����
		b[k][0] = b[k][1] = nn[1];//y����
		c[k][0] = c[k][1] = nn[2];//z����
	}
	a[Np][0] = a[Np][1] = a[Np-1][0];//�I�[��1�O�Ɠ����Ƃ���
	b[Np][0] = b[Np][1] = b[Np-1][0];
	c[Np][0] = c[Np][1] = c[Np-1][0];


	//�\�ʂ̕`��
	for(k = 0 ;k < Np; k++)
	{
		glBegin(GL_POLYGON);
			glNormal3f(a[k][0],b[k][0],c[k][0]);glVertex3fv(p[k][0]);
			glNormal3f(a[k][1],b[k][1],c[k][1]);glVertex3fv(p[k][1]);
			glNormal3f(a[k+1][1],b[k+1][1],c[k+1][1]);glVertex3fv(p[k+1][1]);
			glNormal3f(a[k+1][0],b[k+1][0],c[k+1][0]);glVertex3fv(p[k+1][0]);
		glEnd();
	}
	//���ʂ̕`��
	for(k = 0 ;k < Np; k++)
	{
		glBegin(GL_POLYGON);
			glNormal3f(-a[k][0],-b[k][0],-c[k][0]);glVertex3fv(p[k][0]);
			glNormal3f(-a[k+1][0],-b[k+1][0],-c[k+1][0]);glVertex3fv(p[k+1][0]);
			glNormal3f(-a[k+1][1],-b[k+1][1],-c[k+1][1]);glVertex3fv(p[k+1][1]);
			glNormal3f(-a[k][1],-b[k][1],-c[k][1]);glVertex3fv(p[k][1]);
		glEnd();
	}
}

//-----------------------------------------------------------------------------
//�����`�̃��b�V���ix-y����,���S�����_�j
//���������C���������̕��͌Œ�
void skSolidMesh(int Nx, float* data)
{
	//�S�̂̕�,�����ǂ����1(�����`)
	int i, j, Ny;
	float p[101][101][3]; //���_���W
	float a[101][101], b[101][101], c[101][101];//
	float nn[3];
	float pitchX, pitchY;


	if(Nx > 100) Nx = 100;
	Ny = Nx;
	//�Z���̃T�C�Y
	pitchX = 1.0f / (float)Nx;
	pitchY = 1.0f / (float)Ny;

	//�e���_�̍��W
	for(j = 0; j <= Ny; j++){
		for(i = 0; i <= Nx; i++){
			p[i][j][0] = (float)(i - Nx / 2) * pitchX;
			p[i][j][1] = (float)(j - Ny / 2) * pitchY;
			p[i][j][2] = data[j * (Nx+1) + i];
		}
	}
	//�e�ʂ̖@������
	for(i = 0; i < Nx; i++){
		for(j = 0; j < Ny; j++){
			calcNormal(p[i][j],p[i+1][j],p[i][j+1],nn);
			a[i][j] = nn[0];//x����
			b[i][j] = nn[1];//y����
			c[i][j] = nn[2];//z����
		}
	}
	//j=Ny;
	for(i = 0; i < Nx; i++) {
		a[i][Ny] = a[i][Ny-1];
		b[i][Ny] = b[i][Ny-1];
		c[i][Ny] = c[i][Ny-1];
	}
	//i=Nx
	for(j = 0; j <= Ny; j++) {
		a[Nx][j] = a[Nx-1][j];
		b[Nx][j] = b[Nx-1][j];
		c[Nx][j] = c[Nx-1][j];
	}


	//�\�ʂ̕`��
	for(i = 0 ;i < Nx; i++){
		for(j = 0; j < Ny; j++){
			glBegin(GL_POLYGON);
				glNormal3f(a[i][j],b[i][j],c[i][j]);glVertex3fv(p[i][j]);
				glNormal3f(a[i+1][j],b[i+1][j],c[i+1][j]);glVertex3fv(p[i+1][j]);
				glNormal3f(a[i+1][j+1],b[i+1][j+1],c[i+1][j+1]);glVertex3fv(p[i+1][j+1]);
				glNormal3f(a[i][j+1],b[i][j+1],c[i][j+1]);glVertex3fv(p[i][j+1]);
			glEnd();
		}
	}
	//���ʂ̕`��
	for(i = 0 ;i < Nx; i++){
		for(j = 0; j < Ny; j++){
			glBegin(GL_POLYGON);
				glNormal3f(-a[i][j],-b[i][j],-c[i][j]);glVertex3fv(p[i][j]);
				glNormal3f(-a[i][j+1],-b[i][j+1],-c[i][j+1]);glVertex3fv(p[i][j+1]);
				glNormal3f(-a[i+1][j+1],-b[i+1][j+1],-c[i+1][j+1]);glVertex3fv(p[i+1][j+1]);
				glNormal3f(-a[i+1][j],-b[i+1][j],-c[i+1][j]);glVertex3fv(p[i+1][j]);
			glEnd();
		}
	}
}

