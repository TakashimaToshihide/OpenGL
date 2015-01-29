//mode
#define T_DECAL  0
#define T_MODULATE 1
//texType
#define T_NON 0
#define T_PLANAR1 1
#define T_PLANAR2 2
#define T_SPHERICAL 3
#define T_CYLINDRICAL 4
#define T_SOLID 5
//ø����̍ő�T�C�Y
#define T_MAX 128   

//------------------------------------------------------------------------------
void skTexSquare()//y-z����
{
	static float p[4][3]={{0.0,-0.5,-0.5},{0.0,0.5,-0.5},{0.0,0.5,0.5},
	{0.0,-0.5,0.5} };

	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);
		glNormal3f(1.0f,0.0f,0.0f); //x�����̖@��
		//ø�������W�ƒ��_�ԍ��Ƃ̑Ή��t��
		glTexCoord2f(0.0f, 0.0f); glVertex3fv(p[0]);
		glTexCoord2f(1.0f, 0.0f); glVertex3fv(p[1]);
		glTexCoord2f(1.0f, 1.0f); glVertex3fv(p[2]);
		glTexCoord2f(0.0f, 1.0f); glVertex3fv(p[3]);
	glEnd();
	//���ʂ�����ø����
	glBegin(GL_QUADS);
		glNormal3f(-1.0f,0.0f,0.0f); //x�����̖@��
		//ø�������W�ƒ��_�ԍ��Ƃ̑Ή��t��
		glTexCoord2f(0.0f, 0.0f); glVertex3fv(p[0]);
		glTexCoord2f(1.0f, 0.0f); glVertex3fv(p[3]);
		glTexCoord2f(1.0f, 1.0f); glVertex3fv(p[2]);
		glTexCoord2f(0.0f, 1.0f); glVertex3fv(p[1]);
	glEnd();
	glDisable(GL_TEXTURE_2D);
}

//-----------------------------------------------------------------------
//���ʂ���ϯ��ݸ�
void skTexCube1()
{
	float p[8][3]={{0.5f,0.5f,0.5f},{-0.5f,0.5f,0.5f},{-0.5f,-0.5f,0.5f},
	{0.5f,-0.5f,0.5f},{0.5f,0.5f,-0.5f},{-0.5f,0.5f,-0.5f},{-0.5f,-0.5f,0-0.5f},{0.5f,-0.5f,-0.5f}};
//glTexImage2D(GL_TEXTURE_2D,0,4,T_MAX,T_MAX,0,GL_RGBA,GL_UNSIGNED_BYTE,texImage);

	glBegin(GL_QUADS);
		glNormal3f(0.0f,0.0f,1.0f); //z����
		glVertex3fv(p[0]); glVertex3fv(p[1]);
		glVertex3fv(p[2]); glVertex3fv(p[3]);
	glEnd();

	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);
		glNormal3f(1.0f,0.0f,0.0f); //x����(���ʁj
		glTexCoord2f(1.0f, 1.0f); glVertex3fv(p[0]);
		glTexCoord2f(0.0f, 1.0f); glVertex3fv(p[3]);
		glTexCoord2f(0.0f, 0.0f); glVertex3fv(p[7]);
		glTexCoord2f(1.0f, 0.0f); glVertex3fv(p[4]);
	glEnd();
	glDisable(GL_TEXTURE_2D);

	glBegin(GL_QUADS);
		glNormal3f(0.0f,1.0f,0.0f); //y����
		glVertex3fv(p[0]); glVertex3fv(p[4]);
		glVertex3fv(p[5]); glVertex3fv(p[1]);
	glEnd();

	glBegin(GL_QUADS);
	 	glNormal3f(-1.0f,0.0f,0.0f); //-x����
		glVertex3fv(p[1]); glVertex3fv(p[5]);
		glVertex3fv(p[6]); glVertex3fv(p[2]);
	glEnd();

	glBegin(GL_QUADS);
		glNormal3f(0.0f,-1.0f,0.0f); //-y����
		glVertex3fv(p[2]); glVertex3fv(p[6]);
		glVertex3fv(p[7]); glVertex3fv(p[3]);
	glEnd();

	glBegin(GL_QUADS);
		glNormal3f(0.0f,0.0f,-1.0f); //-z����
		glVertex3fv(p[4]); glVertex3fv(p[7]);
		glVertex3fv(p[6]); glVertex3fv(p[5]);
	glEnd();
}
//-----------------------------------------------------------------------
//�U�ʂɓ����͗l��ϯ��ݸ�
void skTexCube2()
{
	float p[8][3]={{0.5f,0.5f,0.5f},{-0.5f,0.5f,0.5f},{-0.5f,-0.5f,0.5f},
	{0.5f,-0.5f,0.5f},{0.5f,0.5f,-0.5f},{-0.5f,0.5f,-0.5f},{-0.5f,-0.5f,0-0.5f},{0.5f,-0.5f,-0.5f}};

	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);//top
		glNormal3f(0.0f,0.0f,1.0f); //z����
		glTexCoord2f(1.0f, 0.0f); glVertex3fv(p[0]);
		glTexCoord2f(1.0f, 1.0f); glVertex3fv(p[1]);
		glTexCoord2f(0.0f, 1.0f); glVertex3fv(p[2]);
		glTexCoord2f(0.0f, 0.0f); glVertex3fv(p[3]);
	glEnd();

	glBegin(GL_QUADS);
		glNormal3f(1.0f,0.0f,0.0f); //x����(���ʁj
		glTexCoord2f(1.0f, 1.0f); glVertex3fv(p[0]);
		glTexCoord2f(0.0f, 1.0f); glVertex3fv(p[3]);
		glTexCoord2f(0.0f, 0.0f); glVertex3fv(p[7]);
		glTexCoord2f(1.0f, 0.0f); glVertex3fv(p[4]);
	glEnd();

	glBegin(GL_QUADS);
		glNormal3f(0.0f,1.0f,0.0f); //y����
		glTexCoord2f(0.0f, 1.0f); glVertex3fv(p[0]);
		glTexCoord2f(0.0f, 0.0f); glVertex3fv(p[4]);
		glTexCoord2f(1.0f, 0.0f); glVertex3fv(p[5]);
		glTexCoord2f(1.0f, 1.0f); glVertex3fv(p[1]);
	glEnd();

	glBegin(GL_QUADS);
	 	glNormal3f(-1.0f,0.0f,0.0f); //-x����
		glTexCoord2f(0.0f, 1.0f); glVertex3fv(p[1]);
		glTexCoord2f(0.0f, 0.0f); glVertex3fv(p[5]);
		glTexCoord2f(1.0f, 0.0f); glVertex3fv(p[6]);
		glTexCoord2f(1.0f, 1.0f); glVertex3fv(p[2]);
	glEnd();

	glBegin(GL_QUADS);
		glNormal3f(0.0f,-1.0f,0.0f); //-y����
		glTexCoord2f(0.0f, 1.0f); glVertex3fv(p[2]);
		glTexCoord2f(0.0f, 0.0f); glVertex3fv(p[6]);
		glTexCoord2f(1.0f, 0.0f); glVertex3fv(p[7]);
		glTexCoord2f(1.0f, 1.0f); glVertex3fv(p[3]);
	glEnd();

	glBegin(GL_QUADS);
		glNormal3f(0.0f,0.0f,-1.0f); //-z����
		glTexCoord2f(1.0f, 1.0f); glVertex3fv(p[4]);
		glTexCoord2f(0.0f, 1.0f); glVertex3fv(p[7]);
		glTexCoord2f(0.0f, 0.0f); glVertex3fv(p[6]);
		glTexCoord2f(1.0f, 0.0f); glVertex3fv(p[5]);
        glEnd();
	glDisable(GL_TEXTURE_2D);
}

//--------------------------------------------------------------------------
//���ʑ�(x>0)�ɕ��s���e
void skTexSphere1(int Nxy)//�P�ʒ��a�̋�
{
	int i, j, Nz;
	double phai; //x-y���ʂɑ΂���Ίp�i�ܓx�j
	double theta; //x���ɑ΂���Ίp�i�o�x)
	float p[21][21][3]; //�e�_�̍��W
	float *p1,*p2,*p3,*p4;

	if(Nxy > 20) Nxy = 20;
	Nz = Nxy;
	//���_�ԍ��A�ʔԍ�
	//���W
	for(i=0;i<=Nxy;i++)
	{	theta = 2.0 * M_PI * (double)i / (double)Nxy;
		for(j=0;j<=Nz;j++)
		{	//���[��i=0�i��������݂ẲE�[�A-y�������j
			phai = M_PI / 2.0 - M_PI * (double)j / (double)Nz;
			p[i][j][0] = (float)(0.5 * sin(theta) * cos(phai)); //x���W
			p[i][j][1] = -(float)(0.5 * cos(theta) * cos(phai)); //y
			p[i][j][2] = (float)(0.5 * sin(phai));            //z
		}
	}

	//���s���e
	//�\��(x>0)��ϯ��ݸ�
	glEnable(GL_TEXTURE_2D);
	for(i = 0;i < Nxy/2;i++){
		for(j = 0;j < Nz; j++)//-1;j++)
		{
			p1 = p[i][j];	  p2 = p[i][j+1];
			p3 = p[i+1][j+1]; p4 = p[i+1][j];
			glBegin(GL_QUADS);//���[(p1[1]=-0.5)��ø�������W��0�ƂȂ�悤�ɂ���
				glNormal3fv(p1);glTexCoord2f((0.5f+p1[1]),(0.5f+p1[2]));glVertex3fv(p1);
				glNormal3fv(p2);glTexCoord2f((0.5f+p2[1]),(0.5f+p2[2]));glVertex3fv(p2);
				glNormal3fv(p3);glTexCoord2f((0.5f+p3[1]),(0.5f+p3[2]));glVertex3fv(p3);
				glNormal3fv(p4);glTexCoord2f((0.5f+p4[1]),(0.5f+p4[2]));glVertex3fv(p4);
			glEnd();
		}
	}
	glDisable(GL_TEXTURE_2D);
	//���ʂɂ�ϯ��ݸނ��Ȃ�
	for(i=Nxy/2;i<Nxy;i++){
		for(j=0;j<Nz;j++)//-1;j++)
		{
			p1=p[i][j];
			p2=p[i][j+1];
			p3=p[i+1][j+1];
			p4=p[i+1][j];

			glBegin(GL_QUADS);
				glNormal3fv(p1);glVertex3fv(p1);
				glNormal3fv(p2);glVertex3fv(p2);
				glNormal3fv(p3);glVertex3fv(p3);
				glNormal3fv(p4);glVertex3fv(p4);
			glEnd();
		}
	}
}
//--------------------------------------------------------------------------
//�����e
void skTexSphere2(int Nxy)//�P�ʒ��a�̋�
{
	int i, ii, j, Nz;
	double phai; //x-y���ʂɑ΂���Ίp�i�ܓx�j
	double theta; //x���ɑ΂���Ίp�i�o�x)
	float p[21][21][3]; //�e�_�̍��W
	float *p1,*p2,*p3,*p4;
	float th1, th2, ph1, ph2;

	if(Nxy > 20) Nxy = 20;
	Nz = Nxy;
	//���_�ԍ��A�ʔԍ�
	//���W
	for(i=0;i<=Nxy;i++)
	{	theta = 2.0*M_PI*(double)i/(double)Nxy;
		for(j=0;j<=Nz;j++)
		{
			phai = M_PI/2.0-M_PI*(double)j/(double)Nz;
          	//�^����i=0
			p[i][j][0] = -(float)(0.5*cos(theta)*cos(phai)); //x���W
			p[i][j][1] = -(float)(0.5*sin(theta)*cos(phai)); //y
			p[i][j][2] = (float)(0.5*sin(phai));             //z
		}
	}

	//�����e
	//�S�̂�ϯ��ݸ�
	glEnable(GL_TEXTURE_2D);
	for(i = 0;i < Nxy;i++)
	{
		//if(i == Nxy-1) ii = 0; else ii = i+1;
		ii = i + 1;
		th1 = (float)i / (float)Nxy;//2�΂Ő��K�������p�x
		th2 = (float)ii / (float)Nxy;
		if(i == Nxy-1) th2 = 1.0;
		for(j = 0;j < Nz; j++)
		{
			ph1 = 1.0f - (float)j / (float)Nz;//j=0��1�ƂȂ�悤��
			ph2 = 1.0f - (float)(j+1) / (float)Nz;
			p1 = p[i][j];	 p2 = p[i][j+1];
			p3 = p[ii][j+1]; p4 = p[ii][j];
			glBegin(GL_QUADS);//���[(p1[1]=-0.5)��ø�������W��0�ƂȂ�悤�ɂ���
				glNormal3fv(p1);glTexCoord2f(th1, ph1);glVertex3fv(p1);
				glNormal3fv(p2);glTexCoord2f(th1, ph2);glVertex3fv(p2);
				glNormal3fv(p3);glTexCoord2f(th2, ph2);glVertex3fv(p3);
				glNormal3fv(p4);glTexCoord2f(th2, ph1);glVertex3fv(p4);
			glEnd();
		}
	}
	glDisable(GL_TEXTURE_2D);
}

//----------------------------------------------------------------------
//���ʐ��ʑ�(x>0)�ɕ��s���e
void skTexCylinder1(int Nxy)
{	//���ards�A����hgt�̉~��(���S�����_)
	//Nx--���p�`�̎��(Nx<0)
	float p[40][3]; //���_���W
	double theta0,theta;
	double th1,th2;
	int i,ii;
	float rds = 0.5;//���a
	float hgt = 1.0;//����
//	double pp = 2.0 * M_PI;

	if(Nxy > 20) { Nxy = 20;}
	theta0 = M_PI / (double)Nxy;
	for(i = 0;i < Nxy;i++)
	{   theta = 2.0 * theta0 * (double)i;
		//���[��i=0�Ƃ���i��޼ު�Ď��g���猩�ĉE�[�A-y������)
		p[i][0] = (float)(rds*sin(theta)); //����x����
		p[i][1] = -(float)(rds*cos(theta));//������
		p[i][2] = (float)hgt/2.0f;         //������(����)
		p[i+Nxy][0] = p[i][0];             //�����x����
		p[i+Nxy][1] = p[i][1];             //������
		p[i+Nxy][2] = -(float)hgt/2.0f;    //������
	}

	//���0
	glBegin(GL_POLYGON);
		glNormal3f(0.0f,0.0f,1.0f);
		for(i=0;i<Nxy;i++) glVertex3fv(p[i]);
	glEnd();
	//����
	glBegin(GL_POLYGON);
		glNormal3f(0.0f,0.0f,-1.0f);
		for(i=2*Nxy-1;i>=Nxy;i--) glVertex3fv(p[i]);
	glEnd();

	//���s���e
	//����(x>0�̑��ʔ����ɂ���ø����ϯ��ݸ�:���s���e�j
	glEnable(GL_TEXTURE_2D);
	for(i=0;i<Nxy/2;i++)
	{
		if(i == Nxy-1) ii = 0; else ii = i+1;
		th1 = 2.0*M_PI*(double)i/(double)Nxy;
		th2 = 2.0*M_PI*(double)ii/(double)Nxy;
		glBegin(GL_QUADS);
			glNormal3f((float)sin(th1),-(float)cos(th1),0.0f);
			glTexCoord2f((rds+p[i][1])/(2.0f*rds),1.0f);//ø�������W��i=0�̂Ƃ�s=0�ƂȂ�悤��
													//i=Nxy(�E�[)��rTex�ƂȂ�悤�ɒ��a2.0*rds�Ő��K��
			glVertex3fv(p[i]);
			glTexCoord2f((rds+p[i][1])/(2.0f*rds),0.0f);
			glVertex3fv(p[i+Nxy]);

			glNormal3f((float)sin(th2),-(float)cos(th2),0.0f);
			glTexCoord2f((rds+p[ii][1])/(2.0f*rds),0.0f);
			glVertex3fv(p[ii+Nxy]);

			glTexCoord2f((rds+p[ii][1])/(2.0f*rds),1.0f);
			glVertex3fv(p[ii]);
		glEnd();
	}
	glDisable(GL_TEXTURE_2D);
	//����(x < 0)�ɂ�ϯ��ݸނ��Ȃ�
	for(i = Nxy/2;i < Nxy;i++)
	{
		if(i == Nxy-1) ii = 0; else ii = i+1;
		th1 = 2.0*M_PI*(double)i/(double)Nxy;
		th2 = 2.0*M_PI*(double)ii/(double)Nxy;
		glBegin(GL_QUADS);
			glNormal3f((float)sin(th1),-(float)cos(th1),0.0f);
			glVertex3fv(p[i]); glVertex3fv(p[i+Nxy]);
			glNormal3f((float)sin(th2),-(float)cos(th2),0.0f);
			glVertex3fv(p[ii+Nxy]); glVertex3fv(p[ii]);
		glEnd();
	}
}
//----------------------------------------------------------------------
//���ʂ����ɉ~�����e
void skTexCylinder2(int Nxy)
{
	//���ards�A����hgt�̉~��(���S�����_)
	//Nx--���p�`�̎��(Nx<0)
	float p[40][3]; //���_���W
	double theta0,theta;
	double th1,th2;
	int i,ii;
	float rds = 0.5;//���a
	float hgt = 1.0;//����
	double pp = 2.0 * M_PI;

	if(Nxy > 20) { Nxy = 20;}
	theta0 = M_PI / (double)Nxy;
	for(i = 0;i < Nxy;i++)
	{   theta = 2.0 * theta0 * (double)i;
		//���[��i=0�Ƃ���i��޼ު�Ď��g���猩�ĉE�[�A-y������)
//		p[i][0] = (float)(rds*sin(theta)); //����x����
//		p[i][1] = -(float)(rds*cos(theta));//������
		//�w�ʂ�i=0�Ƃ���
		p[i][0] = -(float)(rds*cos(theta)); //����x����
		p[i][1] = -(float)(rds*sin(theta));//������
		p[i][2] = (float)hgt/2.0f;         //������(����)
		p[i+Nxy][0] = p[i][0];             //�����x����
		p[i+Nxy][1] = p[i][1];             //������
		p[i+Nxy][2] = -(float)hgt/2.0f;    //������
	}

	//���0
	glBegin(GL_POLYGON);
		glNormal3f(0.0f,0.0f,1.0f);
		for(i=0;i<Nxy;i++) glVertex3fv(p[i]);
	glEnd();
	//����
	glBegin(GL_POLYGON);
		glNormal3f(0.0f,0.0f,-1.0f);
		for(i=2*Nxy-1;i>=Nxy;i--) glVertex3fv(p[i]);
	glEnd();

	//�~�����e
	//���ʈ����ϯ��ݸ�
	glEnable(GL_TEXTURE_2D);
	for(i = 0;i < Nxy;i++)
	{
		if(i == Nxy-1) ii = 0;
        else ii = i + 1;
		th1 = (double)i/(double)Nxy;//2�΂Ő��K�������p�x
		th2 = (double)ii/(double)Nxy;
		if(i == Nxy-1) th2 = 1.0;
		glBegin(GL_QUADS);
			glNormal3f(-(float)cos(pp * th1),-(float)sin(pp * th1),0.0f);
			glTexCoord2f((float)th1, 1.0f);//ø�������W��i=0�̂Ƃ�s=0�ƂȂ�悤��
			glVertex3fv(p[i]);
			glTexCoord2f((float)th1, 0.0f);
			glVertex3fv(p[i+Nxy]);

			glNormal3f(-(float)cos(pp * th2),-(float)sin(pp * th2),0.0f);
			glTexCoord2f((float)th2, 0.0f);
			glVertex3fv(p[ii+Nxy]);
			glTexCoord2f((float)th2, 1.0f);
			glVertex3fv(p[ii]);
		glEnd();
	}
	glDisable(GL_TEXTURE_2D);
}

//----------------------------------------------------------------------
//�i�q��I�u�W�F�N�g�ɑ΂���e�N�X�`���}�b�s���O�p�T�u���[�`��
void skTexGrid(int N1, int N2, float pd[][3])
{	//���ʂ���`
	int i, j, np;
	float n1[3], n2[3], n3[3], n4[3];
	float a[101][101], b[101][101], c[101][101];

	//�@������
	for(i = 0;i <= N1;i++)
	  for(j = 0;j <= N2;j++)
	  {
		np = j * (N1+1) + i;
		if(j == 0 )
		{
			if(i == 0) {
				calcNormal(pd[0],pd[N1+1],pd[1],n1);
				a[i][j] = n1[0]; b[i][j] = n1[1]; c[i][j] = n1[2]; }
			else if(i == N1) {
				calcNormal(pd[N1-1],pd[2*N1+1],pd[N1],n1);
				a[i][j] = n1[0]; b[i][j] = n1[1]; c[i][j] = n1[2]; }
			else {
				calcNormal(pd[np],pd[np-1],pd[np+N1+1],n1);//����
				calcNormal(pd[np],pd[np+N1+1],pd[np+1],n2);//�E��
				a[i][j] = (n1[0]+n2[0])/2.0f;
				b[i][j] = (n1[1]+n2[1])/2.0f;
				c[i][j] = (n1[2]+n2[2])/2.0f; }
		}
		else if(j == N2)
		{
			if(i == 0) {
				calcNormal(pd[np],pd[np+1],pd[np-N1-1],n1);
				a[i][j] = n1[0]; b[i][j] = n1[1]; c[i][j] = n1[2]; }
			else if(i == N1) {
				calcNormal(pd[np],pd[np-N1-1],pd[np-1],n1);
				a[i][j] = n1[0]; b[i][j] = n1[1]; c[i][j] = n1[2]; }
			else {
				calcNormal(pd[np],pd[np-N1-1],pd[np-1],n1);//����
				calcNormal(pd[np],pd[np+1],pd[np-N1-1],n2);//�E��
				a[i][j] = (n1[0]+n2[0])/2.0f;
				b[i][j] = (n1[1]+n2[1])/2.0f;
				c[i][j] = (n1[2]+n2[2])/2.0f; }
		}
		else
		{
			if(i == 0) {
				calcNormal(pd[np],pd[np+1],pd[np-N1-1],n1);//��
				calcNormal(pd[np],pd[np+N1+1],pd[np+1],n2);//��
				a[i][j] = (n1[0]+n2[0])/2.0f;
				b[i][j] = (n1[1]+n2[1])/2.0f;
				c[i][j] = (n1[2]+n2[2])/2.0f; }
			else if(i == N1) {
				calcNormal(pd[np],pd[np-N1-1],pd[np-1],n1);//��
				calcNormal(pd[np],pd[np-1],pd[np+N1+1],n2);//��
				a[i][j] = (n1[0]+n2[0])/2.0f;
				b[i][j] = (n1[1]+n2[1])/2.0f;
				c[i][j] = (n1[2]+n2[2])/2.0f; }
			else {//�㉺���E�S�̎O�p�`�̕���
				calcNormal(pd[np],pd[np-N1-1],pd[np-1],n1);//����
				calcNormal(pd[np],pd[np+1],pd[np-N1-1],n2);//�E��
				calcNormal(pd[np],pd[np-1],pd[np+N1+1],n3);//����
				calcNormal(pd[np],pd[np+N1+1],pd[np+1],n4);//�E��
				a[i][j] = (n1[0]+n2[0]+n3[0]+n4[0])/4.0f;
				b[i][j] = (n1[1]+n2[1]+n3[1]+n4[1])/4.0f;
				c[i][j] = (n1[2]+n2[2]+n3[2]+n4[2])/4.0f; }
		}
	}

	//Texture�`��
	glEnable(GL_TEXTURE_2D);
	for(j = 0;j < N2;j++)
	  for(i = 0;i < N1;i++){
		np = i + (N1+1) * j;
		glBegin(GL_TRIANGLES);
			//�����̎O�p�`
			//�e���_�̖@������,ø�������W,���_���W��^����B
			glNormal3d(a[i][j],b[i][j],c[i][j]);//�@������
			glTexCoord2f((float)i/(float)N1,(1.0f-(float)j/(float)N2));//ø�������W
			glVertex3fv(pd[np]);//��غ�݂̒��_���W�i�ȉ��������J��Ԃ��j
			glNormal3d(a[i][j+1],b[i][j+1],c[i][j+1]);
			glTexCoord2f((float)i/(float)N1,(1.0f-(float)(j+1)/(float)N2));
			glVertex3fv(pd[np+N1+1]);
			glNormal3d(a[i+1][j+1],b[i+1][j+1],c[i+1][j+1]);
			glTexCoord2f((float)(i+1)/(float)N1,(1.0f-(float)(j+1)/(float)N2));
			glVertex3fv(pd[np+N1+2]);
			//�E��̎O�p�`
			glNormal3d(a[i][j],b[i][j],c[i][j]);
			glTexCoord2f((float)i/(float)N1,(1.0f-(float)j/(float)N2));
			glVertex3fv(pd[np]);
			glNormal3d(a[i+1][j+1],b[i+1][j+1],c[i+1][j+1]);
			glTexCoord2f((float)(i+1)/(float)N1,(1.0f-(float)(j+1)/(float)N2));
			glVertex3fv(pd[np+N1+2]);
			glNormal3d(a[i+1][j],b[i+1][j],c[i+1][j]);
			glTexCoord2f((float)(i+1)/(float)N1,(1.0f-(float)j/(float)N2));
			glVertex3fv(pd[np+1]);
		glEnd();
	}
	//���ʂ��\�ʂɕ`����ø������`��
	for(j = 0;j < N2;j++)
	  for(i = 0;i < N1;i++){
		np = i+(N1+1)*j;
		glBegin(GL_TRIANGLES);
			//�����̎O�p�`
			//�e���_�̖@�������ɕ���,ø�������W,���_���W��^����i�\�Ƃ͋t�̒��_��,2�Ԗڂ�3�Ԗڂ������j�B
			glNormal3d(-a[i][j],-b[i][j],-c[i][j]);//�@������
			glTexCoord2f((float)i/(float)N1,(1.0f-(float)j/(float)N2));//ø�������W
			glVertex3fv(pd[np]);//��غ�݂̒��_���W�i�ȉ��������J��Ԃ��j
			glNormal3d(-a[i+1][j+1],-b[i+1][j+1],-c[i+1][j+1]);
			glTexCoord2f((float)(i+1)/(float)N1,(1.0f-(float)(j+1)/(float)N2));
			glVertex3fv(pd[np+N1+2]);
			glNormal3d(-a[i][j+1],-b[i][j+1],-c[i][j+1]);
			glTexCoord2f((float)i/(float)N1,(1.0f-(float)(j+1)/(float)N2));
			glVertex3fv(pd[np+N1+1]);
			//�E��̎O�p�`
			glNormal3d(-a[i][j],-b[i][j],-c[i][j]);
			glTexCoord2f((float)i/(float)N1,(1.0f-(float)j/(float)N2));
			glVertex3fv(pd[np]);
			glNormal3d(-a[i+1][j],-b[i+1][j],-c[i+1][j]);
			glTexCoord2f((float)(i+1)/(float)N1,(1.0f-(float)j/(float)N2));
			glVertex3fv(pd[np+1]);
			glNormal3d(-a[i+1][j+1],-b[i+1][j+1],-c[i+1][j+1]);
			glTexCoord2f((float)(i+1)/(float)N1,(1.0f-(float)(j+1)/(float)N2));
			glVertex3fv(pd[np+N1+2]);
		glEnd();
	}

	glDisable(GL_TEXTURE_2D);
}
//-----------------------------------------------------------------------------
//�����`�̃��b�V���ix-y����,���S�����_�j
//���ʂ̔g���ȂǂɓK�p�i�U���Ɋ֌W�Ȃ�x-y�������̒����͈��j
void skTexMesh1(int Nx, float* data)
{
	//�S�̂̕�,�����ǂ����1(�����`)
	int i, j, Ny, np;
	float p[10201][3]; //���_���W
	float pitchX, pitchY;


	if(Nx > 100) Nx = 100;
	Ny = Nx;
	//�Z���̃T�C�Y
	pitchX = 1.0f / (float)Nx;
	pitchY = 1.0f / (float)Ny;

	//�e���_�̍��W
	for(j = 0; j <= Ny; j++){
		for(i = 0; i <= Nx; i++){
			np = j * (Ny + 1) + i;
			p[np][0] = (float)(j - Nx / 2) * pitchX ;
			p[np][1] = (float)(i - Ny / 2) * pitchY ;/// length;
			p[np][2] = data[np];
		}
	}
	skTexGrid(Nx, Ny, p);
}

//-----------------------------------------------------------------------------
//�l�p���z�n�ix-y����,���S�����_�j
//Mesh���g�̐L�яk�݂Ȃ��i�U�����傫���Ȃ��x-y�������Z���Ȃ�j
//�U���͂�������
void skTexMesh2(int Nx, double ratio, double amp,
		       double lamdaX, double lamdaY, double alpha)
{
	float pd[10201][3];//���_���W
	double lengthX, lengthY;
	double thetaX, thetaY;
	double cx, cy, ex, ey;
	double x, y, z, dx, dy;
	int i, j, np, Ny;

	if(Nx > 100) Nx = 100;
	Ny = Nx;

	//�S��
	lengthX = 1.0;
	lengthY = lengthX * ratio;//x�����̒���

	//��{������

	//1�v�f�̻���(�O���b�h�Ԋu)
	ex = lengthX / (double)Nx;
	ey = lengthY / (double)Ny;

	cx = 2.0 * M_PI / lamdaX;
	cy = 2.0 * M_PI / lamdaY;

//	alpha = (float)(alpha * M_PI/180.0);//�����ʑ���rad�ɕϊ�
	//���W�l
	x = 0.0; y = 0.0;
	for(j = 0; j <= Nx; j++){//j��-x����
		thetaX = atan(amp * cx * cos(cx * x + cy * y - alpha));
		dx = ex * cos(thetaX);//����������
		x += dx;
		y = 0.0;
		for(i = 0; i <= Ny; i++){//i��y����
			thetaY = atan(amp * cy * cos(cx * x + cy * y - alpha));
			dy = ey * cos(thetaY);//y��������
			y += dy;
			z = amp * sin(cx * x + cy * y - alpha);//*(y+0.5);
			np = j * (Ny + 1) + i;
			pd[np][0] = (float)x;
			pd[np][1] = (float)y;
			pd[np][2] = (float)z;
		}
	}
	skTexGrid(Nx, Ny, pd);
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//�l�p���z�n�ix-z����,Flag�̍��������_�j
//Mesh���g�̐L�яk�݂Ȃ��i�U�����傫���Ȃ��x�������Z�k�j
//���ɂ�鉡�����L�k����
//�ψʂ�y������
void skTexMesh3(int Nx, double ratio, double amp0, double lambda0, double alpha,
                double wind, double z0)
{
	float pd[10201][3];//���_���W
	double lengthX, lengthZ;
	double thetaX, thetaZ, ampl, fct, xShift;
	double cx, cz, ex, ez;
	double x, y, z, dx, dz, y0;
	int i, j, np, Nz;

    double amp = amp0 / (1.0f + 50.0 * wind );//���̋����Ƃ��U����������
    double lambdaX = lambda0 / (1.0f + 5.0 * wind);//���̋����Ƃ�x�����g���͒Z��
    double lambdaZ = lambda0*(1.0 + 5.0 * wind);        //�������g���͒���
    double weight =  1.0f / (1.0 + 20.0 * wind); //�d�͏d�ݔ䗦
	double gamma = weight * (80.0 * M_PI /180.0); //weight=1�őS�̂�80�x�X��

	if(Nx > 100) Nx = 100;
	Nz = Nx;

	//�S��
	lengthZ = 1.0;//����(z����)
	lengthX = lengthZ * ratio;//����(x����)

	//1�v�f�̻���(�O���b�h�Ԋu)
	ex = lengthX / (double)Nx;
	ez = lengthZ / (double)Nz;

	cx = (2.0 * M_PI) / lambdaX;
	cz = (2.0 * M_PI) / lambdaZ;
//cz=0.0;
	//���W�l
	z = 1.0; x = 0.0;
	for(j = 0;j <= Nz;j++){//j��-z����
		x = 0.0;
        xShift = wind * (0.3 * sin(M_PI*z)) ;
//		ampl = amp * (1.0 + x);
		y0 = amp * sin(cz * z - alpha) ;//pole�ɂ����錳�̈ʒu
//		y0 = amp * sin(alpha - cz * z) ;//pole�ɂ����錳�̈ʒu
        fct = 1.0 + (1.0 - sin(M_PI * z) ) * wind * 0.2;//�������T�C�Y�W��
		for(i = 0;i <= Nx;i++){//i��x����
			ampl = amp * (0.f + x );//x���傫���قǐU����傫���Ȃ�悤��
            //�ψ�(x=0�̈ʒu��pole�ɌŒ�ł���悤��y0�ō�������)
			y = ( ampl* sin(  cx * x * fct + cz * z  - alpha) - y0 ) ;
			np = j * (Nx+1) + i;
			pd[np][1] = (float)y;
			//�d�͂̉e��
			pd[np][0] = (float)((x * cos(gamma) + xShift) * fct) ;
			pd[np][2] = (float)(z - x * sin(gamma));
            if(pd[np][2] < -(float)z0 + 0.01f) pd[np][2] = - (float)z0 + 0.01f;//���ʂ�蕂������
            //���̊i�q�_��x���W
			thetaX = atan(ampl * cx * cos( cx * x + cz * z - alpha));
//			thetaX = atan(ampl * cx * cos(alpha - cx * x - cz * z));
			dx = ex * cos(thetaX) ;//x��������
			x += dx;

		}
		thetaZ = atan(amp * cz * cos(cx * x + cz * z - alpha));
		dz = ez * cos(thetaZ);//����������
		z -= dz;
	}
	skTexGrid(Nx, Nz, pd);
}

//-----------------------------------------------------------------
//��2���֐�
void skTexSuper1(int Nxy, int Nz, double eps1, double eps2, double p1, double p2, double p3)
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
	//����(4�p�`�p�b�`�j
	float th1, th2, ph1, ph2;
	glEnable(GL_TEXTURE_2D);
	for(j = 0;j < Nz;j++){
		ph1 = 1.0f - (float)j / (float)Nz;
		ph2 = 1.0f - (float)(j + 1) / (float)Nz;
		for(i = 0;i < Nxy;i++)
		{	ip = i+1;
			if(ip == Nxy) ip = 0;
			th1 = (float)i / (float)Nxy;
			th2 = (float)ip / (float)Nxy;
			if(i == Nxy-1) th2 = 1.0f;
			glBegin(GL_QUADS);
				glNormal3d(a[i][j],b[i][j],c[i][j]);
				glTexCoord2f(th1, ph1);glVertex3fv(pd[j*Nxy+i]);
				glNormal3d(a[i][j+1],b[i][j+1],c[i][j+1]);
				glTexCoord2f(th1, ph2);glVertex3fv(pd[(j+1)*Nxy+i]);
				glNormal3d(a[ip][j+1],b[ip][j+1],c[ip][j+1]);
				glTexCoord2f(th2, ph2);glVertex3fv(pd[(j+1)*Nxy+ip]);
				glNormal3d(a[ip][j],b[ip][j],c[ip][j]);
				glTexCoord2f(th2, ph1);glVertex3fv(pd[j*Nxy+ip]);
			glEnd();
		}
	}
	glDisable(GL_TEXTURE_2D);
}
//----------------------------------------------------------------------
//�Q�����i�q��I�u�W�F�N�g�ɑ΂���e�N�X�`���}�b�s���O�p�T�u���[�`��
void skTexGridSquare(int numRow, int numCol, float pd[][3])
{	//���ʂ���`
	int i, j, np;
	float n1[3], n2[3], n3[3], n4[3];
	float a[20][20], b[20][20], c[20][20];

	//�@������
	for(i = 0;i < numRow;i++) //������
	  for(j = 0;j < numCol;j++)//�E����
	  {
		//np = j * (N1+1) + i;
        np = i * numCol + j;
		if(i == 0 )
		{
			if(j == 0) {
				calcNormal(pd[0],pd[numCol],pd[1],n1);
				a[i][j] = n1[0]; b[i][j] = n1[1]; c[i][j] = n1[2]; }
			else if(j == numCol-1) {
				calcNormal(pd[numCol-2],pd[2*numCol-1],pd[numCol-1],n1);
				a[i][j] = n1[0]; b[i][j] = n1[1]; c[i][j] = n1[2]; }
			else {
				calcNormal(pd[np],pd[np-1],pd[np+numCol],n1);//����
				calcNormal(pd[np],pd[np+numCol],pd[np+1],n2);//�E��
				a[i][j] = (n1[0]+n2[0])/2.0f;
				b[i][j] = (n1[1]+n2[1])/2.0f;
				c[i][j] = (n1[2]+n2[2])/2.0f; }
		}
		else if(i == numRow-1)
		{
			if(j == 0) {
				calcNormal(pd[np],pd[np+1],pd[np-numCol],n1);
				a[i][j] = n1[0]; b[i][j] = n1[1]; c[i][j] = n1[2]; }
			else if(j == numCol-1) {
				calcNormal(pd[np],pd[np-numCol],pd[np-1],n1);
				a[i][j] = n1[0]; b[i][j] = n1[1]; c[i][j] = n1[2]; }
			else {
				calcNormal(pd[np],pd[np-numCol],pd[np-1],n1);//����
				calcNormal(pd[np],pd[np+1],pd[np-numCol],n2);//�E��
				a[i][j] = (n1[0]+n2[0])/2.0f;
				b[i][j] = (n1[1]+n2[1])/2.0f;
				c[i][j] = (n1[2]+n2[2])/2.0f; }
		}
		else
		{
			if(j == 0) {
				calcNormal(pd[np],pd[np+1],pd[np-numCol],n1);//��
				calcNormal(pd[np],pd[np+numCol],pd[np+1],n2);//��
				a[i][j] = (n1[0]+n2[0])/2.0f;
				b[i][j] = (n1[1]+n2[1])/2.0f;
				c[i][j] = (n1[2]+n2[2])/2.0f; }
			else if(j == numCol-1) {
				calcNormal(pd[np],pd[np-numCol],pd[np-1],n1);//��
				calcNormal(pd[np],pd[np-1],pd[np+numCol],n2);//��
				a[i][j] = (n1[0]+n2[0])/2.0f;
				b[i][j] = (n1[1]+n2[1])/2.0f;
				c[i][j] = (n1[2]+n2[2])/2.0f; }
			else {//�㉺���E�S�̎O�p�`�̕���
				calcNormal(pd[np],pd[np-numCol],pd[np-1],n1);//����
				calcNormal(pd[np],pd[np+1],pd[np-numCol],n2);//�E��
				calcNormal(pd[np],pd[np-1],pd[np+numCol],n3);//����
				calcNormal(pd[np],pd[np+numCol],pd[np+1],n4);//�E��
				a[i][j] = (n1[0]+n2[0]+n3[0]+n4[0])/4.0f;
				b[i][j] = (n1[1]+n2[1]+n3[1]+n4[1])/4.0f;
				c[i][j] = (n1[2]+n2[2]+n3[2]+n4[2])/4.0f; }
		}
	}

	//Texture�`��
	glEnable(GL_TEXTURE_2D);
	for(i = 0;i < numRow-1;i++)
	  for(j = 0;j < numCol-1;j++){
		//np = i + (N1+1) * j;
        np = i * numCol + j;
		glBegin(GL_TRIANGLES);
			//�����̎O�p�`
			//�e���_�̖@������,ø�������W,���_���W��^����B
			glNormal3d(a[i][j],b[i][j],c[i][j]);//�@������
			glTexCoord2f((float)j/(float)(numCol-1),(1.0f-(float)i/(float)(numRow-1)));//ø�������W
			glVertex3fv(pd[np]);//��غ�݂̒��_���W�i�ȉ��������J��Ԃ��j
			glNormal3d(a[i+1][j],b[i+1][j],c[i+1][j]);
			glTexCoord2f((float)j/(float)(numCol-1),(1.0f-(float)(i+1)/(float)(numRow-1)));
			glVertex3fv(pd[np+numCol]);
			glNormal3d(a[i+1][j+1],b[i+1][j+1],c[i+1][j+1]);
			glTexCoord2f((float)(j+1)/(float)(numCol-1),(1.0f-(float)(i+1)/(float)(numRow-1)));
			glVertex3fv(pd[np+numCol+1]);
			//�E��̎O�p�`
			glNormal3d(a[i][j],b[i][j],c[i][j]);
			glTexCoord2f((float)j/(float)(numCol-1),(1.0f-(float)i/(float)(numRow-1)));
			glVertex3fv(pd[np]);
			glNormal3d(a[i+1][j+1],b[i+1][j+1],c[i+1][j+1]);
			glTexCoord2f((float)(j+1)/(float)(numCol-1),(1.0f-(float)(i+1)/(float)(numRow-1)));
			glVertex3fv(pd[np+numCol+1]);
			glNormal3d(a[i][j+1],b[i][j+1],c[i][j+1]);
			glTexCoord2f((float)(j+1)/(float)(numCol-1),(1.0f-(float)i/(float)(numRow-1)));
			glVertex3fv(pd[np+1]);
		glEnd();
	}
	//���ʂ��\�ʂɕ`����ø������`��
	for(i = 0;i < numRow-1;i++)
	  for(j = 0;j < numCol-1;j++){
        np = i * numCol + j;
		glBegin(GL_TRIANGLES);
			//�����̎O�p�`
			//�e���_�̖@�������ɕ���,ø�������W,���_���W��^����i�\�Ƃ͋t�̒��_��,2�Ԗڂ�3�Ԗڂ������j�B
			glNormal3d(-a[i][j],-b[i][j],-c[i][j]);//�@������
			glTexCoord2f((float)j/(float)(numCol-1),(1.0f-(float)i/(float)(numRow-1)));//ø�������W
			glVertex3fv(pd[np]);//��غ�݂̒��_���W�i�ȉ��������J��Ԃ��j
			glNormal3d(-a[i+1][j+1],-b[i+1][j+1],-c[i+1][j+1]);
			glTexCoord2f((float)(j+1)/(float)(numCol-1),(1.0f-(float)(i+1)/(float)(numRow-1)));
			glVertex3fv(pd[np+numCol+1]);
			glNormal3d(-a[i+1][j],-b[i+1][j],-c[i+1][j]);
			glTexCoord2f((float)j/(float)(numCol-1),(1.0f-(float)(i+1)/(float)(numRow-1)));
			glVertex3fv(pd[np+numCol]);
			//�E��̎O�p�`
			glNormal3d(-a[i][j],-b[i][j],-c[i][j]);
			glTexCoord2f((float)j/(float)(numCol-1),(1.0f-(float)i/(float)(numRow-1)));
			glVertex3fv(pd[np]);
			glNormal3d(-a[i][j+1],-b[i][j+1],-c[i][j+1]);
			glTexCoord2f((float)(j+1)/(float)(numCol-1),(1.0f-(float)i/(float)(numRow-1)));
			glVertex3fv(pd[np+1]);
			glNormal3d(-a[i+1][j+1],-b[i+1][j+1],-c[i+1][j+1]);
			glTexCoord2f((float)(j+1)/(float)(numCol-1),(1.0f-(float)(i+1)/(float)(numRow-1)));
			glVertex3fv(pd[np+numCol+1]);
		glEnd();
	}

	glDisable(GL_TEXTURE_2D);
}
//------------------------------------------------------------------------
void skTexGrid2(int numRow, int numCol, float pd[][3], bool flagShadow)
{
	int i, j, np;
	float n1[3], n2[3], n3[3], n4[3];
	float a[21][21], b[21][21], c[21][21];
    float shadowDiffuse[] = {0.2f,0.2f,0.2f,0.3f};//�e�̊g�U��

	//�@������
	for(i = 0;i < numRow;i++) //������
	  for(j = 0;j < numCol;j++)//�E����
	  {
        np = i * numCol + j;
		if(i == 0 )
		{
			if(j == 0) {
				calcNormal(pd[0],pd[numCol],pd[1],n1);
				a[i][j] = n1[0]; b[i][j] = n1[1]; c[i][j] = n1[2]; }
			else if(j == numCol-1) {
				calcNormal(pd[numCol-2],pd[2*numCol-1],pd[numCol-1],n1);
				a[i][j] = n1[0]; b[i][j] = n1[1]; c[i][j] = n1[2]; }
			else {
				calcNormal(pd[np],pd[np-1],pd[np+numCol],n1);//����
				calcNormal(pd[np],pd[np+numCol],pd[np+1],n2);//�E��
				a[i][j] = (n1[0]+n2[0])/2.0f;
				b[i][j] = (n1[1]+n2[1])/2.0f;
				c[i][j] = (n1[2]+n2[2])/2.0f; }
		}
		else if(i == numRow-1)
		{
			if(j == 0) {
				calcNormal(pd[np],pd[np+1],pd[np-numCol],n1);
				a[i][j] = n1[0]; b[i][j] = n1[1]; c[i][j] = n1[2]; }
			else if(j == numCol-1) {
				calcNormal(pd[np],pd[np-numCol],pd[np-1],n1);
				a[i][j] = n1[0]; b[i][j] = n1[1]; c[i][j] = n1[2]; }
			else {
				calcNormal(pd[np],pd[np-numCol],pd[np-1],n1);//����
				calcNormal(pd[np],pd[np+1],pd[np-numCol],n2);//�E��
				a[i][j] = (n1[0]+n2[0])/2.0f;
				b[i][j] = (n1[1]+n2[1])/2.0f;
				c[i][j] = (n1[2]+n2[2])/2.0f; }
		}
		else
		{
			if(j == 0) {
				calcNormal(pd[np],pd[np+1],pd[np-numCol],n1);//��
				calcNormal(pd[np],pd[np+numCol],pd[np+1],n2);//��
				a[i][j] = (n1[0]+n2[0])/2.0f;
				b[i][j] = (n1[1]+n2[1])/2.0f;
				c[i][j] = (n1[2]+n2[2])/2.0f; }
			else if(j == numCol-1) {
				calcNormal(pd[np],pd[np-numCol],pd[np-1],n1);//��
				calcNormal(pd[np],pd[np-1],pd[np+numCol],n2);//��
				a[i][j] = (n1[0]+n2[0])/2.0f;
				b[i][j] = (n1[1]+n2[1])/2.0f;
				c[i][j] = (n1[2]+n2[2])/2.0f; }
			else {//�㉺���E�S�̎O�p�`�̕���
				calcNormal(pd[np],pd[np-numCol],pd[np-1],n1);//����
				calcNormal(pd[np],pd[np+1],pd[np-numCol],n2);//�E��
				calcNormal(pd[np],pd[np-1],pd[np+numCol],n3);//����
				calcNormal(pd[np],pd[np+numCol],pd[np+1],n4);//�E��
				a[i][j] = (n1[0]+n2[0]+n3[0]+n4[0])/4.0f;
				b[i][j] = (n1[1]+n2[1]+n3[1]+n4[1])/4.0f;
				c[i][j] = (n1[2]+n2[2]+n3[2]+n4[2])/4.0f; }
		}
	}
	//Texture�`��
    if(flagShadow == true) goto SHADOW_DISP;
	glEnable(GL_TEXTURE_2D);
	for(i = 0;i < numRow-1;i++)
	  for(j = 0;j < numCol-1;j++){
        np = i * numCol + j;
		glBegin(GL_TRIANGLES);
			//�����̎O�p�`
			//�e���_�̖@������,ø�������W,���_���W��^����B
			glNormal3d(a[i][j],b[i][j],c[i][j]);//�@������
			glTexCoord2f((float)j/(float)(numCol-1),(1.0f-(float)i/(float)(numRow-1)));//ø�������W
			glVertex3fv(pd[np]);//��غ�݂̒��_���W�i�ȉ��������J��Ԃ��j
			glNormal3d(a[i+1][j],b[i+1][j],c[i+1][j]);
			glTexCoord2f((float)j/(float)(numCol-1),(1.0f-(float)(i+1)/(float)(numRow-1)));
			glVertex3fv(pd[np+numCol]);
			glNormal3d(a[i+1][j+1],b[i+1][j+1],c[i+1][j+1]);
			glTexCoord2f((float)(j+1)/(float)(numCol-1),(1.0f-(float)(i+1)/(float)(numRow-1)));
			glVertex3fv(pd[np+numCol+1]);
			//�E��̎O�p�`
			glNormal3d(a[i][j],b[i][j],c[i][j]);
			glTexCoord2f((float)j/(float)(numCol-1),(1.0f-(float)i/(float)(numRow-1)));
			glVertex3fv(pd[np]);
			glNormal3d(a[i+1][j+1],b[i+1][j+1],c[i+1][j+1]);
//diffuse[0] = 1.0; diffuse[1] =1.0; diffuse[2] = 0.0;
//glMaterialfv(GL_FRONT,GL_DIFFUSE,diffuse);
    		glTexCoord2f((float)(j+1)/(float)(numCol-1),(1.0f-(float)(i+1)/(float)(numRow-1)));
			glVertex3fv(pd[np+numCol+1]);
			glNormal3d(a[i][j+1],b[i][j+1],c[i][j+1]);
//diffuse[0] = 1.0; diffuse[1] =0.0; diffuse[2] = 1.0;
//glMaterialfv(GL_FRONT,GL_DIFFUSE,diffuse);
			glTexCoord2f((float)(j+1)/(float)(numCol-1),(1.0f-(float)i/(float)(numRow-1)));
			glVertex3fv(pd[np+1]);
		glEnd();
	}
	glDisable(GL_TEXTURE_2D);
    return;

SHADOW_DISP:;
	for(i = 0;i < numRow-1;i++)
	  for(j = 0;j < numCol-1;j++){
        np = i * numCol + j;
        glMaterialfv(GL_FRONT,GL_DIFFUSE,shadowDiffuse);

		glBegin(GL_QUADS);
			glNormal3d(a[i][j],b[i][j],c[i][j]);//�@������
			glVertex3fv(pd[np]);//��غ�݂̒��_���W�i�ȉ��������J��Ԃ��j
			glNormal3d(a[i+1][j],b[i+1][j],c[i+1][j]);
			glVertex3fv(pd[np+numCol]);
			glNormal3d(a[i+1][j+1],b[i+1][j+1],c[i+1][j+1]);
			glVertex3fv(pd[np+numCol+1]);
			glNormal3d(a[i][j+1],b[i][j+1],c[i][j+1]);
			glVertex3fv(pd[np+1]);
		glEnd();
	}

}
//----------------------------------------------------------------------
//�����̊i�q��I�u�W�F�N�g�ɑ΂���e�N�X�`���}�b�s���O�p�T�u���[�`��
void skTexGridCube(int numRow, int numCol, int numStk, float p1[][3],
        float p2[][3], float p3[][3], float p4[][3], float p5[][3],
        float p6[][3], bool flagShadow)
{
    //���ׂē����e�N�X�`���[
    //�㑤
    skTexGrid2(numRow, numCol, p1, flagShadow);
    //����
    skTexGrid2(numStk, numCol, p2, flagShadow);
    //�E��
    skTexGrid2(numStk, numRow, p3, flagShadow);
    //����
    skTexGrid2(numStk, numRow, p4, flagShadow);
    //����
    skTexGrid2(numStk, numCol, p5, flagShadow);
    //����
    skTexGrid2(numRow, numCol, p6, flagShadow);
}
//----------------------------------------------------------------------------
//���i�q��I�u�W�F�N�g�ɑ΂���e�N�X�`���}�b�s���O�p�T�u���[�`��
//2004.07.13����
//�㉺�̒��_���c�����̌��Ɋ܂߂�
void skTexGridSphere(int Nxy, int Nz, float pd[][3], bool flagShadow)
{
	int i,j,ip,im,np,npL,npR,npU,npD;
//	double ct,phai,theta,z,fz;
	float a[10][10],b[10][10],c[10][10];//���_�̖@������
	float n[10][3];//���_�̖@������
    float shadowDiffuse[] = {0.2f,0.2f,0.2f,0.3f};//�e�̊g�U��

    if(Nxy > 10) Nxy = 10;
    if(Nz > 10) Nz = 10;
    //j=0�̂�pd[0]�`pd[Nxy-1]�͓���Top���_
    //j=Nz-1�̂Ƃ�np=(Nz-1)*Nxy�Ƃ���
    //pd[np]�`pd[np+numCol-1]�͓���Bottom���_

    //Top
    //Top��numCol�̎O�p�`�̖@���x�N�g���𕽋ω�
    a[0][0] = 0.0; b[0][0] = 0.0; c[0][0] = 0.0;
	for(i = 0;i < Nxy;i++){
 		ip = i+1;
		if(ip == Nxy) ip = 0;
        calcNormal(pd[i],pd[i+Nxy],pd[ip+Nxy],n[i]);
        a[0][0] += n[i][0];
        b[0][0] += n[i][1];
        c[0][0] += n[i][2];
    }
    a[0][0] = a[0][0] / (float)Nxy;
    b[0][0] = b[0][0] / (float)Nxy;
    c[0][0] = c[0][0] / (float)Nxy;
    for(i = 1;i < Nxy; i++){ //Top���_�͂��ׂē���
        a[i][0] = a[0][0]; b[i][0] = b[0][0]; c[i][0] = c[0][0];
    }

	//���ʂ̖@������
    //�ׂ荇��4�̎O�p�`�̖@���x�N�g���𕽋ω�
	for(i = 0;i < Nxy;i++){
 		ip = i+1;
		if(ip == Nxy) ip = 0;
		im = i-1;
		if(i == 0) im = Nxy-1;
		for(j = 1;j < Nz-1;j++)
		{
			np =  j * Nxy + i;//���ړ_
			npL = j * Nxy + im;//����
			npR = j * Nxy + ip;//�E��
			npU = np - Nxy;//��
			npD = np + Nxy;//��
            calcNormal(pd[np],pd[npU],pd[npL],n[0]);//�O���猩�č���
            calcNormal(pd[np],pd[npR],pd[npU],n[1]);//�E��
            calcNormal(pd[np],pd[npL],pd[npD],n[2]);//�O���猩�č���
            calcNormal(pd[np],pd[npD],pd[npR],n[3]);//�E��

			a[i][j] = (float)((n[0][0]+n[1][0]+n[2][0]+n[3][0])/4.0f);//������
			b[i][j] = (float)((n[0][1]+n[1][1]+n[2][1]+n[3][1])/4.0f);//��
			c[i][j] = (float)((n[0][2]+n[1][2]+n[2][2]+n[3][2])/4.0f);//��
		}
	}
    //Bottom
    //Bottom��numCol�̎O�p�`�̖@���x�N�g���𕽋ω�
    a[0][Nxy-1] = 0.0; b[0][Nxy-1] = 0.0; c[0][Nxy-1] = 0.0;
	for(i = 0;i < Nxy;i++){
        np = (Nz-1) * Nxy +i;
 		ip = i+1;
		if(ip == Nxy) ip = 0;
        calcNormal(pd[np],pd[np-i+ip-Nxy],pd[np-Nxy],n[i]);
        a[0][Nz-1] += n[i][0];
        b[0][Nz-1] += n[i][1];
        c[0][Nz-1] += n[i][2];
    }
    a[0][Nz-1] = a[0][Nz-1] / (float)Nxy;
    b[0][Nz-1] = b[0][Nz-1] / (float)Nxy;
    c[0][Nz-1] = c[0][Nz-1] / (float)Nxy;
    for(i = 1;i < Nxy; i++){ //Top���_�͂��ׂē���
        a[i][Nz-1] = a[0][Nz-1];
        b[i][Nz-1] = b[0][Nz-1];
        c[i][Nz-1] = c[0][Nz-1];
    }

	//�\��
    if(flagShadow == true) goto SHADOW_DISP;
	//����(4�p�`�p�b�`�j
	float th1, th2, ph1, ph2;
	glEnable(GL_TEXTURE_2D);
	for(j = 0;j < Nz-1;j++){
		ph1 = 1.0f - (float)j / (float)(Nz-1);
		ph2 = 1.0f - (float)(j + 1) / (float)(Nz-1);
		for(i = 0;i < Nxy;i++)
		{	ip = i+1;
			if(ip == Nxy) ip = 0;
			th1 = (float)i / (float)Nxy;
			th2 = (float)(i+1) / (float)Nxy;
			//if(i == Nxy-1) th2 = 1.0f;
			glBegin(GL_QUADS);
				glNormal3d(a[i][j],b[i][j],c[i][j]);
				glTexCoord2f(th1, ph1);glVertex3fv(pd[j*Nxy+i]);
				glNormal3d(a[i][j+1],b[i][j+1],c[i][j+1]);
				glTexCoord2f(th1, ph2);glVertex3fv(pd[(j+1)*Nxy+i]);
				glNormal3d(a[ip][j+1],b[ip][j+1],c[ip][j+1]);
				glTexCoord2f(th2, ph2);glVertex3fv(pd[(j+1)*Nxy+ip]);
				glNormal3d(a[ip][j],b[ip][j],c[ip][j]);
				glTexCoord2f(th2, ph1);glVertex3fv(pd[j*Nxy+ip]);
			glEnd();
		}
	}
	glDisable(GL_TEXTURE_2D);
    return;
SHADOW_DISP:;
	for(j = 0;j < Nz-1;j++){
		for(i = 0;i < Nxy;i++)
		{	ip = i+1;
			if(ip == Nxy) ip = 0;
            glMaterialfv(GL_FRONT,GL_DIFFUSE,shadowDiffuse);
			glBegin(GL_QUADS);
				glNormal3d(a[i][j],b[i][j],c[i][j]);
                glVertex3fv(pd[j*Nxy+i]);
				glNormal3d(a[i][j+1],b[i][j+1],c[i][j+1]);
				glVertex3fv(pd[(j+1)*Nxy+i]);
				glNormal3d(a[ip][j+1],b[ip][j+1],c[ip][j+1]);
				glVertex3fv(pd[(j+1)*Nxy+ip]);
				glNormal3d(a[ip][j],b[ip][j],c[ip][j]);
				glVertex3fv(pd[j*Nxy+ip]);
			glEnd();
		}
	}
}

//----------------------------------------------------------------------
//���ʂ����ɉ~�����e
void skTexGridCylinder(int Nxy, int Nz, float p[][3], bool flagShadow)
{
    //p[0]�͏��̒��S�C�Ō��p[]�͉���̒��S
	int i,j,ip,im,np,npL,npR,npU,npD;
//	double ct,theta,z,fz;
	float a[11][10],b[11][10],c[11][10];//���_�̖@������
    float shadowDiffuse[] = {0.2f,0.2f,0.2f,0.3f};//�e�̊g�U��

	float n[11][3];//���_�̖@������

    if(Nxy > 10) Nxy = 10;
    if(Nz > 10) Nz = 10;

    //�S���_��
    int num = Nxy * Nz + 2;

    //Top
    //��ꒆ�S��Top��numCol�̎O�p�`�̖@���x�N�g���𕽋ω�
    a[0][0] = 0.0; b[0][0] = 0.0; c[0][0] = 0.0;
    //���̎��͂̒��_�ׂ͗荇��2�̎O�p�`�̖@���x�N�g���𕽋ω�
    for(i = 1;i <= Nxy; i++){
        a[i][0] = 0.0; b[i][0] = 0.0; c[i][0] = 0.0;
    }
	for(i = 0;i < Nxy;i++){
 		ip = i+1;
		if(ip == Nxy) ip = 0;
        calcNormal(p[0],p[i+1],p[ip+1],n[i]);
        a[0][0] += n[i][0];
        b[0][0] += n[i][1];
        c[0][0] += n[i][2];
        a[i+1][0] += n[i][0];
        b[i+1][0] += n[i][1];
        c[i+1][0] += n[i][2];
        a[ip+1][0] += n[i][0];
        b[ip+1][0] += n[i][1];
        c[ip+1][0] += n[i][2];
    }
    a[0][0] = a[0][0] / (float)Nxy;
    b[0][0] = b[0][0] / (float)Nxy;
    c[0][0] = c[0][0] / (float)Nxy;
    for(i = 1; i < Nxy; i++){
        a[i][0] = a[i][0] / 2.0f;
        b[i][0] = b[i][0] / 2.0f;
        c[i][0] = c[i][0] / 2.0f;
    }
    //Top�̕\��
    for(i = 1; i <= Nxy; i++){
        ip = i + 1;
        if(ip == Nxy + 1) ip = 1;
		glBegin(GL_TRIANGLES);
			glNormal3f(a[0][0],b[0][0],c[0][0]);glVertex3fv(p[0]);
			glNormal3f(a[i][0],b[i][0],c[i][0]);glVertex3fv(p[i]);
			glNormal3f(a[ip][0],b[ip][0],c[ip][0]);glVertex3fv(p[ip]);
		glEnd();
	}

    //Bottom
    //���ꒆ�S��Bottom��numCol�̎O�p�`�̖@���x�N�g���𕽋ω�
    a[0][0] = 0.0; b[0][0] = 0.0; c[0][0] = 0.0;
    //����̎��͂̒��_�ׂ͗荇��2�̎O�p�`�̖@���x�N�g���𕽋ω�
    for(i = 1;i <= Nxy; i++){
        a[i][0] = 0.0; b[i][0] = 0.0; c[i][0] = 0.0;
    }
	for(i = 0;i < Nxy;i++){
 		ip = i+1;
		if(ip == Nxy) ip = 0;
        calcNormal(p[num-1],p[num-2-i],p[num-2-ip],n[i]);
        a[0][0] += n[i][0];
        b[0][0] += n[i][1];
        c[0][0] += n[i][2];
        a[i+1][0] += n[i][0];
        b[i+1][0] += n[i][1];
        c[i+1][0] += n[i][2];
        a[ip+1][0] += n[i][0];
        b[ip+1][0] += n[i][1];
        c[ip+1][0] += n[i][2];
    }
    a[0][0] = a[0][0] / (float)Nxy;
    b[0][0] = b[0][0] / (float)Nxy;
    c[0][0] = c[0][0] / (float)Nxy;
    for(i = 1; i < Nxy; i++){
        a[i][0] = a[i][0] / 2.0f;
        b[i][0] = b[i][0] / 2.0f;
        c[i][0] = c[i][0] / 2.0f;
    }
    //Bottom�̕\��
    for(i = 1; i <= Nxy; i++){
        ip = i + 1;
        if(ip == Nxy + 1) ip = 1;
		glBegin(GL_TRIANGLES);
			glNormal3f(a[0][0],b[0][0],c[0][0]);glVertex3fv(p[num-1]);
			glNormal3f(a[i][0],b[i][0],c[i][0]);glVertex3fv(p[num-1-i]);
			glNormal3f(a[ip][0],b[ip][0],c[ip][0]); glVertex3fv(p[num-1-ip]);
		glEnd();
	}

	//���ʂ̖@������
    //���̎��͒��_�̑��ʂɑ΂���@��������
    //�ׂ荇��2�̎O�p�`�̖@���x�N�g���𕽋ω�
	for(i = 0;i < Nxy;i++){
 		ip = i+1;
		if(ip == Nxy) ip = 0;
		im = i-1;
		if(i == 0) im = Nxy-1;
        //j = 0;
        np = i+1;//���ړ_
        npL = im+1;//����
        npR = ip+1;//�E��
        //npU = np - Nxy;//��
        npD = np + Nxy;//��
        calcNormal(p[np],p[npL],p[npD],n[0]);//�O���猩�č���
        calcNormal(p[np],p[npD],p[npR],n[1]);//�E��

        a[i][0] = (float)((n[0][0]+n[1][0])/2.0f);//������
        b[i][0] = (float)((n[0][1]+n[1][1])/2.0f);//��
        c[i][0] = (float)((n[0][2]+n[1][2])/2.0f);//��
	}

    //�����悤�ɉ���̒��_�ɑ΂���
	for(i = 0;i < Nxy;i++){
 		ip = i+1;
		if(ip == Nxy) ip = 0;
		im = i-1;
		if(i == 0) im = Nxy-1;
        j = Nz - 1;
        np = j*Nxy+i+1;//���ړ_
        npL = j*Nxy+im+1;//����
        npR = j*Nxy+ip+1;//�E��
        npU = np - Nxy;//��
        //npD = np + Nxy;//��
        calcNormal(p[np],p[npU],p[npL],n[0]);//�O���猩�č���
        calcNormal(p[np],p[npR],p[npU],n[1]);//�E��

        a[i][j] = (float)((n[0][0]+n[1][0])/2.0f);//������
        b[i][j] = (float)((n[0][1]+n[1][1])/2.0f);//��
        c[i][j] = (float)((n[0][2]+n[1][2])/2.0f);//��
	}

    //���̑��ʂ̒��_�ɑ΂��Ăׂ͗荇��4�̎O�p�`�̖@���x�N�g���𕽋ω�
	for(i = 0;i < Nxy;i++){
 		ip = i+1;
		if(ip == Nxy) ip = 0;
		im = i-1;
		if(i == 0) im = Nxy-1;
		for(j = 1;j < Nz-1;j++)
		{
			np =  j * Nxy + i + 1;//���ړ_
			npL = j * Nxy + im + 1;//����
			npR = j * Nxy + ip + 1;//�E��
			npU = np - Nxy;//��
			npD = np + Nxy;//��
            calcNormal(p[np],p[npU],p[npL],n[0]);//�O���猩�č���
            calcNormal(p[np],p[npR],p[npU],n[1]);//�E��
            calcNormal(p[np],p[npL],p[npD],n[2]);//�O���猩�č���
            calcNormal(p[np],p[npD],p[npR],n[3]);//�E��

			a[i][j] = (float)((n[0][0]+n[1][0]+n[2][0]+n[3][0])/4.0f);//������
			b[i][j] = (float)((n[0][1]+n[1][1]+n[2][1]+n[3][1])/4.0f);//��
			c[i][j] = (float)((n[0][2]+n[1][2]+n[2][2]+n[3][2])/4.0f);//��
		}
	}

	//����(4�p�`�p�b�`�j
    if(flagShadow == true) goto SHADOW_DISP;
	float th1, th2, ph1, ph2;
	glEnable(GL_TEXTURE_2D);
	for(j = 0;j < Nz-1;j++){
		ph1 = 1.0f - (float)j / (float)(Nz-1);
		ph2 = 1.0f - (float)(j + 1) / (float)(Nz-1);
		for(i = 0;i < Nxy;i++)
		{	ip = i+1;
			if(ip == Nxy) ip = 0;
			th1 = (float)i / (float)Nxy;
			th2 = (float)(i+1) / (float)Nxy;
			//if(i == Nxy-1) th2 = 1.0f;
			glBegin(GL_QUADS);
				glNormal3d(a[i][j],b[i][j],c[i][j]);
				glTexCoord2f(th1, ph1);glVertex3fv(p[j*Nxy+i+1]);
				glNormal3d(a[i][j+1],b[i][j+1],c[i][j+1]);
				glTexCoord2f(th1, ph2);glVertex3fv(p[(j+1)*Nxy+i+1]);
				glNormal3d(a[ip][j+1],b[ip][j+1],c[ip][j+1]);
				glTexCoord2f(th2, ph2);glVertex3fv(p[(j+1)*Nxy+ip+1]);
				glNormal3d(a[ip][j],b[ip][j],c[ip][j]);
				glTexCoord2f(th2, ph1);glVertex3fv(p[j*Nxy+ip+1]);
			glEnd();
		}
	}
	glDisable(GL_TEXTURE_2D);
    return;
SHADOW_DISP:;
	for(j = 0;j < Nz-1;j++){
		for(i = 0;i < Nxy;i++)
		{	ip = i+1;
			if(ip == Nxy) ip = 0;
            glMaterialfv(GL_FRONT,GL_DIFFUSE,shadowDiffuse);
			glBegin(GL_QUADS);
				glNormal3d(a[i][j],b[i][j],c[i][j]);
				glVertex3fv(p[j*Nxy+i+1]);
				glNormal3d(a[i][j+1],b[i][j+1],c[i][j+1]);
				glVertex3fv(p[(j+1)*Nxy+i+1]);
				glNormal3d(a[ip][j+1],b[ip][j+1],c[ip][j+1]);
				glVertex3fv(p[(j+1)*Nxy+ip+1]);
				glNormal3d(a[ip][j],b[ip][j],c[ip][j]);
				glVertex3fv(p[j*Nxy+ip+1]);
			glEnd();
		}
	}
}
//------------------------------------------------------------------------
//2�����i�q��I�u�W�F�N�g�Ɋi�q�Ԋu�ɏ]���Đ�,�V�A���̃`�F�b�N�͗l�𒼐ڕ`��
void skDrawCheck2(int numRow, int numCol, float pd[][3], bool flagShadow)
{
	int i, j, np;
	float n1[3], n2[3], n3[3], n4[3];
	float a[10][10], b[10][10], c[10][10];
    float shadowDiffuse[] = {0.2f,0.2f,0.2f,0.3f};//�e�̊g�U��

	//�@������
	for(i = 0;i < numRow;i++) //������
	  for(j = 0;j < numCol;j++)//�E����
	  {
        np = i * numCol + j;
		if(i == 0 )
		{
			if(j == 0) {
				calcNormal(pd[0],pd[numCol],pd[1],n1);
				a[i][j] = n1[0]; b[i][j] = n1[1]; c[i][j] = n1[2]; }
			else if(j == numCol-1) {
				calcNormal(pd[numCol-2],pd[2*numCol-1],pd[numCol-1],n1);
				a[i][j] = n1[0]; b[i][j] = n1[1]; c[i][j] = n1[2]; }
			else {
				calcNormal(pd[np],pd[np-1],pd[np+numCol],n1);//����
				calcNormal(pd[np],pd[np+numCol],pd[np+1],n2);//�E��
				a[i][j] = (n1[0]+n2[0])/2.0f;
				b[i][j] = (n1[1]+n2[1])/2.0f;
				c[i][j] = (n1[2]+n2[2])/2.0f; }
		}
		else if(i == numRow-1)
		{
			if(j == 0) {
				calcNormal(pd[np],pd[np+1],pd[np-numCol],n1);
				a[i][j] = n1[0]; b[i][j] = n1[1]; c[i][j] = n1[2]; }
			else if(j == numCol-1) {
				calcNormal(pd[np],pd[np-numCol],pd[np-1],n1);
				a[i][j] = n1[0]; b[i][j] = n1[1]; c[i][j] = n1[2]; }
			else {
				calcNormal(pd[np],pd[np-numCol],pd[np-1],n1);//����
				calcNormal(pd[np],pd[np+1],pd[np-numCol],n2);//�E��
				a[i][j] = (n1[0]+n2[0])/2.0f;
				b[i][j] = (n1[1]+n2[1])/2.0f;
				c[i][j] = (n1[2]+n2[2])/2.0f; }
		}
		else
		{
			if(j == 0) {
				calcNormal(pd[np],pd[np+1],pd[np-numCol],n1);//��
				calcNormal(pd[np],pd[np+numCol],pd[np+1],n2);//��
				a[i][j] = (n1[0]+n2[0])/2.0f;
				b[i][j] = (n1[1]+n2[1])/2.0f;
				c[i][j] = (n1[2]+n2[2])/2.0f; }
			else if(j == numCol-1) {
				calcNormal(pd[np],pd[np-numCol],pd[np-1],n1);//��
				calcNormal(pd[np],pd[np-1],pd[np+numCol],n2);//��
				a[i][j] = (n1[0]+n2[0])/2.0f;
				b[i][j] = (n1[1]+n2[1])/2.0f;
				c[i][j] = (n1[2]+n2[2])/2.0f; }
			else {//�㉺���E�S�̎O�p�`�̕���
				calcNormal(pd[np],pd[np-numCol],pd[np-1],n1);//����
				calcNormal(pd[np],pd[np+1],pd[np-numCol],n2);//�E��
				calcNormal(pd[np],pd[np-1],pd[np+numCol],n3);//����
				calcNormal(pd[np],pd[np+numCol],pd[np+1],n4);//�E��
				a[i][j] = (n1[0]+n2[0]+n3[0]+n4[0])/4.0f;
				b[i][j] = (n1[1]+n2[1]+n3[1]+n4[1])/4.0f;
				c[i][j] = (n1[2]+n2[2]+n3[2]+n4[2])/4.0f; }
		}
	}

    float diffuse[4] ;
    int k;
	//�`�F�b�N�͗l�`��
	for(i = 0;i < numRow-1;i++)
	  for(j = 0;j < numCol-1;j++){
        if(fmod(i+j, 2.0) == 0.0)
            { diffuse[0] = 1.0f; diffuse[1] = 0.0f; diffuse[2] = 0.0f; }
        else
            { diffuse[0] = 0.0f; diffuse[1] = 1.0f; diffuse[2] = 1.0f; }
        diffuse[3] = 1.0f;
        if(flagShadow == true){
            for(k = 0; k < 4; k++) diffuse[k] = shadowDiffuse[k];
        }
        glMaterialfv(GL_FRONT,GL_DIFFUSE,diffuse);

        np = i * numCol + j;
		glBegin(GL_TRIANGLES);
			//�����̎O�p�`
			//�e���_�̖@������,ø�������W,���_���W��^����B
			glNormal3d(a[i][j],b[i][j],c[i][j]);//�@������
			glVertex3fv(pd[np]);//��غ�݂̒��_���W�i�ȉ��������J��Ԃ��j
			glNormal3d(a[i+1][j],b[i+1][j],c[i+1][j]);
			glVertex3fv(pd[np+numCol]);
			glNormal3d(a[i+1][j+1],b[i+1][j+1],c[i+1][j+1]);
			glVertex3fv(pd[np+numCol+1]);
			//�E��̎O�p�`
			glNormal3d(a[i][j],b[i][j],c[i][j]);
			glVertex3fv(pd[np]);
			glNormal3d(a[i+1][j+1],b[i+1][j+1],c[i+1][j+1]);
			glVertex3fv(pd[np+numCol+1]);
			glNormal3d(a[i][j+1],b[i][j+1],c[i][j+1]);
			glVertex3fv(pd[np+1]);
		glEnd();
	}
}
//----------------------------------------------------------------------
void skDrawCheckCube(int numRow, int numCol, int numStk, float p1[][3],
        float p2[][3], float p3[][3], float p4[][3], float p5[][3], float p6[][3], bool flagShadow)
{
    //���ׂē����e�N�X�`���[
    //�㑤
    skDrawCheck2(numRow, numCol, p1, flagShadow);
    //����
    skDrawCheck2(numStk, numCol, p2, flagShadow);
    //�E��
    skDrawCheck2(numStk, numRow, p3, flagShadow);
    //����
    skDrawCheck2(numStk, numRow, p4, flagShadow);
    //����
    skDrawCheck2(numStk, numCol, p5, flagShadow);
    //����
    skDrawCheck2(numRow, numCol, p6, flagShadow);
}
//------------------------------------------------------------------------------
//���i�q��I�u�W�F�N�g�ɑ΂���T�u���[�`��
//�i�q�Ԋu�ɏ]���ĐԁC�V�A���̃`�F�b�N�͗l�𒼐ڕ`��
void skDrawCheckSphere(int Nxy, int Nz, float pd[][3], bool flagShadow)
{
	int i,j,ip,im,np,npL,npR,npU,npD;
	float a[10][10],b[10][10],c[10][10];//���_�̖@������
    float shadowDiffuse[] = {0.2f,0.2f,0.2f,0.3f};//�e�̊g�U��

	float n[10][3];//���_�̖@������

    if(Nxy > 10) Nxy = 10;
    if(Nz > 10) Nz = 10;

    //Top
    //Top��numCol�̎O�p�`�̖@���x�N�g���𕽋ω�
    a[0][0] = 0.0; b[0][0] = 0.0; c[0][0] = 0.0;
	for(i = 0;i < Nxy;i++){
 		ip = i+1;
		if(ip == Nxy) ip = 0;
        calcNormal(pd[i],pd[i+Nxy],pd[ip+Nxy],n[i]);
        a[0][0] += n[i][0];
        b[0][0] += n[i][1];
        c[0][0] += n[i][2];
    }
    a[0][0] = a[0][0] / (float)Nxy;
    b[0][0] = b[0][0] / (float)Nxy;
    c[0][0] = c[0][0] / (float)Nxy;
    for(i = 1;i < Nxy; i++){ //Top���_�͂��ׂē���
        a[i][0] = a[0][0]; b[i][0] = b[0][0]; c[i][0] = c[0][0];
    }

	//���ʂ̖@������
    //�ׂ荇��4�̎O�p�`�̖@���x�N�g���𕽋ω�
	for(i = 0;i < Nxy;i++){
 		ip = i+1;
		if(ip == Nxy) ip = 0;
		im = i-1;
		if(i == 0) im = Nxy-1;
		for(j = 1;j < Nz-1;j++)
		{
			np =  j * Nxy + i;//���ړ_
			npL = j * Nxy + im;//����
			npR = j * Nxy + ip;//�E��
			npU = np - Nxy;//��
			npD = np + Nxy;//��
            calcNormal(pd[np],pd[npU],pd[npL],n[0]);//�O���猩�č���
            calcNormal(pd[np],pd[npR],pd[npU],n[1]);//�E��
            calcNormal(pd[np],pd[npL],pd[npD],n[2]);//�O���猩�č���
            calcNormal(pd[np],pd[npD],pd[npR],n[3]);//�E��

			a[i][j] = (float)((n[0][0]+n[1][0]+n[2][0]+n[3][0])/4.0f);//������
			b[i][j] = (float)((n[0][1]+n[1][1]+n[2][1]+n[3][1])/4.0f);//��
			c[i][j] = (float)((n[0][2]+n[1][2]+n[2][2]+n[3][2])/4.0f);//��
		}
	}
    //Bottom
    //Bottom��numCol�̎O�p�`�̖@���x�N�g���𕽋ω�
    a[0][Nxy-1] = 0.0; b[0][Nxy-1] = 0.0; c[0][Nxy-1] = 0.0;
	for(i = 0;i < Nxy;i++){
        np = (Nz-1) * Nxy +i;
 		ip = i+1;
		if(ip == Nxy) ip = 0;
        calcNormal(pd[np],pd[np-i+ip-Nxy],pd[np-Nxy],n[i]);
        a[0][Nz-1] += n[i][0];
        b[0][Nz-1] += n[i][1];
        c[0][Nz-1] += n[i][2];
    }
    a[0][Nz-1] = a[0][Nz-1] / (float)Nxy;
    b[0][Nz-1] = b[0][Nz-1] / (float)Nxy;
    c[0][Nz-1] = c[0][Nz-1] / (float)Nxy;
    for(i = 1;i < Nxy; i++){ //Top���_�͂��ׂē���
        a[i][Nz-1] = a[0][Nz-1];
        b[i][Nz-1] = b[0][Nz-1];
        c[i][Nz-1] = c[0][Nz-1];
    }

	//�\��
    if(flagShadow == true) goto SHADOW_DISP;
    float diffuse[4];
	//����(4�p�`�p�b�`�j
	for(j = 0;j < Nz-1;j++){
		for(i = 0;i < Nxy;i++)
		{
            if(fmod(i+j, 2.0) == 0.0)
                { diffuse[0] = 1.0f; diffuse[1] = 0.0f; diffuse[2] = 0.0f; }
            else
                { diffuse[0] = 0.0f; diffuse[1] = 1.0f; diffuse[2] = 1.0f; }
            diffuse[3] = 1.0f;
            glMaterialfv(GL_FRONT,GL_DIFFUSE,diffuse);

        	ip = i+1;
			if(ip == Nxy) ip = 0;
			glBegin(GL_QUADS);
				glNormal3d(a[i][j],b[i][j],c[i][j]);
				glVertex3fv(pd[j*Nxy+i]);
				glNormal3d(a[i][j+1],b[i][j+1],c[i][j+1]);
				glVertex3fv(pd[(j+1)*Nxy+i]);
				glNormal3d(a[ip][j+1],b[ip][j+1],c[ip][j+1]);
				glVertex3fv(pd[(j+1)*Nxy+ip]);
				glNormal3d(a[ip][j],b[ip][j],c[ip][j]);
				glVertex3fv(pd[j*Nxy+ip]);
			glEnd();
		}
	}
    return;
SHADOW_DISP:;
	for(j = 0;j < Nz-1;j++){
		for(i = 0;i < Nxy;i++)
		{	ip = i+1;
			if(ip == Nxy) ip = 0;
            glMaterialfv(GL_FRONT,GL_DIFFUSE,shadowDiffuse);
			glBegin(GL_QUADS);
				glNormal3d(a[i][j],b[i][j],c[i][j]);
                glVertex3fv(pd[j*Nxy+i]);
				glNormal3d(a[i][j+1],b[i][j+1],c[i][j+1]);
				glVertex3fv(pd[(j+1)*Nxy+i]);
				glNormal3d(a[ip][j+1],b[ip][j+1],c[ip][j+1]);
				glVertex3fv(pd[(j+1)*Nxy+ip]);
				glNormal3d(a[ip][j],b[ip][j],c[ip][j]);
				glVertex3fv(pd[j*Nxy+ip]);
			glEnd();
		}
	}
}
//----------------------------------------------------------------------------
//�~���i�q��I�u�W�F�N�g�ɑ΂���T�u���[�`��
//�i�q�Ԋu�ɏ]���ĐԁC�V�A���̃`�F�b�N�͗l�𒼐ڕ`��
void skDrawCheckCylinder(int Nxy, int Nz, float p[][3], bool flagShadow)
{
	int i,j,ip,im,np,npL,npR,npU,npD;
	float a[11][10],b[11][10],c[11][10];//���_�̖@������
	float n[11][3];//���_�̖@������
    float diffuse[4];
    float shadowDiffuse[] = {0.2f,0.2f,0.2f,0.3f};//�e�̊g�U��

    if(Nxy > 10) Nxy = 10;
    if(Nz > 10) Nz = 10;

    //�S���_��
    int num = Nxy * Nz + 2;

    //Top
    //��ꒆ�S��Top��numCol�̎O�p�`�̖@���x�N�g���𕽋ω�
    a[0][0] = 0.0; b[0][0] = 0.0; c[0][0] = 0.0;
    //���̎��͂̒��_�ׂ͗荇��2�̎O�p�`�̖@���x�N�g���𕽋ω�
    for(i = 1;i <= Nxy; i++){
        a[i][0] = 0.0; b[i][0] = 0.0; c[i][0] = 0.0;
    }
	for(i = 0;i < Nxy;i++){
 		ip = i+1;
		if(ip == Nxy) ip = 0;
        calcNormal(p[0],p[i+1],p[ip+1],n[i]);
        a[0][0] += n[i][0];
        b[0][0] += n[i][1];
        c[0][0] += n[i][2];
        a[i+1][0] += n[i][0];
        b[i+1][0] += n[i][1];
        c[i+1][0] += n[i][2];
        a[ip+1][0] += n[i][0];
        b[ip+1][0] += n[i][1];
        c[ip+1][0] += n[i][2];
    }
    a[0][0] = a[0][0] / (float)Nxy;
    b[0][0] = b[0][0] / (float)Nxy;
    c[0][0] = c[0][0] / (float)Nxy;
    for(i = 1; i < Nxy; i++){
        a[i][0] = a[i][0] / 2.0f;
        b[i][0] = b[i][0] / 2.0f;
        c[i][0] = c[i][0] / 2.0f;
    }
    //Top�̕\��
    for(i = 1; i <= Nxy; i++){
        ip = i + 1;
        if(ip == Nxy + 1) ip = 1;
        if(fmod(i, 2.0) == 0.0)
        { diffuse[0] = 1.0f; diffuse[1] = 0.0f; diffuse[2] = 0.0f; }
        else
        { diffuse[0] = 0.0f; diffuse[1] = 1.0f; diffuse[2] = 1.0f; }
        diffuse[3] = 1.0f;
        if(flagShadow == true) glMaterialfv(GL_FRONT,GL_DIFFUSE,shadowDiffuse);
        else glMaterialfv(GL_FRONT,GL_DIFFUSE,diffuse);
		glBegin(GL_TRIANGLES);
			glNormal3f(a[0][0],b[0][0],c[0][0]);glVertex3fv(p[0]);
			glNormal3f(a[i][0],b[i][0],c[i][0]);glVertex3fv(p[i]);
			glNormal3f(a[ip][0],b[ip][0],c[ip][0]);glVertex3fv(p[ip]);
		glEnd();
	}

    //Bottom
    //���ꒆ�S��Bottom��numCol�̎O�p�`�̖@���x�N�g���𕽋ω�
    a[0][0] = 0.0; b[0][0] = 0.0; c[0][0] = 0.0;
    //����̎��͂̒��_�ׂ͗荇��2�̎O�p�`�̖@���x�N�g���𕽋ω�
    for(i = 1;i <= Nxy; i++){
        a[i][0] = 0.0; b[i][0] = 0.0; c[i][0] = 0.0;
    }
	for(i = 0;i < Nxy;i++){
 		ip = i+1;
		if(ip == Nxy) ip = 0;
        calcNormal(p[num-1],p[num-2-i],p[num-2-ip],n[i]);
        a[0][0] += n[i][0];
        b[0][0] += n[i][1];
        c[0][0] += n[i][2];
        a[i+1][0] += n[i][0];
        b[i+1][0] += n[i][1];
        c[i+1][0] += n[i][2];
        a[ip+1][0] += n[i][0];
        b[ip+1][0] += n[i][1];
        c[ip+1][0] += n[i][2];
    }
    a[0][0] = a[0][0] / (float)Nxy;
    b[0][0] = b[0][0] / (float)Nxy;
    c[0][0] = c[0][0] / (float)Nxy;
    for(i = 1; i < Nxy; i++){
        a[i][0] = a[i][0] / 2.0f;
        b[i][0] = b[i][0] / 2.0f;
        c[i][0] = c[i][0] / 2.0f;
    }
    //Bottom�̕\��
    for(i = 1; i <= Nxy; i++){
        ip = i + 1;
        if(ip == Nxy + 1) ip = 1;
        if(fmod(i, 2.0) == 0.0)
        { diffuse[0] = 1.0f; diffuse[1] = 0.0f; diffuse[2] = 0.0f; }
        else
        { diffuse[0] = 0.0f; diffuse[1] = 1.0f; diffuse[2] = 1.0f; }
        diffuse[3] = 1.0f;
        if(flagShadow == true) glMaterialfv(GL_FRONT,GL_DIFFUSE,shadowDiffuse);
        else glMaterialfv(GL_FRONT,GL_DIFFUSE,diffuse);
		glBegin(GL_TRIANGLES);
			glNormal3f(a[0][0],b[0][0],c[0][0]);glVertex3fv(p[num-1]);
			glNormal3f(a[i][0],b[i][0],c[i][0]);glVertex3fv(p[num-1-i]);
			glNormal3f(a[ip][0],b[ip][0],c[ip][0]); glVertex3fv(p[num-1-ip]);
		glEnd();
	}

	//���ʂ̖@������
    //���̎��͒��_�̑��ʂɑ΂���@��������
    //�ׂ荇��2�̎O�p�`�̖@���x�N�g���𕽋ω�
	for(i = 0;i < Nxy;i++){
 		ip = i+1;
		if(ip == Nxy) ip = 0;
		im = i-1;
		if(i == 0) im = Nxy-1;
        //j = 0;
        np = i+1;//���ړ_
        npL = im+1;//����
        npR = ip+1;//�E��
        //npU = np - Nxy;//��
        npD = np + Nxy;//��
        calcNormal(p[np],p[npL],p[npD],n[0]);//�O���猩�č���
        calcNormal(p[np],p[npD],p[npR],n[1]);//�E��

        a[i][0] = (float)((n[0][0]+n[1][0])/2.0f);//������
        b[i][0] = (float)((n[0][1]+n[1][1])/2.0f);//��
        c[i][0] = (float)((n[0][2]+n[1][2])/2.0f);//��
	}

    //�����悤�ɉ���̒��_�ɑ΂���
	for(i = 0;i < Nxy;i++){
 		ip = i+1;
		if(ip == Nxy) ip = 0;
		im = i-1;
		if(i == 0) im = Nxy-1;
        j = Nz - 1;
        np = j*Nxy+i+1;//���ړ_
        npL = j*Nxy+im+1;//����
        npR = j*Nxy+ip+1;//�E��
        npU = np - Nxy;//��
        //npD = np + Nxy;//��
        calcNormal(p[np],p[npU],p[npL],n[0]);//�O���猩�č���
        calcNormal(p[np],p[npR],p[npU],n[1]);//�E��

        a[i][j] = (float)((n[0][0]+n[1][0])/2.0f);//������
        b[i][j] = (float)((n[0][1]+n[1][1])/2.0f);//��
        c[i][j] = (float)((n[0][2]+n[1][2])/2.0f);//��
	}

    //���̑��ʂ̒��_�ɑ΂��Ăׂ͗荇��4�̎O�p�`�̖@���x�N�g���𕽋ω�
	for(i = 0;i < Nxy;i++){
 		ip = i+1;
		if(ip == Nxy) ip = 0;
		im = i-1;
		if(i == 0) im = Nxy-1;
		for(j = 1;j < Nz-1;j++)
		{
			np =  j * Nxy + i + 1;//���ړ_
			npL = j * Nxy + im + 1;//����
			npR = j * Nxy + ip + 1;//�E��
			npU = np - Nxy;//��
			npD = np + Nxy;//��
            calcNormal(p[np],p[npU],p[npL],n[0]);//�O���猩�č���
            calcNormal(p[np],p[npR],p[npU],n[1]);//�E��
            calcNormal(p[np],p[npL],p[npD],n[2]);//�O���猩�č���
            calcNormal(p[np],p[npD],p[npR],n[3]);//�E��

			a[i][j] = (float)((n[0][0]+n[1][0]+n[2][0]+n[3][0])/4.0f);//������
			b[i][j] = (float)((n[0][1]+n[1][1]+n[2][1]+n[3][1])/4.0f);//��
			c[i][j] = (float)((n[0][2]+n[1][2]+n[2][2]+n[3][2])/4.0f);//��
		}
	}

    //�\��
    if(flagShadow == true) goto SHADOW_DISP;
	//����(4�p�`�p�b�`�j
	for(j = 0;j < Nz-1;j++){
		for(i = 0;i < Nxy;i++)
		{	ip = i+1;
			if(ip == Nxy) ip = 0;
            if(fmod(i+j, 2.0) == 0.0)
                { diffuse[0] = 1.0f; diffuse[1] = 0.0f; diffuse[2] = 0.0f; }
            else
                { diffuse[0] = 0.0f; diffuse[1] = 1.0f; diffuse[2] = 1.0f; }
            diffuse[3] = 1.0f;
            glMaterialfv(GL_FRONT,GL_DIFFUSE,diffuse);
			glBegin(GL_QUADS);
				glNormal3d(a[i][j],b[i][j],c[i][j]);
				glVertex3fv(p[j*Nxy+i+1]);
				glNormal3d(a[i][j+1],b[i][j+1],c[i][j+1]);
				glVertex3fv(p[(j+1)*Nxy+i+1]);
				glNormal3d(a[ip][j+1],b[ip][j+1],c[ip][j+1]);
				glVertex3fv(p[(j+1)*Nxy+ip+1]);
				glNormal3d(a[ip][j],b[ip][j],c[ip][j]);
				glVertex3fv(p[j*Nxy+ip+1]);
			glEnd();
		}
	}
    return;
SHADOW_DISP:;
	for(j = 0;j < Nz-1;j++){
		for(i = 0;i < Nxy;i++)
		{	ip = i+1;
			if(ip == Nxy) ip = 0;
            glMaterialfv(GL_FRONT,GL_DIFFUSE,shadowDiffuse);
			glBegin(GL_QUADS);
				glNormal3d(a[i][j],b[i][j],c[i][j]);
				glVertex3fv(p[j*Nxy+i+1]);
				glNormal3d(a[i][j+1],b[i][j+1],c[i][j+1]);
				glVertex3fv(p[(j+1)*Nxy+i+1]);
				glNormal3d(a[ip][j+1],b[ip][j+1],c[ip][j+1]);
				glVertex3fv(p[(j+1)*Nxy+ip+1]);
				glNormal3d(a[ip][j],b[ip][j],c[ip][j]);
				glVertex3fv(p[j*Nxy+ip+1]);
			glEnd();
		}
	}
}

//------------------------------------------------------------------------------
//�ȉ��̓\���b�h�e�N�X�`���p
//------------------------------------------------------------------------------
//texImage��CTarget�N���X�Ōv�Z
//----------------------------------------------------------------------------
//�����̂̃\���b�h�e�N�X�`��
void skSolidTexCube(unsigned char texImage[6][T_MAX][T_MAX][3])
{
	//���̓����Ńe�N�X�`�����쐬
	float p[8][3]={{0.5f,0.5f,0.5f},{-0.5f,0.5f,0.5f},{-0.5f,-0.5f,0.5f},
	{0.5f,-0.5f,0.5f},{0.5f,0.5f,-0.5f},{-0.5f,0.5f,-0.5f},{-0.5f,-0.5f,0-0.5f},{0.5f,-0.5f,-0.5f}};
	glTexImage2D(GL_TEXTURE_2D,0,3,T_MAX,T_MAX,0,GL_RGB,GL_UNSIGNED_BYTE,texImage[0]);

	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);//top
		glNormal3f(0.0f,0.0f,1.0f); //z����
		glTexCoord2f(1.0f, 0.0f); glVertex3fv(p[0]);
		glTexCoord2f(1.0f, 1.0f); glVertex3fv(p[1]);
		glTexCoord2f(0.0f, 1.0f); glVertex3fv(p[2]);
		glTexCoord2f(0.0f, 0.0f); glVertex3fv(p[3]);
	glEnd();
	glDisable(GL_TEXTURE_2D);

	//Front(x = 0.5)
	glTexImage2D(GL_TEXTURE_2D,0,3,T_MAX,T_MAX,0,GL_RGB,GL_UNSIGNED_BYTE,texImage[1]);

	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);
		glNormal3f(1.0f,0.0f,0.0f);
		glTexCoord2f(1.0f, 1.0f); glVertex3fv(p[0]);
		glTexCoord2f(0.0f, 1.0f); glVertex3fv(p[3]);
		glTexCoord2f(0.0f, 0.0f); glVertex3fv(p[7]);
		glTexCoord2f(1.0f, 0.0f); glVertex3fv(p[4]);
	glEnd();
	glDisable(GL_TEXTURE_2D);

	//Left
	glTexImage2D(GL_TEXTURE_2D,0,3,T_MAX,T_MAX,0,GL_RGB,GL_UNSIGNED_BYTE,texImage[2]);

	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);
		glNormal3f(0.0f,1.0f,0.0f);
		glTexCoord2f(0.0f, 1.0f); glVertex3fv(p[0]);
		glTexCoord2f(0.0f, 0.0f); glVertex3fv(p[4]);
		glTexCoord2f(1.0f, 0.0f); glVertex3fv(p[5]);
		glTexCoord2f(1.0f, 1.0f); glVertex3fv(p[1]);
	glEnd();
	glDisable(GL_TEXTURE_2D);

	//Rear
	glTexImage2D(GL_TEXTURE_2D,0,3,T_MAX,T_MAX,0,GL_RGB,GL_UNSIGNED_BYTE,texImage[3]);

	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);
	 	glNormal3f(-1.0f,0.0f,0.0f); //-x����
		glTexCoord2f(0.0f, 1.0f); glVertex3fv(p[1]);
		glTexCoord2f(0.0f, 0.0f); glVertex3fv(p[5]);
		glTexCoord2f(1.0f, 0.0f); glVertex3fv(p[6]);
		glTexCoord2f(1.0f, 1.0f); glVertex3fv(p[2]);
	glEnd();
	glDisable(GL_TEXTURE_2D);

	//Right
	glTexImage2D(GL_TEXTURE_2D,0,3,T_MAX,T_MAX,0,GL_RGB,GL_UNSIGNED_BYTE,texImage[4]);

	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);
		glNormal3f(0.0f,-1.0f,0.0f); //-y����
		glTexCoord2f(0.0f, 1.0f); glVertex3fv(p[2]);
		glTexCoord2f(0.0f, 0.0f); glVertex3fv(p[6]);
		glTexCoord2f(1.0f, 0.0f); glVertex3fv(p[7]);
		glTexCoord2f(1.0f, 1.0f); glVertex3fv(p[3]);
	glEnd();
	glDisable(GL_TEXTURE_2D);

	glTexImage2D(GL_TEXTURE_2D,0,3,T_MAX,T_MAX,0,GL_RGB,GL_UNSIGNED_BYTE,texImage[5]);

	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);
		glNormal3f(0.0f,0.0f,-1.0f); //-z����
		glTexCoord2f(1.0f, 1.0f); glVertex3fv(p[4]);
		glTexCoord2f(0.0f, 1.0f); glVertex3fv(p[7]);
		glTexCoord2f(0.0f, 0.0f); glVertex3fv(p[6]);
		glTexCoord2f(1.0f, 0.0f); glVertex3fv(p[5]);
	glEnd();
	glDisable(GL_TEXTURE_2D);
}
//----------------------------------------------------------------------------------------
void skSolidTexSphere(int N, byte texImage[6][T_MAX][T_MAX][3])
{
	glTexImage2D(GL_TEXTURE_2D,0,3,T_MAX,T_MAX,0,GL_RGB,GL_UNSIGNED_BYTE,texImage[0]);

	//�`��
	skTexSphere2(N);
}
//-------------------------------------------------------------------------------------
void skSolidTexCylinder(int Nxy, byte texImage[6][T_MAX][T_MAX][3])
{
	int i, ii;

	//���ards�A����hgt�̉~��(���S�����_)
	//Nx--���p�`�̎��(Nx<0)
	float p[40][3]; //���_���W
	double theta0,theta;
	double th1,th2;
	float rds = 0.5;//���a
	float hgt = 1.0;//����
	double pp = 2.0 * M_PI;

	//���_���W
	if(Nxy > 20) { Nxy = 20;}
	theta0 = pp / (double)Nxy;
	for(i = 0;i < Nxy;i++)
	{   theta = theta0 * (double)i;
		//���[��i=0�Ƃ���i��޼ު�Ď��g���猩�ĉE�[�A-y������)
		p[i][0] = (float)(rds*sin(theta)); //����x����
		p[i][1] = -(float)(rds*cos(theta));//������
		p[i][2] = (float)hgt/2.0f;         //������(����)
		p[i+Nxy][0] = p[i][0];             //�����x����
		p[i+Nxy][1] = p[i][1];             //������
		p[i+Nxy][2] = -(float)hgt/2.0f;    //������
	}

	//Top
	glTexImage2D(GL_TEXTURE_2D,0,3,T_MAX,T_MAX,0,GL_RGB,GL_UNSIGNED_BYTE,texImage[0]);
		
	glEnable(GL_TEXTURE_2D);
	glBegin(GL_POLYGON);
		glNormal3f(0.0f,0.0f,1.0f); //z����
		for(i = 0; i < Nxy; i++) {
			theta = theta0 * i;
			glTexCoord2d((1.0 - cos(theta)) / 2.0, (1.0 - sin(theta)) / 2.0);
			glVertex3fv(p[i]);
		}
	glEnd();
	glDisable(GL_TEXTURE_2D);

	//Bottom
	glTexImage2D(GL_TEXTURE_2D,0,3,T_MAX,T_MAX,0,GL_RGB,GL_UNSIGNED_BYTE,texImage[1]);
		
	glEnable(GL_TEXTURE_2D);
	glBegin(GL_POLYGON);
		glNormal3f(0.0f,0.0f,-1.0f); //-z����
		for(i = 0; i < Nxy; i++) {
			theta = theta0 * (i + 1);
			ii = 2 * Nxy - i - 1;
			glTexCoord2d((1.0 + cos(theta)) / 2.0, (1.0 + sin(theta)) / 2.0);
			glVertex3fv(p[ii]);
		}
	glEnd();
	glDisable(GL_TEXTURE_2D);

	//side
	glTexImage2D(GL_TEXTURE_2D,0,3,T_MAX,T_MAX,0,GL_RGB,GL_UNSIGNED_BYTE,texImage[2]);

	glEnable(GL_TEXTURE_2D);	
	for(i = 0;i < Nxy;i++)
	{
		if(i == Nxy-1) ii = 0; else ii = i + 1;
		th1 = (double)i/(double)Nxy;//2�΂Ő��K�������p�x
		th2 = (double)ii/(double)Nxy;
		if(i == Nxy-1) th2 = 1.0;
		glBegin(GL_QUADS);
			glNormal3f((float)sin(pp * th1),-(float)cos(pp * th1), 0.0f);
			glTexCoord2f((float)th1, 1.0f);//ø�������W��i=0�̂Ƃ�s=0�ƂȂ�悤��
			glVertex3fv(p[i]);
			glTexCoord2f((float)th1, 0.0f);
			glVertex3fv(p[i+Nxy]);

			glNormal3f((float)sin(pp * th2),-(float)cos(pp * th2) ,0.0f);
			glTexCoord2f((float)th2, 0.0f);
			glVertex3fv(p[ii+Nxy]); 
			glTexCoord2f((float)th2, 1.0f);
			glVertex3fv(p[ii]);
		glEnd();
	}
	glDisable(GL_TEXTURE_2D);
}

