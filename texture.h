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
//ﾃｸｽﾁｬの最大サイズ
#define T_MAX 128   

//------------------------------------------------------------------------------
void skTexSquare()//y-z平面
{
	static float p[4][3]={{0.0,-0.5,-0.5},{0.0,0.5,-0.5},{0.0,0.5,0.5},
	{0.0,-0.5,0.5} };

	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);
		glNormal3f(1.0f,0.0f,0.0f); //x方向の法線
		//ﾃｸｽﾁｬｰ座標と頂点番号との対応付け
		glTexCoord2f(0.0f, 0.0f); glVertex3fv(p[0]);
		glTexCoord2f(1.0f, 0.0f); glVertex3fv(p[1]);
		glTexCoord2f(1.0f, 1.0f); glVertex3fv(p[2]);
		glTexCoord2f(0.0f, 1.0f); glVertex3fv(p[3]);
	glEnd();
	//裏面も同じﾃｸｽﾁｬｰ
	glBegin(GL_QUADS);
		glNormal3f(-1.0f,0.0f,0.0f); //x方向の法線
		//ﾃｸｽﾁｬｰ座標と頂点番号との対応付け
		glTexCoord2f(0.0f, 0.0f); glVertex3fv(p[0]);
		glTexCoord2f(1.0f, 0.0f); glVertex3fv(p[3]);
		glTexCoord2f(1.0f, 1.0f); glVertex3fv(p[2]);
		glTexCoord2f(0.0f, 1.0f); glVertex3fv(p[1]);
	glEnd();
	glDisable(GL_TEXTURE_2D);
}

//-----------------------------------------------------------------------
//正面だけﾏｯﾋﾟﾝｸﾞ
void skTexCube1()
{
	float p[8][3]={{0.5f,0.5f,0.5f},{-0.5f,0.5f,0.5f},{-0.5f,-0.5f,0.5f},
	{0.5f,-0.5f,0.5f},{0.5f,0.5f,-0.5f},{-0.5f,0.5f,-0.5f},{-0.5f,-0.5f,0-0.5f},{0.5f,-0.5f,-0.5f}};
//glTexImage2D(GL_TEXTURE_2D,0,4,T_MAX,T_MAX,0,GL_RGBA,GL_UNSIGNED_BYTE,texImage);

	glBegin(GL_QUADS);
		glNormal3f(0.0f,0.0f,1.0f); //z方向
		glVertex3fv(p[0]); glVertex3fv(p[1]);
		glVertex3fv(p[2]); glVertex3fv(p[3]);
	glEnd();

	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);
		glNormal3f(1.0f,0.0f,0.0f); //x方向(正面）
		glTexCoord2f(1.0f, 1.0f); glVertex3fv(p[0]);
		glTexCoord2f(0.0f, 1.0f); glVertex3fv(p[3]);
		glTexCoord2f(0.0f, 0.0f); glVertex3fv(p[7]);
		glTexCoord2f(1.0f, 0.0f); glVertex3fv(p[4]);
	glEnd();
	glDisable(GL_TEXTURE_2D);

	glBegin(GL_QUADS);
		glNormal3f(0.0f,1.0f,0.0f); //y方向
		glVertex3fv(p[0]); glVertex3fv(p[4]);
		glVertex3fv(p[5]); glVertex3fv(p[1]);
	glEnd();

	glBegin(GL_QUADS);
	 	glNormal3f(-1.0f,0.0f,0.0f); //-x方向
		glVertex3fv(p[1]); glVertex3fv(p[5]);
		glVertex3fv(p[6]); glVertex3fv(p[2]);
	glEnd();

	glBegin(GL_QUADS);
		glNormal3f(0.0f,-1.0f,0.0f); //-y方向
		glVertex3fv(p[2]); glVertex3fv(p[6]);
		glVertex3fv(p[7]); glVertex3fv(p[3]);
	glEnd();

	glBegin(GL_QUADS);
		glNormal3f(0.0f,0.0f,-1.0f); //-z方向
		glVertex3fv(p[4]); glVertex3fv(p[7]);
		glVertex3fv(p[6]); glVertex3fv(p[5]);
	glEnd();
}
//-----------------------------------------------------------------------
//６面に同じ模様をﾏｯﾋﾟﾝｸﾞ
void skTexCube2()
{
	float p[8][3]={{0.5f,0.5f,0.5f},{-0.5f,0.5f,0.5f},{-0.5f,-0.5f,0.5f},
	{0.5f,-0.5f,0.5f},{0.5f,0.5f,-0.5f},{-0.5f,0.5f,-0.5f},{-0.5f,-0.5f,0-0.5f},{0.5f,-0.5f,-0.5f}};

	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);//top
		glNormal3f(0.0f,0.0f,1.0f); //z方向
		glTexCoord2f(1.0f, 0.0f); glVertex3fv(p[0]);
		glTexCoord2f(1.0f, 1.0f); glVertex3fv(p[1]);
		glTexCoord2f(0.0f, 1.0f); glVertex3fv(p[2]);
		glTexCoord2f(0.0f, 0.0f); glVertex3fv(p[3]);
	glEnd();

	glBegin(GL_QUADS);
		glNormal3f(1.0f,0.0f,0.0f); //x方向(正面）
		glTexCoord2f(1.0f, 1.0f); glVertex3fv(p[0]);
		glTexCoord2f(0.0f, 1.0f); glVertex3fv(p[3]);
		glTexCoord2f(0.0f, 0.0f); glVertex3fv(p[7]);
		glTexCoord2f(1.0f, 0.0f); glVertex3fv(p[4]);
	glEnd();

	glBegin(GL_QUADS);
		glNormal3f(0.0f,1.0f,0.0f); //y方向
		glTexCoord2f(0.0f, 1.0f); glVertex3fv(p[0]);
		glTexCoord2f(0.0f, 0.0f); glVertex3fv(p[4]);
		glTexCoord2f(1.0f, 0.0f); glVertex3fv(p[5]);
		glTexCoord2f(1.0f, 1.0f); glVertex3fv(p[1]);
	glEnd();

	glBegin(GL_QUADS);
	 	glNormal3f(-1.0f,0.0f,0.0f); //-x方向
		glTexCoord2f(0.0f, 1.0f); glVertex3fv(p[1]);
		glTexCoord2f(0.0f, 0.0f); glVertex3fv(p[5]);
		glTexCoord2f(1.0f, 0.0f); glVertex3fv(p[6]);
		glTexCoord2f(1.0f, 1.0f); glVertex3fv(p[2]);
	glEnd();

	glBegin(GL_QUADS);
		glNormal3f(0.0f,-1.0f,0.0f); //-y方向
		glTexCoord2f(0.0f, 1.0f); glVertex3fv(p[2]);
		glTexCoord2f(0.0f, 0.0f); glVertex3fv(p[6]);
		glTexCoord2f(1.0f, 0.0f); glVertex3fv(p[7]);
		glTexCoord2f(1.0f, 1.0f); glVertex3fv(p[3]);
	glEnd();

	glBegin(GL_QUADS);
		glNormal3f(0.0f,0.0f,-1.0f); //-z方向
		glTexCoord2f(1.0f, 1.0f); glVertex3fv(p[4]);
		glTexCoord2f(0.0f, 1.0f); glVertex3fv(p[7]);
		glTexCoord2f(0.0f, 0.0f); glVertex3fv(p[6]);
		glTexCoord2f(1.0f, 0.0f); glVertex3fv(p[5]);
        glEnd();
	glDisable(GL_TEXTURE_2D);
}

//--------------------------------------------------------------------------
//正面側(x>0)に平行投影
void skTexSphere1(int Nxy)//単位直径の球
{
	int i, j, Nz;
	double phai; //x-y平面に対する偏角（緯度）
	double theta; //x軸に対する偏角（経度)
	float p[21][21][3]; //各点の座標
	float *p1,*p2,*p3,*p4;

	if(Nxy > 20) Nxy = 20;
	Nz = Nxy;
	//頂点番号、面番号
	//座標
	for(i=0;i<=Nxy;i++)
	{	theta = 2.0 * M_PI * (double)i / (double)Nxy;
		for(j=0;j<=Nz;j++)
		{	//左端をi=0（自分からみての右端、-y軸方向）
			phai = M_PI / 2.0 - M_PI * (double)j / (double)Nz;
			p[i][j][0] = (float)(0.5 * sin(theta) * cos(phai)); //x座標
			p[i][j][1] = -(float)(0.5 * cos(theta) * cos(phai)); //y
			p[i][j][2] = (float)(0.5 * sin(phai));            //z
		}
	}

	//平行投影
	//表側(x>0)にﾏｯﾋﾟﾝｸﾞ
	glEnable(GL_TEXTURE_2D);
	for(i = 0;i < Nxy/2;i++){
		for(j = 0;j < Nz; j++)//-1;j++)
		{
			p1 = p[i][j];	  p2 = p[i][j+1];
			p3 = p[i+1][j+1]; p4 = p[i+1][j];
			glBegin(GL_QUADS);//左端(p1[1]=-0.5)でﾃｸｽﾁｬｰ座標が0となるようにする
				glNormal3fv(p1);glTexCoord2f((0.5f+p1[1]),(0.5f+p1[2]));glVertex3fv(p1);
				glNormal3fv(p2);glTexCoord2f((0.5f+p2[1]),(0.5f+p2[2]));glVertex3fv(p2);
				glNormal3fv(p3);glTexCoord2f((0.5f+p3[1]),(0.5f+p3[2]));glVertex3fv(p3);
				glNormal3fv(p4);glTexCoord2f((0.5f+p4[1]),(0.5f+p4[2]));glVertex3fv(p4);
			glEnd();
		}
	}
	glDisable(GL_TEXTURE_2D);
	//裏面にはﾏｯﾋﾟﾝｸﾞしない
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
//球投影
void skTexSphere2(int Nxy)//単位直径の球
{
	int i, ii, j, Nz;
	double phai; //x-y平面に対する偏角（緯度）
	double theta; //x軸に対する偏角（経度)
	float p[21][21][3]; //各点の座標
	float *p1,*p2,*p3,*p4;
	float th1, th2, ph1, ph2;

	if(Nxy > 20) Nxy = 20;
	Nz = Nxy;
	//頂点番号、面番号
	//座標
	for(i=0;i<=Nxy;i++)
	{	theta = 2.0*M_PI*(double)i/(double)Nxy;
		for(j=0;j<=Nz;j++)
		{
			phai = M_PI/2.0-M_PI*(double)j/(double)Nz;
          	//真後ろをi=0
			p[i][j][0] = -(float)(0.5*cos(theta)*cos(phai)); //x座標
			p[i][j][1] = -(float)(0.5*sin(theta)*cos(phai)); //y
			p[i][j][2] = (float)(0.5*sin(phai));             //z
		}
	}

	//球投影
	//全体にﾏｯﾋﾟﾝｸﾞ
	glEnable(GL_TEXTURE_2D);
	for(i = 0;i < Nxy;i++)
	{
		//if(i == Nxy-1) ii = 0; else ii = i+1;
		ii = i + 1;
		th1 = (float)i / (float)Nxy;//2πで正規化した角度
		th2 = (float)ii / (float)Nxy;
		if(i == Nxy-1) th2 = 1.0;
		for(j = 0;j < Nz; j++)
		{
			ph1 = 1.0f - (float)j / (float)Nz;//j=0で1となるように
			ph2 = 1.0f - (float)(j+1) / (float)Nz;
			p1 = p[i][j];	 p2 = p[i][j+1];
			p3 = p[ii][j+1]; p4 = p[ii][j];
			glBegin(GL_QUADS);//左端(p1[1]=-0.5)でﾃｸｽﾁｬｰ座標が0となるようにする
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
//側面正面側(x>0)に平行投影
void skTexCylinder1(int Nxy)
{	//半径rds、高さhgtの円柱(中心が原点)
	//Nx--多角形の種類(Nx<0)
	float p[40][3]; //頂点座標
	double theta0,theta;
	double th1,th2;
	int i,ii;
	float rds = 0.5;//半径
	float hgt = 1.0;//高さ
//	double pp = 2.0 * M_PI;

	if(Nxy > 20) { Nxy = 20;}
	theta0 = M_PI / (double)Nxy;
	for(i = 0;i < Nxy;i++)
	{   theta = 2.0 * theta0 * (double)i;
		//左端をi=0とする（ｵﾌﾞｼﾞｪｸﾄ自身から見て右端、-y軸方向)
		p[i][0] = (float)(rds*sin(theta)); //上底のx成分
		p[i][1] = -(float)(rds*cos(theta));//ｙ成分
		p[i][2] = (float)hgt/2.0f;         //ｚ成分(高さ)
		p[i+Nxy][0] = p[i][0];             //下底のx成分
		p[i+Nxy][1] = p[i][1];             //ｙ成分
		p[i+Nxy][2] = -(float)hgt/2.0f;    //ｚ成分
	}

	//上底0
	glBegin(GL_POLYGON);
		glNormal3f(0.0f,0.0f,1.0f);
		for(i=0;i<Nxy;i++) glVertex3fv(p[i]);
	glEnd();
	//下底
	glBegin(GL_POLYGON);
		glNormal3f(0.0f,0.0f,-1.0f);
		for(i=2*Nxy-1;i>=Nxy;i--) glVertex3fv(p[i]);
	glEnd();

	//平行投影
	//正面(x>0の側面半分にだけﾃｸｽﾁｬｰﾏｯﾋﾟﾝｸﾞ:平行投影）
	glEnable(GL_TEXTURE_2D);
	for(i=0;i<Nxy/2;i++)
	{
		if(i == Nxy-1) ii = 0; else ii = i+1;
		th1 = 2.0*M_PI*(double)i/(double)Nxy;
		th2 = 2.0*M_PI*(double)ii/(double)Nxy;
		glBegin(GL_QUADS);
			glNormal3f((float)sin(th1),-(float)cos(th1),0.0f);
			glTexCoord2f((rds+p[i][1])/(2.0f*rds),1.0f);//ﾃｸｽﾁｬｰ座標でi=0のときs=0となるように
													//i=Nxy(右端)でrTexとなるように直径2.0*rdsで正規化
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
	//裏側(x < 0)にはﾏｯﾋﾟﾝｸﾞしない
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
//側面だけに円筒投影
void skTexCylinder2(int Nxy)
{
	//半径rds、高さhgtの円柱(中心が原点)
	//Nx--多角形の種類(Nx<0)
	float p[40][3]; //頂点座標
	double theta0,theta;
	double th1,th2;
	int i,ii;
	float rds = 0.5;//半径
	float hgt = 1.0;//高さ
	double pp = 2.0 * M_PI;

	if(Nxy > 20) { Nxy = 20;}
	theta0 = M_PI / (double)Nxy;
	for(i = 0;i < Nxy;i++)
	{   theta = 2.0 * theta0 * (double)i;
		//左端をi=0とする（ｵﾌﾞｼﾞｪｸﾄ自身から見て右端、-y軸方向)
//		p[i][0] = (float)(rds*sin(theta)); //上底のx成分
//		p[i][1] = -(float)(rds*cos(theta));//ｙ成分
		//背面をi=0とする
		p[i][0] = -(float)(rds*cos(theta)); //上底のx成分
		p[i][1] = -(float)(rds*sin(theta));//ｙ成分
		p[i][2] = (float)hgt/2.0f;         //ｚ成分(高さ)
		p[i+Nxy][0] = p[i][0];             //下底のx成分
		p[i+Nxy][1] = p[i][1];             //ｙ成分
		p[i+Nxy][2] = -(float)hgt/2.0f;    //ｚ成分
	}

	//上底0
	glBegin(GL_POLYGON);
		glNormal3f(0.0f,0.0f,1.0f);
		for(i=0;i<Nxy;i++) glVertex3fv(p[i]);
	glEnd();
	//下底
	glBegin(GL_POLYGON);
		glNormal3f(0.0f,0.0f,-1.0f);
		for(i=2*Nxy-1;i>=Nxy;i--) glVertex3fv(p[i]);
	glEnd();

	//円筒投影
	//側面一周にﾏｯﾋﾟﾝｸﾞ
	glEnable(GL_TEXTURE_2D);
	for(i = 0;i < Nxy;i++)
	{
		if(i == Nxy-1) ii = 0;
        else ii = i + 1;
		th1 = (double)i/(double)Nxy;//2πで正規化した角度
		th2 = (double)ii/(double)Nxy;
		if(i == Nxy-1) th2 = 1.0;
		glBegin(GL_QUADS);
			glNormal3f(-(float)cos(pp * th1),-(float)sin(pp * th1),0.0f);
			glTexCoord2f((float)th1, 1.0f);//ﾃｸｽﾁｬｰ座標でi=0のときs=0となるように
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
//格子状オブジェクトに対するテクスチャマッピング用サブルーチン
void skTexGrid(int N1, int N2, float pd[][3])
{	//裏面も定義
	int i, j, np;
	float n1[3], n2[3], n3[3], n4[3];
	float a[101][101], b[101][101], c[101][101];

	//法線成分
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
				calcNormal(pd[np],pd[np-1],pd[np+N1+1],n1);//左側
				calcNormal(pd[np],pd[np+N1+1],pd[np+1],n2);//右側
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
				calcNormal(pd[np],pd[np-N1-1],pd[np-1],n1);//左側
				calcNormal(pd[np],pd[np+1],pd[np-N1-1],n2);//右側
				a[i][j] = (n1[0]+n2[0])/2.0f;
				b[i][j] = (n1[1]+n2[1])/2.0f;
				c[i][j] = (n1[2]+n2[2])/2.0f; }
		}
		else
		{
			if(i == 0) {
				calcNormal(pd[np],pd[np+1],pd[np-N1-1],n1);//上
				calcNormal(pd[np],pd[np+N1+1],pd[np+1],n2);//下
				a[i][j] = (n1[0]+n2[0])/2.0f;
				b[i][j] = (n1[1]+n2[1])/2.0f;
				c[i][j] = (n1[2]+n2[2])/2.0f; }
			else if(i == N1) {
				calcNormal(pd[np],pd[np-N1-1],pd[np-1],n1);//上
				calcNormal(pd[np],pd[np-1],pd[np+N1+1],n2);//下
				a[i][j] = (n1[0]+n2[0])/2.0f;
				b[i][j] = (n1[1]+n2[1])/2.0f;
				c[i][j] = (n1[2]+n2[2])/2.0f; }
			else {//上下左右４個の三角形の平均
				calcNormal(pd[np],pd[np-N1-1],pd[np-1],n1);//左上
				calcNormal(pd[np],pd[np+1],pd[np-N1-1],n2);//右上
				calcNormal(pd[np],pd[np-1],pd[np+N1+1],n3);//左下
				calcNormal(pd[np],pd[np+N1+1],pd[np+1],n4);//右下
				a[i][j] = (n1[0]+n2[0]+n3[0]+n4[0])/4.0f;
				b[i][j] = (n1[1]+n2[1]+n3[1]+n4[1])/4.0f;
				c[i][j] = (n1[2]+n2[2]+n3[2]+n4[2])/4.0f; }
		}
	}

	//Texture描画
	glEnable(GL_TEXTURE_2D);
	for(j = 0;j < N2;j++)
	  for(i = 0;i < N1;i++){
		np = i + (N1+1) * j;
		glBegin(GL_TRIANGLES);
			//左下の三角形
			//各頂点の法線方向,ﾃｸｽﾁｬｰ座標,頂点座標を与える。
			glNormal3d(a[i][j],b[i][j],c[i][j]);//法線方向
			glTexCoord2f((float)i/(float)N1,(1.0f-(float)j/(float)N2));//ﾃｸｽﾁｬｰ座標
			glVertex3fv(pd[np]);//ﾎﾟﾘｺﾞﾝの頂点座標（以下これらを繰り返す）
			glNormal3d(a[i][j+1],b[i][j+1],c[i][j+1]);
			glTexCoord2f((float)i/(float)N1,(1.0f-(float)(j+1)/(float)N2));
			glVertex3fv(pd[np+N1+1]);
			glNormal3d(a[i+1][j+1],b[i+1][j+1],c[i+1][j+1]);
			glTexCoord2f((float)(i+1)/(float)N1,(1.0f-(float)(j+1)/(float)N2));
			glVertex3fv(pd[np+N1+2]);
			//右上の三角形
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
	//裏面も表面に描いたﾃｸｽﾁｬｰを描画
	for(j = 0;j < N2;j++)
	  for(i = 0;i < N1;i++){
		np = i+(N1+1)*j;
		glBegin(GL_TRIANGLES);
			//左下の三角形
			//各頂点の法線方向に負号,ﾃｸｽﾁｬｰ座標,頂点座標を与える（表とは逆の頂点列,2番目と3番目を交換）。
			glNormal3d(-a[i][j],-b[i][j],-c[i][j]);//法線方向
			glTexCoord2f((float)i/(float)N1,(1.0f-(float)j/(float)N2));//ﾃｸｽﾁｬｰ座標
			glVertex3fv(pd[np]);//ﾎﾟﾘｺﾞﾝの頂点座標（以下これらを繰り返す）
			glNormal3d(-a[i+1][j+1],-b[i+1][j+1],-c[i+1][j+1]);
			glTexCoord2f((float)(i+1)/(float)N1,(1.0f-(float)(j+1)/(float)N2));
			glVertex3fv(pd[np+N1+2]);
			glNormal3d(-a[i][j+1],-b[i][j+1],-c[i][j+1]);
			glTexCoord2f((float)i/(float)N1,(1.0f-(float)(j+1)/(float)N2));
			glVertex3fv(pd[np+N1+1]);
			//右上の三角形
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
//正方形のメッシュ（x-y平面,中心が原点）
//水面の波動などに適用（振幅に関係なくx-y軸方向の長さは一定）
void skTexMesh1(int Nx, float* data)
{
	//全体の幅,長さどちらも1(正方形)
	int i, j, Ny, np;
	float p[10201][3]; //頂点座標
	float pitchX, pitchY;


	if(Nx > 100) Nx = 100;
	Ny = Nx;
	//セルのサイズ
	pitchX = 1.0f / (float)Nx;
	pitchY = 1.0f / (float)Ny;

	//各頂点の座標
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
//四角い布地（x-y平面,中心が原点）
//Mesh自身の伸び縮みなし（振幅が大きくなればx-y方向が短くなる）
//振幅はｚ軸方向
void skTexMesh2(int Nx, double ratio, double amp,
		       double lamdaX, double lamdaY, double alpha)
{
	float pd[10201][3];//頂点座標
	double lengthX, lengthY;
	double thetaX, thetaY;
	double cx, cy, ex, ey;
	double x, y, z, dx, dy;
	int i, j, np, Ny;

	if(Nx > 100) Nx = 100;
	Ny = Nx;

	//全長
	lengthX = 1.0;
	lengthY = lengthX * ratio;//x方向の長さ

	//基本周期長

	//1要素のｻｲｽﾞ(グリッド間隔)
	ex = lengthX / (double)Nx;
	ey = lengthY / (double)Ny;

	cx = 2.0 * M_PI / lamdaX;
	cy = 2.0 * M_PI / lamdaY;

//	alpha = (float)(alpha * M_PI/180.0);//初期位相をradに変換
	//座標値
	x = 0.0; y = 0.0;
	for(j = 0; j <= Nx; j++){//jが-x方向
		thetaX = atan(amp * cx * cos(cx * x + cy * y - alpha));
		dx = ex * cos(thetaX);//ｚ方向増分
		x += dx;
		y = 0.0;
		for(i = 0; i <= Ny; i++){//iはy方向
			thetaY = atan(amp * cy * cos(cx * x + cy * y - alpha));
			dy = ey * cos(thetaY);//y方向増分
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
//四角い布地（x-z平面,Flagの左下が原点）
//Mesh自身の伸び縮みなし（振幅が大きくなればx方向が短縮）
//風による横方向伸縮あり
//変位はy軸方向
void skTexMesh3(int Nx, double ratio, double amp0, double lambda0, double alpha,
                double wind, double z0)
{
	float pd[10201][3];//頂点座標
	double lengthX, lengthZ;
	double thetaX, thetaZ, ampl, fct, xShift;
	double cx, cz, ex, ez;
	double x, y, z, dx, dz, y0;
	int i, j, np, Nz;

    double amp = amp0 / (1.0f + 50.0 * wind );//風の強いとき振幅を小さく
    double lambdaX = lambda0 / (1.0f + 5.0 * wind);//風の強いときx方向波長は短く
    double lambdaZ = lambda0*(1.0 + 5.0 * wind);        //ｚ方向波長は長く
    double weight =  1.0f / (1.0 + 20.0 * wind); //重力重み比率
	double gamma = weight * (80.0 * M_PI /180.0); //weight=1で全体が80度傾く

	if(Nx > 100) Nx = 100;
	Nz = Nx;

	//全長
	lengthZ = 1.0;//高さ(z方向)
	lengthX = lengthZ * ratio;//長さ(x方向)

	//1要素のｻｲｽﾞ(グリッド間隔)
	ex = lengthX / (double)Nx;
	ez = lengthZ / (double)Nz;

	cx = (2.0 * M_PI) / lambdaX;
	cz = (2.0 * M_PI) / lambdaZ;
//cz=0.0;
	//座標値
	z = 1.0; x = 0.0;
	for(j = 0;j <= Nz;j++){//jが-z方向
		x = 0.0;
        xShift = wind * (0.3 * sin(M_PI*z)) ;
//		ampl = amp * (1.0 + x);
		y0 = amp * sin(cz * z - alpha) ;//poleにおける元の位置
//		y0 = amp * sin(alpha - cz * z) ;//poleにおける元の位置
        fct = 1.0 + (1.0 - sin(M_PI * z) ) * wind * 0.2;//横方向サイズ係数
		for(i = 0;i <= Nx;i++){//iはx方向
			ampl = amp * (0.f + x );//xが大きいほど振幅を大きくなるように
            //変位(x=0の位置でpoleに固定できるようにy0で差し引く)
			y = ( ampl* sin(  cx * x * fct + cz * z  - alpha) - y0 ) ;
			np = j * (Nx+1) + i;
			pd[np][1] = (float)y;
			//重力の影響
			pd[np][0] = (float)((x * cos(gamma) + xShift) * fct) ;
			pd[np][2] = (float)(z - x * sin(gamma));
            if(pd[np][2] < -(float)z0 + 0.01f) pd[np][2] = - (float)z0 + 0.01f;//床面より浮かせる
            //次の格子点のx座標
			thetaX = atan(ampl * cx * cos( cx * x + cz * z - alpha));
//			thetaX = atan(ampl * cx * cos(alpha - cx * x - cz * z));
			dx = ex * cos(thetaX) ;//x方向増分
			x += dx;

		}
		thetaZ = atan(amp * cz * cos(cx * x + cz * z - alpha));
		dz = ez * cos(thetaZ);//ｚ方向増分
		z -= dz;
	}
	skTexGrid(Nx, Nz, pd);
}

//-----------------------------------------------------------------
//超2次関数
void skTexSuper1(int Nxy, int Nz, double eps1, double eps2, double p1, double p2, double p3)
{
	//上下の中心が原点
	int i,j,ip,im,np,npL,npR,npU,npD,k1,k2;
	double ct,phai,theta,z,fz;
	float a[31][31],b[31][31],c[31][31];
	float n1[3],n2[3],n3[3],n4[3];
	double cc;
	float pd[961][3];

	if( Nxy > 30 ) Nxy = 30;
	if( Nz > 30 ) Nz = 30;

	//半径
    float r = 0.5f;
	//上半分
    for(j = 0 ;j <= Nz;j++)
	{
        phai = (M_PI/(double)Nz) * ((double)Nz / 2.0 - (double)j);
        if(phai > 0.0) z = (float)(pow(sin(phai),eps1));//z
        else z = - (float)(pow(sin(fabs(phai)), eps1));
        //形状関数
        if(z < 0.0) fz = (p2-p3)*z + p2;
        else fz = (p1-p2)*z + p2;

        for (i = 0 ;i <= Nxy / 2;i++)
        {
			k1 = Nxy * j + i;//自分から見て左側
			k2 = Nxy * j + Nxy - i;//右側
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


	//側面の法線成分
	for(i = 0;i < Nxy;i++){
 		ip = i+1;
		if(ip == Nxy) ip = 0;
		im = i-1;
		if(i == 0) im = Nxy-1;

		//真上(Top)
		a[i][0] = 0.0; b[i][0] = 0.0; c[i][0] = 1.0;
		//真下（Bottom)
		a[i][Nz] = 0.0; b[i][Nz] = 0.0; c[i][Nz] = -1.0;

		for(j=1;j<Nz;j++)//隣り合う4個の三角形の法線ベクトルを平均化
		{
			np = j*Nxy+i;//注目点
			npL = j*Nxy+im;//左側
			npR = j*Nxy+ip;//右側
			npU = np-Nxy;//上
			npD = np+Nxy;//下
			if(j == 1) {
				n1[0]=0.0; n1[1]=0.0; n1[2]=1.0;//Top
				n2[0]=0.0; n2[1]=0.0; n2[2]=1.0;//Top
				calcNormal(pd[np],pd[npL],pd[npD],n3);//外から見て左下
				calcNormal(pd[np],pd[npD],pd[npR],n4);//右下
			}
			if(j == Nz-1){
				calcNormal(pd[np],pd[npU],pd[npL],n1);//外から見て左上
				calcNormal(pd[np],pd[npR],pd[npU],n2);//右上
				n3[0]=0.0; n3[1]=0.0; n3[2]=-1.0;//Bottom
				n4[0]=0.0; n4[1]=0.0; n4[2]=-1.0;//Bottom
			}
			else {
				calcNormal(pd[np],pd[npU],pd[npL],n1);//外から見て左上
				calcNormal(pd[np],pd[npR],pd[npU],n2);//右上
				calcNormal(pd[np],pd[npL],pd[npD],n3);//外から見て左下
				calcNormal(pd[np],pd[npD],pd[npR],n4);//右下
			}
			a[i][j] = (float)((n1[0]+n2[0]+n3[0]+n4[0])/4.0f);//ｘ方向
			b[i][j] = (float)((n1[1]+n2[1]+n3[1]+n4[1])/4.0f);//ｙ
			c[i][j] = (float)((n1[2]+n2[2]+n3[2]+n4[2])/4.0f);//ｚ
		}
	}

	//表示
	//側面(4角形パッチ）
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
//２次元格子状オブジェクトに対するテクスチャマッピング用サブルーチン
void skTexGridSquare(int numRow, int numCol, float pd[][3])
{	//裏面も定義
	int i, j, np;
	float n1[3], n2[3], n3[3], n4[3];
	float a[20][20], b[20][20], c[20][20];

	//法線成分
	for(i = 0;i < numRow;i++) //下方向
	  for(j = 0;j < numCol;j++)//右方向
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
				calcNormal(pd[np],pd[np-1],pd[np+numCol],n1);//左側
				calcNormal(pd[np],pd[np+numCol],pd[np+1],n2);//右側
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
				calcNormal(pd[np],pd[np-numCol],pd[np-1],n1);//左側
				calcNormal(pd[np],pd[np+1],pd[np-numCol],n2);//右側
				a[i][j] = (n1[0]+n2[0])/2.0f;
				b[i][j] = (n1[1]+n2[1])/2.0f;
				c[i][j] = (n1[2]+n2[2])/2.0f; }
		}
		else
		{
			if(j == 0) {
				calcNormal(pd[np],pd[np+1],pd[np-numCol],n1);//上
				calcNormal(pd[np],pd[np+numCol],pd[np+1],n2);//下
				a[i][j] = (n1[0]+n2[0])/2.0f;
				b[i][j] = (n1[1]+n2[1])/2.0f;
				c[i][j] = (n1[2]+n2[2])/2.0f; }
			else if(j == numCol-1) {
				calcNormal(pd[np],pd[np-numCol],pd[np-1],n1);//上
				calcNormal(pd[np],pd[np-1],pd[np+numCol],n2);//下
				a[i][j] = (n1[0]+n2[0])/2.0f;
				b[i][j] = (n1[1]+n2[1])/2.0f;
				c[i][j] = (n1[2]+n2[2])/2.0f; }
			else {//上下左右４個の三角形の平均
				calcNormal(pd[np],pd[np-numCol],pd[np-1],n1);//左上
				calcNormal(pd[np],pd[np+1],pd[np-numCol],n2);//右上
				calcNormal(pd[np],pd[np-1],pd[np+numCol],n3);//左下
				calcNormal(pd[np],pd[np+numCol],pd[np+1],n4);//右下
				a[i][j] = (n1[0]+n2[0]+n3[0]+n4[0])/4.0f;
				b[i][j] = (n1[1]+n2[1]+n3[1]+n4[1])/4.0f;
				c[i][j] = (n1[2]+n2[2]+n3[2]+n4[2])/4.0f; }
		}
	}

	//Texture描画
	glEnable(GL_TEXTURE_2D);
	for(i = 0;i < numRow-1;i++)
	  for(j = 0;j < numCol-1;j++){
		//np = i + (N1+1) * j;
        np = i * numCol + j;
		glBegin(GL_TRIANGLES);
			//左下の三角形
			//各頂点の法線方向,ﾃｸｽﾁｬｰ座標,頂点座標を与える。
			glNormal3d(a[i][j],b[i][j],c[i][j]);//法線方向
			glTexCoord2f((float)j/(float)(numCol-1),(1.0f-(float)i/(float)(numRow-1)));//ﾃｸｽﾁｬｰ座標
			glVertex3fv(pd[np]);//ﾎﾟﾘｺﾞﾝの頂点座標（以下これらを繰り返す）
			glNormal3d(a[i+1][j],b[i+1][j],c[i+1][j]);
			glTexCoord2f((float)j/(float)(numCol-1),(1.0f-(float)(i+1)/(float)(numRow-1)));
			glVertex3fv(pd[np+numCol]);
			glNormal3d(a[i+1][j+1],b[i+1][j+1],c[i+1][j+1]);
			glTexCoord2f((float)(j+1)/(float)(numCol-1),(1.0f-(float)(i+1)/(float)(numRow-1)));
			glVertex3fv(pd[np+numCol+1]);
			//右上の三角形
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
	//裏面も表面に描いたﾃｸｽﾁｬｰを描画
	for(i = 0;i < numRow-1;i++)
	  for(j = 0;j < numCol-1;j++){
        np = i * numCol + j;
		glBegin(GL_TRIANGLES);
			//左下の三角形
			//各頂点の法線方向に負号,ﾃｸｽﾁｬｰ座標,頂点座標を与える（表とは逆の頂点列,2番目と3番目を交換）。
			glNormal3d(-a[i][j],-b[i][j],-c[i][j]);//法線方向
			glTexCoord2f((float)j/(float)(numCol-1),(1.0f-(float)i/(float)(numRow-1)));//ﾃｸｽﾁｬｰ座標
			glVertex3fv(pd[np]);//ﾎﾟﾘｺﾞﾝの頂点座標（以下これらを繰り返す）
			glNormal3d(-a[i+1][j+1],-b[i+1][j+1],-c[i+1][j+1]);
			glTexCoord2f((float)(j+1)/(float)(numCol-1),(1.0f-(float)(i+1)/(float)(numRow-1)));
			glVertex3fv(pd[np+numCol+1]);
			glNormal3d(-a[i+1][j],-b[i+1][j],-c[i+1][j]);
			glTexCoord2f((float)j/(float)(numCol-1),(1.0f-(float)(i+1)/(float)(numRow-1)));
			glVertex3fv(pd[np+numCol]);
			//右上の三角形
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
    float shadowDiffuse[] = {0.2f,0.2f,0.2f,0.3f};//影の拡散光

	//法線成分
	for(i = 0;i < numRow;i++) //下方向
	  for(j = 0;j < numCol;j++)//右方向
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
				calcNormal(pd[np],pd[np-1],pd[np+numCol],n1);//左側
				calcNormal(pd[np],pd[np+numCol],pd[np+1],n2);//右側
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
				calcNormal(pd[np],pd[np-numCol],pd[np-1],n1);//左側
				calcNormal(pd[np],pd[np+1],pd[np-numCol],n2);//右側
				a[i][j] = (n1[0]+n2[0])/2.0f;
				b[i][j] = (n1[1]+n2[1])/2.0f;
				c[i][j] = (n1[2]+n2[2])/2.0f; }
		}
		else
		{
			if(j == 0) {
				calcNormal(pd[np],pd[np+1],pd[np-numCol],n1);//上
				calcNormal(pd[np],pd[np+numCol],pd[np+1],n2);//下
				a[i][j] = (n1[0]+n2[0])/2.0f;
				b[i][j] = (n1[1]+n2[1])/2.0f;
				c[i][j] = (n1[2]+n2[2])/2.0f; }
			else if(j == numCol-1) {
				calcNormal(pd[np],pd[np-numCol],pd[np-1],n1);//上
				calcNormal(pd[np],pd[np-1],pd[np+numCol],n2);//下
				a[i][j] = (n1[0]+n2[0])/2.0f;
				b[i][j] = (n1[1]+n2[1])/2.0f;
				c[i][j] = (n1[2]+n2[2])/2.0f; }
			else {//上下左右４個の三角形の平均
				calcNormal(pd[np],pd[np-numCol],pd[np-1],n1);//左上
				calcNormal(pd[np],pd[np+1],pd[np-numCol],n2);//右上
				calcNormal(pd[np],pd[np-1],pd[np+numCol],n3);//左下
				calcNormal(pd[np],pd[np+numCol],pd[np+1],n4);//右下
				a[i][j] = (n1[0]+n2[0]+n3[0]+n4[0])/4.0f;
				b[i][j] = (n1[1]+n2[1]+n3[1]+n4[1])/4.0f;
				c[i][j] = (n1[2]+n2[2]+n3[2]+n4[2])/4.0f; }
		}
	}
	//Texture描画
    if(flagShadow == true) goto SHADOW_DISP;
	glEnable(GL_TEXTURE_2D);
	for(i = 0;i < numRow-1;i++)
	  for(j = 0;j < numCol-1;j++){
        np = i * numCol + j;
		glBegin(GL_TRIANGLES);
			//左下の三角形
			//各頂点の法線方向,ﾃｸｽﾁｬｰ座標,頂点座標を与える。
			glNormal3d(a[i][j],b[i][j],c[i][j]);//法線方向
			glTexCoord2f((float)j/(float)(numCol-1),(1.0f-(float)i/(float)(numRow-1)));//ﾃｸｽﾁｬｰ座標
			glVertex3fv(pd[np]);//ﾎﾟﾘｺﾞﾝの頂点座標（以下これらを繰り返す）
			glNormal3d(a[i+1][j],b[i+1][j],c[i+1][j]);
			glTexCoord2f((float)j/(float)(numCol-1),(1.0f-(float)(i+1)/(float)(numRow-1)));
			glVertex3fv(pd[np+numCol]);
			glNormal3d(a[i+1][j+1],b[i+1][j+1],c[i+1][j+1]);
			glTexCoord2f((float)(j+1)/(float)(numCol-1),(1.0f-(float)(i+1)/(float)(numRow-1)));
			glVertex3fv(pd[np+numCol+1]);
			//右上の三角形
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
			glNormal3d(a[i][j],b[i][j],c[i][j]);//法線方向
			glVertex3fv(pd[np]);//ﾎﾟﾘｺﾞﾝの頂点座標（以下これらを繰り返す）
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
//立方体格子状オブジェクトに対するテクスチャマッピング用サブルーチン
void skTexGridCube(int numRow, int numCol, int numStk, float p1[][3],
        float p2[][3], float p3[][3], float p4[][3], float p5[][3],
        float p6[][3], bool flagShadow)
{
    //すべて同じテクスチャー
    //上側
    skTexGrid2(numRow, numCol, p1, flagShadow);
    //正面
    skTexGrid2(numStk, numCol, p2, flagShadow);
    //右側
    skTexGrid2(numStk, numRow, p3, flagShadow);
    //左側
    skTexGrid2(numStk, numRow, p4, flagShadow);
    //裏側
    skTexGrid2(numStk, numCol, p5, flagShadow);
    //下側
    skTexGrid2(numRow, numCol, p6, flagShadow);
}
//----------------------------------------------------------------------------
//球格子状オブジェクトに対するテクスチャマッピング用サブルーチン
//2004.07.13改良
//上下の頂点も縦方向の個数に含める
void skTexGridSphere(int Nxy, int Nz, float pd[][3], bool flagShadow)
{
	int i,j,ip,im,np,npL,npR,npU,npD;
//	double ct,phai,theta,z,fz;
	float a[10][10],b[10][10],c[10][10];//頂点の法線成分
	float n[10][3];//頂点の法線成分
    float shadowDiffuse[] = {0.2f,0.2f,0.2f,0.3f};//影の拡散光

    if(Nxy > 10) Nxy = 10;
    if(Nz > 10) Nz = 10;
    //j=0のきpd[0]～pd[Nxy-1]は同じTop頂点
    //j=Nz-1のときnp=(Nz-1)*Nxyとして
    //pd[np]～pd[np+numCol-1]は同じBottom頂点

    //Top
    //TopのnumCol個の三角形の法線ベクトルを平均化
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
    for(i = 1;i < Nxy; i++){ //Top頂点はすべて同じ
        a[i][0] = a[0][0]; b[i][0] = b[0][0]; c[i][0] = c[0][0];
    }

	//側面の法線成分
    //隣り合う4個の三角形の法線ベクトルを平均化
	for(i = 0;i < Nxy;i++){
 		ip = i+1;
		if(ip == Nxy) ip = 0;
		im = i-1;
		if(i == 0) im = Nxy-1;
		for(j = 1;j < Nz-1;j++)
		{
			np =  j * Nxy + i;//注目点
			npL = j * Nxy + im;//左側
			npR = j * Nxy + ip;//右側
			npU = np - Nxy;//上
			npD = np + Nxy;//下
            calcNormal(pd[np],pd[npU],pd[npL],n[0]);//外から見て左上
            calcNormal(pd[np],pd[npR],pd[npU],n[1]);//右上
            calcNormal(pd[np],pd[npL],pd[npD],n[2]);//外から見て左下
            calcNormal(pd[np],pd[npD],pd[npR],n[3]);//右下

			a[i][j] = (float)((n[0][0]+n[1][0]+n[2][0]+n[3][0])/4.0f);//ｘ方向
			b[i][j] = (float)((n[0][1]+n[1][1]+n[2][1]+n[3][1])/4.0f);//ｙ
			c[i][j] = (float)((n[0][2]+n[1][2]+n[2][2]+n[3][2])/4.0f);//ｚ
		}
	}
    //Bottom
    //BottomのnumCol個の三角形の法線ベクトルを平均化
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
    for(i = 1;i < Nxy; i++){ //Top頂点はすべて同じ
        a[i][Nz-1] = a[0][Nz-1];
        b[i][Nz-1] = b[0][Nz-1];
        c[i][Nz-1] = c[0][Nz-1];
    }

	//表示
    if(flagShadow == true) goto SHADOW_DISP;
	//側面(4角形パッチ）
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
//側面だけに円筒投影
void skTexGridCylinder(int Nxy, int Nz, float p[][3], bool flagShadow)
{
    //p[0]は上底の中心，最後のp[]は下底の中心
	int i,j,ip,im,np,npL,npR,npU,npD;
//	double ct,theta,z,fz;
	float a[11][10],b[11][10],c[11][10];//頂点の法線成分
    float shadowDiffuse[] = {0.2f,0.2f,0.2f,0.3f};//影の拡散光

	float n[11][3];//頂点の法線成分

    if(Nxy > 10) Nxy = 10;
    if(Nz > 10) Nz = 10;

    //全頂点数
    int num = Nxy * Nz + 2;

    //Top
    //上底中心はTopのnumCol個の三角形の法線ベクトルを平均化
    a[0][0] = 0.0; b[0][0] = 0.0; c[0][0] = 0.0;
    //上底の周囲の頂点は隣り合う2つの三角形の法線ベクトルを平均化
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
    //Topの表示
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
    //下底中心はBottomのnumCol個の三角形の法線ベクトルを平均化
    a[0][0] = 0.0; b[0][0] = 0.0; c[0][0] = 0.0;
    //下底の周囲の頂点は隣り合う2つの三角形の法線ベクトルを平均化
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
    //Bottomの表示
    for(i = 1; i <= Nxy; i++){
        ip = i + 1;
        if(ip == Nxy + 1) ip = 1;
		glBegin(GL_TRIANGLES);
			glNormal3f(a[0][0],b[0][0],c[0][0]);glVertex3fv(p[num-1]);
			glNormal3f(a[i][0],b[i][0],c[i][0]);glVertex3fv(p[num-1-i]);
			glNormal3f(a[ip][0],b[ip][0],c[ip][0]); glVertex3fv(p[num-1-ip]);
		glEnd();
	}

	//側面の法線成分
    //上底の周囲頂点の側面に対する法線成分は
    //隣り合う2個の三角形の法線ベクトルを平均化
	for(i = 0;i < Nxy;i++){
 		ip = i+1;
		if(ip == Nxy) ip = 0;
		im = i-1;
		if(i == 0) im = Nxy-1;
        //j = 0;
        np = i+1;//注目点
        npL = im+1;//左側
        npR = ip+1;//右側
        //npU = np - Nxy;//上
        npD = np + Nxy;//下
        calcNormal(p[np],p[npL],p[npD],n[0]);//外から見て左下
        calcNormal(p[np],p[npD],p[npR],n[1]);//右下

        a[i][0] = (float)((n[0][0]+n[1][0])/2.0f);//ｘ方向
        b[i][0] = (float)((n[0][1]+n[1][1])/2.0f);//ｙ
        c[i][0] = (float)((n[0][2]+n[1][2])/2.0f);//ｚ
	}

    //同じように下底の頂点に対して
	for(i = 0;i < Nxy;i++){
 		ip = i+1;
		if(ip == Nxy) ip = 0;
		im = i-1;
		if(i == 0) im = Nxy-1;
        j = Nz - 1;
        np = j*Nxy+i+1;//注目点
        npL = j*Nxy+im+1;//左側
        npR = j*Nxy+ip+1;//右側
        npU = np - Nxy;//上
        //npD = np + Nxy;//下
        calcNormal(p[np],p[npU],p[npL],n[0]);//外から見て左上
        calcNormal(p[np],p[npR],p[npU],n[1]);//右上

        a[i][j] = (float)((n[0][0]+n[1][0])/2.0f);//ｘ方向
        b[i][j] = (float)((n[0][1]+n[1][1])/2.0f);//ｙ
        c[i][j] = (float)((n[0][2]+n[1][2])/2.0f);//ｚ
	}

    //他の側面の頂点に対しては隣り合う4個の三角形の法線ベクトルを平均化
	for(i = 0;i < Nxy;i++){
 		ip = i+1;
		if(ip == Nxy) ip = 0;
		im = i-1;
		if(i == 0) im = Nxy-1;
		for(j = 1;j < Nz-1;j++)
		{
			np =  j * Nxy + i + 1;//注目点
			npL = j * Nxy + im + 1;//左側
			npR = j * Nxy + ip + 1;//右側
			npU = np - Nxy;//上
			npD = np + Nxy;//下
            calcNormal(p[np],p[npU],p[npL],n[0]);//外から見て左上
            calcNormal(p[np],p[npR],p[npU],n[1]);//右上
            calcNormal(p[np],p[npL],p[npD],n[2]);//外から見て左下
            calcNormal(p[np],p[npD],p[npR],n[3]);//右下

			a[i][j] = (float)((n[0][0]+n[1][0]+n[2][0]+n[3][0])/4.0f);//ｘ方向
			b[i][j] = (float)((n[0][1]+n[1][1]+n[2][1]+n[3][1])/4.0f);//ｙ
			c[i][j] = (float)((n[0][2]+n[1][2]+n[2][2]+n[3][2])/4.0f);//ｚ
		}
	}

	//側面(4角形パッチ）
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
//2次元格子状オブジェクトに格子間隔に従って赤,シアンのチェック模様を直接描画
void skDrawCheck2(int numRow, int numCol, float pd[][3], bool flagShadow)
{
	int i, j, np;
	float n1[3], n2[3], n3[3], n4[3];
	float a[10][10], b[10][10], c[10][10];
    float shadowDiffuse[] = {0.2f,0.2f,0.2f,0.3f};//影の拡散光

	//法線成分
	for(i = 0;i < numRow;i++) //下方向
	  for(j = 0;j < numCol;j++)//右方向
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
				calcNormal(pd[np],pd[np-1],pd[np+numCol],n1);//左側
				calcNormal(pd[np],pd[np+numCol],pd[np+1],n2);//右側
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
				calcNormal(pd[np],pd[np-numCol],pd[np-1],n1);//左側
				calcNormal(pd[np],pd[np+1],pd[np-numCol],n2);//右側
				a[i][j] = (n1[0]+n2[0])/2.0f;
				b[i][j] = (n1[1]+n2[1])/2.0f;
				c[i][j] = (n1[2]+n2[2])/2.0f; }
		}
		else
		{
			if(j == 0) {
				calcNormal(pd[np],pd[np+1],pd[np-numCol],n1);//上
				calcNormal(pd[np],pd[np+numCol],pd[np+1],n2);//下
				a[i][j] = (n1[0]+n2[0])/2.0f;
				b[i][j] = (n1[1]+n2[1])/2.0f;
				c[i][j] = (n1[2]+n2[2])/2.0f; }
			else if(j == numCol-1) {
				calcNormal(pd[np],pd[np-numCol],pd[np-1],n1);//上
				calcNormal(pd[np],pd[np-1],pd[np+numCol],n2);//下
				a[i][j] = (n1[0]+n2[0])/2.0f;
				b[i][j] = (n1[1]+n2[1])/2.0f;
				c[i][j] = (n1[2]+n2[2])/2.0f; }
			else {//上下左右４個の三角形の平均
				calcNormal(pd[np],pd[np-numCol],pd[np-1],n1);//左上
				calcNormal(pd[np],pd[np+1],pd[np-numCol],n2);//右上
				calcNormal(pd[np],pd[np-1],pd[np+numCol],n3);//左下
				calcNormal(pd[np],pd[np+numCol],pd[np+1],n4);//右下
				a[i][j] = (n1[0]+n2[0]+n3[0]+n4[0])/4.0f;
				b[i][j] = (n1[1]+n2[1]+n3[1]+n4[1])/4.0f;
				c[i][j] = (n1[2]+n2[2]+n3[2]+n4[2])/4.0f; }
		}
	}

    float diffuse[4] ;
    int k;
	//チェック模様描画
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
			//左下の三角形
			//各頂点の法線方向,ﾃｸｽﾁｬｰ座標,頂点座標を与える。
			glNormal3d(a[i][j],b[i][j],c[i][j]);//法線方向
			glVertex3fv(pd[np]);//ﾎﾟﾘｺﾞﾝの頂点座標（以下これらを繰り返す）
			glNormal3d(a[i+1][j],b[i+1][j],c[i+1][j]);
			glVertex3fv(pd[np+numCol]);
			glNormal3d(a[i+1][j+1],b[i+1][j+1],c[i+1][j+1]);
			glVertex3fv(pd[np+numCol+1]);
			//右上の三角形
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
    //すべて同じテクスチャー
    //上側
    skDrawCheck2(numRow, numCol, p1, flagShadow);
    //正面
    skDrawCheck2(numStk, numCol, p2, flagShadow);
    //右側
    skDrawCheck2(numStk, numRow, p3, flagShadow);
    //左側
    skDrawCheck2(numStk, numRow, p4, flagShadow);
    //裏側
    skDrawCheck2(numStk, numCol, p5, flagShadow);
    //下側
    skDrawCheck2(numRow, numCol, p6, flagShadow);
}
//------------------------------------------------------------------------------
//球格子状オブジェクトに対するサブルーチン
//格子間隔に従って赤，シアンのチェック模様を直接描画
void skDrawCheckSphere(int Nxy, int Nz, float pd[][3], bool flagShadow)
{
	int i,j,ip,im,np,npL,npR,npU,npD;
	float a[10][10],b[10][10],c[10][10];//頂点の法線成分
    float shadowDiffuse[] = {0.2f,0.2f,0.2f,0.3f};//影の拡散光

	float n[10][3];//頂点の法線成分

    if(Nxy > 10) Nxy = 10;
    if(Nz > 10) Nz = 10;

    //Top
    //TopのnumCol個の三角形の法線ベクトルを平均化
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
    for(i = 1;i < Nxy; i++){ //Top頂点はすべて同じ
        a[i][0] = a[0][0]; b[i][0] = b[0][0]; c[i][0] = c[0][0];
    }

	//側面の法線成分
    //隣り合う4個の三角形の法線ベクトルを平均化
	for(i = 0;i < Nxy;i++){
 		ip = i+1;
		if(ip == Nxy) ip = 0;
		im = i-1;
		if(i == 0) im = Nxy-1;
		for(j = 1;j < Nz-1;j++)
		{
			np =  j * Nxy + i;//注目点
			npL = j * Nxy + im;//左側
			npR = j * Nxy + ip;//右側
			npU = np - Nxy;//上
			npD = np + Nxy;//下
            calcNormal(pd[np],pd[npU],pd[npL],n[0]);//外から見て左上
            calcNormal(pd[np],pd[npR],pd[npU],n[1]);//右上
            calcNormal(pd[np],pd[npL],pd[npD],n[2]);//外から見て左下
            calcNormal(pd[np],pd[npD],pd[npR],n[3]);//右下

			a[i][j] = (float)((n[0][0]+n[1][0]+n[2][0]+n[3][0])/4.0f);//ｘ方向
			b[i][j] = (float)((n[0][1]+n[1][1]+n[2][1]+n[3][1])/4.0f);//ｙ
			c[i][j] = (float)((n[0][2]+n[1][2]+n[2][2]+n[3][2])/4.0f);//ｚ
		}
	}
    //Bottom
    //BottomのnumCol個の三角形の法線ベクトルを平均化
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
    for(i = 1;i < Nxy; i++){ //Top頂点はすべて同じ
        a[i][Nz-1] = a[0][Nz-1];
        b[i][Nz-1] = b[0][Nz-1];
        c[i][Nz-1] = c[0][Nz-1];
    }

	//表示
    if(flagShadow == true) goto SHADOW_DISP;
    float diffuse[4];
	//側面(4角形パッチ）
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
//円柱格子状オブジェクトに対するサブルーチン
//格子間隔に従って赤，シアンのチェック模様を直接描画
void skDrawCheckCylinder(int Nxy, int Nz, float p[][3], bool flagShadow)
{
	int i,j,ip,im,np,npL,npR,npU,npD;
	float a[11][10],b[11][10],c[11][10];//頂点の法線成分
	float n[11][3];//頂点の法線成分
    float diffuse[4];
    float shadowDiffuse[] = {0.2f,0.2f,0.2f,0.3f};//影の拡散光

    if(Nxy > 10) Nxy = 10;
    if(Nz > 10) Nz = 10;

    //全頂点数
    int num = Nxy * Nz + 2;

    //Top
    //上底中心はTopのnumCol個の三角形の法線ベクトルを平均化
    a[0][0] = 0.0; b[0][0] = 0.0; c[0][0] = 0.0;
    //上底の周囲の頂点は隣り合う2つの三角形の法線ベクトルを平均化
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
    //Topの表示
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
    //下底中心はBottomのnumCol個の三角形の法線ベクトルを平均化
    a[0][0] = 0.0; b[0][0] = 0.0; c[0][0] = 0.0;
    //下底の周囲の頂点は隣り合う2つの三角形の法線ベクトルを平均化
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
    //Bottomの表示
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

	//側面の法線成分
    //上底の周囲頂点の側面に対する法線成分は
    //隣り合う2個の三角形の法線ベクトルを平均化
	for(i = 0;i < Nxy;i++){
 		ip = i+1;
		if(ip == Nxy) ip = 0;
		im = i-1;
		if(i == 0) im = Nxy-1;
        //j = 0;
        np = i+1;//注目点
        npL = im+1;//左側
        npR = ip+1;//右側
        //npU = np - Nxy;//上
        npD = np + Nxy;//下
        calcNormal(p[np],p[npL],p[npD],n[0]);//外から見て左下
        calcNormal(p[np],p[npD],p[npR],n[1]);//右下

        a[i][0] = (float)((n[0][0]+n[1][0])/2.0f);//ｘ方向
        b[i][0] = (float)((n[0][1]+n[1][1])/2.0f);//ｙ
        c[i][0] = (float)((n[0][2]+n[1][2])/2.0f);//ｚ
	}

    //同じように下底の頂点に対して
	for(i = 0;i < Nxy;i++){
 		ip = i+1;
		if(ip == Nxy) ip = 0;
		im = i-1;
		if(i == 0) im = Nxy-1;
        j = Nz - 1;
        np = j*Nxy+i+1;//注目点
        npL = j*Nxy+im+1;//左側
        npR = j*Nxy+ip+1;//右側
        npU = np - Nxy;//上
        //npD = np + Nxy;//下
        calcNormal(p[np],p[npU],p[npL],n[0]);//外から見て左上
        calcNormal(p[np],p[npR],p[npU],n[1]);//右上

        a[i][j] = (float)((n[0][0]+n[1][0])/2.0f);//ｘ方向
        b[i][j] = (float)((n[0][1]+n[1][1])/2.0f);//ｙ
        c[i][j] = (float)((n[0][2]+n[1][2])/2.0f);//ｚ
	}

    //他の側面の頂点に対しては隣り合う4個の三角形の法線ベクトルを平均化
	for(i = 0;i < Nxy;i++){
 		ip = i+1;
		if(ip == Nxy) ip = 0;
		im = i-1;
		if(i == 0) im = Nxy-1;
		for(j = 1;j < Nz-1;j++)
		{
			np =  j * Nxy + i + 1;//注目点
			npL = j * Nxy + im + 1;//左側
			npR = j * Nxy + ip + 1;//右側
			npU = np - Nxy;//上
			npD = np + Nxy;//下
            calcNormal(p[np],p[npU],p[npL],n[0]);//外から見て左上
            calcNormal(p[np],p[npR],p[npU],n[1]);//右上
            calcNormal(p[np],p[npL],p[npD],n[2]);//外から見て左下
            calcNormal(p[np],p[npD],p[npR],n[3]);//右下

			a[i][j] = (float)((n[0][0]+n[1][0]+n[2][0]+n[3][0])/4.0f);//ｘ方向
			b[i][j] = (float)((n[0][1]+n[1][1]+n[2][1]+n[3][1])/4.0f);//ｙ
			c[i][j] = (float)((n[0][2]+n[1][2]+n[2][2]+n[3][2])/4.0f);//ｚ
		}
	}

    //表示
    if(flagShadow == true) goto SHADOW_DISP;
	//側面(4角形パッチ）
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
//以下はソリッドテクスチャ用
//------------------------------------------------------------------------------
//texImageはCTargetクラスで計算
//----------------------------------------------------------------------------
//立方体のソリッドテクスチャ
void skSolidTexCube(unsigned char texImage[6][T_MAX][T_MAX][3])
{
	//この内部でテクスチャを作成
	float p[8][3]={{0.5f,0.5f,0.5f},{-0.5f,0.5f,0.5f},{-0.5f,-0.5f,0.5f},
	{0.5f,-0.5f,0.5f},{0.5f,0.5f,-0.5f},{-0.5f,0.5f,-0.5f},{-0.5f,-0.5f,0-0.5f},{0.5f,-0.5f,-0.5f}};
	glTexImage2D(GL_TEXTURE_2D,0,3,T_MAX,T_MAX,0,GL_RGB,GL_UNSIGNED_BYTE,texImage[0]);

	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);//top
		glNormal3f(0.0f,0.0f,1.0f); //z方向
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
	 	glNormal3f(-1.0f,0.0f,0.0f); //-x方向
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
		glNormal3f(0.0f,-1.0f,0.0f); //-y方向
		glTexCoord2f(0.0f, 1.0f); glVertex3fv(p[2]);
		glTexCoord2f(0.0f, 0.0f); glVertex3fv(p[6]);
		glTexCoord2f(1.0f, 0.0f); glVertex3fv(p[7]);
		glTexCoord2f(1.0f, 1.0f); glVertex3fv(p[3]);
	glEnd();
	glDisable(GL_TEXTURE_2D);

	glTexImage2D(GL_TEXTURE_2D,0,3,T_MAX,T_MAX,0,GL_RGB,GL_UNSIGNED_BYTE,texImage[5]);

	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);
		glNormal3f(0.0f,0.0f,-1.0f); //-z方向
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

	//描画
	skTexSphere2(N);
}
//-------------------------------------------------------------------------------------
void skSolidTexCylinder(int Nxy, byte texImage[6][T_MAX][T_MAX][3])
{
	int i, ii;

	//半径rds、高さhgtの円柱(中心が原点)
	//Nx--多角形の種類(Nx<0)
	float p[40][3]; //頂点座標
	double theta0,theta;
	double th1,th2;
	float rds = 0.5;//半径
	float hgt = 1.0;//高さ
	double pp = 2.0 * M_PI;

	//頂点座標
	if(Nxy > 20) { Nxy = 20;}
	theta0 = pp / (double)Nxy;
	for(i = 0;i < Nxy;i++)
	{   theta = theta0 * (double)i;
		//左端をi=0とする（ｵﾌﾞｼﾞｪｸﾄ自身から見て右端、-y軸方向)
		p[i][0] = (float)(rds*sin(theta)); //上底のx成分
		p[i][1] = -(float)(rds*cos(theta));//ｙ成分
		p[i][2] = (float)hgt/2.0f;         //ｚ成分(高さ)
		p[i+Nxy][0] = p[i][0];             //下底のx成分
		p[i+Nxy][1] = p[i][1];             //ｙ成分
		p[i+Nxy][2] = -(float)hgt/2.0f;    //ｚ成分
	}

	//Top
	glTexImage2D(GL_TEXTURE_2D,0,3,T_MAX,T_MAX,0,GL_RGB,GL_UNSIGNED_BYTE,texImage[0]);
		
	glEnable(GL_TEXTURE_2D);
	glBegin(GL_POLYGON);
		glNormal3f(0.0f,0.0f,1.0f); //z方向
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
		glNormal3f(0.0f,0.0f,-1.0f); //-z方向
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
		th1 = (double)i/(double)Nxy;//2πで正規化した角度
		th2 = (double)ii/(double)Nxy;
		if(i == Nxy-1) th2 = 1.0;
		glBegin(GL_QUADS);
			glNormal3f((float)sin(pp * th1),-(float)cos(pp * th1), 0.0f);
			glTexCoord2f((float)th1, 1.0f);//ﾃｸｽﾁｬｰ座標でi=0のときs=0となるように
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

