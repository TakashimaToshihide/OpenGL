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
	//幅、高さ１の立方体(立方体の中心を原点)
	float p[8][3]={{0.5f,0.5f,0.5f},{-0.5f,0.5f,0.5f},{-0.5f,-0.5f,0.5f},
	{0.5f,-0.5f,0.5f},{0.5f,0.5f,-0.5f},{-0.5f,0.5f,-0.5f},{-0.5f,-0.5f,0-0.5f},
	{0.5f,-0.5f,-0.5f}};

	glBegin(GL_QUADS);
		//z方向
		glVertex3fv(p[0]); glVertex3fv(p[1]);
		glVertex3fv(p[2]); glVertex3fv(p[3]);
		//x方向
		glVertex3fv(p[0]); glVertex3fv(p[3]);
		glVertex3fv(p[7]); glVertex3fv(p[4]);
		//y方向
		glVertex3fv(p[0]); glVertex3fv(p[4]);
		glVertex3fv(p[5]); glVertex3fv(p[1]);
		//-x方向
		glVertex3fv(p[1]); glVertex3fv(p[5]);
		glVertex3fv(p[6]); glVertex3fv(p[2]);
		//-y方向
		glVertex3fv(p[2]); glVertex3fv(p[6]);
		glVertex3fv(p[7]); glVertex3fv(p[3]);
		//-z方向
		glVertex3fv(p[4]); glVertex3fv(p[7]);
		glVertex3fv(p[6]); glVertex3fv(p[5]);
	glEnd();
}
//---------------------------------------------------------------------------
void skSolidCube(void)
{   //幅、高さ１の立方体
	//各点の座標(立方体の中心を原点)
	float p[8][3]={{0.5f,0.5f,0.5f},{-0.5f,0.5f,0.5f},{-0.5f,-0.5f,0.5f},
	{0.5f,-0.5f,0.5f},{0.5f,0.5f,-0.5f},{-0.5f,0.5f,-0.5f},{-0.5f,-0.5f,0-0.5f},{0.5f,-0.5f,-0.5f}};

	glBegin(GL_QUADS);
		glNormal3f(0.0f,0.0f,1.0f); //z方向
		glVertex3fv(p[0]); glVertex3fv(p[1]);
		glVertex3fv(p[2]); glVertex3fv(p[3]);

		glNormal3f(1.0f,0.0f,0.0f); //x方向
		glVertex3fv(p[0]); glVertex3fv(p[3]);
		glVertex3fv(p[7]); glVertex3fv(p[4]);

		glNormal3f(0.0f,1.0f,0.0f); //y方向
		glVertex3fv(p[0]); glVertex3fv(p[4]);
		glVertex3fv(p[5]); glVertex3fv(p[1]);

		glNormal3f(-1.0f,0.0f,0.0f); //-x方向
		glVertex3fv(p[1]); glVertex3fv(p[5]);
		glVertex3fv(p[6]); glVertex3fv(p[2]);

		glNormal3f(0.0f,-1.0f,0.0f); //-y方向
		glVertex3fv(p[2]); glVertex3fv(p[6]);
		glVertex3fv(p[7]); glVertex3fv(p[3]);

		glNormal3f(0.0f,0.0f,-1.0f); //-z方向
		glVertex3fv(p[4]); glVertex3fv(p[7]);
		glVertex3fv(p[6]); glVertex3fv(p[5]);
	glEnd();
}

//--------------------------------------------------------------------------
void skWireSphere(int Nxy)//単位直径の球
{
	int i,j,Nz;
//	int np[21][21];//np[17][17];//i,j点の頂点番号 (i:x-y断面分割点  j:ｚ軸方向分割点)
	double phai; //x-y平面に対する偏角（緯度）
	double theta; //x軸に対する偏角（経度)
	float p[21][21][3];//p[289][3]; //各点の座標
	float p1[3],p2[3],p3[3],p4[3];

	//if(Nxy > 16) Nxy = 16;
	if(Nxy > 20) Nxy = 20;
	Nz = Nxy;
	//頂点番号、座標
	for(i=0;i<=Nxy;i++)
	{	theta = 2.0*M_PI*(double)i/(double)Nxy;
		for(j=0;j<=Nz;j++)
		{	
			phai = M_PI/2.0-M_PI*(double)j/(double)Nz;
			p[i][j][0] = (float)(0.5*cos(theta)*cos(phai)); //x座標
			p[i][j][1] = (float)(0.5*sin(theta)*cos(phai)); //y
			p[i][j][2] = (float)(0.5*sin(phai));            //z
		}
	}

	//頂点列を定義
	for(i=0;i<Nxy;i++)
		for(j=0;j<Nz-1;j++)
		{ //面の頂点
			p1[0] = p[i][j][0]; //x座標
			p1[1] = p[i][j][1]; //y座標
			p1[2] = p[i][j][2]; //z座標
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
	//底
	j = Nz-1;
	for(i=0;i<Nxy;i++)
		{ //面の頂点
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
void skSolidSphere(int Nxy)//単位直径の球
{  
	int i,j,Nz;
//	int np[21][21];//i,j点の頂点番号 (i:x-y断面分割点  j:ｚ軸方向分割点)
	double phai; //x-y平面に対する偏角（緯度）
	double theta; //x軸に対する偏角（経度)
	float p[21][21][3]; //各点の座標
	float p1[3],p2[3],p3[3],p4[3];

	if(Nxy > 21) Nxy = 20;
	Nz = Nxy;
	//頂点番号,座標
	for(i=0;i<=Nxy;i++)
	{	theta = 2.0*M_PI*(double)i/(double)Nxy;
		for(j=0;j<=Nz;j++)
		{
			phai = M_PI/2.0-M_PI*(double)j/(double)Nz;
			p[i][j][0] = (float)(0.5*cos(theta)*cos(phai)); //x座標
			p[i][j][1] = (float)(0.5*sin(theta)*cos(phai)); //y
			p[i][j][2] = (float)(0.5*sin(phai));            //z
		}
	}
    
	//各パッチを描画
	for(i=0;i<Nxy;i++)
		for(j=0;j<Nz-1;j++)
		{ //面の頂点   x座標                       y座標                       z座標
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
	//底
	j = Nz-1;
	for(i=0;i<Nxy;i++)
		{ //面の頂点   x座標                       y座標                       z座標
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
void skSolidSemiSphere(int Nxy)//単位直径の半球
{  
	int i,j,Nz;
	double phai; //x-y平面に対する偏角（緯度）
	double theta; //x軸に対する偏角（経度)
	float p[21][11][3]; //各点の座標
	float p1[3],p2[3],p3[3],p4[3];

	if(Nxy > 21) Nxy = 20;
	Nz = Nxy/2;
	//頂点番号、座標
	for(i=0;i<=Nxy;i++)
	{	theta = 2.0*M_PI*(double)i/(double)Nxy;
		for(j=0;j<=Nz;j++)
		{	
			phai = M_PI*(1.0-(double)j/(double)Nz)/2.0;
			p[i][j][0] = (float)(0.5*cos(theta)*cos(phai)); //x座標
			p[i][j][1] = (float)(0.5*sin(theta)*cos(phai)); //y
			p[i][j][2] = (float)(0.5*sin(phai));            //z
		}
	}
    
	//各パッチの描画
	for(i=0;i<Nxy;i++)
		for(j=0;j<Nz;j++)
		{ //面の頂点   x座標                       y座標                       z座標
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
	//底
	j = Nz;
	glBegin(GL_POLYGON);
		glNormal3f(0.0f,0.0f,-1.0f);//-z方向
		for(i=Nxy;i>0;i--) glVertex3fv(p[i][j]);
	glEnd();
}
//--------------------------------------------------------------------------
void skWireSemiSphere(int Nxy)//単位直径の半球
{  
	int i,j,Nz;
	double phai; //x-y平面に対する偏角（緯度）
	double theta; //x軸に対する偏角（経度)
	float p[21][11][3]; //各点の座標
	float p1[3],p2[3],p3[3],p4[3];

	if(Nxy > 21) Nxy = 20;
	Nz = Nxy/2;
	//頂点番号、座標
	for(i=0;i<=Nxy;i++)
	{	theta = 2.0*M_PI*(double)i/(double)Nxy;
		for(j=0;j<=Nz;j++)
		{	
			phai = M_PI*(1.0-(double)j/(double)Nz)/2.0;
			p[i][j][0] = (float)(0.5*cos(theta)*cos(phai)); //x座標
			p[i][j][1] = (float)(0.5*sin(theta)*cos(phai)); //y
			p[i][j][2] = (float)(0.5*sin(phai));            //z
		}
	}
    
	//各パッチの描画
	for(i=0;i<Nxy;i++)
		for(j=0;j<Nz;j++)
		{ //面の頂点   x座標                       y座標                       z座標
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
	//底
	j = Nz;
	glBegin(GL_POLYGON);
		for(i=Nxy;i>0;i--) glVertex3fv(p[i][j]);
	glEnd();
}

//----------------------------------------------------------------------
void skSolidPrism(int Nxy)//単位直径、単位高さの多角柱
{	
	//Nxy--分割数
	float p[41][3]; //頂点座標
	float nn[3];
	double theta0,theta;
	int i,ii;

	if(Nxy>21) { Nxy = 20;}
	theta0 = 2*M_PI/(double)Nxy;
	for(i=0;i<Nxy;i++)
	{   theta = theta0*(double)i;
		p[i][0] = (float)(0.5*cos(theta)); //上底のx成分
		p[i][1] = (float)(0.5*sin(theta)); //ｙ成分
		p[i][2] = 0.5f;                //ｚ成分(高さ)
		p[i+Nxy][0] = (float)(0.5*cos(theta)); //下底のx成分
		p[i+Nxy][1] = (float)(0.5*sin(theta)); //ｙ成分
		p[i+Nxy][2] = -0.5f;               //ｚ成分
	}

	//上底
	glBegin(GL_POLYGON);
		glNormal3f(0.0f,0.0f,1.0f);
		for(i=0;i<Nxy;i++) glVertex3fv(p[i]); 
	glEnd();
	//下底
	glBegin(GL_POLYGON);
		glNormal3f(0.0f,0.0f,-1.0f); 
		for(i=2*Nxy-1;i>=Nxy;i--) glVertex3fv(p[i]); 
	glEnd();

	//側面
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
void skWirePrism(int Nxy)//単位直径、単位高さの多角柱
{	
	//Nxy--分割数
	float p[41][3]; //頂点座標
	double theta0,theta;
	int i,ii;

	if(Nxy>21) { Nxy = 20;}
	theta0 = 2*M_PI/(double)Nxy;
	for(i=0;i<Nxy;i++)
	{   theta = theta0*(double)i;
		p[i][0] = (float)(0.5*cos(theta)); //上底のx成分
		p[i][1] = (float)(0.5*sin(theta)); //ｙ成分
		p[i][2] = 0.5f;                //ｚ成分(高さ)
		p[i+Nxy][0] = (float)(0.5*cos(theta)); //下底のx成分
		p[i+Nxy][1] = (float)(0.5*sin(theta)); //ｙ成分
		p[i+Nxy][2] = -0.5f;               //ｚ成分
	}

	//上底
	glBegin(GL_POLYGON);
		for(i=0;i<Nxy;i++) glVertex3fv(p[i]);
	glEnd();
	//下底
	glBegin(GL_POLYGON);
		for(i=2*Nxy-1;i>=Nxy;i--) glVertex3fv(p[i]);
	glEnd();

	//側面
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
void skSolidCylinder(int Nxy)//単位直径、単位高さ
{
	//Nxy--分割数
	float p[41][3]; //頂点座標
	double theta0,theta;
	int i,ii;

	if(Nxy>21) { Nxy = 20;}
	theta0 = 2*M_PI/(double)Nxy;
	for(i=0;i<Nxy;i++)
	{   theta = theta0*(double)i;
		p[i][0] = (float)(0.5*cos(theta)); //上底のx成分
		p[i][1] = (float)(0.5*sin(theta)); //ｙ成分
		p[i][2] = 0.5f;                //ｚ成分(高さ)
		p[i+Nxy][0] = (float)(0.5*cos(theta)); //下底のx成分
		p[i+Nxy][1] = (float)(0.5*sin(theta)); //ｙ成分
		p[i+Nxy][2] = -0.5f;               //ｚ成分
	}

	//上底
	glBegin(GL_POLYGON);
		glNormal3f(0.0f,0.0f,1.0f);
		for(i=0;i<Nxy;i++) glVertex3fv(p[i]); 
	glEnd();
	//下底
	glBegin(GL_POLYGON);
		glNormal3f(0.0f,0.0f,-1.0f); 
		for(i = 2*Nxy-1;i >= Nxy;i--) glVertex3fv(p[i]); 
	glEnd();

	//側面
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
void skWireCylinder(int Nxy)//単位直径、単位高さ
{	
	//Nxy--分割数
	float p[41][3]; //頂点座標
	double theta0,theta;
	int i,ii;

	if(Nxy>20) { Nxy = 20;}
	theta0 = 2*M_PI/(double)Nxy;
	for(i=0;i<Nxy;i++)
	{   theta = theta0*(double)i;
		p[i][0] = (float)(0.5*cos(theta)); //上底のx成分
		p[i][1] = (float)(0.5*sin(theta)); //ｙ成分
		p[i][2] = 0.5f;                //ｚ成分(高さ)
		p[i+Nxy][0] = (float)(0.5*cos(theta)); //下底のx成分
		p[i+Nxy][1] = (float)(0.5*sin(theta)); //ｙ成分
		p[i+Nxy][2] = -0.5f;               //ｚ成分
	}

	//上底
	glBegin(GL_POLYGON);
		for(i=0;i<Nxy;i++) glVertex3fv(p[i]); 
	glEnd();
	//下底
	glBegin(GL_POLYGON);
		for(i = 2*Nxy-1;i >= Nxy;i--) glVertex3fv(p[i]); 
	glEnd();

	//側面
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
void skSolidPyramid(int Nxy)//多角錘
{   //下底の中心を物体の中心
	//直径、高さ１の台形
	float nn[3];//法線成分
	float p[21][3]; //頂点座標
	double r; //半径
	double theta0,theta;
	int i,ii;

	if(Nxy >21) { Nxy = 20; }
	theta0 = 2.0*M_PI/(double)Nxy;
	r=0.5;
	//頂点
	p[0][0] =p[0][1] = 0.0;
	p[0][2] = 1.0;
	for(i=1;i<=Nxy;i++)
	{   theta = theta0*(double)i;
		p[i][0] = (float)(r*cos(theta)); //下底のx成分
		p[i][1] = (float)(r*sin(theta)); //y成分
		p[i][2] = 0.0f;                  //ｚ成分
	}

	//下底
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
void skWirePyramid(int Nxy)//多角錘
{   //下底の中心を物体の中心
	//直径、高さ１の台形
	float p[21][3]; //頂点座標
	double r; //半径
	double theta0,theta;
	int i,ii;

	if(Nxy >21) { Nxy = 20; }
	theta0 = 2.0*M_PI/(double)Nxy;
	r=0.5;
	//頂点
	p[0][0] =p[0][1] = 0.0;
	p[0][2] = 1.0;
	for(i=1;i<=Nxy;i++)
	{   theta = theta0*(double)i;
		p[i][0] = (float)(r*cos(theta)); //下底のx成分
		p[i][1] = (float)(r*sin(theta)); //y成分
		p[i][2] = 0.0f;                  //ｚ成分
	}

	//下底
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
void skSolidCone(int Nxy)//円錘
{   //下底の中心を物体の中心
	//直径、高さ１の台形
	float nx1,nx2,ny1,ny2,nxy,nz;//法線成分
	float p[21][3]; //頂点座標
	double r; //半径
	double theta0,theta,alpha;
	int i,ii;

	if(Nxy >21) { Nxy = 21; }
	theta0 = 2.0*M_PI/(double)Nxy;
	r=0.5;
	//頂点座標
	p[0][0] =p[0][1] = 0.0;
	p[0][2] = 1.0;
	//底の座標
	for(i=1;i<=Nxy;i++)
	{   theta = theta0*(double)i;
		p[i][0] = (float)(r*cos(theta)); //底のx成分
		p[i][1] = (float)(r*sin(theta)); //y成分
		p[i][2] = 0.0f;                  //ｚ成分
	}

	//底
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
			//底の頂点1
			theta = theta0*(double)i;
			nx1 = nxy*(float)cos(theta);
			ny1 = nxy*(float)sin(theta);
			glNormal3d(nx1,ny1,nz); glVertex3fv(p[i]);
			//底の頂点2
			theta = theta0*(double)ii;
			nx2 = nxy*(float)cos(theta); 
			ny2 = nxy*(float)sin(theta);
			glNormal3d(nx2,ny2,nz); glVertex3fv(p[ii]);
			//頂点
			nx1 = (nx1+nx2)/2.0f; 
			ny1 = (ny1+ny2)/2.0f;
			glNormal3d(nx1,ny1,nz); glVertex3fv(p[0]);			
		glEnd();
	}
}
//-----------------------------------------------------------------------
void skWireCone(int Nxy)//円錘
{   //下底の中心を物体の中心
	//直径、高さ１の台形
	float p[21][3]; //頂点座標
	double r; //半径
	double theta0,theta;
	int i,ii;

	if(Nxy >21) { Nxy = 21; }
	theta0 = 2.0*M_PI/(double)Nxy;
	r=0.5;
	//頂点座標
	p[0][0] =p[0][1] = 0.0;
	p[0][2] = 1.0;
	//底の座標
	for(i=1;i<=Nxy;i++)
	{   theta = theta0*(double)i;
		p[i][0] = (float)(r*cos(theta)); //底のx成分
		p[i][1] = (float)(r*sin(theta)); //y成分
		p[i][2] = 0.0f;                  //ｚ成分
	}

	//底
	glBegin(GL_POLYGON);
		for(i=Nxy;i>=1;i--) glVertex3fv(p[i]); 
	glEnd();

	for(i=1;i<=Nxy;i++)
	{
		if(i==Nxy) ii=1; else ii=i+1;
		glBegin(GL_TRIANGLES);
			//底の頂点1
			glVertex3fv(p[i]);
			//底の頂点2
			glVertex3fv(p[ii]);
			//頂点
			glVertex3fv(p[0]);
		glEnd();
	}
}
//-----------------------------------------------------------------------
void skSolidFrustum(int Nxy,float ratio)//円錐台
{   //下底の中心を物体の中心
	//r1:下底の半径，r2:上底の半径
	//直径、高さ１の台形,上底の幅ｻｲｽﾞはratio倍,
	double nn[3];//法線成分
	float p[61][3]; //頂点座標
	float r1, r2;
	double theta0, theta, alpha, nxy;
	int i,ii;

	if(Nxy >31) { Nxy = 31; }
	theta0 = 2.0*M_PI/(double)Nxy;
	r1 = 0.5;
	r2 = ratio * r1;
	for(i = 0;i < Nxy;i++)
	{   theta = theta0*(double)i;
		p[i][0] = r2*(float)cos(theta); //上底のx成分
		p[i][1] = r2*(float)sin(theta); //ｙ成分
		p[i][2] = 1.0f;                 //ｚ成分
		p[i+Nxy][0] = r1*(float)cos(theta); //下底のx成分
		p[i+Nxy][1] = r1*(float)sin(theta); //y成分
		p[i+Nxy][2] = 0.0f;                  //ｚ成分
	}

	//上底
	glBegin(GL_POLYGON);
		glNormal3f(0.0f,0.0f,1.0f);
		for(i=0;i<Nxy;i++) glVertex3fv(p[i]);
	glEnd();
	//下底
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
void skWireFrustum(int Nxy, float ratio)//円錐台
{   //下底の中心を物体の中心
	//r1:下底の半径，r2:上底の半径
	//直径、高さ１の台形,上底の幅ｻｲｽﾞはratio倍,
	float p[61][3]; //頂点座標
	float r1, r2;
	double theta0, theta;//, alpha;
	int i,ii;

	if(Nxy >31) { Nxy = 31; }
	theta0 = 2.0*M_PI/(double)Nxy;
	r1 = 0.5;
	r2 = ratio * r1;
	for(i = 0;i < Nxy;i++)
	{   theta = theta0*(double)i;
		p[i][0] = r2*(float)cos(theta); //上底のx成分
		p[i][1] = r2*(float)sin(theta); //ｙ成分
		p[i][2] = 1.0f;                 //ｚ成分
		p[i+Nxy][0] = r1*(float)cos(theta); //下底のx成分
		p[i+Nxy][1] = r1*(float)sin(theta); //y成分
		p[i+Nxy][2] = 0.0f;                  //ｚ成分
	}

	//上底
	glBegin(GL_POLYGON);
		for(i=0;i<Nxy;i++) glVertex3fv(p[i]); 
	glEnd();
	//下底
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
//厚みの無い正方形
void skSolidSquare(void)
{
	float p[4][3];
	p[0][0] = 0.5f;  p[0][1] = 0.5f;  p[0][2] = 0.0f;
	p[1][0] = -0.5f; p[1][1] = 0.5f;  p[1][2] = 0.0f;
	p[2][0] = -0.5f; p[2][1] = -0.5f; p[2][2] = 0.0f;
	p[3][0] = 0.5f;  p[3][1] = -0.5f; p[3][2] = 0.0f;

	glBegin(GL_QUADS);
		glNormal3f(0.0f,0.0f,1.0f);//ｚ方向
		glVertex3fv(p[0]); glVertex3fv(p[1]);
		glVertex3fv(p[2]); glVertex3fv(p[3]);
	glEnd();
}
//--------------------------------------------------------
//厚みの無い正方形
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
//厚みの無い方形枠(外側下部の長さは１、高さ１、法線方向ｘ方向,底はz=0）
void skSolidFrame0(double wTop, double wFrame)
{   //wFrame:frameの幅(0.5以下)
	//wTop:上部の幅（下部1に対する）
	float p[8][3];

	//座標
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
//厚みの無い方形枠(外側下部の長さは１、高さ１、法線方向ｘ方向,底はz=0）
void skWireFrame0(double wTop, double wFrame)
{   //wFrame:frameの幅(0.5以下)
	//wTop:上部の幅（下部1に対する）
	float p[8][3];

	//座標
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
void skSolidTorus(int Nm, int Ns, float ratio)//2003.3.6変更
{	//Nm:主円周分割点数,Ns:断面円周分割点数
	int i,ii,j,jj;
	double phai,ph1,ph2; //主円周分割点のx軸に対する偏角
	double theta,th1,th2; //断面円周分割点のx-z平面に対する偏角
	float rr[21], zz[21];
	float p[21][21][3]; //各点の座標
	float p1[3],p2[3],p3[3],p4[3];
	float r1, r2;

	r1 = 0.5;//主円周軸半径
	r2 = ratio * r1;//断面半径
	//if(r1 < r2) { MessageBox(NULL,"skSolidTorusにおいて r1 < r2 となりました ","半径の値",MB_OK);return;}
	if( Nm > 21 ) Nm = 21;// { MessageBox(NULL,"skSolidTorusにおいてNm>20 となりました ","Nmの値",MB_OK);return;}
	if( Ns > 21 ) Ns = 21;//{ MessageBox(NULL,"skSolidTorusにおいてNc>20 となりました ","Ncの値",MB_OK);return;}

	//基本断面(x-z)の座標
	for(j = 0; j < Ns; j++)
	{	theta = M_PI-2.0*M_PI*(double)j/(double)Ns;
		rr[j] = r1 + r2 * (float)cos(theta); //原点からの距離
		zz[j] = r2 * (float)sin(theta);//ｚ
	}
	//他の断面の座標
	for(i = 0; i < Nm; i++)
	{	phai = 2.0*M_PI*(double)i/(double)Nm;
		for(j = 0; j < Ns; j++)
		{	
			p[i][j][0] = rr[j] * (float)cos(phai); //x座標
			p[i][j][1] = rr[j] * (float)sin(phai); //y
			p[i][j][2] = zz[j];   //z
		}
	}

	//頂点列を定義し描画
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
			//面の頂点   x座標                       y座標                       z座標
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
//自在管(基本姿勢でｚ軸が管の軸方向）
//直径はdata[0],data[1]
//要素の長さは可変,data[2]から
void skSolidTube(int Nxy,int Nz, float *data)
{
	//3軸回転（dataの，角度単位はdeg）
	//Nxy--xy断面分割数
	//Nz---軸方向分割数
	//data---回転角度(最初の2個はx,y方向の直径，2番目からNz個は要素の長さ）
	//xyza---上底の中心軸座標と方向（階層構造とするときの戻り値)
	float p[101][21][3]; //頂点座標
	double x0, y0, z0, x, y, z, xx, yy, zz;
	float xyz[2][3];
	float pitch[101];
	float nn[3];
	float a[101][20], b[101][20], c[101][20];//側面の法線
	double sX, sY;//x,y方向直径
	double theta0,theta;
	double alpha0, beta0, gamma0;//1つ前の中心軸の回転角度
	double alpha, beta, gamma;//現在の節面の回転角度
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
	x0 = 0.0; y0 = 0.0; z0 = 0.0;//下底の中心
	//k = 0 のとき
	theta0 = 2*M_PI/(double)Nxy;
	for(i = 0;i < Nxy; i++)
	{   theta = theta0 * (double)i;
		//各頂点の座標
		p[0][i][0] =  (float)(0.5*cos(theta) * sX); //x成分
		p[0][i][1] =  (float)(0.5*sin(theta) * sY); //ｙ成分
		p[0][i][2] = 0.0f;                          //ｚ成分
		//法線方向
		a[0][i] = p[0][i][0];
		b[0][i] = p[0][i][1];
		c[0][i] = p[0][i][2];
	}

	for( k = 0; k < Nz; k++){
		alpha0 += data[3*k+Nz+3] * pp;//中心軸はdata[3*k]ずつx軸回転（累積される)
		if(k == Nz) alpha = alpha0 ;//+ data[3*k+3] * pp;
		else alpha = alpha0 + data[3*k+Nz+6] * pp / 2.0;
		beta0 += data[3*k+Nz+4] * pp;
		if(k == Nz) beta = beta0 ;//+ data[3*k+4] * pp;
		else beta = beta0 + data[3*k+Nz+7] * pp / 2.0;
		gamma0 += data[3*k+Nz+5] * pp;
		if(k == Nz) gamma = gamma0 ;//+ data[3*k+5] * pp;
		else gamma = gamma0 + data[3*k+Nz+8] * pp / 2.0;

		//中心の平行移動量(元はx=0,y=0)
		x = 0.0; y = 0.0; z = pitch[k]; //pitch[0]=0としておくこと
		//x軸回転
		xx = x;
		yy = -z * (float)sin(alpha0);
		zz = z * (float)cos(alpha0) ;
		//y軸回転
		x = xx; y = yy; z = zz;
		xx = z * (float)sin(beta0);
		yy = y;
		zz = z * (float)cos(beta0);
		//z軸回転
		x = xx; y = yy; z = zz;
		xx = x * (float)cos(gamma0) - y * (float)sin(gamma0);
		yy = x * (float)sin(gamma0) + y * (float)cos(gamma0);
		zz = z;
		//平行移動
		x0 += xx; y0 += yy; z0 += zz;

		if(k == Nz){//上底の中心座標
			xyz[1][0] = (float)x0;
			xyz[1][1] = (float)y0;
			xyz[1][2] = (float)z0;
		}
		for(i = 0;i < Nxy; i++)
		{
			theta = theta0 * (double)i;
			//平行移動前の回転前の各頂点座標（x-y平面)
			x = (float)(0.5 * cos(theta)) * sX;
			y = (float)(0.5 * sin(theta)) * sY;
			z = 0.0f;
			//alphaだけx 軸回転
			xx = x ;
			yy = y * (float)cos(alpha) - z * (float)sin(alpha);
			zz = y * (float)sin(alpha) + z * (float)cos(alpha);
			//betaだけｙ軸回転
			x = xx; y = yy; z = zz;
			xx = x * (float)cos(beta) + z * (float)sin(beta);
			yy = y;
			zz = -x * (float)sin(beta) + z * (float)cos(beta);
			//gammaだけz 軸回転
			x = xx; y = yy; z = zz;
			xx = x * (float)cos(gamma) - y * (float)sin(gamma);
			yy = x * (float)sin(gamma) + y * (float)cos(gamma);
			zz = z;
			//法線方向は回転だけ
			a[k][i] = (float)xx;
			b[k][i] = (float)yy;
			c[k][i] = (float)zz;
			//中心座標だけ平行移動
			p[k][i][0] = (float)(xx + x0);
			p[k][i][1] = (float)(yy + y0);
			p[k][i][2] = (float)(zz + z0);
		}//i
	}//k
/*	//上端の位置，回転角（階層構造にするための戻り値）
	xyza[0] = (float)x0;
	xyza[1] = (float)y0;
	xyza[2] = (float)z0;
	xyza[3] = (float)alpha;
	xyza[4] = (float)beta;
	xyza[5] = (float)gamma; */
	//上底
	glBegin(GL_POLYGON);
		calcNormal(xyz[1], p[Nz][0], p[Nz][1],nn);
		glNormal3f(nn[0],nn[1],nn[2]);
		for(i = 0;i < Nxy; i++) glVertex3fv(p[Nz][i]);
	glEnd();
	//下底
	xyz[0][0] = 0.0; xyz[0][1] = 0.0; xyz[0][2] = 0.0;
	glBegin(GL_POLYGON);
		calcNormal(xyz[0], p[0][1], p[0][0],nn);
		glNormal3f(nn[0],nn[1],nn[2]);
		for(i = Nxy-1;i >= 0; i--) glVertex3fv(p[0][i]);
	glEnd();

	//側面
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
	//Nm:主円周分割点数,Ns:断面円周分割点数,Np:ピッチ数
	//radius:主円周軸半径,ratio:断面半径/主円周軸半径,length:全長
	int i, ii, j, jj, k, k1, k2;
	double phai, ph1, ph2; //主円周分割点のx軸に対する偏角
	double theta, th1, th2; //断面円周分割点のx-z平面に対する偏角
	float p[21][21][101][3];//p[44100][3]; //各点の座標
	float p1[3], p2[3], p3[3], p4[3];
	double rr[21], zz[21];
	double pitch, dp, hh;
	double r1, r2;

	r1 = radius; //主円周軸半径
	r2 = ratio * r1;//断面半径
	if( Nm > 20 ) Nm = 20;
	if( Ns > 20 ) Ns = 20;
	if( Np > 100 ) Np = 100;
	pitch = length / (double)Np;
    //縮んだときの制限
	if(pitch < 2 * r2) pitch = 2.0f * r2 ;

	dp = pitch / (double)Nm;

	//基本断面(x-z)の座標
	for(j = 0; j < Ns; j++)
	{	theta = M_PI - 2.0 * M_PI*(double)j/(double)Ns;
		rr[j] = r1 + r2 * cos(theta); //原点からの距離
		zz[j] = r2 * sin(theta);//ｚ
	}

	//他の断面の座標
	hh = 0;
	for(k = 0; k < Np; k++)
		for(i = 0; i < Nm; i++)
		{	phai = 2.0 * M_PI * (double)i/(double)Nm;
			for(j = 0; j < Ns; j++)
			{
				p[i][j][k][0] = (float)(rr[j] * cos(phai)); //x座標
				p[i][j][k][1] = (float)(rr[j] * sin(phai)); //y
				p[i][j][k][2] = (float)(zz[j] + hh) ;              //z
			}
			hh += dp;//中心軸の高さがdpずつ増加
		}

	//最終端(k=Np-1,i=Nm)
	k = Np - 1; i = Nm;
	for(j = 0; j < Ns; j++){
		phai = 0.0;
		p[i][j][k][0] = (float)(rr[j] * cos(phai)); //x座標
		p[i][j][k][1] = (float)(rr[j] * sin(phai)); //y
		p[i][j][k][2] = (float)(zz[j] + hh) ;            //z
	}

	//頂点列を定義し描画
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
			//面の頂点   x座標                       y座標                       z座標
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
	//始端（k=0)
	glBegin(GL_POLYGON);
		glNormal3d(0.0,-1.0,0.0);
		for(j = Ns-1; j >= 0; j--) glVertex3fv(p[0][j][0]);
	glEnd();
	//終端（k=Np)
	glBegin(GL_POLYGON);
		glNormal3d(0.0,1.0,0.0);
		for(j = 0; j < Ns; j++) glVertex3fv(p[Nm][j][Np-1]);
	glEnd();

}
//-----------------------------------------------------------------
//超2次関数
void skSolidSuper1(int Nxy, int Nz, double eps1, double eps2, double p1, double p2, double p3)
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
	//上底
	for(i = 0;i < Nxy;i++)
	{	ip = i+1;
		if(ip == Nxy) ip = 0;

		glBegin(GL_TRIANGLES);
			glNormal3d(a[i][0],b[i][0],c[i][0]);glVertex3fv(pd[i]);
			glNormal3d(a[i][1],b[i][1],c[i][1]);glVertex3fv(pd[Nxy+i]);
			glNormal3d(a[ip][1],b[ip][1],c[ip][1]);glVertex3fv(pd[Nxy+ip]);
		glEnd();
	}
	//下底
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
	//側面(4角形パッチ）
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
//上半分の超2次関数
void skSolidSuper2(int Nxy, int Nz, double eps1, double eps2, double p1, double p2)
{
	//底の中心が原点
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
	//上半分だけ

    for(j = 0 ;j <= Nz ;j++)
	{
        phai = (M_PI/(double)Nz) * ((double)Nz  - (double)j)/2.0;
		z = (float)(pow(sin(phai),eps1));//z
//		fz = (taper - 1.0) * z + 1.0;
//		fz = 0.5f * (taper + 1.0f + (taper -1.0f) * z);
        //形状関数
        fz = (p1-p2)*z + p2;
        for (i = 0 ;i<= Nxy / 2;i++)
        {
			k1 = Nxy * j + i;//外から見て右側
			k2 = Nxy * j + Nxy - i;//左側
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

	//側面の法線成分
	for(i = 0;i < Nxy; i++){
 	  ip = i + 1;
	  if(ip==Nxy) ip = 0;
	  im=i-1;
	  if(i==0) im = Nxy - 1;

		//真上(Top)
		a[i][0] = 0.0; b[i][0] = 0.0; c[i][0] = 1.0;
	  for(j=1;j<Nz;j++)
	  {
		np = j*Nxy+i;//注目点
		npL = j*Nxy+im;//左側
		npR = j*Nxy+ip;//右側
		npU = np-Nxy;//上
		npD = np+Nxy;//下
		if(j==1) {
			n1[0]=0.0; n1[1]=0.0; n1[2]=1.0;
			n2[0]=0.0; n2[1]=0.0; n2[2]=1.0;
			calcNormal(pd[np],pd[npL],pd[npD],n3);//外から見て左下
			calcNormal(pd[np],pd[npD],pd[npR],n4);//右下
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
	  j = Nz;//一番下の側面(上の2個の三角形の平均）
		calcNormal(pd[np],pd[npU],pd[npL],n1);//外から見て左上
		calcNormal(pd[np],pd[npR],pd[npU],n2);//右上
		a[i][j] = (float)((n1[0]+n2[0])/2.0f);
		b[i][j] = (float)((n1[1]+n2[1])/2.0f);
		c[i][j] = (float)((n1[2]+n2[2])/2.0f);
	}

	//表示
	//上底
	for(i = 0;i < Nxy;i++)
	{	ip = i+1;
		if(ip == Nxy) ip = 0;

		glBegin(GL_TRIANGLES);
			glNormal3d(a[i][0],b[i][0],c[i][0]);glVertex3fv(pd[i]);
			glNormal3d(a[i][1],b[i][1],c[i][1]);glVertex3fv(pd[Nxy+i]);
			glNormal3d(a[ip][1],b[ip][1],c[ip][1]);glVertex3fv(pd[Nxy+ip]);
		glEnd();
	}
	//下底
	glBegin(GL_POLYGON);
		glNormal3f(0.0f,0.0f,-1.0f);
		for(i=(Nz+1)*Nxy-1;i>=Nz*Nxy;i--) glVertex3fv(pd[i]);
	glEnd();
	//側面(4角形パッチ）
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
//１次元のテープ
void skSolidTape1(int Np, float *data)//単位幅、
{	//セルのサイズ幅１,全長1
	//軸方向の長さが一定となるようにセル自身の長さは可変）
	//基本姿勢の軸方向はx軸
	//Np---x軸方向セル数
	float pitch ;
	float p[501][2][3]; //頂点座標
	float a[501][2], b[501][2], c[501][2];//法線方向
	float nn[3];
	int k;

	if( Np > 500 ) Np = 500;
	pitch = 1.0f / (float)Np;

	//各頂点の座標
	for(k = 0; k <= Np; k++)
	{
		//y > 0側の座標
		p[k][0][0] = (float)k * pitch;
		p[k][0][1] = 0.5f;
		p[k][0][2] = data[k];
		//y < 0 側の座標
		p[k][1][0] = (float)k * pitch;
		p[k][1][1] = -0.5f;
		p[k][1][2] = data[k];
	}
	//各面の法線方向
	for(k = 0; k < Np; k++)
	{
		calcNormal(p[k][0],p[k][1],p[k+1][1],nn);
		a[k][0] = a[k][1] = nn[0];//x成分
		b[k][0] = b[k][1] = nn[1];//y成分
		c[k][0] = c[k][1] = nn[2];//z成分
	}
	a[Np][0] = a[Np][1] = a[Np-1][0];//終端は1つ前と同じとする
	b[Np][0] = b[Np][1] = b[Np-1][0];
	c[Np][0] = c[Np][1] = c[Np-1][0];

	//表面の描画
	for(k = 0 ;k < Np; k++)
	{
		glBegin(GL_POLYGON);
			glNormal3f(a[k][0],b[k][0],c[k][0]);glVertex3fv(p[k][0]);
			glNormal3f(a[k][1],b[k][1],c[k][1]);glVertex3fv(p[k][1]);
			glNormal3f(a[k+1][1],b[k+1][1],c[k+1][1]);glVertex3fv(p[k+1][1]);
			glNormal3f(a[k+1][0],b[k+1][0],c[k+1][0]);glVertex3fv(p[k+1][0]);
		glEnd();
	}
	//裏面の描画
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
	//セル(幅１,長さ1/Np,テープ全長1）
	//セル自身の長さは固定（軸方向の長さは変化する）
	//基本姿勢の軸方向はx軸
	//data[k]は正弦波の勾配データ（単位はrad）
	//Np---軸方向ピッチ数
	float pitch ;
	float p[501][2][3]; //頂点座標
	float nn[3];
	float a[501][2], b[501][2], c[501][2];//法線方向
	int k;

	if( Np > 500 ) Np = 500;
	pitch = 1.0f / (float)Np;

	//k = 0(原点）
	p[0][0][0] = p[0][1][0] = 0.0;      //x成分
	p[0][0][1] = 0.5; p[0][1][1] = -0.5;//y成分
	p[0][0][2] = p[0][1][2] = 0.0f;//z成分

	for(k = 1; k <= Np; k++)
	{
		p[k][0][0] = p[k][1][0] = p[k-1][0][0] + (float)cos(data[k-1]) * pitch;//x成分
		p[k][0][1] = 0.5; p[k][1][1] = -0.5;                                   //y成分
		p[k][0][2] = p[k][1][2] = p[k-1][0][2] + (float)sin(data[k-1]) * pitch;//z成分
	}

	//各面の法線方向
	for(k = 0; k < Np; k++)
	{
		calcNormal(p[k][0],p[k][1],p[k+1][1],nn);
		a[k][0] = a[k][1] = nn[0];//x成分
		b[k][0] = b[k][1] = nn[1];//y成分
		c[k][0] = c[k][1] = nn[2];//z成分
	}
	a[Np][0] = a[Np][1] = a[Np-1][0];//終端は1つ前と同じとする
	b[Np][0] = b[Np][1] = b[Np-1][0];
	c[Np][0] = c[Np][1] = c[Np-1][0];


	//表面の描画
	for(k = 0 ;k < Np; k++)
	{
		glBegin(GL_POLYGON);
			glNormal3f(a[k][0],b[k][0],c[k][0]);glVertex3fv(p[k][0]);
			glNormal3f(a[k][1],b[k][1],c[k][1]);glVertex3fv(p[k][1]);
			glNormal3f(a[k+1][1],b[k+1][1],c[k+1][1]);glVertex3fv(p[k+1][1]);
			glNormal3f(a[k+1][0],b[k+1][0],c[k+1][0]);glVertex3fv(p[k+1][0]);
		glEnd();
	}
	//裏面の描画
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
//正方形のメッシュ（x-y平面,中心が原点）
//ｘ軸方向，ｙ軸方向の幅は固定
void skSolidMesh(int Nx, float* data)
{
	//全体の幅,長さどちらも1(正方形)
	int i, j, Ny;
	float p[101][101][3]; //頂点座標
	float a[101][101], b[101][101], c[101][101];//
	float nn[3];
	float pitchX, pitchY;


	if(Nx > 100) Nx = 100;
	Ny = Nx;
	//セルのサイズ
	pitchX = 1.0f / (float)Nx;
	pitchY = 1.0f / (float)Ny;

	//各頂点の座標
	for(j = 0; j <= Ny; j++){
		for(i = 0; i <= Nx; i++){
			p[i][j][0] = (float)(i - Nx / 2) * pitchX;
			p[i][j][1] = (float)(j - Ny / 2) * pitchY;
			p[i][j][2] = data[j * (Nx+1) + i];
		}
	}
	//各面の法線方向
	for(i = 0; i < Nx; i++){
		for(j = 0; j < Ny; j++){
			calcNormal(p[i][j],p[i+1][j],p[i][j+1],nn);
			a[i][j] = nn[0];//x成分
			b[i][j] = nn[1];//y成分
			c[i][j] = nn[2];//z成分
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


	//表面の描画
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
	//裏面の描画
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

