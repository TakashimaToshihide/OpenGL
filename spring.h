//spring.h
enum SKind { S_SPRING, S_CYLINDER};
class CSpring
{
public:
//    String kind;
    SKind kind;
    CVector vSize;
    CVector vPos;
    CVector vEuler;
    int row1, col1, row2, col2, stk1, stk2;
    int n0, np1, np2, np[4];
	float diffuse[4]; //拡散光係数
	float specular0;  //鏡面反射光係数
	float ambient0;   //環境光係数
	float highlight;  //光沢
	int num1; //分割数(Nm:主円周分割点数)
	int num2; //分割数(Ns:断面円周分割点数)
    int num3; //バネの巻き数
    double radius;//バネの主円周軸半径(実寸）
    double length0;//バネの自然長
    double length;//バネの長さ(実寸，vSizeをすべて1）
    double constK;//ばね定数
	double ratio;//Torus,Spring(中心軸半径に対する断面半径比率
    //hinge
    int num;//hingeに隣接するparticleのペア個数
    int row[6];//hingeに隣接するparticleの行番号
    int col[6];//hingeに隣接するparticleの列番号
    int stk[6];//hingeに隣接するparticleの積番号

	//ﾒﾝﾊﾞ関数
	CSpring();       //コンストラクタ
	void draw(bool flagShadow);
};
//-----------------------------------------------------------------
CSpring::CSpring()
{

	kind = S_SPRING;
//	kind = S_CYLINDER;
    vSize = CVector(1.0, 1.0, 1.0);//必ず1
    vPos = CVector();//原点
    vEuler = CVector(0.0, 0.0, 0.0);//基本姿勢ではｚ軸が中心軸
    //ばねの色（質点はtarget.h)
    diffuse[0] = 0.2f;
    diffuse[1] = 0.99f;
    diffuse[2] = 0.99f;
    num1 = 6;
    num2 = 8;
    num3 = 10;//巻き数
    radius = 0.06;//主円周軸半径（実寸)
    ratio = 0.22;
    length0 = 0.5;
}
//-----------------------------------------------------------------
void CSpring::draw(bool flagShadow)
{
    float shadowDiffuse[] = {0.2f,0.2f,0.2f,0.3f};//影の拡散光
    float shadowSpecular[] = {0.0f,0.0f,0.0f,1.0f};//影の鏡面光

	//環境光と反射光はR,G,B同じ係数で指定されている
    float ambient[4];
    float specular[4];
	ambient[0] = ambient0; ambient[1] = ambient0; ambient[2] = ambient0; ambient[3] = 1.0;
	specular[0] = specular0; specular[1] = specular0; specular[2] = specular0; specular[3]=1.0;

	if(flagShadow == false){
		glMaterialfv(GL_FRONT,GL_DIFFUSE,diffuse);
		glMaterialfv(GL_FRONT,GL_AMBIENT,ambient);
		glMaterialfv(GL_FRONT,GL_SPECULAR,specular);
		glMaterialf(GL_FRONT,GL_SHININESS,highlight);
	}
	else{//影
		glMaterialfv(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,shadowDiffuse);
		glMaterialfv(GL_FRONT,GL_SPECULAR,shadowSpecular);
	}

	glPushMatrix();

    if(kind == S_SPRING)
    {
        //現在位置
	    glTranslated(vPos.x, vPos.y, vPos.z);//平行移動

	    //回転
        glRotated(vEuler.z, 0.0, 0.0, 1.0);//z軸回転
        glRotated(vEuler.y, 0.0, 1.0, 0.0);//y軸回転
        glRotated(vEuler.x, 1.0, 0.0, 0.0);//x軸回転
        
        glRotated(90.0, 0.0, 1.0, 0.0);//基本姿勢がｘ軸方向になるように
        skSolidSpring(num1, num2, num3, radius, ratio, length) ;
    }
    else if(kind == S_CYLINDER){//円柱
        vSize.x = length;//軸方向がx軸方向とたあとのサイズ
        vSize.y = 2.0 * radius;
        vSize.z = 2.0 * radius;
        //現在位置
	    glTranslated(vPos.x, vPos.y, vPos.z);//平行移動

	    //回転
        glRotated(vEuler.z, 0.0, 0.0, 1.0);//z軸回転
        glRotated(vEuler.y, 0.0, 1.0, 0.0);//y軸回転
        glRotated(vEuler.x, 1.0, 0.0, 0.0);//x軸回転

	    //スケーリング
	    glScaled(vSize.x, vSize.y, vSize.z);

        glRotated(90.0, 0.0, 1.0, 0.0);//基本姿勢で軸方向がx軸方向になるように
        glTranslated(0.0, 0.0, 0.5);
	   	skSolidCylinder(num1);
    }
	glPopMatrix();
}

