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
	float diffuse[4]; //�g�U���W��
	float specular0;  //���ʔ��ˌ��W��
	float ambient0;   //�����W��
	float highlight;  //����
	int num1; //������(Nm:��~�������_��)
	int num2; //������(Ns:�f�ʉ~�������_��)
    int num3; //�o�l�̊�����
    double radius;//�o�l�̎�~�������a(�����j
    double length0;//�o�l�̎��R��
    double length;//�o�l�̒���(�����CvSize�����ׂ�1�j
    double constK;//�΂˒萔
	double ratio;//Torus,Spring(���S�����a�ɑ΂���f�ʔ��a�䗦
    //hinge
    int num;//hinge�ɗאڂ���particle�̃y�A��
    int row[6];//hinge�ɗאڂ���particle�̍s�ԍ�
    int col[6];//hinge�ɗאڂ���particle�̗�ԍ�
    int stk[6];//hinge�ɗאڂ���particle�̐ϔԍ�

	//���ފ֐�
	CSpring();       //�R���X�g���N�^
	void draw(bool flagShadow);
};
//-----------------------------------------------------------------
CSpring::CSpring()
{

	kind = S_SPRING;
//	kind = S_CYLINDER;
    vSize = CVector(1.0, 1.0, 1.0);//�K��1
    vPos = CVector();//���_
    vEuler = CVector(0.0, 0.0, 0.0);//��{�p���ł͂��������S��
    //�΂˂̐F�i���_��target.h)
    diffuse[0] = 0.2f;
    diffuse[1] = 0.99f;
    diffuse[2] = 0.99f;
    num1 = 6;
    num2 = 8;
    num3 = 10;//������
    radius = 0.06;//��~�������a�i����)
    ratio = 0.22;
    length0 = 0.5;
}
//-----------------------------------------------------------------
void CSpring::draw(bool flagShadow)
{
    float shadowDiffuse[] = {0.2f,0.2f,0.2f,0.3f};//�e�̊g�U��
    float shadowSpecular[] = {0.0f,0.0f,0.0f,1.0f};//�e�̋��ʌ�

	//�����Ɣ��ˌ���R,G,B�����W���Ŏw�肳��Ă���
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
	else{//�e
		glMaterialfv(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,shadowDiffuse);
		glMaterialfv(GL_FRONT,GL_SPECULAR,shadowSpecular);
	}

	glPushMatrix();

    if(kind == S_SPRING)
    {
        //���݈ʒu
	    glTranslated(vPos.x, vPos.y, vPos.z);//���s�ړ�

	    //��]
        glRotated(vEuler.z, 0.0, 0.0, 1.0);//z����]
        glRotated(vEuler.y, 0.0, 1.0, 0.0);//y����]
        glRotated(vEuler.x, 1.0, 0.0, 0.0);//x����]
        
        glRotated(90.0, 0.0, 1.0, 0.0);//��{�p�������������ɂȂ�悤��
        skSolidSpring(num1, num2, num3, radius, ratio, length) ;
    }
    else if(kind == S_CYLINDER){//�~��
        vSize.x = length;//��������x�������Ƃ����Ƃ̃T�C�Y
        vSize.y = 2.0 * radius;
        vSize.z = 2.0 * radius;
        //���݈ʒu
	    glTranslated(vPos.x, vPos.y, vPos.z);//���s�ړ�

	    //��]
        glRotated(vEuler.z, 0.0, 0.0, 1.0);//z����]
        glRotated(vEuler.y, 0.0, 1.0, 0.0);//y����]
        glRotated(vEuler.x, 1.0, 0.0, 0.0);//x����]

	    //�X�P�[�����O
	    glScaled(vSize.x, vSize.y, vSize.z);

        glRotated(90.0, 0.0, 1.0, 0.0);//��{�p���Ŏ�������x�������ɂȂ�悤��
        glTranslated(0.0, 0.0, 0.5);
	   	skSolidCylinder(num1);
    }
	glPopMatrix();
}

