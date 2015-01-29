//target.h

enum Kind {SPHERE, CUBE, CYLINDER, PRISM, PYRAMID, CONE, FRUSTUM, TORUS,
           SEMI_SPHERE, SUPER1, SUPER2, TAPE1, TAPE2, MESH, SQUARE,
           TEX_MESH1, TEX_MESH2, TEX_MESH3, TEX_MESH4};
enum TargetDir { TD_DOWN,TD_UP,TD_FORE,TD_BACK,TD_RIGHT,TD_LEFT};
//�`��I�u�W�F�N�g
class CTarget{ 
public:
	//���ޕϐ�
	Kind kind;     //CTarget�̎��
	CVector vSize; //�T�C�Y
	CVector vPos;  //�ʒu
    CVector vEuler;//�I�C���[�p�i�p���\��)
    CVector vP[40];//������,�~���̒��_���W
    CVector vNormalFacet[6];//�����̖̂ʂ̖@���x�N�g��
    double dir;    //�i�s����
    //�����l��ۑ����Ă��������Ƃ��Ɏg�p
    CVector vPos0;
    CVector vEuler0;
    double dir0;
    //�^���V�~�����[�V�����p
    double mass;        //����
    CVector vForce;    //�O��
    CVector vForce0;   //�O�͏����l�i���ݒl�Ƃ��Ă����p)
    CVector vAccel;    //�����x
    CVector vVelocity; //���x
    CVector vVelocity0;//�����x
    CVector vOmega;    //�p���x
    CVector vOmega0;   //���p���x
    CVector vAlpha;    //�p�����x
    double omega0;      //�C�ӎ����̊p���x
    CVector vAxis;     //���̉�]��
    CVector vGravityToPoint;//�d�S����Փ˓_�܂ł̃x�N�g��
    CMatrix mInertia; //�������[�����g�̃e���\��
    CMatrix mInertiaInverse;//���̋t�s��
    CQuaternion q;   //�l����
    double speed;     //���x�̐�Βl
    double speed0;
	bool flagFixed;
    //�}�e���A������
	float diffuse[4]; //�g�U���W��
	float specular0;  //���ʔ��ˌ��W��
	float ambient0;   //�����W��
	float highlight;  //����
	short num1; //������(Nxy,Nm)
	short num2; //������(Nz,Ns)
	float ratio;//Frustum(����ɑ΂����ꔼ�a�䗦)
				 //Torus(���S�����a�ɑ΂���f�ʔ��a�䗦
	float eps1, eps2, p1, p2, p3;//���Q���֐�
    float amp, lambdaX, lambdaY, lambdaZ, alpha, weight, wind;//TEX_MESH2,TEX_MESH3
    float data[10201]; //�g�̃f�[�^(100*100)
    int texMode;//DECAL,MODULATE
    int texType;//PARALLEL1,PARALLEL2,CYLINDRICAL,SPHERICAL,SOLID
    byte texImage[T_MAX][T_MAX][4];
    byte texImageS[6][T_MAX][T_MAX][3]; //SolidTexture�p
    short numStress;//particle�ɍ�p���鉞�͂̌�
    float stress;//���͂̕��ϒl
    short posType;//�p(0)�C��(1)�C��(2)�C����(3)

	//���ފ֐�
	CTarget();       //�R���X�g���N�^
	void draw(bool flagShadowDrawing);
    void setTexture();
	void makeTexture(System::Drawing::Bitmap* bitmap);
    void calcSolidTexCube(CVector vSize) ;
    void calcSolidTexSphere(CVector vSize);
    void calcSolidTexCylinder(CVector vSize);
    void parabola(double beta, double omega0, double e, double mu, double dt);
    void projection(double e, double mu, double cv, double ci, double dt);
    void getVertexOfCube();
    void getVertexOfCylinder();
    void calcInertia();
    bool collisionCubeToCube2(CTarget& trgt, CVector& vNormal) ;
    bool collisionSphereToCube2(CTarget& trgt,CVector& vNormal);
    bool collisionCubeToSphere2(CTarget& trgt, CVector& vNormal) ;
    bool collisionCubeToCube3(CTarget& trgt, CVector& vNormal) ;
    bool collisionSphereToCube3(CTarget& trgt,CVector& vNormal);
    bool collisionCubeToSphere3(CTarget& trgt, CVector& vNormal) ;
    bool collisionCubeToCylinder3(CTarget& trgt, CVector& vNormal) ;
    bool collisionCylinderToCube3(CTarget& trgt, CVector& vNormal) ;
    bool collisionSphereToCylinder3(CTarget& trgt,CVector& vNormal);
    bool collisionCylinderToSphere3(CTarget& trgt,CVector& vNormal);
    bool collisionCylinderToCylinder3(CTarget& trgt,CVector& vNormal);
    bool collisionSphereToSphere3(CTarget& trgt,CVector& vNormal);
    double getAngularAccelOnFloor(double mu);
    double getTorqueOnFloor(double mu);
    void rolling(double mu, double dt);
    void sliding(double mu, double dt);
    void getDirVector(TargetDir targetDir, CVector& p);
    void calcRotation(double* p);
    double checkCollisionToFloor();
private:
    void getImageGrain(double x, double y, double z, CVector vSize, byte* rgb);
    double getMaxH(void);
    void calcFrictionEffect(double e, double mu, double dt) ;
    void calcFrictionEffectSub(double &vh, double &omega, double e, double mu, double dt) ;
    bool getPointCubeToCube(CTarget& trgt, CVector& vNormal, CVector& vPoint);
    bool getPointCylinderToCube(CTarget& trgt, int& faceNo, CVector& vPoint);
};
//------------------------------------------------------
//�R���X�g���N�^
CTarget::CTarget()
{
	kind = CUBE;
    vSize = CVector(0.1f, 0.1f, 0.1f);
    vPos = CVector();//���_
    vPos.z = vSize.z / 2.0f;//�ʒu�͵�޼ު�Ă̒��S
//    flagFixed = false;
    //���̐F�i�g�U����ԐF�ɐݒ�)
    diffuse[0] = 1.0;//red
    diffuse[1] = 0.0;//green
    diffuse[2] = 0.0;//blue
    diffuse[3] = 1.0;//alpha�l�i�s�����ɐݒ�)
	ambient0 = 0.4f;
	specular0 = 1.0f;
    highlight = 80.0f;
    num1 = 10;
    num2 = 10;
    ratio = 0.5;//Frustum,Torus
    eps1 = 1.0;
    eps2 = 1.0;
    p1 = p2 = p3 = 1.0;//��2���֐��̌`�����Ұ�
}
//---------------------------------------------------------
void CTarget::draw(bool flagShadow)
{
    float shadowDiffuse[] = {0.2f,0.25f,0.25f,0.3f};//�e�̊g�U��
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
    //���݈ʒu
	glTranslated(vPos.x, vPos.y, vPos.z);//���s�ړ�

	//��]
    glRotated(vEuler.z, 0.0f, 0.0f, 1.0f);//z����]
    glRotated(vEuler.y, 0.0f, 1.0f, 0.0f);//y����]
    glRotated(vEuler.x, 1.0f, 0.0f, 0.0f);//x����]
	//�X�P�[�����O
	glScaled(vSize.x, vSize.y, vSize.z);

	//��޼ު�Ă̕`��
	//SolidModel
    if(kind == CUBE && texType == T_NON)
        skSolidCube();
    else if(kind == SPHERE && texType == T_NON)
        skSolidSphere(num1);
    else if(kind == PRISM)
	   	skSolidPrism(num1);
    else if(kind == CYLINDER && texType == T_NON)//�~��
	   	skSolidCylinder(num1);
    else if(kind == PYRAMID)//���p�`
	  	skSolidPyramid(num1);
    else if(kind == CONE)//�~��
	 	skSolidCone(num1);
    else if(kind == FRUSTUM)//�~����
		skSolidFrustum(num1,ratio);
    else if(kind == TORUS)//�g�[���X
   		skSolidTorus(num1,num2,ratio);
    else if(kind == SEMI_SPHERE)
		skSolidSemiSphere(num1);
    else if(kind == SUPER1 && texType == T_NON)//���Q���֐��i�㉺�Ώ́j
	   	skSolidSuper1(num1,num2,eps1,eps2,p1,p2,p3);
    else if(kind == SUPER2)//���Q���֐��i�㕔�����j
	  	skSolidSuper2(num1,num2,eps1,eps2,p1,p2);
    //�g�̃v���~�e�B�u�Ƃ���3��ޒǉ�
    else if(kind == TAPE1)
        skSolidTape1(num1, data);
    else if(kind == TAPE2)
        skSolidTape2(num1, data);
    else if(kind == MESH)
        skSolidMesh(num1, data);
    //Texture
    else if(kind == SQUARE && texType != T_NON)//���ʓ��e
        skTexSquare();
    else if(kind == CUBE && texType == T_PLANAR1)//���ʓ��e
        skTexCube1();
    else if(kind == CUBE && texType == T_PLANAR2)//���ʓ��e
        skTexCube2();
    else if(kind == SPHERE && texType == T_PLANAR1)//���ʓ��e
        skTexSphere1(num1);
    else if(kind == SPHERE && texType == T_SPHERICAL)//���ʓ��e
        skTexSphere2(num1);
    else if(kind == CYLINDER && texType == T_PLANAR1)//���s���e
        skTexCylinder1(num1);
    else if(kind == CYLINDER && texType == T_CYLINDRICAL)//�~�����e
        skTexCylinder2(num1);
    else if(kind == SUPER1 && texType != T_NON)
        skTexSuper1(num1,num2,eps1,eps2,p1,p2,p3);
    else if(kind == TEX_MESH1)
        skTexMesh1(num1, data);
    else if(kind == TEX_MESH2)
        skTexMesh2(num1, ratio, amp, lambdaX, lambdaY, alpha) ;
    else if(kind == TEX_MESH3)
        skTexMesh3(num1, ratio, amp, lambdaX, alpha, wind, vPos.z);
    //Solid Texture
    else if(kind == CUBE && texType == T_SOLID)
        skSolidTexCube(texImageS);
    else if(kind == SPHERE && texType == T_SOLID)
        skSolidTexSphere(num1, texImageS);
    else if(kind == CYLINDER && texType == T_SOLID)
        skSolidTexCylinder(num1, texImageS);

	glPopMatrix();
}
//-----------------------------------------------------------------------------
void CTarget::setTexture()
{
	glTexImage2D(GL_TEXTURE_2D,0,4,T_MAX,T_MAX,0,GL_RGBA,GL_UNSIGNED_BYTE,texImage);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
	if(texMode == T_DECAL)
		glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_DECAL);
	else//(texMode == T_MODULATE)
		glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE);
}
//---------------------------------------------------------------------------
void CTarget::makeTexture(System::Drawing::Bitmap* bitmap)
{//bitmap����e�N�X�`���摜�̍쐬
    int i, j, ii, jj;
    int nx = 128;
    int ny = 128;
    double pWidth = (double)bitmap->Width;
    double pHeight = (double)bitmap->Height;
    if( pHeight == 0.0) {
        MessageBox(NULL,"makeTexture","�摜������܂���",MB_OK);
        return;
    }
//    ratio = pWidth / pHeight;

    for(i = 0; i < nx; i++)
    {
    	for(j = 0; j < ny; j++)
        {
            ii = i * (int)(pWidth / (double)nx);
            jj = j * (int)(pHeight / (double)ny);
            texImage[ny-j-1][i][0] = (bitmap->GetPixel(ii,jj)).R;
            texImage[ny-j-1][i][1] = (bitmap->GetPixel(ii,jj)).G;
            texImage[ny-j-1][i][2] = (bitmap->GetPixel(ii,jj)).B;
            texImage[ny-j-1][i][3] = (byte)(255 * diffuse[3]);
        }
    }
}
//------------------------------------------------------------------------------
double CTarget::getTorqueOnFloor(double mu)
{
    double a, b, a2, b2, a3, b3, c, ss , r, h;
    CVector vFore, vUp, vLeft;
	double g = 9.8f;//m/s^2
    double torque;//���ʖ��C�ɂ���]��

    //�������𒲍�
    getDirVector(TD_FORE, vFore);
    getDirVector(TD_UP, vUp);
    getDirVector(TD_LEFT, vLeft);

    if(kind == SPHERE) torque = 0.01f;//�����Ȓl

    else if(kind == CYLINDER){
        r = vSize.x / 2.0f; h = vSize.z;
        if(fabs(vUp.z) >= 0.99)//fabs(vFore.z) && fabs(vUp.z) >= fabs(vLeft.z) )
        {//��{�p���ŉ�]���Ă���Ƃ�
            torque = 2.0f * mu * mass * g * r / 3.0f;
        }
        else if(fabs(vUp.z) < 0.02){//���ɂȂ��ĉ�]���Ă���Ƃ�
            torque = mu * mass * g * h / 4.0f;
        }
        else torque = 0.01f;//�������l�Ƃ���

    }
    else if(kind == CUBE){
        if(fabs(vUp.z) > 0.99)
        {
            a = vSize.x; b = vSize.y;
        }
        else if(fabs(vFore.z) >0.99)
        {
            a = vSize.y; b = vSize.z;
        }
        else if(fabs(vLeft.z) > 0.99)
        {
            a = vSize.x; b = vSize.z;
        }
        else
        {
            torque = 0.01f;
            if(vOmega.z > 0.0) torque = -0.01f;
            return(torque);

        }
        a2 = (a/2.0f)*(a/2.0f);
        a3 = (a/2.0f)*(a/2.0f)*(a/2.0f);
        b2 = (b/2.0f)*(b/2.0f);
        b3 = (b/2.0f)*(b/2.0f)*(b/2.0f);
        c =  sqrt(a2+b2);
        ss = (a*b*c - b3*log(-a/2.0f + c)
          + b3*log(a/2.0f + c) - a3*log(-b/2.0f + c) + a3*log(b/2.0f + c) ) / (3.0f*a*b);
        torque = mu * mass * g * ss;
    }
    if(vOmega.z > 0.0) torque = - torque;//�����ɂȂ�悤��
    return( torque ); //
}

//------------------------------------------------------------------------------
//Floor�����]����~��,�p���̖��C�ɂ��p�����x���v�Z
double CTarget::getAngularAccelOnFloor(double mu)
{
    double a, b, a2, b2, a3, b3, c, ss , r, h;
    CVector vFore, vUp, vLeft;
	double g = 9.8f;//m/s^2
    double alpha;//�p�����x
    double inertia;

    //�������𒲍�
    getDirVector(TD_FORE, vFore);
    getDirVector(TD_UP, vUp);
    getDirVector(TD_LEFT, vLeft);
//if(kind == SPHERE) alpha = 5.0;
    if(kind == CYLINDER){
        r = vSize.x / 2.0f; h = vSize.z;
        if(fabs(vUp.z) >= fabs(vFore.z) && fabs(vUp.z) >= fabs(vLeft.z) )
        {//��{�p���ŉ�]���Ă���Ƃ�
            alpha = 4.0f * mu * g / (3.0f * r);
        }
        else{//���ɂȂ��ĉ�]���Ă���Ƃ�
            alpha = mu*g*h / (r*r + h*h/3.0f);
        }
    }
    else if(kind == CUBE){

        if(fabs(vFore.z) >= fabs(vLeft.z) && fabs(vFore.z) >= fabs(vUp.z) )
        {
            a = vSize.y; b = vSize.z;
        }
        else if(fabs(vLeft.z) >= fabs(vFore.z) && fabs(vLeft.z) >= fabs(vUp.z) )
        {
            a = vSize.x; b = vSize.z;
        }
        else if(fabs(vLeft.z) >= fabs(vFore.z) && fabs(vLeft.z) >= fabs(vUp.z) )
        {
            a = vSize.x; b = vSize.z;
        }
        else if(fabs(vUp.z) >= fabs(vFore.z) && fabs(vUp.z) >= fabs(vLeft.z) )
        {
            a = vSize.x; b = vSize.y;
        }
        else{
            MessageBox(NULL,"�v�Z�װ" ,"calcAngularAccelOnFloor",MB_OK); return(0);
        }
        a2 = (a/2.0f)*(a/2.0f);
        a3 = (a/2.0f)*(a/2.0f)*(a/2.0f);
        b2 = (b/2.0f)*(b/2.0f);
        b3 = (b/2.0f)*(b/2.0f)*(b/2.0f);
        c =  sqrt(a2+b2);
        ss = (a*b*c - b3*log(-a/2.0f + c)
          + b3*log(a/2.0f + c) - a3*log(-b/2.0f + c) + a3*log(b/2.0f + c) ) / (3.0f*a*b);
        inertia = (a * a + b * b) / 12.0f;
        alpha = mu * g * ss / inertia ;
    }
    return( alpha );
}  
//-------------------------------------------------------------------
//dt��̈ʒu�C��]�p�C���x�����߂�
void CTarget::rolling(double mu, double dt)
{                   //�]���薀�C�W��,���ԍ���
	double dist;
	double g = 9.8f;
	double pp = M_PI / 180.0;
    double ang;
    CQuaternion q;
    //��]���͐����Ői�s�����ɒ��s���鎲
    vAxis.x = cos((dir + 90.0) * pp);
    vAxis.y = sin((dir + 90.0) * pp);
    vAxis.z = 0.0;

    speed -= mu * g * dt; //speed�̍X�V
    dist = speed * dt;    //�ړ�����
    ang = (2.0f * dist / (vSize.z * pp)); //��]�p�i�]����^���ɂ���]�p)
    vPos.x += dist * cos(dir * pp);
    vPos.y += dist * sin(dir * pp);
    q = makeQFromAxis( ang , vAxis); //Quaternion���쐬
    vEuler = makeEulerFromEuler(vEuler, q );//�I�C���[�p�̍X�V
}
//-------------------------------------------------------------------
//dt��̈ʒu�C��]�p�C���x�����߂�
//�I�u�W�F�N�g����{�p���̂Ƃ�����������
void CTarget::sliding(double mu, double dt)
{                   //�����C�W��,���ԍ���
	double a, aa, dist;
	double g = 9.8f;
	double pp = M_PI / 180.0;
    //��]��
    vAxis = CVector(0.0, 0.0, 1.0);

    //�����x
    a = - mu * g;
    //�p�����x(���C�ɂ��p�����x�̌����j
    aa = getAngularAccelOnFloor(mu);//rad/s^2�Ŏ擾
    if(vOmega.z > 0.0) aa = -aa;

    //�V���x
    speed += a * dt;
    //�V�p���x
    vOmega.z += aa * dt;

    //�ړ�����
    dist = speed * dt;
    //�ʒu�̍X�V
    vPos.x += dist * cos(dir * pp);
    vPos.y += dist * sin(dir * pp);
    //��]�p�x�̍X�V
    vEuler.z += vOmega.z * dt / pp; //�I�C���[�p��deg�P��
}
//-----------------------------------------------------------------------------
//�p���璆�S�܂ł̍����̍ő�l�i���ƒ����̂����j
double CTarget::getMaxH(void)
{
	short i;
    CVector xyz;
	double h[4],maxH;

	if(kind == SPHERE)
        return(vSize.z / 2.0f);//���̔��a��Ԃ�

	for(i = 0;i < 4;i++){//������
		if(i == 0){
            xyz = CVector(vSize.x / 2.0f, vSize.y / 2.0f, vSize.z / 2.0f);
        }
		else if(i == 1){
            xyz = CVector(-vSize.x / 2.0f, vSize.y / 2.0f, vSize.z / 2.0f);
        }
		else if(i == 2){
            xyz = CVector(-vSize.x / 2.0f, -vSize.y / 2.0f, vSize.z / 2.0f);
        }
		else{
            xyz = CVector(vSize.x / 2.0f, -vSize.y / 2.0f, vSize.z / 2.0f);
        }

		//x����]
        xyz.rotX(vEuler.x);
		//y����]
        xyz.rotY(vEuler.y);
		h[i] = fabs(xyz.z);//���͋��߂�K�v�Ȃ�
	}

	//�ő�l�����߂�
	maxH = h[0];
	for(i=1;i<4;i++){
		if(maxH < h[i]) maxH = h[i];
	}
	return(maxH);
}
//------------------------------------------------------------------------------

void CTarget::projection(double e, double mu, double cv, double ci, double dt)
                //�����W��, �����C��R,�S����R�W��,������R�W��,���ݎ���
{
    double minH;
    double rikiseki ;
    CVector vNormal;//���ڕ��̂��猩���Փ˖ʂ̖@������
    CVector vp;//�Փ˖ʂ̍������x�i���`���x�{�p���x�ɂ��ڐ����x)
    CVector gtp;//�d�S����Փ˓_�܂ł̃x�N�g��
    CVector vTorque = CVector();
    double ang;
    CQuaternion q;
	double pp = M_PI / 180.0;

    //�͂̏����l
    vForce = vForce0;//vForce0��main���Őݒ�

    //���݂̑��x
    speed = vVelocity.magnitude();
    //��C�̒�R
    vForce -= cv * vVelocity;//cv�͑��x�ɔ�Ⴗ��S����R�W��
    vForce -= ci * speed * vVelocity;//ci�͑��x��2��ɔ�Ⴗ�銵����R�W��

    //�����x
    vAccel = vForce / mass;
    //���x
    vVelocity += vAccel * dt;//���ʂ����炵�Ă�-z�����̑��x��������
    //�ʒu
    vPos += vVelocity * dt ;
    //��]�p�x���X�V����]��̃I�C���[�p�����߂�
    ang = omega0 * dt / pp; //deg
    if(fabs(vEuler.y) == 90.0 && ang < 0.1) ang = 0.1f;//�}90�x�ɂ�����lock�h�~
    q = makeQFromAxis( ang , vAxis); //Quaternion���쐬
    vEuler = makeEulerFromEuler(vEuler, q);//��]��̃I�C���[�p


    //���݂̊p���x
    vOmega = omega0 * vAxis; //rad
    //���ʂƂ̏Փ˔���i�Œ�_��Ԃ�)
    minH = checkCollisionToFloor();//�Փ˓_�Ȃǂ�����ٰ�݂ŋ��߂Ă���
    if(minH < 0.0) vPos.z += fabs(minH) ;//�߂荞�݁C���ݍ��ݖh�~
    if( minH <= 0.0 && vVelocity.z <= 0.0)  //Floor�ɏՓ�
    {
        if(kind == SPHERE ){
            //�V�d�S���x�̂������i�ēx���̂��̑��x�ŕ����^���j
            vVelocity.z = - e * vVelocity.z;
        }
        else if(kind == CUBE || kind == CYLINDER){
            //�܂����C�Ȃ��ōl����
            vNormal = CVector(0.0, 0.0, -1.0);//�Փ˖ʂ̖@������(target������̖@���j
            //�d�S����Փ˓_�܂ł̃x�N�g���icheckCollisionToFloor()���[�`���ŋ��߂Ă���)
            gtp = vGravityToPoint;
            //�Փ˓_�̍������x
            vp = vVelocity + (vOmega ^ gtp); //�Փ˓_���x(�x�N�g���j

            //�͐�
            rikiseki = - (e + 1.0f) * (vNormal * vp) / (1.0f/mass + vNormal*((mInertia.inverse()*(gtp^vNormal))^gtp));
            //�V�p���x
            vOmega += mInertiaInverse * (gtp ^ vNormal) * rikiseki ;
            //�V�d�S���x
            vVelocity += (vNormal * rikiseki) / mass ;
        }
        else{
            MessageBox(NULL,"���̌`�󕨑̂ɂ͖��Ή�","projection",MB_OK);
        }
        //���ʖ��C������ꍇ
        calcFrictionEffect(e, mu, dt);
        //���ʖ��C�ɂ�錸���g���N
        vTorque.z = getTorqueOnFloor(mu);
	}
    //��]���C�g���N�ɂ��p�����x
    vAlpha = mInertiaInverse * vTorque;
//    vAlpha = mInertiaInverse *( vTorque - (vOmega ^ (mInertia * vOmega)));
    //�p���x�̍X�V
    vOmega += vAlpha * dt;
    //�V�p���x�̑傫���Ɖ�]��
    omega0 = vOmega.magnitude();
    vAxis = vOmega.normalize2();
}

//----------------------------------------------------------------
//�Փ˂̂Ƃ����C������Ή�]���x���ω�����
void CTarget::calcFrictionEffect(double e, double mu, double dt)
{
    double omega;

    calcFrictionEffectSub(vVelocity.x, vOmega.y, e, mu, dt);
    omega = -vOmega.x;
    calcFrictionEffectSub(vVelocity.y, omega, e, mu, dt);
    vOmega.x = -omega;
}

//----------------------------------------------------------------
void CTarget::calcFrictionEffectSub(double &vh, double &omega, double e, double mu, double dt)
{
    double vhc, muc, rr;
    double g = 9.8f;

//    if(fabs(vVelocity.z) < 0.001 && fabs(vh) < 0.001) return;
    if(kind == SPHERE )
    {
        rr = vSize.x / 2.0f;

        //���݂�omega�Ŋ��炸�ɓ]����Ƃ��̐������x����
        vhc = rr * omega ;
        //�ՊE���C�W��
        if(fabs(vVelocity.z) <= 0.01) muc = 1000.0;
        else muc = ((2.0/7.0) * fabs((vh-vhc) / (1.0+e) / vVelocity.z));
//wsprintf(buf , "%d, %d, %d",(int)(vh*1000), (int)(muc * 1000),(int)(vVelocity.z*1000));
//MessageBox(NULL,buf ,"aaa",MB_OK);
        if(mu >= muc) {
            vh = (5.0f * vh + 2.0f * vhc) / 7.0f;
            omega = vh / rr;
//wsprintf(buf , "%d, %d, %d",(int)(vh*1000.0),(int)(vhc*1000.0),(int)(omega*1000.0));
//MessageBox(NULL,buf ,"aaa",MB_OK);
        }
        else{//vVelocity.z =0�Ȃ�Ώ�ɂ��̏��
            if(vh > vhc){
                vh -= mu * (1.0f+e) * fabs(vVelocity.z);
                omega += 2.5f * mu * (1.0f+e) * fabs(vVelocity.z) / rr;
            }
            else if(vh < vhc){
                vh += mu * (1.0f+e) * fabs(vVelocity.z);
                omega -= 2.5f * mu * (1.0f+e) * fabs(vVelocity.z) / rr;
            }
        }
    }
    else if(kind == CUBE )
    {
        //�����̂𔼌arr�̋��ŋߎ�
        rr = sqrt(vSize.x*vSize.x + vSize.y*vSize.y + vSize.z*vSize.z) / 2.0f;

        //���݂�omega�Ŋ��炸�ɓ]����Ƃ��̐������x����
        vhc = rr * omega; //�����̂𔼌arr�̋��ŋߎ�
        //�ՊE���C�W��
//char buf[30];
//if(vVelocity.z ==0.0){
//wsprintf(buf , "%d, %d, %d",(int)(vh*1000.0),(int)(omega*1000.0),(int)(vVelocity.z*1000.0));
//MessageBox(NULL,buf ,"aaa",MB_OK);  }
        if(fabs(vVelocity.z) <= 0.005) {
            //�]���炸�Ɋ���
            if(vh > 0.0) vh -= mu * g * dt;
            else         vh += mu * g * dt;
        }
        else
        {
            if(fabs(vVelocity.z) <= 0.01) muc = 1000.0;
            else muc = (0.25 * fabs((vh-vhc) / (1.0+e) / vVelocity.z));
            if(mu >= muc) {
                vh = (3.0f*vh + vhc)/4.0f;
                omega = vh / rr;
            }
            else{ //vVelocity.z =0�Ȃ�Ώ�ɂ��̏��
                if(vh > vhc){
                    vh -= mu * ((1.0+e) * fabs(vVelocity.z));
                    omega += 3.0f * mu * ((1.0f+e) * fabs(vVelocity.z)) / rr;
                }
                else if(vh < vhc){
                    vh += mu * ((1.0+e) * fabs(vVelocity.z));
                    omega -= 3.0f * mu * ((1.0+e) * fabs(vVelocity.z)) / rr;
                }
            }
        }
    }
    else//�~��
    {
        rr = vSize.x / 2.0f;

        //���݂�omega�Ŋ��炸�ɓ]����Ƃ��̐������x����
        vhc = rr * omega ;
        //�ՊE���C�W��
//        if(fabs(vVelocity.z) <= 0.01) muc = 1000.0;
        if(fabs(vVelocity.z) <= 0.005) {//muc = 1000.0;
            //�]���炸�Ɋ���
            if(vh > 0.0) vh -= mu * g * dt;
            else         vh += mu * g * dt;
        }
        else {
            muc = ((1.0/3.0) * fabs((vh-vhc) / ((1.0+e) * vVelocity.z)));
            if(mu >= muc) {
                vh = (2.0f * vh + vhc) / 3.0f;
                omega = vh / rr;
            }
            else{
                if(vh > vhc){
                    vh -= mu * ((1.0+e) * fabs(vVelocity.z));
                    omega += 2.0f * mu * ((1.0+e) * fabs(vVelocity.z)) / rr;
                }
                else if(vh < vhc){
                    vh += mu * ((1.0+e) * fabs(vVelocity.z));
                    omega -= 2.0f * (mu * (1.0+e) * fabs(vVelocity.z)) / rr;
                }
            }
        }
    }
}
//-------------------------------------------------------------------------
//Floor�Ƃ̏Փ˂𒲂�Floor����Œ�ʒu�̍�����Ԃ�
//���̒l�����ł���ΏՓ�
double CTarget::checkCollisionToFloor()
{
    double minH;
    int i, cnt, numVertex;// maxNum;
//    int faceNo, numCollision[6], collisionP[8][4];//�Փ˖ʔԍ��C�Փ˓_���C�Փ˓_�ԍ�
    double eps = 0.001f;//1mm�ȉ��Ȃ�Փ˂Ƃ���
    double ss;//���ł���Δ��a�C�����̂ł���ΑΊp���̔���
    CVector vPointCollision;
    //�ʂ̒��_�ԍ�(CUBE)
/*    int nFace[6][4] = { {0, 1, 2, 3}, {0, 3, 7, 4}, {0, 4, 5, 1},
                       {1, 5, 6, 2}, {2, 6, 7, 3}, {4, 7, 6, 5} }; */
    CVector vCollision[20];//�Փ˓_�̒��_�ԍ�(�~��)
//char buf[40];
    //���̂Ƃ݂Ȃ����Ƃ��̔��a
    if(kind == CUBE ) {//������
        ss = sqrt(vSize.x * vSize.x + vSize.y * vSize.y + vSize.z * vSize.z) / 2.0f;
    }
    else if(kind == SPHERE){ //��
        ss = vSize.x / 2.0f;
    }
    else{//�~��
        ss = sqrt(vSize.x*vSize.x + vSize.z*vSize.z);
    }

    minH = 1000.0;
    if( vPos.z <= ss)
    {//�Փ˂̉\������
        if(kind == SPHERE)
        {
            minH = vPos.z - ss;
            vGravityToPoint.x = 0.0;
            vGravityToPoint.y = 0.0;
            vGravityToPoint.z = -vSize.z / 2.0f;

            return minH;
        }
        else
        {
            if(kind == CUBE)
            {
                getVertexOfCube();
                numVertex = 8;
            }
            else//�~��
            {
                getVertexOfCylinder();
                numVertex = 2 * num1;
            }
            cnt = 0;//�Փ˓_��
            minH = 1000.0;
            vPointCollision = CVector(0.0, 0.0, 0.0);
            for(i = 0 ; i < numVertex; i++) {
                if(vP[i].z <= eps) { vPointCollision += vP[i]; cnt++; }
                if(vP[i].z < minH )  minH = vP[i].z;
            }

            if(minH <= 0.0)  //Floor�ƏՓ�
            {
                vPointCollision /= (double)cnt;
                vGravityToPoint = vPointCollision - vPos;
                return minH;
            }
        }
    }
    return minH;
//    return 1000.0;//�Փ˂Ȃ�
}
//------------------------------------------------------------------------------
void CTarget::calcInertia()
{
    double Ixx, Iyy, Izz;
    if(kind == SPHERE)
    {
        Ixx = Iyy = Izz = vSize.x * vSize.y * mass / 10.0f; //vSize�͒��a
    }
    else if(kind == CYLINDER)//
    {
        Ixx = vSize.x * vSize.y * mass / 16.0f + vSize.z * vSize.z * mass / 12.0f;
        Iyy = Ixx;
        Izz = vSize.z * vSize.z * mass / 8.0f;
    }
    else if(kind == CUBE)//������(�p��)
    {
        Ixx = mass * (vSize.y*vSize.y + vSize.z*vSize.z) / 12.0f;
        Iyy = mass * (vSize.x*vSize.x + vSize.z*vSize.z) / 12.0f;
        Izz = mass * (vSize.x*vSize.x + vSize.y*vSize.y) / 12.0f;
    }
    //�������[�����g �̃e���\��
    mInertia = CMatrix(Ixx, 0.0, 0.0,
                       0.0, Iyy, 0.0,
                       0.0, 0.0, Izz );
    mInertiaInverse = mInertia.inverse();

}


//-------------------------------------------------------------------
//�c������Ƃ��ȂǁC��̕��������߂�Ƃ��ɗp����
void CTarget::getDirVector(TargetDir targetDir, CVector& p)
{
	if(targetDir == TD_FORE)       p = CVector(1.0, 0.0, 0.0);
	else if(targetDir == TD_BACK)  p = CVector(-1.0, 0.0, 0.0);
	else if(targetDir == TD_RIGHT) p = CVector(0.0, -1.0, 0.0);
	else if(targetDir == TD_LEFT)  p = CVector(0.0, 1.0, 0.0);
	else if(targetDir == TD_UP)    p = CVector(0.0, 0.0, 1.0);
	else if(targetDir == TD_DOWN)  p = CVector(0.0, 0.0, -1.0);//{ vP[0] = 0.0; vP[1] = 0.0; vP[2] = -1.0; }
	else {MessageBox(NULL," ���̕����ɂ͑Ή����Ă��Ȃ�","getDirVector",MB_OK);return;}

	p = rotate(p, vEuler);
}
//-----------------------------------------------------------------------------
void CTarget::getVertexOfCube()
{
    int i;
    vP[0].x = vSize.x / 2.0f; vP[0].y = vSize.y / 2.0f; vP[0].z = vSize.z / 2.0f;
    vP[1].x =-vSize.x / 2.0f; vP[1].y = vSize.y / 2.0f; vP[1].z = vSize.z / 2.0f;
    vP[2].x =-vSize.x / 2.0f; vP[2].y =-vSize.y / 2.0f; vP[2].z = vSize.z / 2.0f;
    vP[3].x = vSize.x / 2.0f; vP[3].y =-vSize.y / 2.0f; vP[3].z = vSize.z / 2.0f;
    vP[4].x = vSize.x / 2.0f; vP[4].y = vSize.y / 2.0f; vP[4].z =-vSize.z / 2.0f;
    vP[5].x =-vSize.x / 2.0f; vP[5].y = vSize.y / 2.0f; vP[5].z =-vSize.z / 2.0f;
    vP[6].x =-vSize.x / 2.0f; vP[6].y =-vSize.y / 2.0f; vP[6].z =-vSize.z / 2.0f;
    vP[7].x = vSize.x / 2.0f; vP[7].y =-vSize.y / 2.0f; vP[7].z =-vSize.z / 2.0f;

    //World���W�ɕϊ�
    for(i = 0; i < 8; i++){
        //��]
        vP[i] = rotate(vP[i], vEuler);
        //���S���W�������s�ړ�
        vP[i] += vPos;
    }
}
//-----------------------------------------------------------------------------
void CTarget::getVertexOfCylinder()
{
    int i;
    double theta, theta0;

    if(num1 > 20) { num1 = 20;}
	theta0 = (2.0*M_PI/(double)num1);
	for(i = 0;i < num1;i++)
	{   theta = theta0*(double)i;
		vP[i].x = (0.5*cos(theta)*vSize.x); //����x����
		vP[i].y = (0.5*sin(theta)*vSize.y); //������
		vP[i].z = 0.5f * vSize.z;                   //������(����)
		vP[i+num1].x = vP[i].x;                     //�����x����
		vP[i+num1].y = vP[i].y;                     //������
		vP[i+num1].z = - 0.5f * vSize.z;             //������
	}
    //���̒��S
    vP[2*num1].x = 0.0;
    vP[2*num1].y = 0.0;
    vP[2*num1].z = 0.5f * vSize.z;
    //����̒��S
    vP[2*num1+1].x = 0.0;
    vP[2*num1+1].y = 0.0;
    vP[2*num1+1].z = -0.5f * vSize.z;

    //World���W�ɕϊ�
    for(i = 0; i < 2 * num1 + 2; i++){
        //��]
        vP[i] = rotate(vP[i], vEuler);
        //���S���W�������s�ړ�
        vP[i] += vPos;
    }
}

//------------------------------------------------------------------------------
//�����̂ƒ����̂̏Փ˔���i�Q����)
//�����̂͊�{�p���iz���������)�ł��邱��
//���ڍ��̂̒��_���Ώۍ��̂̓����i���E����܂�)�ɂ���Ƃ������Փ˂Ƃ���
bool CTarget::collisionCubeToCube2(CTarget& trgt, CVector& vNormal)
{
    int i, j, ii, jj;
    double a[4], b[4], c[4];
    double len[4];//�ӂ̒���
    double dist[4];//���_����ӂ܂ł̋���
    double f;//���莮
    double min;
    int minNo, VertexNo, kaisu;
    CVector vCollision;//�Փ˓_
    CVector vn;

    getVertexOfCube();//���ڍ���(#1)�̒��_�̈ʒu���擾

    trgt.getVertexOfCube();//�ՓˑΏ�(#2)�̒��_�̈ʒu���擾

    //#1�̒��_��#2�̈ʒu�֌W
    for(j = 0; j < 4; j++)//j��#2�̒��_�ԍ�
    {
        if(j == 0 || j == 2) len[j] = trgt.vSize.x;
        else                 len[j] = trgt.vSize.y;
        jj = j + 1;
        if(jj == 4) jj = 0;
        a[j] = trgt.vP[j].y - trgt.vP[jj].y;
        b[j] = trgt.vP[jj].x - trgt.vP[j].x;
        c[j] = trgt.vP[j].x * trgt.vP[jj].y - trgt.vP[jj].x * trgt.vP[j].y;
    }
    for(i = 0; i < 4; i++) //#1�̒��_
    {
        kaisu = 0;//���莮�����ƂȂ��
        for(j = 0; j < 4; j++) //#2�̕Ӕԍ�
        {
            f = a[j] * vP[i].x + b[j] * vP[i].y + c[j];
            if( f < 0.0 ) break;  //
            //f���S�Đ��̂Ƃ��Փ�
            dist[j] = fabs(f) / len[j];//�ӂ܂ł̋���
            kaisu ++;
        }
        if( kaisu == 4) //#1�̒��_��#2�̓���
        {
            VertexNo = i;//#2�ɏՓ˂��Ă���#1�̒��_�ԍ�
            //�����P�̒��_���Փ˂��Ă��邩����
            ii = i+1;
            if(ii == 4) ii = 0;
            kaisu = 0;//���莮�����ƂȂ��
            for(j = 0; j < 4; j++)
            { //#2�̕Ӕԍ�j
                f = a[j] * vP[ii].x + b[j] * vP[ii].y + c[j];
                if( f < 0.0 ) break;
                kaisu ++;
            }
            if(kaisu == 4) goto calcNormalVector2;
            kaisu = 0;
            ii = i - 1;
            if(ii == -1) ii = 3;
            for(j = 0; j < 4; j++)
            { //#2�̕Ӕԍ�
                f = a[j] * vP[ii].x + b[j] * vP[ii].y + c[j];
                if( f < 0.0 ) goto calcNormalVector1;
                kaisu ++;
            }
            if(kaisu == 4) goto calcNormalVector2;
        }
    }

    return false;//���Փ�

calcNormalVector1:;
    //�ł��߂��Ӕԍ�minNo�����߂��̖@������vNormal������
    //#1�̒��_��#2�̕ӂɏՓ�
    min = dist[0];
    minNo = 0;
    for(j = 1; j < 4; j++)
    {
        if( min > dist[j] )
        {
            min = dist[j];
            minNo = j;
        }
    }
    //���ڍ��̂��猩���@�������itrgt�̖ʕ����̔��Ε����j
    if(minNo == 0)      vNormal = CVector(0.0, -1.0, 0.0);//��{�p���̖@������
    else if(minNo == 1) vNormal = CVector(1.0, 0.0, 0.0);
    else if(minNo == 2) vNormal = CVector(0.0, 1.0, 0.0);
    else                vNormal = CVector(-1.0, 0.0, 0.0);
    vNormal.rotZ(trgt.vEuler.z);//���̖��ł͂�����]����  */
    //�Փ˓_��#1�̒��_�ԍ�VertexNo��VertexNo+4�̒��_(0~3�͒����̂̏��j
    vCollision = (vP[VertexNo] + vP[VertexNo+4]) / 2.0;
    //�d�S����Փ˓_�܂ł̃x�N�g��
    vGravityToPoint = vCollision - vPos ;
    trgt.vGravityToPoint = vCollision - trgt.vPos ;
    return true;

calcNormalVector2:;
    vCollision = ( (vP[i] + vP[i+4]) + (vP[ii] + vP[ii+4]) )/ 4.0;
    vNormal = vPos >> vCollision;
    //�d�S����Փ˓_�܂ł̃x�N�g��
    vGravityToPoint = vCollision - vPos ;
    trgt.vGravityToPoint = vCollision - trgt.vPos ;
    return true;
}
//------------------------------------------------------------------------------
//���i�~��)�ƒ����̂̏Փ˔���i�Q����)
//���̏������S�Đ�������Ƃ��Փ˂Ƃ���
//��(�~��)�̒��S����ӂ܂ł̋��� <= ��(�~��)�̔��a
//�Q�̍��̂̒��S�������̂̕ӂ̔��Α��ɑ���
//��(�~��)�̒��S�������̂̕ӂɑ΂��ĕ��̈悪�P��
bool CTarget::collisionSphereToCube2(CTarget& trgt,CVector& vNormal)
{
    int j, jj;
    double a[4], b[4], c[4];
    double len[4];//�ӂ̒���
    double dist[4];//���_����ӂ܂ł̋���
    double f, f2;//���莮
    int kaisu;
    double rr =vSize.x / 2.0f;//���܂��͉~���̔��a
    CVector vCollision;//�Փ˓_
    trgt.getVertexOfCube();//�ՓˑΏ�(#2)�̒��_�̈ʒu���擾

    kaisu = 0;
    for(j = 0; j < 4; j++) //#2
    {
        if(j == 0 || j == 2) len[j] = trgt.vSize.x;
        else                 len[j] = trgt.vSize.y;
        jj = j + 1;
        if(jj == 4) jj = 0;
        a[j] = trgt.vP[j].y - trgt.vP[jj].y;
        b[j] = trgt.vP[jj].x - trgt.vP[j].x;
        c[j] = trgt.vP[j].x * trgt.vP[jj].y - trgt.vP[jj].x * trgt.vP[j].y;
        f = a[j] * vPos.x + b[j] * vPos.y + c[j];
        if(f < 0.0) kaisu ++;
    }
    for(j = 0; j < 4; j++)
    {   //#1�̒��S����#2�̕ӂ܂ł̋���
        f = a[j] * vPos.x + b[j] * vPos.y + c[j];
        dist[j] = fabs(f) / len[j];//�ӂ܂ł̋���
        if(dist[j] < rr && kaisu == 1)//�Փ˂̉\������
        {
//wsprintf(buf , "%d, %d, %d, %d",j,(int)(dist[j]*1000.0),(int)(f*1000.0),(int)(f2*1000.0));
//MessageBox(NULL,buf ,"collision5",MB_OK);
            f2 = a[j] * trgt.vPos.x + b[j] * trgt.vPos.y + c[j];
            if(f * f2 < 0.0)//2�̒��S���ӂ�����Ŕ��Α��ɂ���Ƃ��Փ�
            {
                if(j == 0)      vNormal = CVector( 0.0, -1.0, 0.0);//#2�̓�����#1���猩�Ė@������
                else if(j == 1) vNormal = CVector( 1.0,  0.0, 0.0);
                else if(j == 2) vNormal = CVector( 0.0,  1.0, 0.0);
                else            vNormal = CVector(-1.0,  0.0, 0.0);
                vNormal.rotZ(trgt.vEuler.z);
                //�Փ˓_
                vCollision.x = vPos.x + (vSize.x/2.0f) * (vNormal.x);//
                vCollision.y = vPos.y + (vSize.y/2.0f) * (vNormal.y);//
                vCollision.z = vPos.z;
                //�d�S����Փ˓_�܂ł̃x�N�g��
                vGravityToPoint = vCollision - vPos ;
                trgt.vGravityToPoint = vCollision - trgt.vPos ;
                return true;
            }
        }
    }
    return false;
}
//------------------------------------------------------------------------------
//�����̂Ƌ�(�~��)�̏Փ˔���i�Q����)
//�����̂̒��_���~���̓����ɑ��݂���Ƃ��Փ˂Ƃ���D
bool CTarget::collisionCubeToSphere2(CTarget& trgt, CVector& vNormal)
{
    int i;
//    double len[4];//�ӂ̒���
    double dist[4];//���_����ӂ܂ł̋���
//    double f2;//���莮
    double x, y;
    CVector vCollision;//�Փ˓_

    getVertexOfCube();//���_�̈ʒu���擾(vP[j]�Ɋi�[
    for(i = 0; i < 4; i++)
    {
        x = vP[i].x - trgt.vPos.x;
        y = vP[i].y - trgt.vPos.y;
        dist[i] = sqrt( x * x + y * y); //���C�~���̒��S���璼���̂̒��_�܂ł̋���
        if(dist[i] < trgt.vSize.x / 2.0)
        {//�Փ�            vP[i]�͒����̂̏�ӁCvP[i+4]�͂��̉���
            vNormal = (vP[i]+vP[i+4])/2.0 >> trgt.vPos;//�@������
            //�d�S����Փ˓_�܂ł̃x�N�g��
            vGravityToPoint = (vP[i]+vP[i+4])/2.0 - vPos; //������(���S����Փ˓_�֌������x�N�g��)
            trgt.vGravityToPoint = (vP[i]+vP[i+4])/2.0 - trgt.vPos;
            return true;
        }
    }
    return false;
}

//------------------------------------------------------------------------------
//�����̂ƒ����̂̏Փ˔���i�R����)
//���ڍ��̂̒��_���Ώۍ��̂̓����i���E����܂�)�ɂ���Ƃ��Փ�
//�ӂƕӂ���������Ƃ��Փ�
bool CTarget::collisionCubeToCube3(CTarget& trgt, CVector& vNormal)
{
    //�ʂ��\�����钸�_�ԍ�
    int vs[6][4] = { {0,1,2,3}, {0,3,7,4}, {0,4,5,1},
                     {1,5,6,2}, {2,6,7,3}, {4,7,6,5} };

    int i, j, k;
//    double d[6];
    CVector nfvP[4];//��_�ƕӂō��ʂ̖@���x�N�g��
    double f;//���莮
    double min, dd;
    int minNo, VertexNo[8], kaisu, cnt;
    CVector vCollision;//�Փ˓_
    CVector vPoint;
//    CVector vn, aa, bb;
    CVector ps;//�ʂƂ̌�_

    getVertexOfCube();//���ڍ���(#1)�̒��_�̈ʒu���擾
    trgt.getVertexOfCube();//�ՓˑΏ�(#2)�̒��_�̈ʒu���擾

    //#2�̖ʂ̖@���x�N�g��
    for(j = 0; j < 6; j++)//j��#2�̖ʔԍ�
    {
        trgt.vNormalFacet[j] = (trgt.vP[vs[j][1]] - trgt.vP[vs[j][0]]) ^ (trgt.vP[vs[j][2]] - trgt.vP[vs[j][1]]) ;
        trgt.vNormalFacet[j].normalize();
    }

    //���ڊp���̑S�Ă̒��_�ɂ��đΏۍ��̂̓����ɂ��邩�ǂ����𒲍�
    cnt = 0;//�Ώۍ��̓����ɂ��钍�ڍ��̂̒��_��
    for(i = 0; i < 8; i++) //#1�̒��_
    {
        kaisu = 0;//���莮�����ƂȂ��
        for(j = 0; j < 6; j++) //#2�̖ʔԍ�
        {
            f = trgt.vNormalFacet[j] * (vP[i] - trgt.vP[vs[j][0]]);
            if( f > 0.001 ) break;  //
            //f���S�ĕ��̂Ƃ��Փ�
//            dist[i][j] = fabs(f) / trgt.vNormalFacet[j].magnitude() ;//�ʂ܂ł̋���
            kaisu ++;
        }
        if( kaisu == 6) //#1�̒��_��#2�̓���
        {
            VertexNo[cnt] = i;//#2�ɏՓ˂��Ă���#1�̒��_�ԍ�
            cnt++;
        }
    }
    if(cnt == 0) goto EdgeAndEdge;

    //#2�̖ʂɏՓ˂��Ă���#1�̒��_�����ׂċ��ߕ��ϒl���Փ˓_�Ƃ���
//    vNormal = CVector(0.0, 0.0, 0.0);//�@���x�N�g���N���A
    vCollision = CVector(0.0, 0.0, 0.0);//�Փ˓_�̃N���A
    for(k = 0; k < cnt; k++) 
    {
        vCollision += vP[VertexNo[k]];//�Փ˓_��ǉ�
    }
    //�Փ˓_
    vCollision /= (double)cnt;//���ϒl
    //�ӂƕӂ��������Ă���΂��̕����̌�_��ǉ�������
    if(getPointCubeToCube(trgt, vNormal, vPoint) == true)
    {
        vCollision += vPoint;
        vCollision /= 2.0;
    }

    //�ł��߂��ʔԍ�minNo������
    f = trgt.vNormalFacet[0] * (vCollision - trgt.vP[vs[0][0]]);
    min = fabs(f) / trgt.vNormalFacet[0].magnitude() ;//�ʂ܂ł̋���
    minNo = 0;
    for(j = 1; j < 6; j++)//j�͑Ώۍ��̖̂ʔԍ�
    {
        //if(fabs(vSpeed * trgt.vNormalFacet[j]) < 0.5) continue;//���x�ɕ��s�Ȗʂ�����
        f = trgt.vNormalFacet[j] * (vCollision - trgt.vP[vs[j][0]]);
        dd = fabs(f) / trgt.vNormalFacet[j].magnitude() ;//�ʂ܂ł̋���
        if( min > dd )
        {
            min = dd;
            minNo = j;
        }
    }
    //���̖ʂ̖@���x�N�g���𔽓]
    vNormal = (-1.0) * trgt.vNormalFacet[minNo];

    //�d�S����Փ˓_�܂ł̃x�N�g��
    vGravityToPoint = vCollision - vPos ; //���ڍ��̑��ivPos�͒��S���W)
    trgt.vGravityToPoint = vCollision - trgt.vPos ;//�Ώۍ��̑�
    return true;

EdgeAndEdge:;//�ӂƕӂ̌���
    if(getPointCubeToCube(trgt, vNormal, vCollision) == false) return false;

    if(vNormal.x == 0.0 && vNormal.y == 0.0 && vNormal.z == 0.0 )
    {   //�㉺�i�܂��͍��E�j�̖ʂƌ���
        //�ł��߂��ʔԍ�minNo������
        f = trgt.vNormalFacet[0] * (vCollision - trgt.vP[vs[0][0]]);
        min = fabs(f) / trgt.vNormalFacet[0].magnitude() ;//�ʂ܂ł̋���
        minNo = 0;
        for(j = 1; j < 6; j++)//j�͑Ώۍ��̖̂ʔԍ�
        {
            f = trgt.vNormalFacet[j] * (vCollision - trgt.vP[vs[j][0]]);
            dd = fabs(f) / trgt.vNormalFacet[j].magnitude() ;//�ʂ܂ł̋���
            if( min > dd )
            {
                min = dd;
                minNo = j;
            }
        }
        vNormal = trgt.vNormalFacet[minNo];
    }
    vNormal.normalize();//�ׂ荇���ʂ̂Ƃ��͘a�̐��K���ŋߎ�
    //���̖ʂ̖@���x�N�g���𔽓]
    vNormal.reverse();

    //�d�S����Փ˓_�܂ł̃x�N�g��
    vGravityToPoint = vCollision - vPos ; //���ڍ��̑��ivPos�͒��S���W)
    trgt.vGravityToPoint = vCollision - trgt.vPos ;//�Ώۍ��̑�

    return true;
}
//---------------------------------------------------------------------
//���ڍ��̂̕ӂ��Ώۍ��̖̂ʂƌ������Ă���_���������ς��Ƃ�Փ˓_�Ƃ���
//�ӂ̗��[���P�̖ʂ̊O���ł���΂��̕ӂ͌������Ȃ�
//�ʂƖʂ��\����ɃN���X���Ă���Ƃ�������
bool CTarget::getPointCubeToCube(CTarget& trgt, CVector& vNormal, CVector& vPoint)
{
    //�ʂ��\�����钸�_�ԍ�
    int vs[6][4] = { {0,1,2,3}, {0,3,7,4}, {0,4,5,1},
                     {1,5,6,2}, {2,6,7,3}, {4,7,6,5} };
    //�ӂ��\�����钸�_�ԍ�
    int ve[12][2] = { {0,1}, {1,2}, {2,3}, {3,0},
                      {0,4}, {1,5}, {2,6}, {3,7},
                      {4,5}, {5,6}, {6,7}, {7,4} };

    int i, j, k, kp, kaisu, cnt;
    double fa, fb, tt;
    CVector vNormal0[4];//��_�ƕӂō��ʂ̖@���x�N�g��
    CVector vPoint0;    //��_

    kaisu = 0; //���ڍ��̂̕ӂ��Ώۍ��̖̂ʂƌ��������
    vPoint = CVector(0.0, 0.0, 0.0); //�����_
    vNormal = CVector(0.0, 0.0, 0.0); //�ʂ̖@�������̘a
    for(i = 0; i < 12; i++) //���ڊp���̕�
    {
        for(j = 0; j < 6; j++)//�Ώۊp���̖�
        {   //�ӂ̒��_���Ώۍ��̖̂ʂ̐��̈悩���̈悩
            fa = trgt.vNormalFacet[j] * (vP[ve[i][0]] - trgt.vP[vs[j][0]]);
            fb = trgt.vNormalFacet[j] * (vP[ve[i][1]] - trgt.vP[vs[j][0]]);
            if(fa * fb >= 0.0) continue;//�����̈�ɂ���Ό������Ȃ�
            tt = fa / (fa - fb);
            vPoint0 = (vP[ve[i][1]] - vP[ve[i][0]]) * tt + vP[ve[i][0]];//�ʂ��܂ޕ��ʂƂ̌�_
            cnt = 0;
            for(k = 0; k < 4; k++)//��_�����ʂ̕�
            {
                kp = k+1;
                if(kp == 4) kp = 0;
                vNormal0[k] = (trgt.vP[vs[j][k]] - vPoint0) ^ (trgt.vP[vs[j][kp]] - trgt.vP[vs[j][k]]) ;
                if(trgt.vNormalFacet[j] * vNormal0[k] < 0.0) break;//�P�ł����Ȃ�΂��̖ʂƂ͌������Ȃ�
                cnt++;
            }
            if(cnt == 4){//����
                kaisu++;
                vPoint += vPoint0;
                vNormal += trgt.vNormalFacet[j];
            }
        }
    }

    if(kaisu != 2 && kaisu != 4) return false;//�����Ȃ�
    //kaisu=4�͖ʂƖʂ��N���X
    vPoint /= (double)kaisu;
    return true;
}

//------------------------------------------------------------------------------
//���ƒ����̂̏Փ˔���i�R����)
//���̒��S����ӂ܂ł̋��� <= ���̔��a�ł���
//���̒��S�������̂̕ӂɑ΂��Đ��̈悪�P�̏ꍇ
//�����̂̕ӂƌ�������Ƃ����Փˁi���̈悪2�j
bool CTarget::collisionSphereToCube3(CTarget& trgt,CVector& vNormal)
{
    //�ʂ��\�����钸�_�ԍ�
    int vs[6][4] = { {0,1,2,3}, {0,3,7,4}, {0,4,5,1},
                     {1,5,6,2}, {2,6,7,3}, {4,7,6,5} };
    //�ӂ��\�����钸�_�ԍ�(�ŏ���2�j�Ɩʔԍ�
    int ve[12][4] = { {0,1,0,2}, {1,2,0,3}, {2,3,0,4}, {3,0,0,1},
                      {0,4,1,2}, {1,5,2,3}, {2,6,3,4}, {3,7,1,4},
                      {4,5,2,5}, {5,6,3,5}, {6,7,4,5}, {7,4,1,5} };
    int j, k, no1, no2;//, minJ;
    double dist[6];//���̒��S����ʂ܂ł̋���
    double f, d;//���莮
    int faceNo[6], cnt   ;
    double rr =vSize.x / 2.0f;//���̔��a
    CVector vCollision;//�Փ˓_
    CVector vDir;

    trgt.getVertexOfCube();//�ՓˑΏ�(#2)�̒��_�̈ʒu���擾
    //�Ώۊp���̖ʂ̖@���x�N�g��
    for(j = 0; j < 6; j++)//j��#2�̖ʔԍ�
    {
        trgt.vNormalFacet[j] = (trgt.vP[vs[j][1]] - trgt.vP[vs[j][0]]) ^ (trgt.vP[vs[j][2]] - trgt.vP[vs[j][1]]) ;
        trgt.vNormalFacet[j].normalize();
    }
    //���̒��S����Ώۊp���̖ʂ܂ł̋����𒲍�
    cnt = 0;//���̒��S���p���̖ʂ̐��̈�ɂ����
    for(j = 0; j < 6; j++) //�Ώۊp���̖ʔԍ�
    {
        f = trgt.vNormalFacet[j] * (vPos - trgt.vP[vs[j][0]]);
        if( f >= 0.0 ) {
            faceNo[cnt] = j;
            dist[faceNo[cnt]] = f;// / trgt.vNormalFacet[j].magnitude() ;//�ʂ܂ł̋���
            cnt++;
            vCollision = vPos - (trgt.vNormalFacet[j] * f); //cnt=1�̂Ƃ��̌�_���
        }
    }
    if(cnt == 1){ //�ʂƏՓ˂̉\������
        if(dist[faceNo[0]] > rr) return false;
        vNormal = (-1.0) * trgt.vNormalFacet[faceNo[0]] ;
        //�d�S����Փ˓_�܂ł̃x�N�g��
        vGravityToPoint = vCollision - vPos ; //���ڍ��̑��ivPos�͒��S���W)
        trgt.vGravityToPoint = vCollision - trgt.vPos ;//�Ώۍ��̑�
        return true;
    }
    else if(cnt == 2) { //�ӂƌ����̉\������
        //�ʂ����L����Ӕԍ�
        for(k = 0; k < 12; k++){
            if( faceNo[0] == ve[k][2] && faceNo[1] == ve[k][3] ){
                no1 = ve[k][0]; no2 = ve[k][1]; //��������ӂ̒��_�ԍ�
                break;
            }
        }
        //�ӂ̕����x�N�g��
        vDir = trgt.vP[no1] >> trgt.vP[no2]; //
        f = vDir * (trgt.vP[no1] - vPos) ;// (vP[no1] - trgt.vPos).magnitude();
        vCollision = trgt.vP[no1] - vDir * f;//���̒��S����ӂ֐��������낵�����̌�_
        d = (vCollision - vPos).magnitude();

        if(d > rr) return false;
        
        vNormal = vPos >> vCollision ;//���̒��S����Փ˓_�ւ̒P�ʃx�N�g��
        //�d�S����Փ˓_�܂ł̃x�N�g��
        vGravityToPoint = vCollision - vPos ; //���ڍ��̑��ivPos�͋��̒��S���W)
        trgt.vGravityToPoint = vCollision - trgt.vPos ;//�Ώۍ��̑�
        return true;
    }
    else return false;
}
//------------------------------------------------------------------------------
//������(�p��)�Ƌ��̏Փ˔���i�R����)
//�����̂̒��_�����̓����ɑ��݂���Ƃ��Փ�
bool CTarget::collisionCubeToSphere3(CTarget& trgt, CVector& vNormal)
{
    int i;
    double dist;
    int cnt;
    CVector vCollision;//�Փ˓_
//    CVector vDist;

    getVertexOfCube();//���_�̈ʒu���擾(vP[i]�Ɋi�[)
    vCollision = CVector(0.0, 0.0, 0.0);
    cnt = 0;
    for(i = 0; i < 8; i++)
    {
        //vDist = trgt.vPos - vP[i];
        dist = (trgt.vPos - vP[i]).magnitude(); //�����S���璼���̂̒��_�܂ł̋���
        if(dist < trgt.vSize.x / 2.0) //���̔��a�ȉ��Ȃ�
        {//�Փ�
            cnt++;
            vCollision += vP[i];
        }
    }
    if(cnt == 0) return false;

    vCollision /= (double)cnt;
    vNormal = vCollision >> trgt.vPos ; //�@������
    //�d�S����Փ˓_�܂ł̃x�N�g��
    vGravityToPoint = vCollision - vPos; //������(���S����Փ˓_�֌������x�N�g��)
    trgt.vGravityToPoint = vCollision - trgt.vPos;
    return true;
}
//------------------------------------------------------------------------------
//������(�p��)�Ɖ~���̏Փ˔���i�R����)
//�����̂̒��_���~���̓����ɑ��݂���Ƃ��Փ�
bool CTarget::collisionCubeToCylinder3(CTarget& trgt, CVector& vNormal)
{
    int i;
    double dist;//���_����~�����S���܂ł̋���
    double h1, h2;//���C����܂ł̋���
    double aveH1, aveH2;//���̕��ϒl
//    double x, y, min;
    int cnt;
    CVector vCollision;//�Փ˓_
    CVector vCenter;//�~�����S�������x�N�g��
    CVector vKoten;//�����̂̒��_����~���̒��S���։��낵�������̌�_
    CVector vAveKoten;//���̕��ϒl

    getVertexOfCube();//�����̂̒��_�̈ʒu���擾(vP[i]�Ɋi�[)
    trgt.getVertexOfCylinder();//�Ώۉ~���̒��_
    vCenter = trgt.vP[2*num1+1] >> trgt.vP[2*num1] ;//���S���P�ʃx�N�g��
    vCollision = CVector(0.0, 0.0, 0.0);
    cnt = 0;
    vAveKoten = CVector(0.0,0.0,0.0);
    aveH1 = 0.0; aveH2 = 0.0;
    for(i = 0; i < 8; i++)
    {
        vKoten = trgt.vPos-vCenter*(vCenter*(trgt.vPos - vP[i]));
        dist = (vKoten - vP[i]).magnitude(); //�����̂̒��_����~�����S���܂ł̋���
        if(dist < trgt.vSize.x / 2.0) //�~���̔��a�ȉ��Ȃ�Փ˂̉\������
        {
            h1 = fabs((vCenter * (vP[i] - trgt.vP[2*trgt.num1])) / vCenter.magnitude()); //���܂ł̋���
            h2 = fabs((vCenter * (vP[i] - trgt.vP[2*trgt.num1+1])) / vCenter.magnitude());//����܂ł̋���
            if(h1 <= trgt.vSize.z && h2 <= trgt.vSize.z){//�Փ�
                cnt++;
                vCollision += vP[i];
                vAveKoten += vKoten;
                aveH1 += h1;
                aveH2 += h2;
//wsprintf(buf , "%d, %d, %d, %d, %d",cnt, i,(int)(h1*1000.0), (int)(h2*1000.0),(int)(dist * 1000.0));
//MessageBox(NULL,buf ,"collisionA",MB_OK);
            }
        }
    }

    if(cnt == 0) return false;
    vCollision /= (double)cnt;
    vAveKoten /= (double)cnt;
    //�@�����������߂�
    if(aveH1 < 0.005 && aveH2 > trgt.vSize.z-0.005){//���ƏՓ�
        vNormal = (-1.0) * vCenter;
    }
    else if(aveH2 < 0.005 && aveH1 > trgt.vSize.z-0.005){//����ƏՓ�
        vNormal = vCenter;
    }
    else{
        vNormal = vCollision - vAveKoten;
        vNormal.normalize();
    }
    //�d�S����Փ˓_�܂ł̃x�N�g��
    vGravityToPoint = vCollision - vPos; //������(���S����Փ˓_�֌������x�N�g��)
    trgt.vGravityToPoint = vCollision - trgt.vPos;
    return true;
}
//------------------------------------------------------------------------------
//�~���ƒ����̂̏Փ˔���i�R����)
//�~���̒��_�������̂̓����i���E����܂�)�ɂ���Ƃ��Փ˂Ƃ���
//�ӂƕӂ̌������l��
bool CTarget::collisionCylinderToCube3(CTarget& trgt, CVector& vNormal)
{
    //�ʂ��\�����钸�_�ԍ��i�Ώے����́j
    int vs[6][4] = { {0,1,2,3}, {0,3,7,4}, {0,4,5,1},
                     {1,5,6,2}, {2,6,7,3}, {4,7,6,5} };
    int i, j, k;
    double dist[40][6];//���_����ʂ܂ł̋���
    double f;//���莮
    double min, dd;
    int minNo, VertexNo[8], kaisu, cnt, fNo[2], faceNo;
    CVector vCollision;//�Փ˓_
    CVector vPoint;
    CVector vCenter;//�~���̒��S��
    CVector vLine;//���̂̒��S�Ԃ����Ԓ���

    getVertexOfCylinder();//���ڍ���(#1,�~��)�̒��_�̈ʒu���擾
    trgt.getVertexOfCube();//�ՓˑΏ�(#2,������)�̒��_�̈ʒu���擾

    //#2�̖ʂ̖@���x�N�g��
    for(j = 0; j < 6; j++)//j��#2�̖ʔԍ�
    {
        trgt.vNormalFacet[j] = (trgt.vP[vs[j][1]] - trgt.vP[vs[j][0]]) ^ (trgt.vP[vs[j][2]] - trgt.vP[vs[j][1]]) ;
        trgt.vNormalFacet[j].normalize();
    }

    //���ډ~���̑S�Ă̒��_�ɂ��đΏۍ��̂̓����ɂ��邩�ǂ����𒲍�
    cnt = 0;//�Ώۍ��̓����ɂ��钍�ڍ��̂̒��_��
    for(i = 0; i < 2*num1; i++) //#1�̒��_
    {
        kaisu = 0;//���莮�����ƂȂ��
        for(j = 0; j < 6; j++) //#2�̖ʔԍ�
        {
            f = trgt.vNormalFacet[j] * (vP[i] - trgt.vP[vs[j][0]]);
            if( f > 0.001 ) break;  //
            //f���S�ĕ��̂Ƃ��Փ�
            dist[i][j] = fabs(f) / trgt.vNormalFacet[j].magnitude() ;//�ʂ܂ł̋���
            kaisu ++;
        }
        if( kaisu == 6) //#1�̒��_��#2�̓���
        {
            VertexNo[cnt] = i;//#2�ɏՓ˂��Ă���#1�̒��_�ԍ�
            cnt++;
        }
    }
    if(cnt == 0) goto EdgeAndEdge;
    
    vCollision = CVector(0.0, 0.0, 0.0);//�Փ˓_�̃N���A
    for(k = 0; k < cnt; k++)
    {
        vCollision += vP[VertexNo[k]];//�Փ˓_��ǉ�
    }
    //�Փ˓_
    vCollision /= (double)cnt;//���ϒl

//wsprintf(buf , "%d",cnt);
//MessageBox(NULL,buf ,"collisionA",MB_OK);

    //�����̖̂ʂƉ~���̑��ʂ��Փ�
    //�~�����S���̃x�N�g��
    vCenter = vP[2*num1+1] >> vP[2*num1];//���ꂩ����̒��S
    kaisu = 0;//�~���̒��S���������̖̂ʂɕ��s�Ŗʂ̊O���ɂ����
    for(j = 0; j < 6; j++)//�e�ʂɂ���
    {
        if(fabs(trgt.vNormalFacet[j] * vCenter) < 0.001){//���s(�����̖̂ʂ̖@���ƒ����j
            dd = trgt.vNormalFacet[j] * (vPos - trgt.vP[vs[j][0]]) / trgt.vNormalFacet[j].magnitude();
            if(dd > 0.0 && dd <= vSize.x / 2.0 + 0.001){
                fNo[kaisu] = j; kaisu ++;
            }
        }
    }

    if(kaisu == 1 || kaisu == 2)
    {
        if(kaisu == 1){//�����̖̂ʂƉ~���̑��ʂ����s�ɏՓ�
            vNormal = (-1.0) * trgt.vNormalFacet[fNo[0]];
        }
        else{//�����̂̕ӂƉ~���̑��ʂ����s�ɏՓ�(kaisu=2)
            vNormal = (-1.0) * (trgt.vNormalFacet[fNo[0]] + trgt.vNormalFacet[fNo[1]]) / (-2.0);
            vNormal.normalize();
        }
        //��_(2�̒��S����ʂ֍~�낵�������̌�_�̕���)
        if(getPointCylinderToCube(trgt, faceNo, vPoint) == true)
        {
            vCollision += vPoint;
            vCollision /= 2.0;
        }
        //�d�S����Փ˓_�܂ł̃x�N�g��
        vGravityToPoint = vCollision - vPos ; //���ڍ��̑��ivPos�͒��S���W)
        trgt.vGravityToPoint = vCollision - trgt.vPos ;//�Ώۍ��̑�
//wsprintf(buf , "%d, %d, %d, %d",kaisu, (int)(vNormal.x*1000.0), (int)(vNormal.y*1000.0),(int)(vNormal.z*1000.0));
//MessageBox(NULL,buf ,"collisionA",MB_OK);
        return true;
    }

    min = 1000.0;
    for(j = 0; j < 6; j++)//�Փ˓_����ł��߂��ʔԍ�
    {

        dd = fabs(trgt.vNormalFacet[j] * (vCollision - trgt.vP[vs[j][0]]))/trgt.vNormalFacet[j].magnitude() ;//�ʂ܂ł̋���
        if(min > dd){
            min = dd;
            minNo = j;
        }
    }

    //���̖ʂ̖@���x�N�g���𔽓]
    vNormal = -1.0 * trgt.vNormalFacet[minNo];

//    vNormal.normalize();
    //�d�S����Փ˓_�܂ł̃x�N�g��
    vGravityToPoint = vCollision - vPos ; //���ڍ��̑��ivPos�͒��S���W)
    trgt.vGravityToPoint = vCollision - trgt.vPos ;//�Ώۍ��̑�
    return true;

EdgeAndEdge:;//�ӂƕӂ̌���
    if(getPointCylinderToCube(trgt, faceNo, vCollision) == false) return false;
    //�~�����S���̃x�N�g��
    vCenter = vP[2*num1+1] >> vP[2*num1];//���ꂩ����̒��S
    //�Փ˓_����~�����S���ɉ��낵�������̌�_
    CVector vKoten = vP[2*num1+1] - vCenter * ( vCenter * (vP[2*num1+1] - vCollision));
    //�~�����猩���Փ˓_�̖@���x�N�g��
    vNormal = vKoten >> vCollision;
    //�~���̏��܂��͉��ꂪ�Փ˂��Ă���ꍇ
    //�~�������猩�����Α��x
    CVector vSpeed = (vVelocity - trgt.vVelocity).normalize2();

    if(vCenter * vSpeed > 0.5){ //��ꑤ���ՓˑΏ�
        if((trgt.vNormalFacet[faceNo] + vCenter).magnitude() <= 0.1) vNormal = vCenter; //���ɕ��s
    }
    else if(vCenter * vSpeed < -0.5){ //���ꑤ���ՓˑΏ�
        if((trgt.vNormalFacet[faceNo] - vCenter).magnitude() <= 0.1) vNormal = (-1.0) * vCenter;//����ɕ��s
    }
    //�d�S����Փ˓_�܂ł̃x�N�g��
    vGravityToPoint = vCollision - vPos ; //���ڍ��̑��ivPos�͒��S���W)
    trgt.vGravityToPoint = vCollision - trgt.vPos ;//�Ώۍ��̑�

    return true;
}
//---------------------------------------------------------------------
bool CTarget::getPointCylinderToCube(CTarget& trgt, int& faceNo, CVector& vPoint)
{
    //���ڍ��̂�edge���Ώۍ��̖̂ʂƌ������Ă���_���������ς��Ƃ�Փ˓_�Ƃ���
    //�����̗��[���P�̖ʂ̊O���ł���΂��̐����͌������Ȃ�
    //���ׂĂ̐����ɑ΂��ď�̏��������������΁C�Q�̍��͓̂Ɨ�
    //�ʂ��\�����钸�_�ԍ��i�Ώے����́j
    int vs[6][4] = { {0,1,2,3}, {0,3,7,4}, {0,4,5,1},
                     {1,5,6,2}, {2,6,7,3}, {4,7,6,5} };
    int i, j, k, kp, kaisu, cnt;
    double fa, fb, tt;
    CVector vNormal0[4];//��_�ƕӂō��ʂ̖@���x�N�g��
    CVector vPoint0;    //��_

    //�ӂ��\�����钸�_�ԍ�(���ډ~��)
    int ve[20][2];
    for(i = 0; i < num1; i++){
        ve[i][0] = i; ve[i][1] = i + num1;
    }

    kaisu = 0;
    vPoint = CVector(0.0, 0.0, 0.0);
    for(i = 0; i < num1; i++) //���ډ~���̐���
    {
        for(j = 0; j < 6; j++)//�Ώۊp���̖�
        {
            fa = trgt.vNormalFacet[j] * (vP[ve[i][0]] - trgt.vP[vs[j][0]]);
            fb = trgt.vNormalFacet[j] * (vP[ve[i][1]] - trgt.vP[vs[j][0]]);
            if(fa * fb >= 0.0) continue;//���̖ʂƂ͌������Ȃ�
            tt = fa / (fa - fb);
            vPoint0 = (vP[ve[i][1]] - vP[ve[i][0]]) * tt + vP[ve[i][0]];//�ʂ��܂ޕ��ʂƂ̌�_
            cnt = 0;
            for(k = 0; k < 4; k++)//��_�����ʂ̕�
            {
                kp = k+1;
                if(kp == 4) kp = 0;
                vNormal0[k] = (trgt.vP[vs[j][k]] - vPoint0) ^ (trgt.vP[vs[j][kp]] - trgt.vP[vs[j][k]]) ;
                if(trgt.vNormalFacet[j] * vNormal0[k] < 0.0) break;//�P�ł����Ȃ�΂��̖ʂƂ͌������Ȃ�
                cnt++;
            }
            if(cnt == 4){//����
                kaisu++; //�p���̖ʂƌ������Ă����
                vPoint += vPoint0;
                faceNo = j;
//wsprintf(buf , "%d, %d, %d, %d, %d",kaisu, i,j,(int)(ps.x*1000.0), (int)(ps.y*1000.0), (int)(ps.z*1000.0));
//MessageBox(NULL,buf ,"collisionA",MB_OK);
            }
        }
    }

    if(kaisu == 0) return false;//�Q�̍��͓̂Ɨ�

    vPoint /= (double)kaisu;
    return true;
}
//------------------------------------------------------------------
//���Ɖ~���̏Փ˔���i�R����)
//���̒��S����~���̒��S���܂ł̋��� <= ���̔��a�{�~���̔��a(���ʏՓˁj
//���܂��͉���ƏՓˁi���̒��S�͏��܂��͉���̊O��)
//�ǂ�����Փ˓_���~���̏��Ɖ���̒��Ԃɑ���
bool CTarget::collisionSphereToCylinder3(CTarget& trgt,CVector& vNormal)
{
    double dist;//���̒��S����~�����S���܂ł̋���
    double h1, h2;
    double rr =vSize.x / 2.0f;//���̔��a
    CVector vCollision;//�Փ˓_
    CVector vKoten;//���̒��S����~�����S���։��낵�������̌�_

    trgt.getVertexOfCylinder();//�Ώۉ~���̒��_
    //�~�����S���x�N�g��
    CVector vCenter = trgt.vP[2*trgt.num1 + 1] >> trgt.vP[2*trgt.num1];
    //���̒��S����~�����S���։��낵�������̌�_
    vKoten = trgt.vP[2*trgt.num1+1]-vCenter*(vCenter*(trgt.vP[2*trgt.num1+1] - vPos));//vPos�͋��̒��S
//char buf[40];
    //���̒��S�����_�܂ł̋���
    dist = (vKoten - vPos).magnitude();
    //�Փ˂̍Œ����
    if(dist > rr + trgt.vSize.x / 2.0) return false;

    //���̒��S������A����܂ł̋���
    h1 = fabs((vCenter * (vPos - trgt.vP[2*trgt.num1])));// / vCenter.magnitude()); //���܂ł̋���
    h2 = fabs((vCenter * (vPos - trgt.vP[2*trgt.num1+1])));// / vCenter.magnitude());//����܂ł̋���
    //���ʏՓ�
    if(h1 < trgt.vSize.z && h2 < trgt.vSize.z){
        vNormal = vKoten >> vPos;
        vCollision = vKoten + vNormal * (trgt.vSize.x / 2.0f) ;
    }
    //���܂��͉���ƏՓ�
    else{
        if(h1 <= rr && dist <= trgt.vSize.x / 2.0) {//��ꑤ�ŏՓ�
//            vNormal = (-1.0) * vCenter;
            vNormal = vCenter;
            vCollision = vPos - vNormal * (vNormal * (vPos - trgt.vP[2*num1]));
//vPos = vCollision + vNormal * vSize.x / 2.0;
//            vNormal.reverse();
        }
        else if(h2 <= rr && dist <= trgt.vSize.x / 2.0){ //���ꑤ�ŏՓ�
            vNormal = (-1.0) * vCenter;
            vCollision = vPos - vNormal * (vNormal * (vPos - trgt.vP[2*num1+1]));
//vPos = vCollision + vNormal * vSize.x / 2.0;
//            vNormal.reverse();
        }
        else return false;
    }
    vNormal.reverse();//���ڍ��̑����猩���@������
    //�Փ˓_
//    vCollision = vPos + rr * vNormal;
    //�d�S����Փ˓_�܂ł̃x�N�g��
    vGravityToPoint = vCollision - vPos ; //���ڍ��̑��ivPos�͋��̒��S���W)
    trgt.vGravityToPoint = vCollision - trgt.vPos ;//�Ώۍ��̑�
    return true;
}
//------------------------------------------------------------------------------
//�~���Ƌ��̏Փ˔���i�R����)
//�~���̒��_�����̓����ɑ��݂���Ƃ��Փ�
bool CTarget::collisionCylinderToSphere3(CTarget& trgt, CVector& vNormal)
{
    int i, cnt;
    double dist;
    CVector vCollision;//�Փ˓_
//    CVector vDist;

    getVertexOfCylinder();//���_�̈ʒu���擾(vP[i]�Ɋi�[)
    vCollision = CVector(0.0, 0.0, 0.0);
    cnt = 0;
    for(i = 0; i < 2*num1; i++)
    {
        dist = (trgt.vPos - vP[i]).magnitude(); //�����S����~���̒��_�܂ł̋���
        if(dist < trgt.vSize.x / 2.0) //���̔��a�ȉ��Ȃ�
        {//�Փ�
            cnt++;
            vCollision += vP[i];
        }
    }
    if(cnt == 0) return false;

    vCollision /= (double)cnt;
    vNormal = vCollision >> trgt.vPos;//�@������
    //�d�S����Փ˓_�܂ł̃x�N�g��
    vGravityToPoint = vCollision - vPos; //������(���S����Փ˓_�֌������x�N�g��)
    trgt.vGravityToPoint = vCollision - trgt.vPos;
    return true;
}
//------------------------------------------------------------------------------
//�~�����m�̏Փ˔���i�R����)
//�~���̒��_�������̉~���̓����i���E����܂�)�ɂ���Ƃ��Փ�
//�ӂƕӂ̌����ɂ��Փ�
bool CTarget::collisionCylinderToCylinder3(CTarget& trgt, CVector& vNormal)
{
//    double l1,m1,n1,l2,m2,n2;//���S���̕����x�N�g��
    double h1, h2, aveH1, aveH2, dist0, dist;
    CVector vCollision;//�Փ˓_
    CVector vKoten;//�Փ˓_����Ώۉ~���̒��S���։��낵�������̌�_
//    CVector vNormal0[3];
    CVector vp, vq, vKosaP, vKosaQ;
    double s, t, p2, q2, d;
    int ve[20][2]; //�ӂ��\�����钸�_�ԍ�(���ډ~��)
    int i, cnt;

    getVertexOfCylinder();//���ډ~��(#1,�~��)�̒��_�̈ʒu���擾
    trgt.getVertexOfCylinder();//�ՓˑΏ�(#2,�~��)�̒��_�̈ʒu���擾
    CVector vCenter1 = vP[2*num1+1] >> vP[2*num1];//���ډ~���̒��S���x�N�g��
    CVector vCenter2 = trgt.vP[2*num1+1] >> trgt.vP[2*num1];//�Ώۉ~���̒��S���x�N�g��
    //���S���Ԃ̍ŒZ����
    if( (vCenter1 - vCenter2).magnitude() == 0.0
        || (vCenter1 + vCenter2).magnitude() == 0.0){  //�Q���������s
        //���ډ~���̒��S����Ώے��S���։��낵�������̌�_
//        vKoten = trgt.vP[2*trgt.num1+1]-vCenter2*(vCenter2*(trgt.vP[2*trgt.num1+1] - vPos));
        vKoten = trgt.vPos - vCenter2 * (vCenter2 * (trgt.vPos - vPos));
        dist0 = (vKoten - vPos).magnitude();
    }
    else
        dist0 = fabs((vPos - trgt.vPos) * (vCenter1 ^ vCenter2));

    //�ՓˍŒ����
    if(dist0 > vSize.x/2.0 + trgt.vSize.x/2.0) return false;

    //���ډ~���̑S�Ă̒��_�ɂ��đΏۍ��̂̓����ɂ��邩�ǂ����𒲍�
    cnt = 0;//�Ώۍ��̓����ɂ��钍�ڍ��̂̒��_��
    aveH1 = 0.0; aveH2 = 0.0;
    for(i = 0; i < 2*num1; i++) //#1�̒��_
    {
        vKoten = trgt.vP[2*trgt.num1+1]-vCenter2*(vCenter2*(trgt.vP[2*trgt.num1+1] - vP[i]));
        dist = (vKoten - vP[i]).magnitude(); //���ډ~���̒��_����Ώۉ~�����S���܂ł̋���
        if(dist < trgt.vSize.x / 2.0) //�~���̔��a�ȉ��Ȃ�Փ˂̉\������
        {
            h1 = fabs((vCenter2 * (vKoten - trgt.vP[2*trgt.num1]))); //���܂ł̋���
            h2 = fabs((vCenter2 * (vKoten - trgt.vP[2*trgt.num1+1])));//����܂ł̋���
            if(h1 <= trgt.vSize.z+0.0001 && h2 <= trgt.vSize.z+0.0001){//�Փ�
                cnt++;
                aveH1 += h1;
                aveH2 += h2;
                vCollision += vP[i];
            }
        }
    }

    if(cnt == 0) goto EdgeAndEdge;
    //#2�̖ʂɏՓ˂��Ă���#1�̒��_�̕��ϒl���Փ˓_�Ƃ���
    vCollision /= (double)cnt;
    //���S���։��낵�������̌�_
    vKoten = trgt.vP[2*trgt.num1+1]-vCenter2*(vCenter2*(trgt.vP[2*trgt.num1+1] - vCollision));
    //�@�����������߂�
    if(aveH1 < 0.001 && aveH2 > trgt.vSize.z-0.001){//���ƏՓ�
        vNormal = (-1.0) * vCenter2;
    }
    else if(aveH2 < 0.001 && aveH1 > trgt.vSize.z-0.001){//����ƏՓ�
        vNormal = vCenter2;
    }
    else{
        vNormal = vCollision - vKoten;
        vNormal.normalize();
    }
    //�d�S����Փ˓_�܂ł̃x�N�g��
    vGravityToPoint = vCollision - vPos ; //���ڍ��̑��ivPos�͒��S���W)
    trgt.vGravityToPoint = vCollision - trgt.vPos ;//�Ώۍ��̑�
    return true;

EdgeAndEdge:;//�ӂƕӂ̌���
    //���ډ~����edge���Ώۉ~���ƌ������Ă���_���������ς��Ƃ�Փ˓_�Ƃ���
    //edge�ƑΏۉ~���̒��S���ƍŒZ�ƂȂ�ʒu�����߂��̓_���ǂ���̉~���ɑ΂��Ă�
    //���Ɖ���̊Ԃɑ��݂���Ό���

    //���ډ~����edge�̒��_�ԍ�
    for(i = 0; i < num1; i++){
        ve[i][0] = i; ve[i][1] = i + num1;
    }

    vq = vCenter2;
    q2 = vq * vq;
    cnt = 0;
    vCollision = CVector(0.0, 0.0, 0.0);
    for(i = 0; i < num1; i++)
    {
        vp = (vP[ve[i][0]] - vP[ve[i][1]]).normalize2();
        p2 = vp * vp;
        d = p2*q2-(vp*vq)*(vp*vq);
        if(d == 0){
            //MessageBox(NULL,"collisionCylinderToCylinder3���[�`���ɂ�����d=0" ,"Error",MB_OK);
            return false;
        }
        s = ( q2*((trgt.vP[2*trgt.num1+1]-vP[ve[i][1]])*vp)
            + (vp*vq)*((vP[ve[i][1]]-trgt.vP[2*trgt.num1+1])*vq) ) / d;
        t = ( vp*vq*((trgt.vP[2*trgt.num1+1]-vP[ve[i][1]])*vp)
            + (p2)*((vP[ve[i][1]]-trgt.vP[2*trgt.num1+1])*vq) ) / d;
        //�����_
        vKosaP = vP[ve[i][1]] + s * vp;
        vKosaQ = trgt.vP[2*num1+1] + t * vq;
        if( (vKosaP-vKosaQ).magnitude() > trgt.vSize.x / 2.0) continue;
        //���ډ~�����̔���
        h1 = fabs((vp * (vKosaP - vP[ve[i][0]]))); //���܂ł̋���
        h2 = fabs((vp * (vKosaP - vP[ve[i][1]])));//����܂ł̋���
        if(h1 > vSize.z || h2 > vSize.z) return false;
        //�Ώۉ~�����̔���
        h1 = fabs((vq * (vKosaQ - trgt.vP[2*trgt.num1]))); //���܂ł̋���
        h2 = fabs((vq * (vKosaQ - trgt.vP[2*trgt.num1+1])));//����܂ł̋���
        if(h1 > trgt.vSize.z || h2 > trgt.vSize.z) return false;
        vCollision += vKosaP;
        cnt++;
    }
    if(cnt == 0) return false;
    vCollision /= (double)cnt;
    //���ډ~�����S���Ƃ̌�_
    vKoten = vP[2*num1+1]-vCenter1*(vCenter1*(vP[2*num1+1] - vCollision));
    vNormal = vKoten >> vCollision;
    //�d�S����Փ˓_�܂ł̃x�N�g��
    vGravityToPoint = vCollision - vPos ; //���ڍ��̑��ivPos�͒��S���W)
    trgt.vGravityToPoint = vCollision - trgt.vPos ;//�Ώۍ��̑�

    return true;
}
//------------------------------------------------------------------------------
//�����m�̏Փ�(�K���Փˁj
bool CTarget::collisionSphereToSphere3(CTarget& trgt, CVector& vNormal)
{
    vNormal = vPos >> trgt.vPos;//��0���狅1�֌������P�ʖ@���x�N�g��
    //�d�S����Փ˓_�܂ł̃x�N�g��
    vGravityToPoint = (vSize.x/2.0f) * vNormal ;
    trgt.vGravityToPoint = (-trgt.vSize.x/2.0f) * vNormal ;
    return true;
}
//------------------------------------------------------------------------------
//�ȉ��̓\���b�h�e�N�X�`���p
//--------------------------------------------------------------------------------
//texImage��O�����č쐬
void CTarget::calcSolidTexCube(CVector vSize)
{
	int i, j, k;
	double x, y, z;
    byte rgb[3];

	//Top
	z = vSize.z / 2.0f;
	for(j = 0; j < T_MAX; j++)
	{
		x = vSize.x * (0.5f - (double)j / T_MAX);//j��0��texImage�̉���
		for(i = 0; i < T_MAX; i++)
		{
			y = vSize.y * ( -0.5f + (double)i / T_MAX);
			//getImageGrain(i, j, x, y, z, 0);
			getImageGrain(x, y, z, vSize, rgb);
            for(k = 0; k < 3; k++) texImageS[0][j][i][k] = rgb[k];
            //texImageS[0][j][i][3] = 255;//�s����
		}
	}

	//Front(x = 0.5)
	x = vSize.x / 2.0f;
	for(j = 0; j < T_MAX; j++)
	{
		z = vSize.z * (-0.5f + (double)j / T_MAX);//texImage�͉�������
		for(i = 0; i < T_MAX; i++)
		{
			y = vSize.y * (-0.5f + (double)i / T_MAX);
			//getImageGrain(i, j, x, y, z, 1);
			getImageGrain(x, y, z, vSize, rgb);
            for(k = 0; k < 3; k++) texImageS[1][j][i][k] = rgb[k];
            //texImageS[1][j][i][3] = 255;//�s����
		}
	}

	//Left
	y = vSize.y / 2.0f;
	for(j = 0; j < T_MAX; j++)
	{
		z = vSize.z * (-0.5f + (double)j / T_MAX);//texImage�͉�������
		for(i = 0; i < T_MAX; i++)
		{
			x = vSize.x * (0.5f - (double)i / T_MAX);
			getImageGrain(x, y, z, vSize, rgb);
            for(k = 0; k < 3; k++) texImageS[2][j][i][k] = rgb[k];
            //texImageS[2][j][i][3] = 255;//�s����
		}
	}

	//Rear
	x = - vSize.x / 2.0f;
	for(j = 0; j < T_MAX; j++)
	{
		z = vSize.z * (-0.5f + (double)j / T_MAX);//texImage�͉�������
		for(i = 0; i < T_MAX; i++)
		{
			y = vSize.y * (0.5f - (double)i / T_MAX);
			getImageGrain(x, y, z, vSize, rgb);
            for(k = 0; k < 3; k++) texImageS[3][j][i][k] = rgb[k];
            //texImageS[3][j][i][3] = 255;//�s����
		}
	}

	//Right
	y = - vSize.y / 2.0f;
	for(j = 0; j < T_MAX; j++)
	{
		z = vSize.z * (-0.5f + (double)j / T_MAX);//texImage�͉�������
		for(i = 0; i < T_MAX; i++)
		{
			x = vSize.x * (-0.5f + (double)i / T_MAX);
//			getImageGrain(i, j, x, y, z, 4);
			getImageGrain(x, y, z, vSize, rgb);
            for(k = 0; k < 3; k++) texImageS[4][j][i][k] = rgb[k];
            //texImageS[4][j][i][3] = 255;//�s����
		}
	}

	//Bottom(z = -0.5)
	z = - vSize.z / 2.0f;
	for(j = 0; j < T_MAX; j++)
	{
		x = vSize.x * (-0.5f + (double)j / T_MAX);//j��0��texImage�̉���
		for(i = 0; i < T_MAX; i++)
		{
			y = vSize.y * (-0.5f + (double)i / T_MAX);
//			getImageGrain(i, j, x, y, z, 5);
			getImageGrain(x, y, z, vSize, rgb);
            for(k = 0; k < 3; k++) texImageS[5][j][i][k] = rgb[k];
            //texImageS[5][j][i][3] = 255;//�s����
		}
	}
}
//---------------------------------------------------------------------------------------
void CTarget::calcSolidTexSphere(CVector vSize)
{
	int i, j, k;
	double phi; //x-y���ʂɑ΂���Ίp�i�ܓx�j
	double theta; //x���ɑ΂���Ίp�i�o�x)
	double x, y, z, rr, r1, r2;
    byte rgb[3];

	//�e�N�`���[�̍쐬(������쐬)
	for(j = 0; j < T_MAX; j++)
	{
		phi = (M_PI * (double)j / (double)T_MAX);//phi�͓�ɂ���̈ܓx
		z = - vSize.z * cos(phi) / 2.0f;//z=0�͓��
		rr =   sin(phi) ;
		r1 = vSize.x * rr; r2 = vSize.y * rr;
		for(i = 0; i < T_MAX; i++)
		{
			theta = (2.0f * M_PI * (double)i / (double)T_MAX);//x=0����̌o�x
			x = r1 * cos(theta) ;
			y = r2 * sin(theta) ;
			getImageGrain(x, y, z, vSize, rgb);
            for(k = 0; k < 3; k++) texImageS[0][j][i][k] = rgb[k];
            //texImageS[0][j][i][3] = 255;//�s����
		}
	}
}
//----------------------------------------------------------------------------------------
void CTarget::calcSolidTexCylinder(CVector vSize)
{
	int i, j, k;
	double x, y, z;
	double th;
	double pp = 2.0 * M_PI;
    byte rgb[3];

	//Top��ø����
	z = (vSize.z / 2.0);
	for(j = 0; j < T_MAX; j++)//
	{
		x = vSize.x * (0.5f - (double)j / T_MAX);//j��0��texImage�̉���
		for(i = 0; i < T_MAX; i++)
		{
			y = vSize.y * ( - 0.5f + (double)i / T_MAX);
			getImageGrain(x, y, z, vSize, rgb);
            for(k = 0; k < 3; k++) texImageS[0][j][i][k] = rgb[k];
		}
	}

	//Bottom
	z = - vSize.z / 2.0f;
	for(j = 0; j < T_MAX; j++)//
	{
		x = vSize.x * ( 0.5f - (double)j / T_MAX);//j��0��texImage�̉���
		for(i = 0; i < T_MAX; i++)
		{
			y = vSize.y * ( 0.5f - (double)i / T_MAX);
			getImageGrain(x, y, z, vSize, rgb);
            for(k = 0; k < 3; k++) texImageS[1][j][i][k] = rgb[k];
		}
	}

	//side
	for(j = 0; j < T_MAX; j++)//
	{
		z = vSize.z * (- 0.5f + (double)j / T_MAX);//j��0��texImage�̉���
		for(i = 0; i < T_MAX; i++)
		{
			th = (pp * (double)i / T_MAX - M_PI / 2.0);
			x = vSize.x * cos(th) / 2.0f;
			y = vSize.y * sin(th) / 2.0f;
			getImageGrain(x, y, z, vSize, rgb);
            for(k = 0; k < 3; k++) texImageS[2][j][i][k] = rgb[k];
		}
	}
}

//--------------------------------------------------------------------------------
void CTarget::getImageGrain(double x, double y, double z, CVector vSize, byte* rgb)
{
	double period = 2.0;//vSize=10���
	double dist;
	double a1, a2;
	double x0 = 0.08;
	double y0 = 0.5;//random(100)/100.0;//  0.1;
	double r, g, b;

	//�U��
	a1 = 30.0;
	a2 = 20.0;
	//�ؖڂ̎���
	period = 0.05f * (1.0f - 0.005f * z);
	//���S����̋���
	dist = sqrt((x-x0) * (x-x0) + (y-y0) * (y-y0));
	//�F����
	r = a1 * (sin(2.0 * M_PI * dist / period));
	g = a2 * (sin(2.0 * M_PI * dist / period));
	b = 50.0;
	//red
	if( r > -0.8 * a1 )
        rgb[0] = 220 + (byte)r;
	else
        rgb[0] = 150 + (byte)r;
	//green
	if( g > - 0.8 * a2)
        rgb[1] = 120 + (byte)g;
	else
        rgb[1] = 100 + (byte)g;
	//blue
    rgb[2] = (byte)b;
 }
//----------------------------------------------------------------------------------------

