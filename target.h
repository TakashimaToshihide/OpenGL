//target.h

enum Kind {SPHERE, CUBE, CYLINDER, PRISM, PYRAMID, CONE, FRUSTUM, TORUS,
           SEMI_SPHERE, SUPER1, SUPER2, TAPE1, TAPE2, MESH, SQUARE,
           TEX_MESH1, TEX_MESH2, TEX_MESH3, TEX_MESH4};
enum TargetDir { TD_DOWN,TD_UP,TD_FORE,TD_BACK,TD_RIGHT,TD_LEFT};
//描画オブジェクト
class CTarget{ 
public:
	//ﾒﾝﾊﾞ変数
	Kind kind;     //CTargetの種類
	CVector vSize; //サイズ
	CVector vPos;  //位置
    CVector vEuler;//オイラー角（姿勢表現)
    CVector vP[40];//直方体,円柱の頂点座標
    CVector vNormalFacet[6];//直方体の面の法線ベクトル
    double dir;    //進行方向
    //初期値を保存しておきたいときに使用
    CVector vPos0;
    CVector vEuler0;
    double dir0;
    //運動シミュレーション用
    double mass;        //質量
    CVector vForce;    //外力
    CVector vForce0;   //外力初期値（現在値としても利用)
    CVector vAccel;    //加速度
    CVector vVelocity; //速度
    CVector vVelocity0;//初速度
    CVector vOmega;    //角速度
    CVector vOmega0;   //初角速度
    CVector vAlpha;    //角加速度
    double omega0;      //任意軸回りの角速度
    CVector vAxis;     //その回転軸
    CVector vGravityToPoint;//重心から衝突点までのベクトル
    CMatrix mInertia; //慣性モーメントのテンソル
    CMatrix mInertiaInverse;//その逆行列
    CQuaternion q;   //四元数
    double speed;     //速度の絶対値
    double speed0;
	bool flagFixed;
    //マテリアル特性
	float diffuse[4]; //拡散光係数
	float specular0;  //鏡面反射光係数
	float ambient0;   //環境光係数
	float highlight;  //光沢
	short num1; //分割数(Nxy,Nm)
	short num2; //分割数(Nz,Ns)
	float ratio;//Frustum(下底に対する上底半径比率)
				 //Torus(中心軸半径に対する断面半径比率
	float eps1, eps2, p1, p2, p3;//超２次関数
    float amp, lambdaX, lambdaY, lambdaZ, alpha, weight, wind;//TEX_MESH2,TEX_MESH3
    float data[10201]; //波のデータ(100*100)
    int texMode;//DECAL,MODULATE
    int texType;//PARALLEL1,PARALLEL2,CYLINDRICAL,SPHERICAL,SOLID
    byte texImage[T_MAX][T_MAX][4];
    byte texImageS[6][T_MAX][T_MAX][3]; //SolidTexture用
    short numStress;//particleに作用する応力の個数
    float stress;//応力の平均値
    short posType;//角(0)，辺(1)，面(2)，内部(3)

	//ﾒﾝﾊﾞ関数
	CTarget();       //コンストラクタ
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
//コンストラクタ
CTarget::CTarget()
{
	kind = CUBE;
    vSize = CVector(0.1f, 0.1f, 0.1f);
    vPos = CVector();//原点
    vPos.z = vSize.z / 2.0f;//位置はｵﾌﾞｼﾞｪｸﾄの中心
//    flagFixed = false;
    //物体色（拡散光を赤色に設定)
    diffuse[0] = 1.0;//red
    diffuse[1] = 0.0;//green
    diffuse[2] = 0.0;//blue
    diffuse[3] = 1.0;//alpha値（不透明に設定)
	ambient0 = 0.4f;
	specular0 = 1.0f;
    highlight = 80.0f;
    num1 = 10;
    num2 = 10;
    ratio = 0.5;//Frustum,Torus
    eps1 = 1.0;
    eps2 = 1.0;
    p1 = p2 = p3 = 1.0;//超2次関数の形状ﾊﾟﾗﾒｰﾀ
}
//---------------------------------------------------------
void CTarget::draw(bool flagShadow)
{
    float shadowDiffuse[] = {0.2f,0.25f,0.25f,0.3f};//影の拡散光
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
    //現在位置
	glTranslated(vPos.x, vPos.y, vPos.z);//平行移動

	//回転
    glRotated(vEuler.z, 0.0f, 0.0f, 1.0f);//z軸回転
    glRotated(vEuler.y, 0.0f, 1.0f, 0.0f);//y軸回転
    glRotated(vEuler.x, 1.0f, 0.0f, 0.0f);//x軸回転
	//スケーリング
	glScaled(vSize.x, vSize.y, vSize.z);

	//ｵﾌﾞｼﾞｪｸﾄの描画
	//SolidModel
    if(kind == CUBE && texType == T_NON)
        skSolidCube();
    else if(kind == SPHERE && texType == T_NON)
        skSolidSphere(num1);
    else if(kind == PRISM)
	   	skSolidPrism(num1);
    else if(kind == CYLINDER && texType == T_NON)//円柱
	   	skSolidCylinder(num1);
    else if(kind == PYRAMID)//多角形
	  	skSolidPyramid(num1);
    else if(kind == CONE)//円錐
	 	skSolidCone(num1);
    else if(kind == FRUSTUM)//円錐台
		skSolidFrustum(num1,ratio);
    else if(kind == TORUS)//トーラス
   		skSolidTorus(num1,num2,ratio);
    else if(kind == SEMI_SPHERE)
		skSolidSemiSphere(num1);
    else if(kind == SUPER1 && texType == T_NON)//超２次関数（上下対称）
	   	skSolidSuper1(num1,num2,eps1,eps2,p1,p2,p3);
    else if(kind == SUPER2)//超２次関数（上部だけ）
	  	skSolidSuper2(num1,num2,eps1,eps2,p1,p2);
    //波のプリミティブとして3種類追加
    else if(kind == TAPE1)
        skSolidTape1(num1, data);
    else if(kind == TAPE2)
        skSolidTape2(num1, data);
    else if(kind == MESH)
        skSolidMesh(num1, data);
    //Texture
    else if(kind == SQUARE && texType != T_NON)//平面投影
        skTexSquare();
    else if(kind == CUBE && texType == T_PLANAR1)//平面投影
        skTexCube1();
    else if(kind == CUBE && texType == T_PLANAR2)//平面投影
        skTexCube2();
    else if(kind == SPHERE && texType == T_PLANAR1)//平面投影
        skTexSphere1(num1);
    else if(kind == SPHERE && texType == T_SPHERICAL)//球面投影
        skTexSphere2(num1);
    else if(kind == CYLINDER && texType == T_PLANAR1)//平行投影
        skTexCylinder1(num1);
    else if(kind == CYLINDER && texType == T_CYLINDRICAL)//円筒投影
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
{//bitmapからテクスチャ画像の作成
    int i, j, ii, jj;
    int nx = 128;
    int ny = 128;
    double pWidth = (double)bitmap->Width;
    double pHeight = (double)bitmap->Height;
    if( pHeight == 0.0) {
        MessageBox(NULL,"makeTexture","画像がありません",MB_OK);
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
    double torque;//床面摩擦による回転力

    //軸方向を調査
    getDirVector(TD_FORE, vFore);
    getDirVector(TD_UP, vUp);
    getDirVector(TD_LEFT, vLeft);

    if(kind == SPHERE) torque = 0.01f;//小さな値

    else if(kind == CYLINDER){
        r = vSize.x / 2.0f; h = vSize.z;
        if(fabs(vUp.z) >= 0.99)//fabs(vFore.z) && fabs(vUp.z) >= fabs(vLeft.z) )
        {//基本姿勢で回転しているとき
            torque = 2.0f * mu * mass * g * r / 3.0f;
        }
        else if(fabs(vUp.z) < 0.02){//横になって回転しているとき
            torque = mu * mass * g * h / 4.0f;
        }
        else torque = 0.01f;//小さい値とする

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
    if(vOmega.z > 0.0) torque = - torque;//減速になるように
    return( torque ); //
}

//------------------------------------------------------------------------------
//Floor上を回転する円柱,角柱の摩擦による角加速度を計算
double CTarget::getAngularAccelOnFloor(double mu)
{
    double a, b, a2, b2, a3, b3, c, ss , r, h;
    CVector vFore, vUp, vLeft;
	double g = 9.8f;//m/s^2
    double alpha;//角加速度
    double inertia;

    //軸方向を調査
    getDirVector(TD_FORE, vFore);
    getDirVector(TD_UP, vUp);
    getDirVector(TD_LEFT, vLeft);
//if(kind == SPHERE) alpha = 5.0;
    if(kind == CYLINDER){
        r = vSize.x / 2.0f; h = vSize.z;
        if(fabs(vUp.z) >= fabs(vFore.z) && fabs(vUp.z) >= fabs(vLeft.z) )
        {//基本姿勢で回転しているとき
            alpha = 4.0f * mu * g / (3.0f * r);
        }
        else{//横になって回転しているとき
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
            MessageBox(NULL,"計算ｴﾗｰ" ,"calcAngularAccelOnFloor",MB_OK); return(0);
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
//dt後の位置，回転角，速度を求める
void CTarget::rolling(double mu, double dt)
{                   //転がり摩擦係数,時間刻み
	double dist;
	double g = 9.8f;
	double pp = M_PI / 180.0;
    double ang;
    CQuaternion q;
    //回転軸は水平で進行方向に直行する軸
    vAxis.x = cos((dir + 90.0) * pp);
    vAxis.y = sin((dir + 90.0) * pp);
    vAxis.z = 0.0;

    speed -= mu * g * dt; //speedの更新
    dist = speed * dt;    //移動距離
    ang = (2.0f * dist / (vSize.z * pp)); //回転角（転がり運動による回転角)
    vPos.x += dist * cos(dir * pp);
    vPos.y += dist * sin(dir * pp);
    q = makeQFromAxis( ang , vAxis); //Quaternionを作成
    vEuler = makeEulerFromEuler(vEuler, q );//オイラー角の更新
}
//-------------------------------------------------------------------
//dt後の位置，回転角，速度を求める
//オブジェクトが基本姿勢のときだけ正しい
void CTarget::sliding(double mu, double dt)
{                   //動摩擦係数,時間刻み
	double a, aa, dist;
	double g = 9.8f;
	double pp = M_PI / 180.0;
    //回転軸
    vAxis = CVector(0.0, 0.0, 1.0);

    //加速度
    a = - mu * g;
    //角加速度(摩擦による角加速度の減速）
    aa = getAngularAccelOnFloor(mu);//rad/s^2で取得
    if(vOmega.z > 0.0) aa = -aa;

    //新速度
    speed += a * dt;
    //新角速度
    vOmega.z += aa * dt;

    //移動距離
    dist = speed * dt;
    //位置の更新
    vPos.x += dist * cos(dir * pp);
    vPos.y += dist * sin(dir * pp);
    //回転角度の更新
    vEuler.z += vOmega.z * dt / pp; //オイラー角はdeg単位
}
//-----------------------------------------------------------------------------
//角から中心までの高さの最大値（球と直方体だけ）
double CTarget::getMaxH(void)
{
	short i;
    CVector xyz;
	double h[4],maxH;

	if(kind == SPHERE)
        return(vSize.z / 2.0f);//球の半径を返す

	for(i = 0;i < 4;i++){//立方体
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

		//x軸回転
        xyz.rotX(vEuler.x);
		//y軸回転
        xyz.rotY(vEuler.y);
		h[i] = fabs(xyz.z);//他は求める必要なし
	}

	//最大値を求める
	maxH = h[0];
	for(i=1;i<4;i++){
		if(maxH < h[i]) maxH = h[i];
	}
	return(maxH);
}
//------------------------------------------------------------------------------

void CTarget::projection(double e, double mu, double cv, double ci, double dt)
                //反発係数, 動摩擦抵抗,粘性抵抗係数,慣性抵抗係数,刻み時間
{
    double minH;
    double rikiseki ;
    CVector vNormal;//注目物体から見た衝突面の法線方向
    CVector vp;//衝突面の合成速度（線形速度＋角速度による接線速度)
    CVector gtp;//重心から衝突点までのベクトル
    CVector vTorque = CVector();
    double ang;
    CQuaternion q;
	double pp = M_PI / 180.0;

    //力の初期値
    vForce = vForce0;//vForce0はmain側で設定

    //現在の速度
    speed = vVelocity.magnitude();
    //空気の抵抗
    vForce -= cv * vVelocity;//cvは速度に比例する粘性抵抗係数
    vForce -= ci * speed * vVelocity;//ciは速度の2乗に比例する慣性抵抗係数

    //加速度
    vAccel = vForce / mass;
    //速度
    vVelocity += vAccel * dt;//床面を滑らしても-z方向の速度が生じる
    //位置
    vPos += vVelocity * dt ;
    //回転角度を更新し回転後のオイラー角を求める
    ang = omega0 * dt / pp; //deg
    if(fabs(vEuler.y) == 90.0 && ang < 0.1) ang = 0.1f;//±90度におけるlock防止
    q = makeQFromAxis( ang , vAxis); //Quaternionを作成
    vEuler = makeEulerFromEuler(vEuler, q);//回転後のオイラー角


    //現在の角速度
    vOmega = omega0 * vAxis; //rad
    //床面との衝突判定（最低点を返す)
    minH = checkCollisionToFloor();//衝突点などもこのﾙｰﾁﾝで求めている
    if(minH < 0.0) vPos.z += fabs(minH) ;//めり込み，沈み込み防止
    if( minH <= 0.0 && vVelocity.z <= 0.0)  //Floorに衝突
    {
        if(kind == SPHERE ){
            //新重心速度のｚ成分（再度このこの速度で放物運動）
            vVelocity.z = - e * vVelocity.z;
        }
        else if(kind == CUBE || kind == CYLINDER){
            //まず摩擦なしで考える
            vNormal = CVector(0.0, 0.0, -1.0);//衝突面の法線方向(target側からの法線）
            //重心から衝突点までのベクトル（checkCollisionToFloor()ルーチンで求めている)
            gtp = vGravityToPoint;
            //衝突点の合成速度
            vp = vVelocity + (vOmega ^ gtp); //衝突点速度(ベクトル）

            //力積
            rikiseki = - (e + 1.0f) * (vNormal * vp) / (1.0f/mass + vNormal*((mInertia.inverse()*(gtp^vNormal))^gtp));
            //新角速度
            vOmega += mInertiaInverse * (gtp ^ vNormal) * rikiseki ;
            //新重心速度
            vVelocity += (vNormal * rikiseki) / mass ;
        }
        else{
            MessageBox(NULL,"この形状物体には未対応","projection",MB_OK);
        }
        //床面摩擦がある場合
        calcFrictionEffect(e, mu, dt);
        //床面摩擦による減速トルク
        vTorque.z = getTorqueOnFloor(mu);
	}
    //回転摩擦トルクによる角加速度
    vAlpha = mInertiaInverse * vTorque;
//    vAlpha = mInertiaInverse *( vTorque - (vOmega ^ (mInertia * vOmega)));
    //角速度の更新
    vOmega += vAlpha * dt;
    //新角速度の大きさと回転軸
    omega0 = vOmega.magnitude();
    vAxis = vOmega.normalize2();
}

//----------------------------------------------------------------
//衝突のとき摩擦があれば回転速度が変化する
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

        //現在のomegaで滑らずに転がるときの水平速度成分
        vhc = rr * omega ;
        //臨界摩擦係数
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
        else{//vVelocity.z =0ならば常にこの状態
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
        //直方体を半径rrの球で近似
        rr = sqrt(vSize.x*vSize.x + vSize.y*vSize.y + vSize.z*vSize.z) / 2.0f;

        //現在のomegaで滑らずに転がるときの水平速度成分
        vhc = rr * omega; //直方体を半径rrの球で近似
        //臨界摩擦係数
//char buf[30];
//if(vVelocity.z ==0.0){
//wsprintf(buf , "%d, %d, %d",(int)(vh*1000.0),(int)(omega*1000.0),(int)(vVelocity.z*1000.0));
//MessageBox(NULL,buf ,"aaa",MB_OK);  }
        if(fabs(vVelocity.z) <= 0.005) {
            //転がらずに滑る
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
            else{ //vVelocity.z =0ならば常にこの状態
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
    else//円柱
    {
        rr = vSize.x / 2.0f;

        //現在のomegaで滑らずに転がるときの水平速度成分
        vhc = rr * omega ;
        //臨界摩擦係数
//        if(fabs(vVelocity.z) <= 0.01) muc = 1000.0;
        if(fabs(vVelocity.z) <= 0.005) {//muc = 1000.0;
            //転がらずに滑る
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
//Floorとの衝突を調べFloorから最低位置の高さを返す
//その値が負であれば衝突
double CTarget::checkCollisionToFloor()
{
    double minH;
    int i, cnt, numVertex;// maxNum;
//    int faceNo, numCollision[6], collisionP[8][4];//衝突面番号，衝突点個数，衝突点番号
    double eps = 0.001f;//1mm以下なら衝突とする
    double ss;//球であれば半径，直方体であれば対角線の半分
    CVector vPointCollision;
    //面の頂点番号(CUBE)
/*    int nFace[6][4] = { {0, 1, 2, 3}, {0, 3, 7, 4}, {0, 4, 5, 1},
                       {1, 5, 6, 2}, {2, 6, 7, 3}, {4, 7, 6, 5} }; */
    CVector vCollision[20];//衝突点の頂点番号(円柱)
//char buf[40];
    //球体とみなしたときの半径
    if(kind == CUBE ) {//直方体
        ss = sqrt(vSize.x * vSize.x + vSize.y * vSize.y + vSize.z * vSize.z) / 2.0f;
    }
    else if(kind == SPHERE){ //球
        ss = vSize.x / 2.0f;
    }
    else{//円柱
        ss = sqrt(vSize.x*vSize.x + vSize.z*vSize.z);
    }

    minH = 1000.0;
    if( vPos.z <= ss)
    {//衝突の可能性あり
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
            else//円柱
            {
                getVertexOfCylinder();
                numVertex = 2 * num1;
            }
            cnt = 0;//衝突点個数
            minH = 1000.0;
            vPointCollision = CVector(0.0, 0.0, 0.0);
            for(i = 0 ; i < numVertex; i++) {
                if(vP[i].z <= eps) { vPointCollision += vP[i]; cnt++; }
                if(vP[i].z < minH )  minH = vP[i].z;
            }

            if(minH <= 0.0)  //Floorと衝突
            {
                vPointCollision /= (double)cnt;
                vGravityToPoint = vPointCollision - vPos;
                return minH;
            }
        }
    }
    return minH;
//    return 1000.0;//衝突なし
}
//------------------------------------------------------------------------------
void CTarget::calcInertia()
{
    double Ixx, Iyy, Izz;
    if(kind == SPHERE)
    {
        Ixx = Iyy = Izz = vSize.x * vSize.y * mass / 10.0f; //vSizeは直径
    }
    else if(kind == CYLINDER)//
    {
        Ixx = vSize.x * vSize.y * mass / 16.0f + vSize.z * vSize.z * mass / 12.0f;
        Iyy = Ixx;
        Izz = vSize.z * vSize.z * mass / 8.0f;
    }
    else if(kind == CUBE)//直方体(角柱)
    {
        Ixx = mass * (vSize.y*vSize.y + vSize.z*vSize.z) / 12.0f;
        Iyy = mass * (vSize.x*vSize.x + vSize.z*vSize.z) / 12.0f;
        Izz = mass * (vSize.x*vSize.x + vSize.y*vSize.y) / 12.0f;
    }
    //慣性モーメント のテンソル
    mInertia = CMatrix(Ixx, 0.0, 0.0,
                       0.0, Iyy, 0.0,
                       0.0, 0.0, Izz );
    mInertiaInverse = mInertia.inverse();

}


//-------------------------------------------------------------------
//把握するときなど，手の方向を決めるときに用いる
void CTarget::getDirVector(TargetDir targetDir, CVector& p)
{
	if(targetDir == TD_FORE)       p = CVector(1.0, 0.0, 0.0);
	else if(targetDir == TD_BACK)  p = CVector(-1.0, 0.0, 0.0);
	else if(targetDir == TD_RIGHT) p = CVector(0.0, -1.0, 0.0);
	else if(targetDir == TD_LEFT)  p = CVector(0.0, 1.0, 0.0);
	else if(targetDir == TD_UP)    p = CVector(0.0, 0.0, 1.0);
	else if(targetDir == TD_DOWN)  p = CVector(0.0, 0.0, -1.0);//{ vP[0] = 0.0; vP[1] = 0.0; vP[2] = -1.0; }
	else {MessageBox(NULL," この方向には対応していない","getDirVector",MB_OK);return;}

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

    //World座標に変換
    for(i = 0; i < 8; i++){
        //回転
        vP[i] = rotate(vP[i], vEuler);
        //中心座標だけ平行移動
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
		vP[i].x = (0.5*cos(theta)*vSize.x); //上底のx成分
		vP[i].y = (0.5*sin(theta)*vSize.y); //ｙ成分
		vP[i].z = 0.5f * vSize.z;                   //ｚ成分(高さ)
		vP[i+num1].x = vP[i].x;                     //下底のx成分
		vP[i+num1].y = vP[i].y;                     //ｙ成分
		vP[i+num1].z = - 0.5f * vSize.z;             //ｚ成分
	}
    //上底の中心
    vP[2*num1].x = 0.0;
    vP[2*num1].y = 0.0;
    vP[2*num1].z = 0.5f * vSize.z;
    //下底の中心
    vP[2*num1+1].x = 0.0;
    vP[2*num1+1].y = 0.0;
    vP[2*num1+1].z = -0.5f * vSize.z;

    //World座標に変換
    for(i = 0; i < 2 * num1 + 2; i++){
        //回転
        vP[i] = rotate(vP[i], vEuler);
        //中心座標だけ平行移動
        vP[i] += vPos;
    }
}

//------------------------------------------------------------------------------
//直方体と直方体の衝突判定（２次元)
//直方体は基本姿勢（z軸が上方向)であること
//注目剛体の頂点が対象剛体の内部（境界上を含む)にあるときだけ衝突とする
bool CTarget::collisionCubeToCube2(CTarget& trgt, CVector& vNormal)
{
    int i, j, ii, jj;
    double a[4], b[4], c[4];
    double len[4];//辺の長さ
    double dist[4];//頂点から辺までの距離
    double f;//判定式
    double min;
    int minNo, VertexNo, kaisu;
    CVector vCollision;//衝突点
    CVector vn;

    getVertexOfCube();//注目剛体(#1)の頂点の位置を取得

    trgt.getVertexOfCube();//衝突対象(#2)の頂点の位置を取得

    //#1の頂点と#2の位置関係
    for(j = 0; j < 4; j++)//jは#2の頂点番号
    {
        if(j == 0 || j == 2) len[j] = trgt.vSize.x;
        else                 len[j] = trgt.vSize.y;
        jj = j + 1;
        if(jj == 4) jj = 0;
        a[j] = trgt.vP[j].y - trgt.vP[jj].y;
        b[j] = trgt.vP[jj].x - trgt.vP[j].x;
        c[j] = trgt.vP[j].x * trgt.vP[jj].y - trgt.vP[jj].x * trgt.vP[j].y;
    }
    for(i = 0; i < 4; i++) //#1の頂点
    {
        kaisu = 0;//判定式が正となる回数
        for(j = 0; j < 4; j++) //#2の辺番号
        {
            f = a[j] * vP[i].x + b[j] * vP[i].y + c[j];
            if( f < 0.0 ) break;  //
            //fが全て正のとき衝突
            dist[j] = fabs(f) / len[j];//辺までの距離
            kaisu ++;
        }
        if( kaisu == 4) //#1の頂点が#2の内部
        {
            VertexNo = i;//#2に衝突している#1の頂点番号
            //もう１つの頂点が衝突しているか調査
            ii = i+1;
            if(ii == 4) ii = 0;
            kaisu = 0;//判定式が正となる回数
            for(j = 0; j < 4; j++)
            { //#2の辺番号j
                f = a[j] * vP[ii].x + b[j] * vP[ii].y + c[j];
                if( f < 0.0 ) break;
                kaisu ++;
            }
            if(kaisu == 4) goto calcNormalVector2;
            kaisu = 0;
            ii = i - 1;
            if(ii == -1) ii = 3;
            for(j = 0; j < 4; j++)
            { //#2の辺番号
                f = a[j] * vP[ii].x + b[j] * vP[ii].y + c[j];
                if( f < 0.0 ) goto calcNormalVector1;
                kaisu ++;
            }
            if(kaisu == 4) goto calcNormalVector2;
        }
    }

    return false;//未衝突

calcNormalVector1:;
    //最も近い辺番号minNoを求めその法線方向vNormalを決定
    //#1の頂点が#2の辺に衝突
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
    //注目剛体から見た法線方向（trgtの面方向の反対符号）
    if(minNo == 0)      vNormal = CVector(0.0, -1.0, 0.0);//基本姿勢の法線方向
    else if(minNo == 1) vNormal = CVector(1.0, 0.0, 0.0);
    else if(minNo == 2) vNormal = CVector(0.0, 1.0, 0.0);
    else                vNormal = CVector(-1.0, 0.0, 0.0);
    vNormal.rotZ(trgt.vEuler.z);//この問題ではｚ軸回転だけ  */
    //衝突点は#1の頂点番号VertexNoとVertexNo+4の中点(0~3は直方体の上底）
    vCollision = (vP[VertexNo] + vP[VertexNo+4]) / 2.0;
    //重心から衝突点までのベクトル
    vGravityToPoint = vCollision - vPos ;
    trgt.vGravityToPoint = vCollision - trgt.vPos ;
    return true;

calcNormalVector2:;
    vCollision = ( (vP[i] + vP[i+4]) + (vP[ii] + vP[ii+4]) )/ 4.0;
    vNormal = vPos >> vCollision;
    //重心から衝突点までのベクトル
    vGravityToPoint = vCollision - vPos ;
    trgt.vGravityToPoint = vCollision - trgt.vPos ;
    return true;
}
//------------------------------------------------------------------------------
//球（円柱)と直方体の衝突判定（２次元)
//次の条件が全て成立するとき衝突とする
//球(円柱)の中心から辺までの距離 <= 球(円柱)の半径
//２つの剛体の中心が直方体の辺の反対側に存在
//球(円柱)の中心が直方体の辺に対して負領域が１つ
bool CTarget::collisionSphereToCube2(CTarget& trgt,CVector& vNormal)
{
    int j, jj;
    double a[4], b[4], c[4];
    double len[4];//辺の長さ
    double dist[4];//頂点から辺までの距離
    double f, f2;//判定式
    int kaisu;
    double rr =vSize.x / 2.0f;//球または円柱の半径
    CVector vCollision;//衝突点
    trgt.getVertexOfCube();//衝突対象(#2)の頂点の位置を取得

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
    {   //#1の中心から#2の辺までの距離
        f = a[j] * vPos.x + b[j] * vPos.y + c[j];
        dist[j] = fabs(f) / len[j];//辺までの距離
        if(dist[j] < rr && kaisu == 1)//衝突の可能性あり
        {
//wsprintf(buf , "%d, %d, %d, %d",j,(int)(dist[j]*1000.0),(int)(f*1000.0),(int)(f2*1000.0));
//MessageBox(NULL,buf ,"collision5",MB_OK);
            f2 = a[j] * trgt.vPos.x + b[j] * trgt.vPos.y + c[j];
            if(f * f2 < 0.0)//2つの中心が辺を挟んで反対側にあるとき衝突
            {
                if(j == 0)      vNormal = CVector( 0.0, -1.0, 0.0);//#2の内側が#1から見て法線方向
                else if(j == 1) vNormal = CVector( 1.0,  0.0, 0.0);
                else if(j == 2) vNormal = CVector( 0.0,  1.0, 0.0);
                else            vNormal = CVector(-1.0,  0.0, 0.0);
                vNormal.rotZ(trgt.vEuler.z);
                //衝突点
                vCollision.x = vPos.x + (vSize.x/2.0f) * (vNormal.x);//
                vCollision.y = vPos.y + (vSize.y/2.0f) * (vNormal.y);//
                vCollision.z = vPos.z;
                //重心から衝突点までのベクトル
                vGravityToPoint = vCollision - vPos ;
                trgt.vGravityToPoint = vCollision - trgt.vPos ;
                return true;
            }
        }
    }
    return false;
}
//------------------------------------------------------------------------------
//直方体と球(円柱)の衝突判定（２次元)
//直方体の頂点が円柱の内部に存在するとき衝突とする．
bool CTarget::collisionCubeToSphere2(CTarget& trgt, CVector& vNormal)
{
    int i;
//    double len[4];//辺の長さ
    double dist[4];//頂点から辺までの距離
//    double f2;//判定式
    double x, y;
    CVector vCollision;//衝突点

    getVertexOfCube();//頂点の位置を取得(vP[j]に格納
    for(i = 0; i < 4; i++)
    {
        x = vP[i].x - trgt.vPos.x;
        y = vP[i].y - trgt.vPos.y;
        dist[i] = sqrt( x * x + y * y); //球，円柱の中心から直方体の頂点までの距離
        if(dist[i] < trgt.vSize.x / 2.0)
        {//衝突            vP[i]は直方体の上辺，vP[i+4]はその下辺
            vNormal = (vP[i]+vP[i+4])/2.0 >> trgt.vPos;//法線方向
            //重心から衝突点までのベクトル
            vGravityToPoint = (vP[i]+vP[i+4])/2.0 - vPos; //直方体(中心から衝突点へ向かうベクトル)
            trgt.vGravityToPoint = (vP[i]+vP[i+4])/2.0 - trgt.vPos;
            return true;
        }
    }
    return false;
}

//------------------------------------------------------------------------------
//直方体と直方体の衝突判定（３次元)
//注目剛体の頂点が対象剛体の内部（境界上を含む)にあるとき衝突
//辺と辺が交差するとき衝突
bool CTarget::collisionCubeToCube3(CTarget& trgt, CVector& vNormal)
{
    //面を構成する頂点番号
    int vs[6][4] = { {0,1,2,3}, {0,3,7,4}, {0,4,5,1},
                     {1,5,6,2}, {2,6,7,3}, {4,7,6,5} };

    int i, j, k;
//    double d[6];
    CVector nfvP[4];//交点と辺で作る面の法線ベクトル
    double f;//判定式
    double min, dd;
    int minNo, VertexNo[8], kaisu, cnt;
    CVector vCollision;//衝突点
    CVector vPoint;
//    CVector vn, aa, bb;
    CVector ps;//面との交点

    getVertexOfCube();//注目剛体(#1)の頂点の位置を取得
    trgt.getVertexOfCube();//衝突対象(#2)の頂点の位置を取得

    //#2の面の法線ベクトル
    for(j = 0; j < 6; j++)//jは#2の面番号
    {
        trgt.vNormalFacet[j] = (trgt.vP[vs[j][1]] - trgt.vP[vs[j][0]]) ^ (trgt.vP[vs[j][2]] - trgt.vP[vs[j][1]]) ;
        trgt.vNormalFacet[j].normalize();
    }

    //注目角柱の全ての頂点について対象剛体の内部にあるかどうかを調査
    cnt = 0;//対象剛体内部にある注目剛体の頂点個数
    for(i = 0; i < 8; i++) //#1の頂点
    {
        kaisu = 0;//判定式が負となる回数
        for(j = 0; j < 6; j++) //#2の面番号
        {
            f = trgt.vNormalFacet[j] * (vP[i] - trgt.vP[vs[j][0]]);
            if( f > 0.001 ) break;  //
            //fが全て負のとき衝突
//            dist[i][j] = fabs(f) / trgt.vNormalFacet[j].magnitude() ;//面までの距離
            kaisu ++;
        }
        if( kaisu == 6) //#1の頂点が#2の内部
        {
            VertexNo[cnt] = i;//#2に衝突している#1の頂点番号
            cnt++;
        }
    }
    if(cnt == 0) goto EdgeAndEdge;

    //#2の面に衝突している#1の頂点をすべて求め平均値を衝突点とする
//    vNormal = CVector(0.0, 0.0, 0.0);//法線ベクトルクリア
    vCollision = CVector(0.0, 0.0, 0.0);//衝突点のクリア
    for(k = 0; k < cnt; k++) 
    {
        vCollision += vP[VertexNo[k]];//衝突点を追加
    }
    //衝突点
    vCollision /= (double)cnt;//平均値
    //辺と辺が交差していればその部分の交点を追加し平均
    if(getPointCubeToCube(trgt, vNormal, vPoint) == true)
    {
        vCollision += vPoint;
        vCollision /= 2.0;
    }

    //最も近い面番号minNoを決定
    f = trgt.vNormalFacet[0] * (vCollision - trgt.vP[vs[0][0]]);
    min = fabs(f) / trgt.vNormalFacet[0].magnitude() ;//面までの距離
    minNo = 0;
    for(j = 1; j < 6; j++)//jは対象剛体の面番号
    {
        //if(fabs(vSpeed * trgt.vNormalFacet[j]) < 0.5) continue;//速度に平行な面を除く
        f = trgt.vNormalFacet[j] * (vCollision - trgt.vP[vs[j][0]]);
        dd = fabs(f) / trgt.vNormalFacet[j].magnitude() ;//面までの距離
        if( min > dd )
        {
            min = dd;
            minNo = j;
        }
    }
    //その面の法線ベクトルを反転
    vNormal = (-1.0) * trgt.vNormalFacet[minNo];

    //重心から衝突点までのベクトル
    vGravityToPoint = vCollision - vPos ; //注目剛体側（vPosは中心座標)
    trgt.vGravityToPoint = vCollision - trgt.vPos ;//対象剛体側
    return true;

EdgeAndEdge:;//辺と辺の交差
    if(getPointCubeToCube(trgt, vNormal, vCollision) == false) return false;

    if(vNormal.x == 0.0 && vNormal.y == 0.0 && vNormal.z == 0.0 )
    {   //上下（または左右）の面と交差
        //最も近い面番号minNoを決定
        f = trgt.vNormalFacet[0] * (vCollision - trgt.vP[vs[0][0]]);
        min = fabs(f) / trgt.vNormalFacet[0].magnitude() ;//面までの距離
        minNo = 0;
        for(j = 1; j < 6; j++)//jは対象剛体の面番号
        {
            f = trgt.vNormalFacet[j] * (vCollision - trgt.vP[vs[j][0]]);
            dd = fabs(f) / trgt.vNormalFacet[j].magnitude() ;//面までの距離
            if( min > dd )
            {
                min = dd;
                minNo = j;
            }
        }
        vNormal = trgt.vNormalFacet[minNo];
    }
    vNormal.normalize();//隣り合う面のときは和の正規化で近似
    //その面の法線ベクトルを反転
    vNormal.reverse();

    //重心から衝突点までのベクトル
    vGravityToPoint = vCollision - vPos ; //注目剛体側（vPosは中心座標)
    trgt.vGravityToPoint = vCollision - trgt.vPos ;//対象剛体側

    return true;
}
//---------------------------------------------------------------------
//注目剛体の辺が対象剛体の面と交差している点を加え平均をとり衝突点とする
//辺の両端が１つの面の外側であればその辺は交差しない
//面と面が十字状にクロスしているときも判定
bool CTarget::getPointCubeToCube(CTarget& trgt, CVector& vNormal, CVector& vPoint)
{
    //面を構成する頂点番号
    int vs[6][4] = { {0,1,2,3}, {0,3,7,4}, {0,4,5,1},
                     {1,5,6,2}, {2,6,7,3}, {4,7,6,5} };
    //辺を構成する頂点番号
    int ve[12][2] = { {0,1}, {1,2}, {2,3}, {3,0},
                      {0,4}, {1,5}, {2,6}, {3,7},
                      {4,5}, {5,6}, {6,7}, {7,4} };

    int i, j, k, kp, kaisu, cnt;
    double fa, fb, tt;
    CVector vNormal0[4];//交点と辺で作る面の法線ベクトル
    CVector vPoint0;    //交点

    kaisu = 0; //注目剛体の辺が対象剛体の面と交差する回数
    vPoint = CVector(0.0, 0.0, 0.0); //交差点
    vNormal = CVector(0.0, 0.0, 0.0); //面の法線方向の和
    for(i = 0; i < 12; i++) //注目角柱の辺
    {
        for(j = 0; j < 6; j++)//対象角柱の面
        {   //辺の頂点が対象剛体の面の正領域か負領域か
            fa = trgt.vNormalFacet[j] * (vP[ve[i][0]] - trgt.vP[vs[j][0]]);
            fb = trgt.vNormalFacet[j] * (vP[ve[i][1]] - trgt.vP[vs[j][0]]);
            if(fa * fb >= 0.0) continue;//同じ領域にあれば交差しない
            tt = fa / (fa - fb);
            vPoint0 = (vP[ve[i][1]] - vP[ve[i][0]]) * tt + vP[ve[i][0]];//面を含む平面との交点
            cnt = 0;
            for(k = 0; k < 4; k++)//交点を持つ面の辺
            {
                kp = k+1;
                if(kp == 4) kp = 0;
                vNormal0[k] = (trgt.vP[vs[j][k]] - vPoint0) ^ (trgt.vP[vs[j][kp]] - trgt.vP[vs[j][k]]) ;
                if(trgt.vNormalFacet[j] * vNormal0[k] < 0.0) break;//１つでも負ならばこの面とは交差しない
                cnt++;
            }
            if(cnt == 4){//交差
                kaisu++;
                vPoint += vPoint0;
                vNormal += trgt.vNormalFacet[j];
            }
        }
    }

    if(kaisu != 2 && kaisu != 4) return false;//交差なし
    //kaisu=4は面と面がクロス
    vPoint /= (double)kaisu;
    return true;
}

//------------------------------------------------------------------------------
//球と直方体の衝突判定（３次元)
//球の中心から辺までの距離 <= 球の半径でかつ
//球の中心が直方体の辺に対して正領域が１つの場合
//直方体の辺と交差するときも衝突（正領域が2つ）
bool CTarget::collisionSphereToCube3(CTarget& trgt,CVector& vNormal)
{
    //面を構成する頂点番号
    int vs[6][4] = { {0,1,2,3}, {0,3,7,4}, {0,4,5,1},
                     {1,5,6,2}, {2,6,7,3}, {4,7,6,5} };
    //辺を構成する頂点番号(最初の2個）と面番号
    int ve[12][4] = { {0,1,0,2}, {1,2,0,3}, {2,3,0,4}, {3,0,0,1},
                      {0,4,1,2}, {1,5,2,3}, {2,6,3,4}, {3,7,1,4},
                      {4,5,2,5}, {5,6,3,5}, {6,7,4,5}, {7,4,1,5} };
    int j, k, no1, no2;//, minJ;
    double dist[6];//球の中心から面までの距離
    double f, d;//判定式
    int faceNo[6], cnt   ;
    double rr =vSize.x / 2.0f;//球の半径
    CVector vCollision;//衝突点
    CVector vDir;

    trgt.getVertexOfCube();//衝突対象(#2)の頂点の位置を取得
    //対象角柱の面の法線ベクトル
    for(j = 0; j < 6; j++)//jは#2の面番号
    {
        trgt.vNormalFacet[j] = (trgt.vP[vs[j][1]] - trgt.vP[vs[j][0]]) ^ (trgt.vP[vs[j][2]] - trgt.vP[vs[j][1]]) ;
        trgt.vNormalFacet[j].normalize();
    }
    //球の中心から対象角柱の面までの距離を調査
    cnt = 0;//球の中心が角柱の面の正領域にある回数
    for(j = 0; j < 6; j++) //対象角柱の面番号
    {
        f = trgt.vNormalFacet[j] * (vPos - trgt.vP[vs[j][0]]);
        if( f >= 0.0 ) {
            faceNo[cnt] = j;
            dist[faceNo[cnt]] = f;// / trgt.vNormalFacet[j].magnitude() ;//面までの距離
            cnt++;
            vCollision = vPos - (trgt.vNormalFacet[j] * f); //cnt=1のときの交点候補
        }
    }
    if(cnt == 1){ //面と衝突の可能性あり
        if(dist[faceNo[0]] > rr) return false;
        vNormal = (-1.0) * trgt.vNormalFacet[faceNo[0]] ;
        //重心から衝突点までのベクトル
        vGravityToPoint = vCollision - vPos ; //注目剛体側（vPosは中心座標)
        trgt.vGravityToPoint = vCollision - trgt.vPos ;//対象剛体側
        return true;
    }
    else if(cnt == 2) { //辺と交差の可能性あり
        //面を共有する辺番号
        for(k = 0; k < 12; k++){
            if( faceNo[0] == ve[k][2] && faceNo[1] == ve[k][3] ){
                no1 = ve[k][0]; no2 = ve[k][1]; //交差する辺の頂点番号
                break;
            }
        }
        //辺の方向ベクトル
        vDir = trgt.vP[no1] >> trgt.vP[no2]; //
        f = vDir * (trgt.vP[no1] - vPos) ;// (vP[no1] - trgt.vPos).magnitude();
        vCollision = trgt.vP[no1] - vDir * f;//球の中心から辺へ垂線を下ろした時の交点
        d = (vCollision - vPos).magnitude();

        if(d > rr) return false;
        
        vNormal = vPos >> vCollision ;//球の中心から衝突点への単位ベクトル
        //重心から衝突点までのベクトル
        vGravityToPoint = vCollision - vPos ; //注目剛体側（vPosは球の中心座標)
        trgt.vGravityToPoint = vCollision - trgt.vPos ;//対象剛体側
        return true;
    }
    else return false;
}
//------------------------------------------------------------------------------
//直方体(角柱)と球の衝突判定（３次元)
//直方体の頂点が球の内部に存在するとき衝突
bool CTarget::collisionCubeToSphere3(CTarget& trgt, CVector& vNormal)
{
    int i;
    double dist;
    int cnt;
    CVector vCollision;//衝突点
//    CVector vDist;

    getVertexOfCube();//頂点の位置を取得(vP[i]に格納)
    vCollision = CVector(0.0, 0.0, 0.0);
    cnt = 0;
    for(i = 0; i < 8; i++)
    {
        //vDist = trgt.vPos - vP[i];
        dist = (trgt.vPos - vP[i]).magnitude(); //球中心から直方体の頂点までの距離
        if(dist < trgt.vSize.x / 2.0) //球の半径以下なら
        {//衝突
            cnt++;
            vCollision += vP[i];
        }
    }
    if(cnt == 0) return false;

    vCollision /= (double)cnt;
    vNormal = vCollision >> trgt.vPos ; //法線方向
    //重心から衝突点までのベクトル
    vGravityToPoint = vCollision - vPos; //直方体(中心から衝突点へ向かうベクトル)
    trgt.vGravityToPoint = vCollision - trgt.vPos;
    return true;
}
//------------------------------------------------------------------------------
//直方体(角柱)と円柱の衝突判定（３次元)
//直方体の頂点が円柱の内部に存在するとき衝突
bool CTarget::collisionCubeToCylinder3(CTarget& trgt, CVector& vNormal)
{
    int i;
    double dist;//頂点から円柱中心軸までの距離
    double h1, h2;//上底，下底までの距離
    double aveH1, aveH2;//その平均値
//    double x, y, min;
    int cnt;
    CVector vCollision;//衝突点
    CVector vCenter;//円柱中心軸方向ベクトル
    CVector vKoten;//直方体の頂点から円柱の中心軸へ下ろした垂線の交点
    CVector vAveKoten;//その平均値

    getVertexOfCube();//直方体の頂点の位置を取得(vP[i]に格納)
    trgt.getVertexOfCylinder();//対象円柱の頂点
    vCenter = trgt.vP[2*num1+1] >> trgt.vP[2*num1] ;//中心軸単位ベクトル
    vCollision = CVector(0.0, 0.0, 0.0);
    cnt = 0;
    vAveKoten = CVector(0.0,0.0,0.0);
    aveH1 = 0.0; aveH2 = 0.0;
    for(i = 0; i < 8; i++)
    {
        vKoten = trgt.vPos-vCenter*(vCenter*(trgt.vPos - vP[i]));
        dist = (vKoten - vP[i]).magnitude(); //直方体の頂点から円柱中心軸までの距離
        if(dist < trgt.vSize.x / 2.0) //円柱の半径以下なら衝突の可能性あり
        {
            h1 = fabs((vCenter * (vP[i] - trgt.vP[2*trgt.num1])) / vCenter.magnitude()); //上底までの距離
            h2 = fabs((vCenter * (vP[i] - trgt.vP[2*trgt.num1+1])) / vCenter.magnitude());//下底までの距離
            if(h1 <= trgt.vSize.z && h2 <= trgt.vSize.z){//衝突
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
    //法線方向を求める
    if(aveH1 < 0.005 && aveH2 > trgt.vSize.z-0.005){//上底と衝突
        vNormal = (-1.0) * vCenter;
    }
    else if(aveH2 < 0.005 && aveH1 > trgt.vSize.z-0.005){//下底と衝突
        vNormal = vCenter;
    }
    else{
        vNormal = vCollision - vAveKoten;
        vNormal.normalize();
    }
    //重心から衝突点までのベクトル
    vGravityToPoint = vCollision - vPos; //直方体(中心から衝突点へ向かうベクトル)
    trgt.vGravityToPoint = vCollision - trgt.vPos;
    return true;
}
//------------------------------------------------------------------------------
//円柱と直方体の衝突判定（３次元)
//円柱の頂点が直方体の内部（境界上を含む)にあるとき衝突とする
//辺と辺の交差も考慮
bool CTarget::collisionCylinderToCube3(CTarget& trgt, CVector& vNormal)
{
    //面を構成する頂点番号（対象直方体）
    int vs[6][4] = { {0,1,2,3}, {0,3,7,4}, {0,4,5,1},
                     {1,5,6,2}, {2,6,7,3}, {4,7,6,5} };
    int i, j, k;
    double dist[40][6];//頂点から面までの距離
    double f;//判定式
    double min, dd;
    int minNo, VertexNo[8], kaisu, cnt, fNo[2], faceNo;
    CVector vCollision;//衝突点
    CVector vPoint;
    CVector vCenter;//円柱の中心軸
    CVector vLine;//剛体の中心間を結ぶ直線

    getVertexOfCylinder();//注目剛体(#1,円柱)の頂点の位置を取得
    trgt.getVertexOfCube();//衝突対象(#2,直方体)の頂点の位置を取得

    //#2の面の法線ベクトル
    for(j = 0; j < 6; j++)//jは#2の面番号
    {
        trgt.vNormalFacet[j] = (trgt.vP[vs[j][1]] - trgt.vP[vs[j][0]]) ^ (trgt.vP[vs[j][2]] - trgt.vP[vs[j][1]]) ;
        trgt.vNormalFacet[j].normalize();
    }

    //注目円柱の全ての頂点について対象剛体の内部にあるかどうかを調査
    cnt = 0;//対象剛体内部にある注目剛体の頂点個数
    for(i = 0; i < 2*num1; i++) //#1の頂点
    {
        kaisu = 0;//判定式が負となる回数
        for(j = 0; j < 6; j++) //#2の面番号
        {
            f = trgt.vNormalFacet[j] * (vP[i] - trgt.vP[vs[j][0]]);
            if( f > 0.001 ) break;  //
            //fが全て負のとき衝突
            dist[i][j] = fabs(f) / trgt.vNormalFacet[j].magnitude() ;//面までの距離
            kaisu ++;
        }
        if( kaisu == 6) //#1の頂点が#2の内部
        {
            VertexNo[cnt] = i;//#2に衝突している#1の頂点番号
            cnt++;
        }
    }
    if(cnt == 0) goto EdgeAndEdge;
    
    vCollision = CVector(0.0, 0.0, 0.0);//衝突点のクリア
    for(k = 0; k < cnt; k++)
    {
        vCollision += vP[VertexNo[k]];//衝突点を追加
    }
    //衝突点
    vCollision /= (double)cnt;//平均値

//wsprintf(buf , "%d",cnt);
//MessageBox(NULL,buf ,"collisionA",MB_OK);

    //直方体の面と円柱の側面が衝突
    //円柱中心軸のベクトル
    vCenter = vP[2*num1+1] >> vP[2*num1];//下底から上底の中心
    kaisu = 0;//円柱の中心軸が直方体の面に平行で面の外側にある回数
    for(j = 0; j < 6; j++)//各面について
    {
        if(fabs(trgt.vNormalFacet[j] * vCenter) < 0.001){//平行(直方体の面の法線と直交）
            dd = trgt.vNormalFacet[j] * (vPos - trgt.vP[vs[j][0]]) / trgt.vNormalFacet[j].magnitude();
            if(dd > 0.0 && dd <= vSize.x / 2.0 + 0.001){
                fNo[kaisu] = j; kaisu ++;
            }
        }
    }

    if(kaisu == 1 || kaisu == 2)
    {
        if(kaisu == 1){//直方体の面と円柱の側面が平行に衝突
            vNormal = (-1.0) * trgt.vNormalFacet[fNo[0]];
        }
        else{//直方体の辺と円柱の側面が平行に衝突(kaisu=2)
            vNormal = (-1.0) * (trgt.vNormalFacet[fNo[0]] + trgt.vNormalFacet[fNo[1]]) / (-2.0);
            vNormal.normalize();
        }
        //交点(2つの中心から面へ降ろした垂線の交点の平均)
        if(getPointCylinderToCube(trgt, faceNo, vPoint) == true)
        {
            vCollision += vPoint;
            vCollision /= 2.0;
        }
        //重心から衝突点までのベクトル
        vGravityToPoint = vCollision - vPos ; //注目剛体側（vPosは中心座標)
        trgt.vGravityToPoint = vCollision - trgt.vPos ;//対象剛体側
//wsprintf(buf , "%d, %d, %d, %d",kaisu, (int)(vNormal.x*1000.0), (int)(vNormal.y*1000.0),(int)(vNormal.z*1000.0));
//MessageBox(NULL,buf ,"collisionA",MB_OK);
        return true;
    }

    min = 1000.0;
    for(j = 0; j < 6; j++)//衝突点から最も近い面番号
    {

        dd = fabs(trgt.vNormalFacet[j] * (vCollision - trgt.vP[vs[j][0]]))/trgt.vNormalFacet[j].magnitude() ;//面までの距離
        if(min > dd){
            min = dd;
            minNo = j;
        }
    }

    //その面の法線ベクトルを反転
    vNormal = -1.0 * trgt.vNormalFacet[minNo];

//    vNormal.normalize();
    //重心から衝突点までのベクトル
    vGravityToPoint = vCollision - vPos ; //注目剛体側（vPosは中心座標)
    trgt.vGravityToPoint = vCollision - trgt.vPos ;//対象剛体側
    return true;

EdgeAndEdge:;//辺と辺の交差
    if(getPointCylinderToCube(trgt, faceNo, vCollision) == false) return false;
    //円柱中心軸のベクトル
    vCenter = vP[2*num1+1] >> vP[2*num1];//下底から上底の中心
    //衝突点から円柱中心軸に下ろした垂線の交点
    CVector vKoten = vP[2*num1+1] - vCenter * ( vCenter * (vP[2*num1+1] - vCollision));
    //円柱から見た衝突点の法線ベクトル
    vNormal = vKoten >> vCollision;
    //円柱の上底または下底が衝突している場合
    //円柱側から見た相対速度
    CVector vSpeed = (vVelocity - trgt.vVelocity).normalize2();

    if(vCenter * vSpeed > 0.5){ //上底側が衝突対象
        if((trgt.vNormalFacet[faceNo] + vCenter).magnitude() <= 0.1) vNormal = vCenter; //上底に平行
    }
    else if(vCenter * vSpeed < -0.5){ //下底側が衝突対象
        if((trgt.vNormalFacet[faceNo] - vCenter).magnitude() <= 0.1) vNormal = (-1.0) * vCenter;//下底に平行
    }
    //重心から衝突点までのベクトル
    vGravityToPoint = vCollision - vPos ; //注目剛体側（vPosは中心座標)
    trgt.vGravityToPoint = vCollision - trgt.vPos ;//対象剛体側

    return true;
}
//---------------------------------------------------------------------
bool CTarget::getPointCylinderToCube(CTarget& trgt, int& faceNo, CVector& vPoint)
{
    //注目剛体のedgeが対象剛体の面と交差している点を加え平均をとり衝突点とする
    //線分の両端が１つの面の外側であればその線分は交差しない
    //すべての線分に対して上の条件が満たされれば，２つの剛体は独立
    //面を構成する頂点番号（対象直方体）
    int vs[6][4] = { {0,1,2,3}, {0,3,7,4}, {0,4,5,1},
                     {1,5,6,2}, {2,6,7,3}, {4,7,6,5} };
    int i, j, k, kp, kaisu, cnt;
    double fa, fb, tt;
    CVector vNormal0[4];//交点と辺で作る面の法線ベクトル
    CVector vPoint0;    //交点

    //辺を構成する頂点番号(注目円柱)
    int ve[20][2];
    for(i = 0; i < num1; i++){
        ve[i][0] = i; ve[i][1] = i + num1;
    }

    kaisu = 0;
    vPoint = CVector(0.0, 0.0, 0.0);
    for(i = 0; i < num1; i++) //注目円柱の線分
    {
        for(j = 0; j < 6; j++)//対象角柱の面
        {
            fa = trgt.vNormalFacet[j] * (vP[ve[i][0]] - trgt.vP[vs[j][0]]);
            fb = trgt.vNormalFacet[j] * (vP[ve[i][1]] - trgt.vP[vs[j][0]]);
            if(fa * fb >= 0.0) continue;//この面とは交差しない
            tt = fa / (fa - fb);
            vPoint0 = (vP[ve[i][1]] - vP[ve[i][0]]) * tt + vP[ve[i][0]];//面を含む平面との交点
            cnt = 0;
            for(k = 0; k < 4; k++)//交点を持つ面の辺
            {
                kp = k+1;
                if(kp == 4) kp = 0;
                vNormal0[k] = (trgt.vP[vs[j][k]] - vPoint0) ^ (trgt.vP[vs[j][kp]] - trgt.vP[vs[j][k]]) ;
                if(trgt.vNormalFacet[j] * vNormal0[k] < 0.0) break;//１つでも負ならばこの面とは交差しない
                cnt++;
            }
            if(cnt == 4){//交差
                kaisu++; //角柱の面と交差している回数
                vPoint += vPoint0;
                faceNo = j;
//wsprintf(buf , "%d, %d, %d, %d, %d",kaisu, i,j,(int)(ps.x*1000.0), (int)(ps.y*1000.0), (int)(ps.z*1000.0));
//MessageBox(NULL,buf ,"collisionA",MB_OK);
            }
        }
    }

    if(kaisu == 0) return false;//２つの剛体は独立

    vPoint /= (double)kaisu;
    return true;
}
//------------------------------------------------------------------
//球と円柱の衝突判定（３次元)
//球の中心から円柱の中心軸までの距離 <= 球の半径＋円柱の半径(側面衝突）
//上底または下底と衝突（球の中心は上底または下底の外部)
//どちらも衝突点が円柱の上底と下底の中間に存在
bool CTarget::collisionSphereToCylinder3(CTarget& trgt,CVector& vNormal)
{
    double dist;//球の中心から円柱中心軸までの距離
    double h1, h2;
    double rr =vSize.x / 2.0f;//球の半径
    CVector vCollision;//衝突点
    CVector vKoten;//球の中心から円柱中心軸へ下ろした垂線の交点

    trgt.getVertexOfCylinder();//対象円柱の頂点
    //円柱中心軸ベクトル
    CVector vCenter = trgt.vP[2*trgt.num1 + 1] >> trgt.vP[2*trgt.num1];
    //球の中心から円柱中心軸へ下ろした垂線の交点
    vKoten = trgt.vP[2*trgt.num1+1]-vCenter*(vCenter*(trgt.vP[2*trgt.num1+1] - vPos));//vPosは球の中心
//char buf[40];
    //球の中心から交点までの距離
    dist = (vKoten - vPos).magnitude();
    //衝突の最低条件
    if(dist > rr + trgt.vSize.x / 2.0) return false;

    //球の中心から上底、下底までの距離
    h1 = fabs((vCenter * (vPos - trgt.vP[2*trgt.num1])));// / vCenter.magnitude()); //上底までの距離
    h2 = fabs((vCenter * (vPos - trgt.vP[2*trgt.num1+1])));// / vCenter.magnitude());//下底までの距離
    //側面衝突
    if(h1 < trgt.vSize.z && h2 < trgt.vSize.z){
        vNormal = vKoten >> vPos;
        vCollision = vKoten + vNormal * (trgt.vSize.x / 2.0f) ;
    }
    //上底または下底と衝突
    else{
        if(h1 <= rr && dist <= trgt.vSize.x / 2.0) {//上底側で衝突
//            vNormal = (-1.0) * vCenter;
            vNormal = vCenter;
            vCollision = vPos - vNormal * (vNormal * (vPos - trgt.vP[2*num1]));
//vPos = vCollision + vNormal * vSize.x / 2.0;
//            vNormal.reverse();
        }
        else if(h2 <= rr && dist <= trgt.vSize.x / 2.0){ //下底側で衝突
            vNormal = (-1.0) * vCenter;
            vCollision = vPos - vNormal * (vNormal * (vPos - trgt.vP[2*num1+1]));
//vPos = vCollision + vNormal * vSize.x / 2.0;
//            vNormal.reverse();
        }
        else return false;
    }
    vNormal.reverse();//注目剛体側から見た法線方向
    //衝突点
//    vCollision = vPos + rr * vNormal;
    //重心から衝突点までのベクトル
    vGravityToPoint = vCollision - vPos ; //注目剛体側（vPosは球の中心座標)
    trgt.vGravityToPoint = vCollision - trgt.vPos ;//対象剛体側
    return true;
}
//------------------------------------------------------------------------------
//円柱と球の衝突判定（３次元)
//円柱の頂点が球の内部に存在するとき衝突
bool CTarget::collisionCylinderToSphere3(CTarget& trgt, CVector& vNormal)
{
    int i, cnt;
    double dist;
    CVector vCollision;//衝突点
//    CVector vDist;

    getVertexOfCylinder();//頂点の位置を取得(vP[i]に格納)
    vCollision = CVector(0.0, 0.0, 0.0);
    cnt = 0;
    for(i = 0; i < 2*num1; i++)
    {
        dist = (trgt.vPos - vP[i]).magnitude(); //球中心から円柱の頂点までの距離
        if(dist < trgt.vSize.x / 2.0) //球の半径以下なら
        {//衝突
            cnt++;
            vCollision += vP[i];
        }
    }
    if(cnt == 0) return false;

    vCollision /= (double)cnt;
    vNormal = vCollision >> trgt.vPos;//法線方向
    //重心から衝突点までのベクトル
    vGravityToPoint = vCollision - vPos; //直方体(中心から衝突点へ向かうベクトル)
    trgt.vGravityToPoint = vCollision - trgt.vPos;
    return true;
}
//------------------------------------------------------------------------------
//円柱同士の衝突判定（３次元)
//円柱の頂点が他方の円柱の内部（境界上を含む)にあるとき衝突
//辺と辺の交差による衝突
bool CTarget::collisionCylinderToCylinder3(CTarget& trgt, CVector& vNormal)
{
//    double l1,m1,n1,l2,m2,n2;//中心軸の方向ベクトル
    double h1, h2, aveH1, aveH2, dist0, dist;
    CVector vCollision;//衝突点
    CVector vKoten;//衝突点から対象円柱の中心軸へ下ろした垂線の交点
//    CVector vNormal0[3];
    CVector vp, vq, vKosaP, vKosaQ;
    double s, t, p2, q2, d;
    int ve[20][2]; //辺を構成する頂点番号(注目円柱)
    int i, cnt;

    getVertexOfCylinder();//注目円柱(#1,円柱)の頂点の位置を取得
    trgt.getVertexOfCylinder();//衝突対象(#2,円柱)の頂点の位置を取得
    CVector vCenter1 = vP[2*num1+1] >> vP[2*num1];//注目円柱の中心軸ベクトル
    CVector vCenter2 = trgt.vP[2*num1+1] >> trgt.vP[2*num1];//対象円柱の中心軸ベクトル
    //中心軸間の最短距離
    if( (vCenter1 - vCenter2).magnitude() == 0.0
        || (vCenter1 + vCenter2).magnitude() == 0.0){  //２直線が平行
        //注目円柱の中心から対象中心軸へ下ろした垂線の交点
//        vKoten = trgt.vP[2*trgt.num1+1]-vCenter2*(vCenter2*(trgt.vP[2*trgt.num1+1] - vPos));
        vKoten = trgt.vPos - vCenter2 * (vCenter2 * (trgt.vPos - vPos));
        dist0 = (vKoten - vPos).magnitude();
    }
    else
        dist0 = fabs((vPos - trgt.vPos) * (vCenter1 ^ vCenter2));

    //衝突最低条件
    if(dist0 > vSize.x/2.0 + trgt.vSize.x/2.0) return false;

    //注目円柱の全ての頂点について対象剛体の内部にあるかどうかを調査
    cnt = 0;//対象剛体内部にある注目剛体の頂点個数
    aveH1 = 0.0; aveH2 = 0.0;
    for(i = 0; i < 2*num1; i++) //#1の頂点
    {
        vKoten = trgt.vP[2*trgt.num1+1]-vCenter2*(vCenter2*(trgt.vP[2*trgt.num1+1] - vP[i]));
        dist = (vKoten - vP[i]).magnitude(); //注目円柱の頂点から対象円柱中心軸までの距離
        if(dist < trgt.vSize.x / 2.0) //円柱の半径以下なら衝突の可能性あり
        {
            h1 = fabs((vCenter2 * (vKoten - trgt.vP[2*trgt.num1]))); //上底までの距離
            h2 = fabs((vCenter2 * (vKoten - trgt.vP[2*trgt.num1+1])));//下底までの距離
            if(h1 <= trgt.vSize.z+0.0001 && h2 <= trgt.vSize.z+0.0001){//衝突
                cnt++;
                aveH1 += h1;
                aveH2 += h2;
                vCollision += vP[i];
            }
        }
    }

    if(cnt == 0) goto EdgeAndEdge;
    //#2の面に衝突している#1の頂点の平均値を衝突点とする
    vCollision /= (double)cnt;
    //中心軸へ下ろした垂線の交点
    vKoten = trgt.vP[2*trgt.num1+1]-vCenter2*(vCenter2*(trgt.vP[2*trgt.num1+1] - vCollision));
    //法線方向を求める
    if(aveH1 < 0.001 && aveH2 > trgt.vSize.z-0.001){//上底と衝突
        vNormal = (-1.0) * vCenter2;
    }
    else if(aveH2 < 0.001 && aveH1 > trgt.vSize.z-0.001){//下底と衝突
        vNormal = vCenter2;
    }
    else{
        vNormal = vCollision - vKoten;
        vNormal.normalize();
    }
    //重心から衝突点までのベクトル
    vGravityToPoint = vCollision - vPos ; //注目剛体側（vPosは中心座標)
    trgt.vGravityToPoint = vCollision - trgt.vPos ;//対象剛体側
    return true;

EdgeAndEdge:;//辺と辺の交差
    //注目円柱のedgeが対象円柱と交差している点を加え平均をとり衝突点とする
    //edgeと対象円柱の中心軸と最短となる位置を求めその点がどちらの円柱に対しても
    //上底と下底の間に存在すれば交差

    //注目円柱のedgeの頂点番号
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
            //MessageBox(NULL,"collisionCylinderToCylinder3ルーチンにおいてd=0" ,"Error",MB_OK);
            return false;
        }
        s = ( q2*((trgt.vP[2*trgt.num1+1]-vP[ve[i][1]])*vp)
            + (vp*vq)*((vP[ve[i][1]]-trgt.vP[2*trgt.num1+1])*vq) ) / d;
        t = ( vp*vq*((trgt.vP[2*trgt.num1+1]-vP[ve[i][1]])*vp)
            + (p2)*((vP[ve[i][1]]-trgt.vP[2*trgt.num1+1])*vq) ) / d;
        //交差点
        vKosaP = vP[ve[i][1]] + s * vp;
        vKosaQ = trgt.vP[2*num1+1] + t * vq;
        if( (vKosaP-vKosaQ).magnitude() > trgt.vSize.x / 2.0) continue;
        //注目円柱側の判定
        h1 = fabs((vp * (vKosaP - vP[ve[i][0]]))); //上底までの距離
        h2 = fabs((vp * (vKosaP - vP[ve[i][1]])));//下底までの距離
        if(h1 > vSize.z || h2 > vSize.z) return false;
        //対象円柱側の判定
        h1 = fabs((vq * (vKosaQ - trgt.vP[2*trgt.num1]))); //上底までの距離
        h2 = fabs((vq * (vKosaQ - trgt.vP[2*trgt.num1+1])));//下底までの距離
        if(h1 > trgt.vSize.z || h2 > trgt.vSize.z) return false;
        vCollision += vKosaP;
        cnt++;
    }
    if(cnt == 0) return false;
    vCollision /= (double)cnt;
    //注目円柱中心軸との交点
    vKoten = vP[2*num1+1]-vCenter1*(vCenter1*(vP[2*num1+1] - vCollision));
    vNormal = vKoten >> vCollision;
    //重心から衝突点までのベクトル
    vGravityToPoint = vCollision - vPos ; //注目剛体側（vPosは中心座標)
    trgt.vGravityToPoint = vCollision - trgt.vPos ;//対象剛体側

    return true;
}
//------------------------------------------------------------------------------
//球同士の衝突(必ず衝突）
bool CTarget::collisionSphereToSphere3(CTarget& trgt, CVector& vNormal)
{
    vNormal = vPos >> trgt.vPos;//球0から球1へ向かう単位法線ベクトル
    //重心から衝突点までのベクトル
    vGravityToPoint = (vSize.x/2.0f) * vNormal ;
    trgt.vGravityToPoint = (-trgt.vSize.x/2.0f) * vNormal ;
    return true;
}
//------------------------------------------------------------------------------
//以下はソリッドテクスチャ用
//--------------------------------------------------------------------------------
//texImageを前もって作成
void CTarget::calcSolidTexCube(CVector vSize)
{
	int i, j, k;
	double x, y, z;
    byte rgb[3];

	//Top
	z = vSize.z / 2.0f;
	for(j = 0; j < T_MAX; j++)
	{
		x = vSize.x * (0.5f - (double)j / T_MAX);//j＝0はtexImageの下部
		for(i = 0; i < T_MAX; i++)
		{
			y = vSize.y * ( -0.5f + (double)i / T_MAX);
			//getImageGrain(i, j, x, y, z, 0);
			getImageGrain(x, y, z, vSize, rgb);
            for(k = 0; k < 3; k++) texImageS[0][j][i][k] = rgb[k];
            //texImageS[0][j][i][3] = 255;//不透明
		}
	}

	//Front(x = 0.5)
	x = vSize.x / 2.0f;
	for(j = 0; j < T_MAX; j++)
	{
		z = vSize.z * (-0.5f + (double)j / T_MAX);//texImageは下から上へ
		for(i = 0; i < T_MAX; i++)
		{
			y = vSize.y * (-0.5f + (double)i / T_MAX);
			//getImageGrain(i, j, x, y, z, 1);
			getImageGrain(x, y, z, vSize, rgb);
            for(k = 0; k < 3; k++) texImageS[1][j][i][k] = rgb[k];
            //texImageS[1][j][i][3] = 255;//不透明
		}
	}

	//Left
	y = vSize.y / 2.0f;
	for(j = 0; j < T_MAX; j++)
	{
		z = vSize.z * (-0.5f + (double)j / T_MAX);//texImageは下から上へ
		for(i = 0; i < T_MAX; i++)
		{
			x = vSize.x * (0.5f - (double)i / T_MAX);
			getImageGrain(x, y, z, vSize, rgb);
            for(k = 0; k < 3; k++) texImageS[2][j][i][k] = rgb[k];
            //texImageS[2][j][i][3] = 255;//不透明
		}
	}

	//Rear
	x = - vSize.x / 2.0f;
	for(j = 0; j < T_MAX; j++)
	{
		z = vSize.z * (-0.5f + (double)j / T_MAX);//texImageは下から上へ
		for(i = 0; i < T_MAX; i++)
		{
			y = vSize.y * (0.5f - (double)i / T_MAX);
			getImageGrain(x, y, z, vSize, rgb);
            for(k = 0; k < 3; k++) texImageS[3][j][i][k] = rgb[k];
            //texImageS[3][j][i][3] = 255;//不透明
		}
	}

	//Right
	y = - vSize.y / 2.0f;
	for(j = 0; j < T_MAX; j++)
	{
		z = vSize.z * (-0.5f + (double)j / T_MAX);//texImageは下から上へ
		for(i = 0; i < T_MAX; i++)
		{
			x = vSize.x * (-0.5f + (double)i / T_MAX);
//			getImageGrain(i, j, x, y, z, 4);
			getImageGrain(x, y, z, vSize, rgb);
            for(k = 0; k < 3; k++) texImageS[4][j][i][k] = rgb[k];
            //texImageS[4][j][i][3] = 255;//不透明
		}
	}

	//Bottom(z = -0.5)
	z = - vSize.z / 2.0f;
	for(j = 0; j < T_MAX; j++)
	{
		x = vSize.x * (-0.5f + (double)j / T_MAX);//j＝0はtexImageの下部
		for(i = 0; i < T_MAX; i++)
		{
			y = vSize.y * (-0.5f + (double)i / T_MAX);
//			getImageGrain(i, j, x, y, z, 5);
			getImageGrain(x, y, z, vSize, rgb);
            for(k = 0; k < 3; k++) texImageS[5][j][i][k] = rgb[k];
            //texImageS[5][j][i][3] = 255;//不透明
		}
	}
}
//---------------------------------------------------------------------------------------
void CTarget::calcSolidTexSphere(CVector vSize)
{
	int i, j, k;
	double phi; //x-y平面に対する偏角（緯度）
	double theta; //x軸に対する偏角（経度)
	double x, y, z, rr, r1, r2;
    byte rgb[3];

	//テクチャーの作成(下から作成)
	for(j = 0; j < T_MAX; j++)
	{
		phi = (M_PI * (double)j / (double)T_MAX);//phiは南極からの緯度
		z = - vSize.z * cos(phi) / 2.0f;//z=0は南極
		rr =   sin(phi) ;
		r1 = vSize.x * rr; r2 = vSize.y * rr;
		for(i = 0; i < T_MAX; i++)
		{
			theta = (2.0f * M_PI * (double)i / (double)T_MAX);//x=0からの経度
			x = r1 * cos(theta) ;
			y = r2 * sin(theta) ;
			getImageGrain(x, y, z, vSize, rgb);
            for(k = 0; k < 3; k++) texImageS[0][j][i][k] = rgb[k];
            //texImageS[0][j][i][3] = 255;//不透明
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

	//Topのﾃｸｽﾁｬｰ
	z = (vSize.z / 2.0);
	for(j = 0; j < T_MAX; j++)//
	{
		x = vSize.x * (0.5f - (double)j / T_MAX);//j＝0はtexImageの下部
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
		x = vSize.x * ( 0.5f - (double)j / T_MAX);//j＝0はtexImageの下部
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
		z = vSize.z * (- 0.5f + (double)j / T_MAX);//j＝0はtexImageの下部
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
	double period = 2.0;//vSize=10を基準
	double dist;
	double a1, a2;
	double x0 = 0.08;
	double y0 = 0.5;//random(100)/100.0;//  0.1;
	double r, g, b;

	//振幅
	a1 = 30.0;
	a2 = 20.0;
	//木目の周期
	period = 0.05f * (1.0f - 0.005f * z);
	//中心からの距離
	dist = sqrt((x-x0) * (x-x0) + (y-y0) * (y-y0));
	//色成分
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

