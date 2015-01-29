//skMath.h

//3次元ベクトル演算
class CVector{

public:
    //メンバ変数
    double x, y, z;
    //コンストラクタ
    CVector();
    CVector(double x0, double y0, double z0);
    //演算子
    void operator+=(CVector a);//加算
    void operator-=(CVector a);//減算
    void operator*=(double s);//スカラ乗算
    void operator/=(double s);//スカラ除算
    friend CVector operator*(CVector a, double s);//スカラ乗算
    friend CVector operator*(double s, CVector a);//スカラ乗算
    friend CVector operator/(CVector a, double s);//スカラ除算
    friend CVector operator+(CVector a, CVector b);//加算
    friend CVector operator-(CVector a, CVector b);//減算
    friend double operator*(CVector a, CVector b);//内積
    friend CVector operator^(CVector a, CVector b);//外積
    friend CVector operator>>(CVector a, CVector b);//aからbへ向かう単位ベクトル
    //メンバ関数
    double magnitude();  //大きさ
    void normalize(void);//大きさ1のベクトルに変換
    void reverse(void);//反転
    CVector normalize2();//正規化したベクトルを返す
    CVector rotate(CVector v, CVector elr);//オイラー角による回転
    void rotX(double angle); //x軸回転
    void rotY(double angle); //y軸回転
    void rotZ(double angle); //z軸回転
    double dist(CVector a, CVector b);//距離
    double getAngle(CVector a, CVector b);//２つのベクトルのなす角度
    CVector getEulerX(CVector a, CVector b);//aからｂへ向かうベクトルのオイラー角
    CVector getEulerZ(CVector a, CVector b);//aからｂへ向かうベクトルのオイラー角
};
//----------------------------------------------------------------------------
//コンストラクタ
CVector::CVector()
{
    x = y = z = 0.0;
}
//----------------------------------------------------------------------------
//コンストラクタ
CVector::CVector(double x0, double y0, double z0)
{
    x = x0;
    y = y0;
    z = z0;
}
//------------------------------------------------------------------------------
void CVector::operator+=(CVector a)
{
    x += a.x;
    y += a.y;
    z += a.z;
}
//------------------------------------------------------------------------------
void CVector::operator-=(CVector a)
{
    x -= a.x;
    y -= a.y;
    z -= a.z;
}
//------------------------------------------------------------------------------
void CVector::operator*=(double s)
{
    x *= s;
    y *= s;
    z *= s;
}
//------------------------------------------------------------------------------
void CVector::operator/=(double s)
{
    x /= s;
    y /= s;
    z /= s;
}
//------------------------------------------------------------------------------
CVector operator+(CVector a, CVector b)
{
    return CVector(a.x+b.x, a.y+b.y, a.z+b.z);
}
//------------------------------------------------------------------------------
CVector operator-(CVector a, CVector b)
{
    return CVector(a.x - b.x, a.y - b.y, a.z - b.z);
}
//------------------------------------------------------------------------------
//スカラ乗算
CVector operator*(CVector a, double s)
{
    return CVector(a.x * s, a.y * s, a.z * s);
}
//------------------------------------------------------------------------------
//スカラ乗算
CVector operator*(double s, CVector a)
{
    return CVector(a.x * s, a.y * s, a.z * s);
}
//------------------------------------------------------------------------------
//スカラ除算
CVector operator/(CVector a, double s)
{
    if(s == 0.0)
        MessageBox(NULL,"スカラ除算の分母が０です！" ,"CVector",MB_OK);

    return CVector(a.x / s, a.y / s, a.z / s);
}

//------------------------------------------------------------------------------
//内積
double operator*(CVector a, CVector b)
{
    return( a.x*b.x + a.y*b.y + a.z*b.z );
}
//------------------------------------------------------------------------------
//外積
CVector operator^(CVector a, CVector b)
{
    double x = a.y * b.z - a.z * b.y;
    double y = a.z * b.x - a.x * b.z;
    double z = a.x * b.y - a.y * b.x;
    return( CVector(x, y, z) );
}
//------------------------------------------------------------------------------
//aからbへ向かう単位ベクトル
CVector operator>>(CVector a, CVector b)
{
    CVector c = b - a;
    c.normalize();
    return c ;
}

//-----------------------------------------------------------------------------
double CVector::magnitude()
{
    return( sqrt(x * x + y * y + z * z) );
}
//-----------------------------------------------------------------------------
void CVector::normalize()
{
    double eps = 0.000001;
    double mag = sqrt(x * x + y * y + z * z) ;
    if(mag <= eps) mag = 1.0;
    x /= mag;
    y /= mag;
    z /= mag;
    if(fabs(x) < eps) x = 0.0;
    if(fabs(y) < eps) y = 0.0;
    if(fabs(z) < eps) z = 0.0;
}
//-----------------------------------------------------------------------------
CVector CVector::normalize2()
{
    double eps = 0.000001;
    double xx, yy, zz;
    double mag = sqrt(x * x + y * y + z * z) ;
    if(mag < eps) mag = 1.0;
    xx = x / mag;
    yy = y / mag;
    zz = z / mag;
    if(fabs(x) < eps) xx = 0.0;
    if(fabs(y) < eps) yy = 0.0;
    if(fabs(z) < eps) zz = 0.0;
    return(CVector(xx, yy, zz));
}

//-----------------------------------------------------------------------------
void CVector::reverse()
{
    x = -x;
    y = -y;
    z = -z;
}

//------------------------------------------------------------------------------------
void CVector::rotX(double angle)
{
	double xx, yy, zz;
	double pp = M_PI / 180.0;

	xx = x; yy = y; zz = z;
	x = xx;
	y = yy * cos(pp * angle) - zz * sin(pp * angle);
	z = yy * sin(pp * angle) + zz * cos(pp * angle);
}
//------------------------------------------------------------------------------------
void CVector::rotY(double angle)
{
	double xx, yy, zz;
	double pp = M_PI / 180.0;

	xx = x; yy = y; zz = z;
	x = xx * cos(pp * angle) + zz * sin(pp * angle);
	y = yy;
	z = - xx * sin(pp * angle) + zz * cos(pp * angle);
}
//------------------------------------------------------------------------------------
void CVector::rotZ(double angle)
{
	double xx, yy, zz;
	double pp = M_PI / 180.0;

	xx = x; yy = y; zz = z;
	x = xx * cos(pp * angle) - yy * sin(pp * angle);
	y = xx * sin(pp * angle) + yy * cos(pp * angle);
	z = zz;
}
//-----------------------------------------------------------------------------
CVector rotate(CVector v, CVector elr)
{
    v.rotX(elr.x);
    v.rotY(elr.y);
    v.rotZ(elr.z);
    return v;
}
//------------------------------------------------------------------------------
//ベクトル間距離
double dist(CVector a, CVector b)
{
    CVector c = b - a;
    return( c.magnitude() );
}

//-----------------------------------------------------------------------------
double getAngle(CVector a, CVector b)
{
    double ang = (a.x*b.x+a.y*b.y+a.z*b.z)/(a.magnitude()*b.magnitude());
    if(ang >= 1.0) ang = 0.0;
    else if (ang <= -1.0) ang = M_PI;
    else ang = acos(ang);
    return ang;//rad単位で返す
}
//-----------------------------------------------------------------------------
CVector getEulerX(CVector a, CVector b)
{//基本姿勢で中心軸がｘ軸方向である棒状オブジェクトのオイラー角
    double cx, cy, cz, len;
    CVector e;
    cx = b.x - a.x;
    cy = b.y - a.y;
    cz = b.z - a.z;
    len = dist(a, b);
    e.x = 0.0;
    if(cz >= len) e.y = -90.0;
    else if(cz <= -len) e.y = 90.0;
    else e.y = - (asin(cz / len) * 180.0 / M_PI);
    if(fabs(cx) <= 0.0001 && fabs(cy) <= 0.0001) e.z = 0.0;
    else e.z = (atan2(cy, cx) * 180.0 / M_PI);
    return e;
}
//------------------------------------------------------------------------------
CVector getEulerZ(CVector a, CVector b)
{//基本姿勢で中心軸がz軸方向である棒状オブジェクトのオイラー角
    double cx, cy, cz, len;
    CVector e;
    cx = b.x - a.x;
    cy = b.y - a.y;
    cz = b.z - a.z;
    len = dist(a, b);
    e.x = 0.0;
    if(cz >= len) e.y = 0.0;
    else if(cz <= -len) e.y = 180.0;
    else e.y = (acos(cz / len) * 180.0 / M_PI);
    if(fabs(cx) <= 0.0001 && fabs(cy) <= 0.0001) e.z = 0.0;
    else e.z = (atan2(cy, cx) * 180.0 / M_PI) ;
    return e;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//3*3のマトリクス
//------------------------------------------------------------------------------
class CMatrix{

public:
    //メンバ変数
    double e11, e12, e13, e21, e22, e23, e31, e32, e33;
    //コンストラクタ
    CMatrix();
    CMatrix(double a11, double a12, double a13,
            double a21, double a22, double a23,
            double a31, double a32, double a33);
    //演算子
    friend CVector operator*(CMatrix m, CVector v);
    friend CVector operator*(CVector v, CMatrix m);
    //メンバ関数
    double det();
    CMatrix inverse();
};
//----------------------------------------------------------------------------
//コンストラクタ
CMatrix::CMatrix()
{
    e11 = e12 = e13 = 0.0;
    e21 = e22 = e23 = 0.0;
    e31 = e32 = e33 = 0.0;
}
//-----------------------------------------------------------------------------
//コンストラクタ
CMatrix::CMatrix(double a11, double a12, double a13,
                 double a21, double a22, double a23,
                 double a31, double a32, double a33)
{
    e11 = a11; e12 = a12; e13 = a13;
    e21 = a21; e22 = a22; e23 = a23;
    e31 = a31; e32 = a32; e33 = a33;
}
//----------------------------------------------------------------------------
//行列式
double CMatrix::det()
{
    return( e11*e22*e33 + e12*e23*e31 + e13*e32*e21
          - e11*e32*e23 - e12*e21*e33 - e31*e22*e13 );
}
//----------------------------------------------------------------------------
CMatrix CMatrix::inverse()
{
    double d = e11*e22*e33 + e12*e23*e31 + e13*e32*e21
             - e11*e32*e23 - e12*e21*e33 - e31*e22*e13 ;
    if( d == 0.0 ) {
        MessageBox(NULL,"逆行列を求めることができません！" ,"CMatrix",MB_OK);
        d = 1.0;
    }
    return CMatrix(
        (e22*e33-e23*e32)/d, -(e12*e33-e13*e32)/d,  (e12*e23-e13*e22)/d,
       -(e21*e33-e23*e31)/d,  (e11*e33-e13*e31)/d, -(e11*e23-e13*e21)/d,
        (e21*e32-e22*e31)/d, -(e11*e32-e12*e31)/d,  (e11*e22-e12*e21)/d );
}
//----------------------------------------------------------------------------
//ベクトル乗算
CVector operator*(CMatrix m, CVector v)
{
    return CVector(
        m.e11*v.x + m.e12*v.y + m.e13*v.z,    //x成分
        m.e21*v.x + m.e22*v.y + m.e23*v.z,    //y成分
        m.e31*v.x + m.e32*v.y + m.e33*v.z );  //z成分
}
//----------------------------------------------------------------------------
//ベクトル乗算
CVector operator*(CVector v, CMatrix m)
{
    return CVector(
        v.x*m.e11 + v.y*m.e21 + v.z*m.e31,
        v.x*m.e12 + v.y*m.e22 + v.z*m.e32,
        v.x*m.e13 + v.y*m.e23 + v.z*m.e33 );
}


//----------------------------------------------------------------------------
//四元数：quaternion
// 定義式 q = s + xi + yj + zk;
//----------------------------------------------------------------------------
class CQuaternion{

public:
    //メンバ変数
    double s;//スカラ値
    CVector v;//ベクトル
    //コンストラクタ
    CQuaternion();
    CQuaternion(double a0, double a1, double a2, double a3);
    //演算子
    CQuaternion operator~();
    void operator+=(CQuaternion q);
    void operator-=(CQuaternion q);
    void operator*=(double a);
    void operator/=(double a);
    friend CQuaternion operator+(CQuaternion p, CQuaternion q);
    friend CQuaternion operator-(CQuaternion p, CQuaternion q);
    friend CQuaternion operator*(CQuaternion p, CQuaternion q);
    friend CQuaternion operator*(CQuaternion p, double s);
    friend CQuaternion operator*(double s, CQuaternion p);
    friend CQuaternion operator*(CVector v, CQuaternion q);
    friend CQuaternion operator*(CQuaternion q, CVector v);
    //メンバ関数
    double magnitude();
    CVector getVector();
    double getScalar();
    CQuaternion qvRotate(CVector v, CQuaternion q);
    CVector qvRotate(CQuaternion q, CVector v);
//    CQuaternion makeQFromEuler(CVector elr);
//    CVector makeEulerFromQ(CQuaternion q);
    CQuaternion makeQFromAxis(double angle, CVector axis);
    CVector makeEulerFromEuler(CVector v, CQuaternion);
};
//-----------------------------------------------------------------------------
//コンストラクタ
CQuaternion::CQuaternion()
{
    s = 0;
    v.x = v.y = v.z = 0.0;
}
//-----------------------------------------------------------------------------
//コンストラクタ
CQuaternion::CQuaternion(double a0, double a1, double a2, double a3)
{
    s = a0;
    v.x = a1;
    v.y = a2;
    v.z = a3;
}
//-----------------------------------------------------------------------------
//共役四元数
CQuaternion CQuaternion::operator~()
{
    return CQuaternion( s, -v.x, -v.y, -v.z);
}
//-----------------------------------------------------------------------------
void CQuaternion::operator+=(CQuaternion q)
{
    s += q.s;
    v.x += q.v.x;
    v.y += q.v.y;
    v.z += q.v.z;
}
//-----------------------------------------------------------------------------
void CQuaternion::operator-=(CQuaternion q)
{
    s -= q.s;
    v.x -= q.v.x;
    v.y -= q.v.y;
    v.z -= q.v.z;
}
//-----------------------------------------------------------------------------
void CQuaternion::operator*=(double a)
{
    s *= a;
    v.x *= a;
    v.y *= a;
    v.z *= a;
}
//-----------------------------------------------------------------------------
void CQuaternion::operator/=(double a)
{
    s /= a;
    v.x /= a;
    v.y /= a;
    v.z /= a;
}

//------------------------------------------------------------------------------
//四元数どうしの和
CQuaternion operator+(CQuaternion p, CQuaternion q)
{
    return CQuaternion( p.s + q.s,
           p.v.x + q.v.x, p.v.y + q.v.y, p.v.z + q.v.z);
}
//------------------------------------------------------------------------------
//四元数どうしの差
CQuaternion operator-(CQuaternion p, CQuaternion q)
{
    return CQuaternion( p.s - q.s,
           p.v.x - q.v.x, p.v.y - q.v.y, p.v.z - q.v.z);
}
//------------------------------------------------------------------------------
//四元数どうしの乗算
CQuaternion operator*(CQuaternion p, CQuaternion q)
{
    CQuaternion qq;
     qq.s = p.s * q.s - p.v * q.v ; //２項目はベクトル内積
     qq.v = p.s * q.v + q.s * p.v + (p.v^q.v); //３項目はベクトル外積
    return qq;
}
//-----------------------------------------------------------------------------
//スカラ乗算
CQuaternion operator*(CQuaternion p, double a)
{
    return CQuaternion(p.s * a, p.v.x * a, p.v.y * a, p.v.z * a);
}
//-----------------------------------------------------------------------------
//スカラ乗算
CQuaternion operator*(double a, CQuaternion p)
{
    return CQuaternion(p.s * a, p.v.x * a, p.v.y * a, p.v.z * a);
}
//------------------------------------------------------------------------------
//ベクトルと四元数の乗算
CQuaternion operator*(CVector v, CQuaternion q)
{
    CQuaternion p, pp;
    p.s = 0.0; p.v = v; //スカラ部が0のクオータニオン
    pp = p * q;
    return pp; 
}
//------------------------------------------------------------------------------
//ベクトルと四元数の乗算
CQuaternion operator*(CQuaternion p, CVector v)
{
    CQuaternion q, pp;
    q.s = 0.0; q.v = v;
    pp = p * q;
    return pp;
}

//-----------------------------------------------------------------------------
//大きさを返す
double CQuaternion::magnitude()
{
    return sqrt( s * s + v.x * v.x + v.y * v.y + v.z * v.z);
}
//-----------------------------------------------------------------------------
//ベクトル部分を返す
CVector CQuaternion::getVector()
{
    return CVector(v.x, v.y, v.z);
}
//-----------------------------------------------------------------------------
//スカラ部を返す
double CQuaternion::getScalar()
{
    return s;
}
//-----------------------------------------------------------------------------
//ベクトルvを４元数qで回転
CVector qvRotate(CQuaternion q, CVector v)
{
    CQuaternion p = q * v * (~q);
    return CVector(p.v.x, p.v.y, p.v.z);
}
/*
//-----------------------------------------------------------------------------
//オイラー角から４元数
CQuaternion makeQFromEuler(CVector elr)
{
    CQuaternion q;
    double cx, cy, cz, sx, sy, sz;

    elr = elr * M_PI / 180.0;
    cx = cos(0.5 * elr.x); sx = sin(0.5 * elr.x); //yaw
    cy = cos(0.5 * elr.y); sy = sin(0.5 * elr.y); //pitch
    cz = cos(0.5 * elr.z); sz = sin(0.5 * elr.z); //roll
    q.s = cx*cy*cz + sx*sy*sz;
    q.v.x = sx*cy*cz - cx*sy*sz;
    q.v.y = cx*sy*cz + sx*cy*sz;
    q.v.z = cx*cy*sz - sx*sy*cz;

    return q;
}
//-----------------------------------------------------------------------------
//４元数からオイラー角
CVector makeEulerFromQ(CQuaternion q)
{   //基本姿勢からq回転
    double eps = 0.999999;
    CVector elr; //返されるオイラー角（ベクトル形式)
    CVector vx = CVector(1.0, 0.0, 0.0);//x軸
    CVector vy = CVector(0.0, 1.0, 0.0);//y軸
    CVector vz = CVector(0.0, 0.0, 1.0);//z軸
    vx = qvRotate(q, vx);
    vy = qvRotate(q, vy);
    vz = qvRotate(q, vz);
    
    if(fabs(vx.z) > eps)//elr.yが±90度に近いときは注意
    {
        elr.x = 0.0;
        elr.y = - (vx.z / fabs(vx.z)) * 90.0;
        if(vy.x > eps) elr.z = -90.0;
        else if( vy.x < -eps) elr.z = 90.0;
        else {
            if(vy.y >= 0.0) elr.z = - asin(vy.x) * 180.0 / M_PI;
            else  elr.z = 180.0 + asin(vy.x) * 180.0 / M_PI;
        }
    }
    else
    {
        //alpha
       if(vy.z == 0.0 && vz.z == 0.0) {
            MessageBox(NULL,"atan2ｴﾗｰ" ,"MakeEulerFromQ",MB_OK);
            elr.x = 0.0;
        }
        else elr.x = atan2(vy.z , vz.z) * 180.0 / M_PI;
        //beta
        if( vx.z > eps) elr.y = -90.0;
        else if( vx.z < -eps) elr.y = 90.0;
        else elr.y = -asin(vx.z) * 180.0 / M_PI;
        //gamma
        elr.z = atan2(vx.y , vx.x) * 180.0 / M_PI;
    }
    return elr;

}    */
//------------------------------------------------------------------------------
//任意の回転軸から４元数を取得(axisはworld座標）
CQuaternion makeQFromAxis(double angle, CVector axis)
{
    axis.normalize();
    double a2 = (angle * M_PI / 360.0);//rad単位でθ/2
    return CQuaternion( cos(a2) ,
         sin(a2) * axis.x, sin(a2) * axis.y, sin(a2) * axis.z );
}

//-----------------------------------------------------------------------------
//基本ベクトルを現在のEuler角で回転し，
//さらにクオータニオンｑによる回転を行い
//回転後の新しいオイラー角を求める
CVector makeEulerFromEuler(CVector elr0, CQuaternion q)
{   //現在姿勢からq回転
    double eps = 0.9999999;
    CVector elr; //返されるオイラー角（ベクトル形式)
    //物体の基本ベクトル
    CVector v1 = CVector(1.0, 0.0, 0.0);//x軸(TD_FORE)
    CVector v2 = CVector(0.0, 1.0, 0.0);//y軸(TD_LEFT)
    CVector v3 = CVector(0.0, 0.0, 1.0);//z軸(TD_UP)
//char buf[20];
    //現在姿勢(elr0)による基本ベクトルの回転
    v1 = rotate(v1, elr0);
    v2 = rotate(v2, elr0);
    v3 = rotate(v3, elr0);
    //クオータニオンによる回転
    v1 = qvRotate(q, v1);
    v2 = qvRotate(q, v2);
    v3 = qvRotate(q, v3);
//wsprintf(buf , "%d, %d, %d",(int)(v2.z*10000000.0),(int)(v3.z*10000000.0),(int)(v1.z*1000.0));
//MessageBox(NULL,buf ,"quaternion1",MB_OK);
    //回転後の基本ベクトルからオイラー角
    if(fabs(v1.z) > eps)//elr.yが±90度に近いときは注意
    {   //x軸回転を0，y軸を90°または-90°として
        //z軸だけを回転すればよい
        elr.x = 0.0;
        elr.y = - (v1.z / fabs(v1.z)) * 90.0;
        if(v2.x > eps) elr.z = -90.0;
        else if( v2.x < -eps) elr.z = 90.0;
        else{
            if(v2.y >= 0.0) elr.z = - (asin(v2.x) * 180.0 / M_PI);
            else  elr.z = (180.0 + asin(v2.x) * 180.0 / M_PI);
        }
    }
    else
    {
        //alpha
       if(v2.z == 0.0 && v3.z == 0.0) {
            MessageBox(NULL,"atan2ｴﾗｰ" ,"MakeEulerFromEuler",MB_OK);
            elr.x = 0.0;
        }
        else elr.x = (atan2(v2.z , v3.z) * 180.0 / M_PI);
        //beta
        if( v1.z > eps) elr.y = -90.0;
        else if( v1.z < -eps) elr.y = 90.0;
        else elr.y = - (asin(v1.z) * 180.0 / M_PI);
        //gamma
        elr.z = (atan2(v1.y , v1.x) * 180.0 / M_PI);
    }
    return elr;
}
//------------------------------------------------------------------------------

//---------------------------------------------------------------------------
//LU分解による連立1次方程式の解法(一般の場合)
//LU分解ルーチン(LU分解し元の係数配列に格納)
void LUDecomposition(int n, double a[][20])
{
    int i, j, k;
    double w;

    //最初の行(Upper)
    for(j = 1; j < n; j++) a[0][j] = a[0][j] / a[0][0];

    for(k = 1; k < n-1; k++)
    {
        //Lower
        for(i = k; i < n; i++)
        {
            w = a[i][k];
            for(j = 0; j < k; j++) w -= a[i][j] * a[j][k];
            a[i][k] = w;
        }
        //Upper
        for(j = k+1; j < n; j++)
        {
            w = a[k][j];
            for(i = 0; i < k; i++) w -= a[k][i] * a[i][j];
            a[k][j] = w / a[k][k];
        }
    }
    //a[n-1][n-1]
    w = a[n-1][n-1];
    for(j = 0; j < n-1; j++) w -= a[n-1][j] * a[j][n-1];
    a[n-1][n-1] = w;
}
//-----------------------------------------------------------------
//LU分解後代入・逆代入により連立1次方程式を解く
void SolutionByLU(int n, double a[][20], double* b)
{
    int i, j;
    double w;

    //代入
    b[0] = b[0] / a[0][0];
    for(i = 1; i < n; i++)
    {
        w = b[i];
        for(j = 0; j < i; j++) w -= a[i][j] * b[j];
        b[i] = w / a[i][i];
    }
    //逆代入
    for(i = n-2; i >= 0; i--)
    {
        w = b[i];
        for(j = i+1; j < n; j++) w -= a[i][j] * b[j];
        b[i] = w;
    }
}
//----------------------------------------------------------------------------
//３次スプライン曲線による補間
void spline3( double p[][3], double q[][3], int numControl, int numDivision)
{
    int i, j, k, ip;
    int n = numControl;
    int m = numControl - 1;//segment個数
    double s[20];//segmentの長さ
    double a[20][20], b[20][3];
    double u, u0;
    double w0[3], w1[3], w2[3], w3[3];
    double bb[20];

//char buf[50];
    //segment長
    //係数マトリクス初期化
    for(j = 0; j < n; j++){
        for(i = 0; i < n; i++) a[i][j] = 0.0;
        b[i][0] = b[i][1] = b[i][2] = 0.0;
    }
    for(i = 0; i < m; i++)
        s[i] = sqrt( (p[i+1][0]-p[i][0]) * (p[i+1][0]-p[i][0])
                   + (p[i+1][1]-p[i][1]) * (p[i+1][1]-p[i][1]) );

    if(p[n-1][0] == p[0][0] && p[n-1][1] == p[0][1]
        && p[n-1][2] == p[0][2]) goto syuuki;//始点と終点が同じ

    //自由境界条件
    a[0][0] = 2.0f * s[0];
    a[0][1] = s[0];
    for(i = 1; i < n-1; i++)
    {
        a[i][i-1] = s[i];
        a[i][i] = 2.0f * (s[i-1] + s[i]);
        a[i][i+1] = s[i-1];
    }
    a[n-1][n-2] = s[n-2];
    a[n-1][n-1] = 2.0f * s[n-2];

    //LU分解
    LUDecomposition(n, a);

    //定数ベクトル
    for(k = 0; k < 3; k++)//k=0:x, k=1:y, k=2:z
    {
        b[0][k] = 3.0f * (p[1][k] - p[0][k]);
        for(i = 1; i < n-1; i++)
        {
            b[i][k] = 3.0f * ((p[i][k]-p[i-1][k]) * s[i] * s[i]
              + (p[i+1][k]-p[i][k]) * s[i-1] * s[i-1]) / (s[i-1]*s[i]);
        }
        b[n-1][k] = 3.0f * (p[n-1][k] - p[n-2][k]);
    }

    //x座標の解
    for(i = 0; i < n; i++) bb[i] = b[i][0];
    SolutionByLU(n, a, bb);
    for(i = 0; i < n; i++) b[i][0] = bb[i];
    //ｙ座標の解
    for(i = 0; i < n; i++) bb[i] = b[i][1];
    SolutionByLU(n, a, bb);
    for(i = 0; i < n; i++) b[i][1] = bb[i];
    //ｚ座標の解
    for(i = 0; i < n; i++) bb[i] = b[i][2];
    SolutionByLU(n, a, bb);
    for(i = 0; i < n; i++) b[i][2] = bb[i];
    //1次微係数はb[i][0],b[i][1],b[i][2]に格納されている
    //補間
    for(k = 0; k < 3; k++)
    {
        for(i = 0; i < m; i++)
        {
            w0[k] = p[i][k];
            w1[k] = b[i][k];//1次微係数
            w2[k] = 3.0f * (p[i+1][k]-p[i][k])/(s[i]*s[i]) - (2.0f*b[i][k]+b[i+1][k])/s[i];
            w3[k] = 2.0f * (p[i][k]-p[i+1][k])/(s[i]*s[i]*s[i]) + (b[i][k]+b[i+1][k])/(s[i]*s[i]);
            u0 = s[i] / numDivision;
            for(j = 0; j < numDivision; j++)
            {
                u = u0 * (double)j;
                q[numDivision * i + j][k] = w0[k] + w1[k]*u + w2[k]*u*u + w3[k]*u*u*u;
            }
        }
    }
    //最終点
    q[m*numDivision][0] = p[n-1][0];
    q[m*numDivision][1] = p[n-1][1];
    q[m*numDivision][2] = p[n-1][2];
    return;

syuuki:;
    //周期境界条件
    a[0][0] = 2.0f * (s[0] + s[n-2]);
    a[0][1] = s[n-2];
    a[0][n-2] = s[0];
    for(i = 1; i < n-1; i++)
    {
        a[i][i-1] = s[i];
        a[i][i] = 2.0f * (s[i-1] + s[i]);
        a[i][i+1] = s[i-1];
    }
//    a[n-2][0] = s[n-3];

    //LU分解
    LUDecomposition(n-1, a);

    //定数ベクトル
    for(k = 0; k < 3; k++)//k=0:x, k=1:y, k=2:z
    {
        b[0][k] = 3.0f * ((p[1][k] - p[0][k])*s[n-2]*s[n-2]
                        - (p[n-2][k]-p[0][k]) *s[0]*s[0]) / (s[0]*s[n-2]);
        for(i = 1; i < n-1; i++)
        {
            b[i][k] = 3.0f * ((p[i][k]-p[i-1][k]) * s[i] * s[i]
              + (p[i+1][k]-p[i][k]) * s[i-1] * s[i-1]) / (s[i-1]*s[i]) ;
        }
//        b[n-2][k] = 3.0 * ((p[n-2][k]-p[n-3][k]) * s[n-2] * s[n-2]
//              + (p[0][k]-p[n-2][k]) * s[n-3] * s[n-3]) / (s[n-3]*s[n-2]) ;
    }

    //x座標の解
    for(i = 0; i < n-1; i++) bb[i] = b[i][0];
    SolutionByLU(n-1, a, bb);
    for(i = 0; i < n-1; i++) b[i][0] = bb[i];
    //ｙ座標の解
    for(i = 0; i < n-1; i++) bb[i] = b[i][1];
    SolutionByLU(n-1, a, bb);
    for(i = 0; i < n-1; i++) b[i][1] = bb[i];
    //z座標の解
    for(i = 0; i < n-1; i++) bb[i] = b[i][2];
    SolutionByLU(n-1, a, bb);
    for(i = 0; i < n-1; i++) b[i][2] = bb[i];
    //1次微係数はb[i][0],b[i][1],b[i][2]に格納されている
    //補間
    for(k = 0; k < 3; k++)
    {
        for(i = 0; i < m; i++)
        {   ip = i+1;
            if(ip == m) ip = 0;
            w0[k] = p[i][k];
            w1[k] = b[i][k];//1次微係数
            w2[k] = 3.0f * (p[ip][k]-p[i][k])/(s[i]*s[i]) - (2.0f*b[i][k]+b[ip][k])/s[i];
            w3[k] = 2.0f * (p[i][k]-p[ip][k])/(s[i]*s[i]*s[i]) + (b[i][k]+b[ip][k])/(s[i]*s[i]);
            u0 = s[i] / numDivision;
            for(j = 0; j < numDivision; j++)
            {
                u = u0 * (double)j;
                q[numDivision * i + j][k] = w0[k] + w1[k]*u + w2[k]*u*u + w3[k]*u*u*u;
            }
        }
    }

}
