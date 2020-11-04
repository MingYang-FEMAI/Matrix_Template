//   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//                                                                   Liner Matrix by Ming YANG 
//                                                               young_94@126.com / 2020/11/02
//                                                     In Southern University of Science and Technology
//                           
//   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#pragma once
#include <iostream>
#include<fstream> 
#include<ostream> 
#include<ctime>
#include<Windows.h>
#include<cstdlib>
#include<math.h>
#include<exception>

using namespace std;

template<class T>
class Matrix {

private:
    int _rows;
    int _cols;
    T* L;                                                   //不可从外部访问
    T* U;                                              //不可从外部访问
    T detm;                                            //行列式

public:
    Matrix();
    Matrix(int N, int M);                                                     //构造函数  
    Matrix(int N, int M, T sig);                                    //带默认值的初始化矩阵
    Matrix(const Matrix<T>& mat);                                //拷贝构造函数  
    ~Matrix();                                                                  //析构函数
    Matrix<T> operator+(Matrix<T>& mat);                //重载+  
    Matrix<T> operator+(const T val);                         //重载 ，矩阵+常数  
    Matrix<T> operator-(Matrix<T>& mat);                //重载+  
    Matrix<T> operator - (const T val);                         //重载 ，矩阵+常数  
    Matrix<T> operator * (Matrix<T>& mat);                   //重载*  
    Matrix<T> operator * (const T mat);                          //重载*  ，矩阵乘以常数
    //Matrix<T> divs (const T& mat1, const T& mat2);                          //定义  ，矩阵除法
    Matrix<T> operator / (const T& mat);               // 重载，矩阵/ 除法
    Matrix<T>& operator=(Matrix<T>& mat);     //重载=  
    Matrix<T>& operator=(const T value);          //重载=   矩阵赋常值
    Matrix<T> ones();  // 生成全为1的矩阵  
    Matrix<T> zeros();    // 生成全为0的矩阵  
    Matrix<T> rands(T min, T max);    // 生成随机数填充的矩阵 
    Matrix<T> rands();   // 生成随机数填充的矩阵   重载
    T* _p;                                                     // 一维数组，按列存储，从0开始,
    unsigned int getRow() const;
    unsigned int getCol() const;
    T get(unsigned int x, unsigned int y);          //
    T getdet();          //返回矩阵行列式DET
    void exchangerow(const unsigned int x, const unsigned int y); 
    void exchangecol(const unsigned int x, const unsigned int y);
    void SetSingElem(const unsigned int x, const unsigned int y, T value);    //矩阵 赋值
    Matrix<T> Expand(const unsigned adx, const unsigned ady,const T val);             //矩阵扩充 扩充部分赋值为  val
    void InsertRow(const unsigned int x, const Matrix<T>& Vect);    //矩阵 插入行
    void InsertCol(const unsigned int y, const Matrix<T>& Vect);    //矩阵 插入列
    void PartInsert(const unsigned int xidx, const unsigned int yidx,const Matrix<T>& mtx);    //矩阵 插入区块 赋值为 mtx
    void PartDelete(const unsigned int xidx, const unsigned int yidx, const unsigned int xlen, const unsigned int ylen);    //矩阵 删除区块 赋值为0
    void DeleteRow(const unsigned int x);                                         //矩阵  删除行
    void DeleteCol(const unsigned int y);                                           //矩阵  删除列
    void SetElements();
    void SetRows(const unsigned int x, const Matrix<T>& Vect);       //单行赋值
    void SetCols(const unsigned int y, const Matrix<T>& Vect);         //单列赋值 
    void Display();
    Matrix<T> TR();  //转置
    Matrix<T> Remaind(const unsigned int x, const unsigned int y);    // 生成除去x，y点的余子式 
    Matrix<T> inv();            //求逆矩阵 
    T iteratdet();                     //迭代 求行列式
    Matrix<T>& LUdet();                        // LU 分解 
    void deleteLU();           //重置LU，释放空间

    friend ostream& operator<<(ostream& out, Matrix<T>& mat)    //重载<<  // friend 函数仅允许类内重载
    {
        for (int i = 0; i < mat._rows; i++)
        {
            for (int j = 0; j < mat._cols; j++)
            {
                out << mat._p[j * mat._rows + i] << "  ";
            }
            out << endl; 
        }
        return out;
    }

    friend istream& operator>>(istream& in, Matrix<T>& mat)     //重载>>  // friend 函数仅允许类内重载
    {
        for (int i = 0; i < mat._rows * mat._cols; i++) {
            in >> mat._p[i];
        }
        return in;
    }
};

template<class T>
Matrix<T>::Matrix() {
    _rows = 0;
    _cols = 0;
    _p = 0;
    L = 0;
    U = 0;
    detm = 0;
}

template<class T>
Matrix<T>::Matrix(int N, int M)                                             //  M N 是代表实际存在的数目，从1开始
{                               
    if (N <= 0 || M <= 0) {
        cout << "N= " << N << "   M=  " << M << endl;
        cout<<"Wrong! initial input invalid "<<endl;
        exit(1);
    }
    _rows = N;
    _cols = M;
    _p = new T[_rows*_cols];
    for (int ii = 0; ii < M * N; ii++) 
    {
        _p[ii]= static_cast<T>(1);       //特别注意
    }
    L = 0;
    U = 0;
    detm = 0;
}

template<class T>
Matrix<T>::Matrix(int N, int M, T sig)                                             //  带初值的初始化
{
    if (N <= 0 || M <= 0) {
        cout << "N= " << N << "   M=  " << M << endl;
        cout << "Wrong! initial input invalid " << endl;
        exit(1);
    }
    _rows = N;
    _cols = M;
    _p = new T[_rows * _cols];
    for (int ii = 0; ii < M * N; ii++)
    {
        _p[ii] = static_cast<T>(sig);
    }
    L = 0;
    U = 0;
    detm = sig;
}

template<class T>
Matrix<T>::Matrix(const Matrix<T>& mat) {          // LU不会被初始化
    _rows = mat._rows;
    _cols = mat._cols;
    detm = mat.detm;    _p = new T[_rows*_cols];
    for (int i = 0; i < _rows*_cols; i++) {
        _p[i] = mat._p[i];
    }
    L = 0;
    U = 0;
}

template<class T>
Matrix<T>::~Matrix() {
    delete[]_p;
}

template<class T>
Matrix<T> Matrix<T>::operator+ (Matrix<T>& mat) {
    if (_rows != mat._rows || _cols != mat._cols) {
        throw "Wrong! expect Matrixes in Same Size";   //抛出异常  
        exit(1);
    }
    else {
        Matrix<T> matSum(_rows, _cols);
        for (int i = 0; i < _rows* _cols; i++) {
                matSum._p[i] = _p[i] + mat._p[i];      // 定义为一维数组，计算简单
        }
        return matSum;
    }
}

template<class T>
Matrix<T> Matrix<T>::operator+ (const T val) {             //重载，矩阵+常数
    Matrix<T> matSum(_rows, _cols);
    for (int i = 0; i < _rows * _cols; i++) {
        thmatSum._p[i] = _p[i] + val;
    }
    return matSum;
}

template<class T>
Matrix<T> Matrix<T>::operator- (Matrix<T>& mat) {
    if (_rows != mat._rows || _cols != mat._cols) {
        throw "Wrong! expect Matrixes in Same Size";   //抛出异常  
        exit(1);
    }
    else {
        Matrix<T> matSub(_rows, _cols);
        for (int i = 0; i < _rows * _cols; i++) {
            matSub._p[i] = _p[i] - mat._p[i];      // 定义为一维数组，计算简单
        }
        return matSub;
    }
}

template<class T>
Matrix<T> Matrix<T>::operator- (const T val) {             //重载，矩阵+常数
    Matrix<T> matSub(_rows, _cols);
    for (int i = 0; i < _rows * _cols; i++) {
        matSub._p[i] = _p[i] - val;
    }
    return matSub;
}

template<class T>
Matrix<T> Matrix<T>::operator* (Matrix<T>& mat)    // 重载，矩阵*矩阵
{
    Matrix<T> matProduct(_rows, mat._cols,0);   // colA=rowB   会被初始化为0；
    if (_cols != mat._rows) 
    {
        throw "Incoordinate Matrix: in Rows and Cols！";
        exit(1);
    }
    else 
    {
        for (int idxO = 0; idxO < _cols; idxO++)
        {
            for (int idxI = 0; idxI < _rows * mat._cols; idxI++)
            {
                //cout << "p      IDX   ( " <<  idxI % _rows+1<<", "<<idxO+1<<" )" << endl;
                //cout << "mat   IDX   ( " <<  idxO+1<< ", " << idxI / _rows+1<<" ) " << endl;
                //cout << "result IDX    ( " << idxI%_rows+1<<", "<<idxI/_rows+1<<" )" <<"    Series = "<<idxO+1<< endl;
                //T temp = _p[idxO * _rows + idxI % _rows] * mat._p[(idxI / _rows) * mat._rows + idxO];
                //cout << temp << endl;
                matProduct._p[idxI] = matProduct._p[idxI] + _p[idxO * _rows + idxI % _rows] * mat._p[(idxI / _rows) * mat._rows + idxO];
                //     A×B=R  ↓↓↓
                //   idxO(idx-Outloop) : 0 -> colA;  idxI(idx -Innerloop) : 0->rowA*colB
                //   R[idxI]+=A[idxO*rowA+idxI%rowA]*B[(idxI/rwoA)*rowB+idxO];
            }
        }
 /*       for (int i = 0; i < _rows; i++) 
        {
            for (int j = 0; j < mat._cols; j++)
            {
                matProduct._p[j * _rows + i] += _p[j * _cols + i] * mat._p[i * mat._rows + j];
            }
        }*/
    }
    return matProduct;
}

template<class T>
Matrix<T> Matrix<T>::operator/ (const T& mat)   // 定义，矩阵  /除法    //  Failed Unreasonablly!!
{
    Matrix<T> matDiv(_rows, _cols);
    matDiv = mat.inv();
    matDiv = (*this) * matDiv;
    return matDiv;
}

template<class T>
Matrix<T> Matrix<T>::operator* (const T val)   // 重载，矩阵*常数
{
    Matrix<T> matProduct(_rows, mat._cols);   
    for (int i = 0; i < _rows * _cols; i++) {
        matProduct._p[i] = _p[i] * val;
    }
    return matProduct;
}

template<class T>
Matrix<T>& Matrix<T>::operator=(Matrix<T>& mat)    //  understand deeply who`s "=" here! 
{
    this->_cols = mat._cols;
    this->_rows = mat._rows;
    if (_cols < 0 || _rows < 0) {  cout << "failed to conduct value! " << endl; exit(1); }
    for (int i = 0; i < mat._rows * mat._cols; i++)
    {
        this->_p[i] = mat._p[i];                      //不考虑传 LU 和 det
    }
    //cout << "<<<   = \n" << endl;
    //this->Display();
    return (*this);
}

template<class T>
Matrix<T>& Matrix<T>::operator=(const T value) 
{
    for (int i = 0; i < _rows * _cols; i++) {
        this->_p[i] = value;
    }
    return (*this);
}
template<class T>
Matrix<T> Matrix<T>:: ones()   // 生成全为1的矩阵  
{
    Matrix<T> temp(_rows, _cols);
    return temp;
}
template<class T>
Matrix<T> Matrix<T>::zeros()    // 生成全为0的矩阵  
{
    Matrix<T> temp(_rows, _cols);
    temp = temp - static_cast<T>(1);
    return temp;
}

template<class T>
Matrix<T> Matrix<T>::rands(T min, T max)    // 生成随机数填充的矩阵 
{
    srand(time(0));
    Matrix<T> temp(_rows, _cols);
    for (int ii = 0; ii < _rows * _cols; ii++)
    {
        T output = static_cast<T>(rand() % int(ceil(max)));
        temp._p[ii] = output / ceil(max) * max + min;
    }
    //temp.Display();
    return temp;
}
template<class T>
Matrix<T> Matrix<T>::rands()    // 生成随机数填充的矩阵   重载
{
    T max = static_cast<T>(10);
    T min = static_cast<T>(0);
    srand(time(0));
    Matrix<T> temp(_rows, _cols);
    for (unsigned int i = 0; i < _rows * _cols; i++)
    {
        T output = static_cast<T>(rand() % int(ceil(max)));
        temp._p[i] = output / ceil(max) * max + min;
    }
    return temp;
}

template <class T>
unsigned int Matrix<T>::getRow() const
{
    return _rows;
}
template <class T>
unsigned int Matrix<T>::getCol() const
{
    return _cols;
}

template <class T>
T Matrix<T>::get(unsigned int x, unsigned int y) 
{ 
    return this->_p[y* _rows + x]; 
}

template <class T>
T Matrix < T>::getdet() 
{
    return detm;
}

template <class T>
void Matrix<T>::exchangerow(const unsigned int x, const unsigned int y)     // 行互换
{
    if (x < _rows && x >= 0 && y >= 0 && y < _rows)
    {
        unsigned int i;
        T temp;
        for (i = 0; i < _cols; i++)
        {
            temp = _p[x+i*_rows];
            _p[x + i * _rows] = _p[y + i * _rows];
            _p[y + i * _rows] = temp;
        }
    }
    else
    {
        cout << "The row(s) you specified is not existed!" << endl;
        exit(1);
    }
}

template <class T>
void Matrix<T>::exchangecol(const unsigned int x, const unsigned int y)   // 列互换
{
    if (x <= _cols && x >= 1 && y >= 1 && y <= _cols)
    {
        unsigned int i;
        T temp;
        for (i = 0; i < _rows; i++)
        {
            temp = _p[_rows*x + i];
            _p[_rows * x + i]= _p[_rows * y + i];
            _p[_rows * y + i] = temp;
        }
    }
    else
    {
        cout << "The col(s) you specified is not existed!" << endl;
        exit(1);
    }
}
template <class T>
void Matrix<T>::SetSingElem(const unsigned int x, const unsigned int y, T value) 
{
    _p[x + y * _rows] = value;
}

template <class T>
Matrix<T> Matrix<T>::TR()     // 转置 
{
    unsigned int i, j;
    Matrix<T> temp(_rows,_cols);
    for (i = 0; i < _rows * _cols; i++)
    {
        int ii = i % _rows;
        int ij = i / _rows;
        temp._p[ii*_rows+ij] = _p[i];
    }
    //cout << "<<<    TR   \n" << endl;
    //temp.Display();
    return temp;
}

template<class T>
void Matrix<T>::SetElements() {
    cout << "请输入矩阵的元素，共" << _rows * _cols << "个：" << endl;
    for (int i = 0; i < _rows* _cols; i++) {
        cin >> _p[i];
    }
}

template<class T>
void Matrix<T>::Display() {
    cout << "Display the current Matrix:" << endl;
    for (int i = 0; i < _rows; i++) 
    {
        for (int j =0 ; j < _cols; j++)
        {
            cout << _p[j * _rows + i] << "  ";
        }
        cout << endl;
    }
    cout << "Next<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
}

template<class T>
void Matrix<T>::SetRows(const unsigned int x, const Matrix<T>& Vect)       //单行赋值
{
    int lenV = 0;
    if (Vect._cols == 1)
    {
        lenV = Vect._rows > _rows ? _rows : Vect._rows;
        for (int ii = 0; ii <lenV ;ii++)
        {
            _p[ii * _rows + x] = Vect._p[ii];
        }
    }
    else if (Vect._rows==1)
    {
        lenV = Vect._cols > _rows ? _rows : Vect._cols;
        for (int ii = 0; ii < lenV; ii++)
        {
            _p[ii * _rows + x] = Vect._p[ii];
        }
    }
    else
    {
        cout << "wrong in input a vector for set rows! " << endl;
        exit(1);
    }
}

template<class T>
void Matrix<T>::SetCols(const unsigned int y, const Matrix<T>& Vect)         //单列赋值
{
    int lenV = 0;
    if (Vect._cols == 1)
    {
        lenV = Vect._rows > _rows ? _rows : Vect._rows;
        //cout << "Len V  " << lenV << endl;
        for (int ii = 0; ii < lenV; ii++)
        {
            _p[y * _rows + ii] = Vect._p[ii];
        }
    }
    else if (Vect._rows == 1)
    {
        lenV = Vect._cols > _rows ? _rows : Vect._cols;
        //cout << "Len V  " << lenV << endl;
        for (int ii = 0; ii < lenV; ii++)
        {
            _p[y * _rows + ii] = Vect._p[ii];
        }
    }
    else
    {
        cout << "wrong in input a vector for set rows! " << endl;
        exit(1);
    }
}

template<class T>
Matrix<T> Matrix<T>::Expand(const unsigned adx, const unsigned ady, const T val)             //矩阵扩充 扩充部分赋值为  val
{
    Matrix<T> tempmtx(_rows + adx, _cols + ady);
    tempmtx = tempmtx - static_cast<T>(1) + val;
    //tempmtx.Display();
    //cout << "<<< aftrt insert \n" << endl;
    tempmtx.PartInsert(0, 0, (*this));
    tempmtx.Display();
    return tempmtx;                                             //??
}

template<class T>
void Matrix<T>::InsertRow(const unsigned int xidx, const Matrix<T>& Vect)    //指定位置 矩阵 插入行，不为实现而实现，应该注意速度，减少交叉调用
{
    if((Vect._cols!=1 && Vect._rows!=1) || xidx>_rows)
    {
    cout << "failed!  invalid insert " << endl;
    exit(1);
    }

    Matrix<T> tempM(_rows + 1, _cols);
    T lenv;
    lenv = Vect._cols > Vect._rows ? Vect._cols : Vect._rows;
    for (int idx = 0; idx < _cols * (_rows + 1); idx++)
    {
        if (idx % (_rows+1) < xidx)
        {
            tempM._p[idx] = _p[idx - idx / (_rows + 1)];
        }
        else if (idx % (_rows+1) == xidx)
        {
            T vp;
            if (idx / (tempM._rows) >= lenv)
            {
                cout << "Warning! Vect is not long enough, add zero at the end ! " << endl;
                vp = static_cast<T>(0);
            }
            else { vp = Vect._p[idx / (_rows + 1)]; }
            tempM._p[idx] = vp;
        }
        else {                                                        // idx%_rows>=xidx;
            tempM._p[idx] = _p[idx-idx / (_rows+1) -1];
        }
    }
    (*this) = tempM;   // ??
}

template<class T>
void Matrix<T>::InsertCol(const unsigned int yidx, const Matrix<T>& Vect)    //指定位置 矩阵 插入列
{
    if ((Vect._cols != 1 && Vect._rows != 1) || yidx>_cols)
    {
        cout << "failed!  invalid insert " << endl;
        exit(1);
    }
    T lenv;
    lenv = Vect._cols > Vect._rows ? Vect._cols : Vect._rows;
    Matrix<T> tempM(_rows, _cols + 1);
    for (int idx = 0; idx < (_cols +1)* _rows ; idx++)
    {
        if (idx / _rows < yidx)
        {
            tempM._p[idx] = _p[idx];
        }
        else if (idx / _rows == yidx)
        {
            T vp;
            if (idx % _rows >= lenv)
            {
                cout << "Warning! Vect is not long enough, add zero at the end ! " << endl;
                vp = static_cast<T>(0);
            }
            else { vp = Vect._p[idx % _rows]; }
            tempM._p[idx] = vp;
        }
        else {                                                        // idx /_rows>=yidx;
            tempM._p[idx] = _p[idx  - _rows];
        }
    }
    (*this) = tempM;   // yes
}

template<class T>
void Matrix<T>::DeleteRow(const unsigned int xidx)    //矩阵  删除行
{
    if (xidx >= _rows)
    {
        cout << "X index overteps the Matrix size " << endl;
        exit(1);
    }
    for (unsigned int ii = 0; ii < _cols * _rows; ii++)                   //此处速度可继续优化提升
    {
        _p[ii - ii / _rows - (ii % _rows > xidx)] = _p[ii];
    }
    _rows = _rows - 1;
}

template<class T>
void Matrix<T>::DeleteCol(const unsigned int yidx)      //矩阵  删除列
{
    if (yidx >= _cols)
    {
        cout << "Y index overteps the Matrix size " << endl;
        exit(1);
    }
    for (unsigned int ii = (yidx+1)*_rows; ii < _cols * _rows; ii++)
    {
        _p[ii - _rows] = _p[ii];      //其他的元素不管了
    }
    _cols = _cols - 1;
}

template<class T>
void Matrix<T>::PartInsert(const unsigned int xidx, const unsigned int yidx, const Matrix<T>& mtx)    //矩阵 插入区块, 大小不变
{
    for (unsigned int idx = 0; idx < mtx._cols * mtx._rows; idx++)
    {
        unsigned int row = idx % mtx._rows + xidx;                 //快速建立联系的关键公式
        unsigned int col = idx / mtx._rows + yidx;
        if (row<=_rows && col<=_cols)
        {
            _p[col*_rows+row] = mtx._p[idx];                         //精巧的设计，达到快速计算的效果
        }
    }
}

template<class T>
void Matrix<T>::PartDelete(const unsigned int xidx, const unsigned int yidx, const unsigned int xlen, const unsigned int ylen)   //矩阵 删除区块值， 赋值0
{
    for (unsigned int idx = 0; idx < xlen * ylen; idx++)
    {
        unsigned int row = idx % xlen + xidx;
        unsigned int col = idx / xlen + yidx;
        if (row <= _rows && col <= _cols)
        {
            _p[col * _rows + row] = 0;                        
        }
    }
}


template <class T>
Matrix<T> Matrix<T>::Remaind(const unsigned int xidx, const unsigned int yidx)    // 生成除去x，y点的余子式 
{
    Matrix<T> temp(_rows-1, _cols-1);
    int nidx = 0;
    for (int idx = 0; idx < _rows * _cols; idx++)
    {
        if (idx % _rows != xidx && idx / _rows != yidx)
        {
            temp._p[nidx] = _p[idx];
            nidx++;
        }
    }
    //temp.Display();
    return temp;
}

template <class T>
Matrix<T> Matrix<T>::inv()            //求矩阵的逆
{
    if (_rows != _cols || _rows < 2)
    {
        cout << "Matrix Should in Qualified Size!" << endl;
        exit(1);
    }
    this->LUdet();             //激活L U
    if (this->detm < static_cast<T>(1e-16) && this->detm > static_cast<T>(-1e-16)) { cout << "DET=0,  inverse don`t exist! " << endl; exit(1); }
    T* invu=new T(_rows* _cols);
    T* invl=new T(_rows* _cols);
    for (int i = 0; i < _rows * _cols; i++)
    {
        invu[i] = 0;
        invl[i] = 0;
    }
    Matrix<T> inverse(_rows, _cols, 0);

    for (int i = 0; i < _rows; i++)                                                          //求矩阵U的逆 
    {
        invu[i * _rows + i] = 1 / U[i*_rows+i];                                     //对角元素的值，直接取倒数
        //invl._p[i * _rows + i] = 1;
        for (int k = i - 1; k >= 0; k--)
        {
            T s = 0;
            for (int j = k + 1; j <= i; j++){ s = s + U[j * _rows + k] * invu[i * _rows + j]; }
            invu[i*_rows+k] = -s / U[k*_rows+k];                         //迭代计算，按列倒序依次得到每一个值，
        }
    }

    for (int i = 0; i < _rows; i++)                                     //求矩阵L的逆
    {
        invl[i*_rows+i] = 1;                                     //对角元素的值，直接取倒数，这里为1
        for (int k = i + 1; k < _rows; k++)
        {
            for (int j = i; j <= k - 1; j++)
            {
                invl[i * _rows + k] = invl[k + i * _rows] - L[k + j * _rows] * invl[j + i * _rows];   //迭代计算，按列顺序依次得到每一个值
            }
        }
    } 

    for (int sidx = 0; sidx < _rows; sidx++)
    {
        for (int bidx = 0; bidx < _rows*_rows; bidx++)
        {
            inverse._p[(bidx/_rows) * _rows + sidx] = inverse._p[(bidx / _rows) * _rows + sidx]+ invu[(bidx %_rows) * _cols + sidx] * invl[(bidx/_rows)*_rows + bidx%_rows];
            //     A×B=R  ↓↓↓
            //  R[(bidx/rowB)*rowA+sidx]+=A[(bidx%rowB)*rowA+sidx]*B[(bidx/rowB)*rowB+bidx%rowB];
        }
    }
    deleteLU();                           //及时释放空间
    return inverse;
}

template <typename T>
Matrix<T>& Matrix<T>::LUdet()   // LU分解求行列式值
{
    if (_rows < 2 || _rows != _cols ) { cout << "Wrong! this matrix is not square! " << endl; exit(1); }
    this->L = new T[_rows * _cols];
    this->U = new T[_rows * _cols];
    this->detm = static_cast<T>(1);
    if(_rows == 2)
    {
        if (_p[0] == 0) { cout << "Wrong! this matrix don`t have det!" << endl; exit(1); }
        this->detm = (_p[0] * _p[3] - _p[1] * _p[2]);
        this->L[0] = 1; this->L[1] = _p[1] / _p[0]; this->L[2] = 0; this->L[3] = 1;
        this->U[0] = _p[0]; this->U[1] = 0; this->U[2] = _p[2]; this->U[3] = _p[3] - (_p[1] * _p[2]) / _p[0];
        return *this;
    }

    for (int i = 0; i < _cols; i++)                   // 第一步，计算L的第一列和U的第一行：U1i=A1i,Li1=Ai1/U11
    {
        U[i*_rows] = _p[i*_rows];
        L[i] = _p[i] / _p[0];
        L[i * _rows + i] = 1;
        //for (int j = i+1; j < _rows; j++)              //Useless Actually
        //{
        //    U[i*_rows+j] = 0;
        //    L[(_rows-i-1)*_rows+ j-1-i] = 0; 
        //}
    }

    for (int r = 1; r < _rows; r++) {         //外层大循环(认为是对角循环)
        for (int i = r; i < _rows; i++) {        //行列循环
            T re = static_cast<T>(0);
            for (int k = 0; k < r; k++) {
                re = re + L[r + k*_rows] * U[k + i*_rows];               // U[k][i]*L[r][k]
            }
            U[r + i*_rows] = _p[r + i*_rows] - re;
            if(i>r)                                                                    //同时循环无误，检查不会出现预先引用
            {
                T re2 = static_cast<T>(0);
                for (int k = 0; k < r; k++) {
                    re2 = re2 + L[i + k*_rows] * U[k + r*_rows];                  // U[k][r]*L[i][k]
                }
                L[i + r*_rows] = (_p[i + r*_rows] - re2) / U[r + r*_rows];
            }
        }
    }

    T temp = 1;
    for(int ii=0;ii<_rows;ii++){ temp = temp * this->U[ii + ii*_rows]; }
    this->detm = temp;
    return (*this);
}

template <typename T>
void Matrix<T>::deleteLU()   //释放 LU 的空间
{
    L=0;
    U = 0;
}

template <typename T>
T Matrix<T>::iteratdet()                                           // 迭代法分解求矩阵行列式的值  错误
{
    if (_rows == _cols)
    {
        T temp = 0;
        if (_rows == 1) { return _p[0]; }
        else if (_rows == 2) {  return _p[0] * _p[3] - _p[1] * _p[2];  }
        else {
            Matrix<T> M0(*this);    // 余子式
            for (int i = 0; i < M0._cols; i++)
            {
                T mpi = M0._p[i*M0._rows];
                M0 = M0.Remaind(0, i);
                temp = temp+pow(-1, (i+1)) * mpi * M0.iteratdet();
            }
        }
        return temp;
    }
    else
    {
        cout << "This matrix is not a square matrix!" << endl;
        exit(1);
    }
}


//  @@@@@@@@@@@@@@@@@   FOR TEST  MATRIX   @@@@@@@@@@@@@@@
// 1-1  Basic  Matrix construction functions

//int main()
//{
//    int m = 4;
//    int n = 4;
//    Matrix<int> mat1(m, n);
//    mat1 = mat1.rands(1, 5);
//    cout << "INIT MATRIX   " << endl;
//    mat1.Display();
//    cout << "<<<\n" << endl;
//    mat1 = mat1.Remaind(1, 1);
//    cout << "Remaind MATRIX   " << endl;
//    mat1.Display();
//    cout << "<<<\n" << endl;
//    Matrix<int> mat1(m, n);
//    Matrix<int> mat2(m, n);
//    Matrix<int> vect(1, 4);
//    Matrix<int> vect2(1, 4);
//    mat1 = mat1.rands(1, 5);
//    cout << "SET VECT1    " << endl;
//    vect.SetElements();
//    cout << "<<<\n" << endl;
//
//    cout << "SET ROWS   " << endl;
//    mat1.SetCols(1, vect);
//    cout << endl;
//    cout << "mat 1  After set rows " << endl;
//    mat1.Display();
//    cout << "<<<\n" << endl;
//
//    cout << "SET VECT 2    " << endl;
//    vect2.SetElements();
//    cout << "<<<\n" << endl;
//
//    cout << "Insert rows " << endl;
//    mat1.InsertRow(1, vect2);
//    cout << endl;
//    cout << "mat 1  After insert rows  " << endl;
//    mat1.Display();
//    cout << "<<<\n" << endl;
//
//    cout << "Insert cols " << endl;
//    mat1.InsertCol(3, vect);
//    cout << endl;
//    cout << "mat 1  After insert cols  " << endl;
//    mat1.Display();
//    cout << "<<<\n" << endl;
//
//    cout << "Delete Rows " << endl;
//    mat1.DeleteRow(2);
//    cout << endl;
//    cout << "mat 1  After Delete Rows  " << 2 << endl;
//    mat1.Display();
//    cout << "<<<\n" << endl;
//
//    cout << "Delete cols " << endl;
//    mat1.DeleteCol(0);
//    cout << endl;
//    cout << "mat 1  After Delete cols  " << 0 << endl;
//    mat1.Display();
//    cout << "<<<\n" << endl;
//
//
//    cout << endl;
//    cout << "Vector 2 " << endl;
//    vect2.Display();
//    mat1.SetCols(2, vect2);
//    cout << "mat 2   " << endl;
//    cout << mat1;
//
//    return 0;
//}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//1-2       LUdet & Remainder

//int main()
//{
//    int m = 3;
//    int n = 3;
//    double tdet = 0;
//    Matrix<double> mat1(m, n);
//    Matrix<double> mat2(m, n);
//    mat1 = mat1.rands(0.1, 5.4);
//    cout << "INIT MATRIX   " << endl;
//    mat1.Display();
//    cout << "<<<\n" << endl;
//    mat2 = mat1.Remaind(0, 0);
//    cout << "Remaind MATRIX 0 - 0  " << endl;
//    mat2.Display();
//    cout << "<<<\n" << endl;
//    tdet = tdet + mat1.get(0, 0) * (mat2.get(0, 0) * mat2.get(1, 1) - mat2.get(1, 0) * mat2.get(0, 1));
//
//    mat2 = mat1.Remaind(0, 1);
//    cout << "Remaind MATRIX 0 - 1  " << endl;
//    mat2.Display();
//    cout << "<<<\n" << endl;
//    tdet = tdet - mat1.get(0, 1) * (mat2.get(0, 0) * mat2.get(1, 1) - mat2.get(1, 0) * mat2.get(0, 1));
//
//    mat2 = mat1.Remaind(0, 2);
//    cout << "Remaind MATRIX 0 - 2  " << endl;
//    mat2.Display();
//    cout << "<<<\n" << endl;
//    tdet = tdet + mat1.get(0, 2) * (mat2.get(0, 0) * mat2.get(1, 1) - mat2.get(1, 0) * mat2.get(0, 1));
//
//    mat1.LUdet();
//    //tdet = mat1.get(0, 0) * mat1.get(1, 1) - mat1.get(1, 0) * mat1.get(0, 1);
//    cout << "DET   is    " << mat1.getdet() << endl;
//    cout << "the true DET   is    " << tdet << endl;
//
//    return 0;
//}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//1-23      Inverse matrix

//int main()
//{
//    int m = 5;
//    int n = 5;
//    double tdet = 0;
//    Matrix<double> mat1(m, n);
//    Matrix<double> mat2(m, n);
//    mat1 = mat1.rands(1, 10);
//    cout << "INIT MATRIX   " << endl;
//    mat1.Display();
//
//    mat1.LUdet();
//    cout << "DET   is    " << mat1.getdet() << endl;
//    mat2 = mat1.inv();
//    cout << "Inv MATRIX   " << endl;
//    mat2.Display();
//    return 0;
//}