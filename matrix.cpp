#define DEBUG_FUN_MATH 1
#include<iostream>
#include<vector>
#include<string>
#include<iomanip>
//#include<cstdlib>
using namespace std;
int gcd(int a,int b){
	
	int c;
	if(a<0)a=-a;
	if(b<0)b=-b;
	//取消递归算法，提升性能
	if(a<b){
		c = a;
		a = b;
		b = c;
	}
	if(b == 0)return 1;
	while(1){
		a = a % b;
		if(a == 0)return b;
		b = b % a;
		if(b == 0)return a;
	}
}
int str2int(string s){
	int q = 0;
	bool minus = false;
	if(s[0]=='-')minus = true;
	for(int i=minus;i<s.length();i++){
		q*=10;
		q += s[i]-'0';
	}
	return minus?-q:q;
}
//将来支持高精度
class Fraction{
private:
	int up,down;
public:
	Fraction(){
		up = 0;
		down = 1;
	}
	Fraction(int a){
		up = a;
		down = 1;
	}
	Fraction(string str){ 
		string temp;
		int i = 0;
		down = 1;
		for(;i<str.length( );i++){
			if(str[i]=='/'){
				string t2;
				for(++i;i<str.length();i++){
					t2 += str[i];
				}
				down = str2int(t2);
				break;
			}
			temp += str[i];
		}
		up = str2int(temp);
	}
	
	Fraction(const char  ch[]){
		string temp;
		temp += ch;
		*this = Fraction(temp);
	}
	
	void Simplify(){
		if(down == 0)down = 1;
		if(up == 0)return;
		bool minus = false;
		if(up < 0){
			minus = true;
			up = -up;
		}
		if(down < 0){
			minus = !minus;
			down = -down;
		}
		int c = gcd(up,down);
		if(c!=0){
			up /= c;
			down /= c;
		}
		if(minus)up = -up;
	}
	friend ostream & operator << (ostream &r,Fraction &s){
		s.Simplify();
		r <<"\t"<< s.up;
		if(s.down != 1 && s.up!= 0)r <<"/"<< s.down;
		return r;
	}
	Fraction operator+(const Fraction &right){
		Fraction temp;
		temp.up = up * right.down + right.up * down;
		temp.down = down * right.down;
		Simplify();
		return temp;
	}
	Fraction operator+(const int &i){
		Fraction temp;
		temp.up = up + down * i;
		Simplify();
		return temp;
	}

	Fraction operator-(const Fraction &right){ 
		Fraction temp;
		temp.up = up * right.down - right.up * down;
		temp.down = down * right.down;
		Simplify();
		return temp;
	}
	Fraction operator-(const int &i){ 
		Fraction temp;
		temp.up = up - down * i; 
		Simplify();
		return temp;
	}

	Fraction operator*(const Fraction &right){
		Fraction temp;
		temp.up = up * right.up;
		temp.down = down * right.down;
		Simplify();
		return temp;
	}
	Fraction operator*(const int &i){
		Fraction temp;
		temp.up = up * i;
		Simplify();
		return temp;
	}

	Fraction operator/(const Fraction &right){ 
		Fraction temp;
		temp.up = up * right.down;
		temp.down = down * right.up;
		Simplify();
		return temp;
	}
	Fraction operator/(const int &i){
		Fraction temp;
		temp.down *= i;
		Simplify();
		return temp;
	}

	Fraction operator-(){
		Fraction temp = (*this);
		temp.up = -temp.up;
		Simplify();
		return temp;
	}

	Fraction& operator+=(const Fraction& right){
		up = up * right.down + right.up * down;
		down *= right.down;
		Simplify();
		return *this;
	} 

	Fraction& operator-=(const Fraction& right){
		up = up * right.down - right.up * down;
		down *= right.down;
		Simplify();
		return *this;
	}

	Fraction& operator*=(const Fraction& right){
		up = up * right.up;
		down *= right.down;
		Simplify();
		return *this;
	}

 	Fraction& operator/=(const Fraction& right){
  		up = up * right.down;
		down = down * right.up;
		Simplify();
		return *this;
  	}

	bool operator<(const Fraction& right){
		return (up*right.down<right.up*down);
	}
	bool operator<=(const Fraction& right){
		return (up*right.down<=right.up*down);
	}
	bool operator==(const Fraction& right){
		return (up*right.down==right.up*down);
	}
	bool operator==(const int &num){
		return (up==down*num);
	}
	bool operator!=(const Fraction& right){
		return (up*right.down!=right.up*down);
	}
	
	Fraction(int u,int d){
		up = u;
		down = d;
	}
	void Set(int u,int d){ 
		up = u;
		down = d;
	}
	void Set(string str){
		Fraction temp(str);
		*this = temp;
	} 

	int GetInt(){
		return up/down;
	}
	double GetDouble(){
		return up*1.0/down;
	}
};


typedef Fraction fra;

template<class T>
class Matrix{
private:
public:
	int row;
	int col;
	T *datas;
	
	void SetSize(int r,int c){
		row = r; col = c;
		int msize = r*c;
		if(datas != NULL)delete [] datas;
		if(msize == 0)return;
		datas = new T[msize];
		for(int i = 0;i<msize;i++){
			*(datas+i) = 0;
		}
	}

	Matrix(int r,int c){
		//先行后列
		//A(r,c) datas[r][c]
		row = r; col = c;
		int msize = r*c;
		if(msize == 0)return;
		datas = new T[msize];
		for(int i = 0;i<msize;i++){
			*(datas+i) = 0;
		}
	} 
	Matrix(){
		row = col = 0;
		datas = NULL;
	}
	~Matrix(){
		//if(row>0)delete [] datas;
	} 
	void CopyTo(Matrix<T> *target){
		if(target == this)return;
		int msize = row*col;
		target->datas = new T[msize];
		target->row = row;
		target->col = col;
		for(int i = 0;i<msize;i++)
			target->datas[i] = datas[i];
	}  

	T * const operator[](const int k){
		return &datas[k*col];
	}
	friend ostream & operator << (ostream &os,Matrix &s){
		for(int y = 0;y<s.row;y++){
			for (int x = 0;x<s.col;x++){
				os<<s[y][x]<<" ";
			}
			os<<endl;
		}
		return os;
	}
	void Set(int r,int c,T num){
		//出于数组储存需要,这里有(0,0)
		(*this)[r][c]=num;
	}  

	//暂时不作矩阵检查
	Matrix operator+(Matrix &right){
		Matrix temp;// = *this;
		CopyTo(&temp);
		for (int x = 0;x<col;x++){
			for(int y = 0;y<row;y++){
				temp[y][x]+=right[y][x];
			}
		}
		return temp;
	}
	Matrix operator-(Matrix &right){
		Matrix temp;// = *this;
		CopyTo(&temp);
		for (int x = 0;x<col;x++){
			for(int y = 0;y<row;y++){
				temp[y][x]-=right[y][x];
			}
		}
		return temp;
	}
	Matrix operator*(Matrix &right){
		Matrix temp(row,right.col);
		for(int x = 0;x<temp.col;x++){
			for(int y=0;y<temp.row;y++){
				T num = 0;//元素[y][x]
				for (int z = 0;z<col;z++){
					//cout<<(*this)[y][z]<<" * "<<right[z][x]<<" + ";
					num += (*this)[y][z]*right[z][x];
				}
				temp[y][x] = num;
				//cout<<endl;
			}
		}
		return temp; 
  	}
	
	Matrix operator*(const T &num){
		Matrix temp(row,col);
		CopyTo(&temp);
		int msize = col * row;
		for(int i=0;i<msize;i++)
			*(temp.datas+i) = num;

		return temp;
	}
	Matrix& operator*=(const T &num){
		int msize = col * row;
		for(int i=0;i<msize;i++)
			*(datas+i) *= num;

		return *this;
	}
	Matrix operator/(const T &num){
		Matrix temp(row,col);
		CopyTo(&temp);
		int msize = col * row;
		for(int i=0;i<msize;i++)
			*(temp.datas+i) /= num;

		return temp;
	}

	Matrix& operator/=(const T &num){
		int msize = col * row;
		for(int i=0;i<msize;i++)
			*(datas+i) /= num;

		return *this;
	}
	// ELEMENTARY ROW OPERATIONS	
	void Identity(){
		for(int x = 0;x<col;x++){
			for(int y = 0;y<row;y++){
				(*this)[y][x] = (x==y)?T(1):T(0);
			}
		}
	}
	void Interchange(int r1,int r2){
		if(r1 == r2)return;
		T temp;
		for(int i = 0;i<col;i++){
			temp = (*this)[r1][i];
			(*this)[r1][i] = (*this)[r2][i];
			(*this)[r2][i] = temp;
		}
	}
	void Scaling(int r,T m){
		for(int i = 0;i<col;i++){
			(*this)[r][i] *= m;
		}
	}
	void Replacement(int r1,int r2,T m){
		for(int i = 0;i<col;i++){ 
			(*this)[r1][i] += (*this)[r2][i] * m;
		}
	}

	//无参数输入时，为自我进行REF
	 void REF(Matrix *result = NULL){
		//Row Echelon Form
		int x,y,z;
		z = 0;
		if(result == NULL)result = this;

		CopyTo(result);

	 	for(x = 0;x<col;x++){
	 		for(y = z;y<row;y++){
				if((*result)[y][x]!=0)break;
			}
			if(y == row)continue;
			result->Interchange(z,y);
			for(int q = y+1;q<row;q++){
				if((*result)[q][x]!=0){
					T sc = -(*result)[q][x]/(*result)[z][x];
					result->Replacement(q,z,sc);
				}
			}
			z++;
		}
	}

	void RREF(Matrix *result = NULL){
		//Reduced Row Echelon Form
		int x,y,z;
		z = 0;
		if(result == NULL)result = this;
		CopyTo(result);

		for(x = 0;x<col;x++){
			for(y = z;y<row;y++){
				if((*result)[y][x]!=0)break;
			}
			if(y == row)continue;
			result->Interchange(z,y);
			result->Scaling(z,T(1)/(*result)[z][x]);
			for(int q = 0;q<row;q++){
				if((*result)[q][x]!=0&&q!=z){
					T sc = -(*result)[q][x];
					result->Replacement(q,z,sc);
				}
			}
			z++;
		}
	} 

	void Rotate180(Matrix<T> *result = NULL){
		Matrix<T> *source;
		if(result == NULL){
			result = this;
			source = new Matrix(row,col);
			CopyTo(source);
		}
		else{
			source = this;
		}
		for(int y=0;y<row;y++){
			for(int x=0;x<col;x++){
				(*result)[y][x] = (*source)[row-1-y][col-1-x];
			}
		}
	}
	
	void Transpose(Matrix<T> *result = NULL){
		Matrix<T> *source;
		if(result == NULL){
			result = this;
			source = new Matrix(col,row);
			CopyTo(source);
		}
		else{
			source = this;
		}
		for(int y=0;y<row;y++){
			for(int x=0;x<col;x++){
				result[x][y] = (*source)[y][x];
			}
		}
	}

	void Convolution(Matrix<T> &mould,Matrix<T> *result = NULL){
		Matrix<T> *source;
		if(result == NULL){
			result = this;
			source = new Matrix(col,row);
			CopyTo(source);
		}
		else{
			source = this;
		}
		int newRow = row - mould.row + 1;
		int newCol = col - mould.col + 1;
		if(newCol < 1||newRow < 1){
			result->SetSize(0,0);
			//(*result)[0][0] = T(0);
			return;
		}
		result->SetSize(newRow,newCol);
		T sum;
		for(int x = 0;x<newCol;x++){
			for(int y = 0;y<newRow;y++){
				sum = T(0);
				for(int mx = 0;mx<mould.col;mx++){
					for(int my = 0;my<mould.row;my++){
						sum += mould[mould.row-1-my][mould.col-1-mx] * (*source)[y+my][x+mx];
					}
				}
				(*result)[y][x] = sum;
			}
		}
	}

	bool Invert(Matrix *ra = NULL){
		if(row!=col)return false;
		Matrix &result = *ra;
		int newcol = col * 2;
		Matrix temp(row,newcol);
		for(int y=0;y<row;y++){
			for(int x=0;x<col;x++){
				temp[y][x] = (*this)[y][x];
			}
		}
		for(int y=0;y<row;y++){
			for(int x=col;x<newcol;x++){
				temp[y][x] = (y==x-col)?T(1):T(0);
			}
		}
		//cout<<temp<<endl;
		temp.RREF();
		//cout<<temp<<endl;
		for(int y=0;y<row;y++){
			for(int x=0;x<col;x++){
				if (temp[y][x] != ((y==x)?T(1):T(0)))return false;
			}
		}
		if(ra !=NULL){
			*ra = Matrix(row,col);
			for(int y=0;y<row;y++){
				for(int x=col;x<=newcol;x++){
					(*ra)[y][x-col] = temp[y][x];
				}
			}
		}
		return true;
	}
}; 

	void SetMatrix(Matrix<Fraction> * mat,string d){
	 	int i = 0;
	 	for(int y = 0;y<mat->row;y++){
			for(int x = 0;x<mat->col;x++){
				string temp;
				while(i<d.length()){
					if(d[i]!=' ')
						temp += d[i];
					else
						if(temp.length()>0)break;
					i++;
	 			}
				(*mat)[y][x] = temp;
	 		}
	 	}
	}
#if DEBUG_FUN_MATH
typedef Matrix<int> mat;
int main(){ 
	/*
	mat a(2,2),b(2,2),cz(2,2);
	a.Set("8 5 -7 -5");
	b.Set("1 1 -7/5 -8/5");
	cout<<endl<<endl;
	//a.Set("1 2 3 4");
	//b.Set("5 6 7 8");
	cz = a*b;
	fra fa("-7/5");
	fra fb("5");
	//cout<<cz<<endl<<"--------------"<<endl;
	//cout<<fa<<endl;
	//cout<<fb<<endl;
	fra fc = fa*fb;

	mat te(2,2);
	//te.Set("1 0 3 0 2 0 1 0 -3 3 0 -2 3 2 1 3 0 0 7 -5");
	//te.Set("3 -6 0 0 0 0 0 3 -6 0 -1 2 0 0 0 0 0 -1 2 0");
	te.Set("2 0 0 2");

	mat dac(4,5);

	//cout<<cz<<endl<<"--------------"<<endl;
	//te.RREF(&dac);
	//cout<<dac<<endl<<"--------------------------"<<endl;
	//te.REF(&dac);
	//cout<<dac<<endl<<"--------------------------"<<endl;
	te.Invert(&dac);
	cout<<dac<<endl;
	//cout<<cz<<endl<<"--------------"<<endl;
	//cz = a*b;
	//cout<<"eee"<<fc.up<<" "<<fc.down<<endl;
	//cout<<fc<<endl<<"d"<<endl;
	//cout<<a<<endl<<" X "<<endl<<b<<endl<<" = "<<endl<<cz<<endl;
	fra att = "21/14";
	cout<<att;
	*/
	mat a(5,5);
	a.Identity();
	//a.Set("1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25");
	//mat b(3,3);
	//b.Set("0 1/5 0 1/5 1/5 1/5 0 1/5 0");
	cout<<a<<endl;
	//cout<<b<<endl;
	//a.Convolution(b);
	//a*= 10;
	cout<<"--------------\n"<<a<<endl;
	mat d(8,8);
	d.Identity();
	d/=8;
	cout<<d;
	return 0;
}
#endif
