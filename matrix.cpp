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
	if(a<b){
		c = a;
		a = b;
		b = c;
	}
	if(b == 0)return a;
	c = a%b;
	if(c == 0)return b;
	return gcd(b,c);
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
		Simplify();
	}
	void Simplify(){
		//if(down == 0)down = 1;
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
	 friend ostream & operator << (ostream &r,const Fraction &s){
		r <<"\t"<< s.up;
		if(s.down != 1 && s.up!= 0)r <<"/"<< s.down;
		return r;
	}
	Fraction operator+(const Fraction &right){
		Fraction temp;
		temp.up = up * right.down + right.up * down;
		temp.down = down * right.down;
		temp.Simplify();
		return temp;
	}
	Fraction operator+(const int &i){
		Fraction temp;
		temp.up = up + down * i;
		temp.Simplify();
		return temp;
	}

	Fraction operator-(const Fraction &right){ 
		Fraction temp;
		temp.up = up * right.down - right.up * down;
		temp.down = down * right.down;
		temp.Simplify();
		return temp;
	}
	Fraction operator-(const int &i){ 
		Fraction temp;
		temp.up = up - down * i; 
		temp.Simplify();
		return temp;
	}

	Fraction operator*(const Fraction &right){
		Fraction temp;
		temp.up = up * right.up;
		temp.down = down * right.down;
		temp.Simplify();
		return temp;
	}
	Fraction operator*(const int &i){
		Fraction temp;
		temp.up = up * i;
		temp.Simplify();
		return temp;
	}

	Fraction operator/(const Fraction &right){ 
		Fraction temp;
		temp.up = up * right.down;
		temp.down = down * right.up;
		temp.Simplify();
		return temp;
	}
	Fraction operator/(const int &i){
		Fraction temp;
		temp.down *= i;
		temp.Simplify();
		return temp;
	}

	Fraction operator-(){
		Fraction temp = (*this);
		temp.up = -temp.up;
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
};


typedef Fraction fra;

template<class T>
class Matrix{
private:
public:
	int row;
	int col;
	T *datas;

	Matrix(int r,int c){
		//先行后列
		//A(r,c) datas[r][c]
		row = r; col = c;
		int msize = r*c;
		datas = new T[msize];
		for(int i = 0;i<msize;i++){
			*(datas+i) = 0;
		}
	}
	Matrix(){
		row = col = 0;
	}
	~Matrix(){
		//if(row>0)delete [] datas;
	} 
/*
	Matrix<T>& operator= (const Matrix<T> & right){
		if(this == &right)return *this;
		int msize = right.row*right.col;
		//if(row>0)delete [] datas;
		datas = new T(msize);
		row = right.row;
		col = right.col;
		for(int i = 0;i<msize;i++)
			datas[i] = right.datas[i];
		return *this;
	}
	*/
	void CopyTo(Matrix<T> *target){
		int msize = row*col;
		target->datas = new T(msize);
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
	void Set(string d){
	 	int i = 0;
	 	for(int y = 0;y<row;y++){
			for(int x = 0;x<col;x++){
				string temp;
				while(i<d.length()){
					if(d[i]!=' ')
						temp += d[i];
					else
						if(temp.length()>0)break;
					i++;
	 			}
				(*this)[y][x] = T(temp);
	 		}
	 	}
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

	// ELEMENTARY ROW OPERATIONS	
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
	void REF(Matrix *ra){
		//Row Echelon Form
		int x,y,z;
		z = 0;
		Matrix &result = *ra;
		CopyTo(ra);
		for(x = 0;x<col;x++){
			for(y = z;y<row;y++){
				if(result[y][x]!=0)break;
			}
			if(y == row)continue;
			result.Interchange(z,y);
			for(int q = y+1;q<row;q++){
				if(result[q][x]!=0){
					T sc = -result[q][x]/result[z][x];
					result.Replacement(q,z,sc);
				}
			}
			z++;
		}
	}
	void RREF(Matrix *ra){
		//Reduced Row Echelon Form
		int x,y,z;
		z = 0;
		Matrix &result = *ra;
		CopyTo(ra);
		for(x = 0;x<col;x++){
			for(y = z;y<row;y++){
				if(result[y][x]!=0)break;
			}
			if(y == row)continue;
			result.Interchange(z,y);
			result.Scaling(z,T(1)/result[z][x]);
			for(int q = 0;q<row;q++){
				if(result[q][x]!=0&&q!=z){
					T sc = -result[q][x];
					result.Replacement(q,z,sc);
				}
			}
			z++;
		}
	}
}; 


typedef Matrix<Fraction> mat;
int main(){ 
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

	mat te(4,5);
	//te.Set("1 0 3 0 2 0 1 0 -3 3 0 -2 3 2 1 3 0 0 7 -5");
	te.Set("3 -6 0 0 0 0 0 3 -6 0 -1 2 0 0 0 0 0 -1 2 0");

	mat dac(4,5);

	//cout<<cz<<endl<<"--------------"<<endl;
	te.RREF(&dac);
	cout<<dac<<endl<<"--------------------------"<<endl;
	te.REF(&dac);
	cout<<dac<<endl<<"--------------------------"<<endl;
	//cout<<cz<<endl<<"--------------"<<endl;
	//cz = a*b;
	//cout<<"eee"<<fc.up<<" "<<fc.down<<endl;
	//cout<<fc<<endl<<"d"<<endl;
	cout<<a<<endl<<" X "<<endl<<b<<endl<<" = "<<endl<<cz<<endl;
	return 0;
}
