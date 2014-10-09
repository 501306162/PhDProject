#ifndef AUXCLASSES_H
#define AUXCLASSES_H

#include <iostream>
#include <math.h>
#include <complex>
#include <vector>

using namespace std;

#define VAL_INVALID 1000000
#define VECT_INVALID Vector3(1000000, 1000000, 1000000)
#define VECT4_INVALID Vector4(1000000, 1000000, 1000000, 1000000)
#define VECT_NULL Vector3(0, 0, 0)
#define VECT_NULL_4 Vector4(0, 0, 0, 0)
#define MAT_NULL Matrice4x4(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

class VectorN
{
	private:
		int dimension;
		double* val;
		
	public:
		VectorN();
		VectorN(int dim, double* v);
		VectorN(const VectorN& vect);
		~VectorN();
		
		int dim() const;
		double& operator[](unsigned int i);
		double operator[](unsigned int i) const;
		VectorN& operator=(const VectorN& vect);
		
		VectorN& operator+=(const VectorN& vect);
		VectorN& operator-=(const VectorN& vect);
		VectorN& operator*=(double fac);
		VectorN& operator/=(double div);
		
		VectorN operator-() const;
		VectorN operator-(const VectorN& vect) const;
		VectorN operator+(const VectorN& vect) const;
		
		VectorN operator*(double cst) const;
		double operator*(const VectorN& vect) const;
		VectorN operator/(double div) const;
		
		bool operator==(const VectorN& vect) const;
		bool operator!=(const VectorN& vect) const;
		
		double norm2() const;
		double norm() const;
};

VectorN vectProduct(const VectorN& V1, const VectorN& V2);

class Vector3
{
	private:
		double a;
		double b;
		double c;
		
	public:
		Vector3();
		Vector3(double A, double B, double C);
		Vector3(double* v);
		Vector3(const Vector3& vect);
		
		//bool isNull() const;
		
		double& operator[](int i);
		double operator[](int i) const;
		Vector3& operator=(const Vector3& vect);
		
		Vector3& operator+=(const Vector3& vect);
		Vector3& operator-=(const Vector3& vect);
		Vector3& operator*=(double fac);
		Vector3& operator/=(double div);
		
		Vector3 operator-() const;
		Vector3 operator-(const Vector3& vect) const;
		Vector3 operator+(const Vector3& vect) const;
		
		Vector3 operator*(double cst) const;
		double operator*(const Vector3& vect) const;
		Vector3 operator/(double div) const;
		Vector3 operator/(Vector3 div) const;
		
		bool operator==(const Vector3& vect) const;
		bool operator!=(const Vector3& vect) const;
		
		double norm2() const;
		double norm() const;
};

Vector3 vectProduct(const Vector3& V1, const Vector3& V2);
Vector3 projectionProduct(const Vector3& V1, const Vector3& V2);

ostream& operator<<(ostream& stream, const Vector3& vect);

class Vector4
{
	private:
		double a;
		double b;
		double c;
		double d;
		
	public:
		Vector4();
		Vector4(double A, double B, double C, double D);
		Vector4(double* v);
		Vector4(const Vector4& vect);
		
		double& operator[](int i);
		double operator[](int i) const;
		Vector4& operator=(const Vector4& vect);
		
		Vector4& operator+=(const Vector4& vect);
		Vector4& operator-=(const Vector4& vect);
		Vector4& operator*=(double fac);
		Vector4& operator/=(double div);
		
		Vector4 operator-() const;
		Vector4 operator-(const Vector4& vect) const;
		Vector4 operator+(const Vector4& vect) const;
		
		Vector4 operator*(double cst) const;
		double operator*(const Vector4& vect) const;
		Vector4 operator/(double div) const;
		
		bool operator<(const Vector4& vect) const;
		
		bool operator==(const Vector4& vect) const;
		bool operator!=(const Vector4& vect) const;
		
		double norm2() const;
		double norm() const;
};

ostream& operator<<(ostream& stream, const Vector4& vect);

class Matrice3x3
{
	private:
		double data[9];
		
	public:
		Matrice3x3();
		Matrice3x3(double a00, double a01, double a02,
			   double a10, double a11, double a12,
			   double a20, double a21, double a22);
		Matrice3x3(Vector3 col0, Vector3 col1, Vector3 col2);
		Matrice3x3(double* d);
		Matrice3x3(const Matrice3x3& mat);
		
		Matrice3x3& operator=(const Matrice3x3& mat);
		
		double& operator() (unsigned row, unsigned col);
		double  operator() (unsigned row, unsigned col) const;
		
		double det();
		Matrice3x3 inv();
		
		Matrice3x3& operator+=(const Matrice3x3& mat);
		
		Matrice3x3 operator*(double cst);
		Vector3 operator*(const Vector3& vect);
		
		friend ostream& operator<<(ostream& stream, const Matrice3x3& mat);
};

Matrice3x3 product_vectvectT(Vector3& v1, Vector3& v2);

class Matrice4x4
{
	private:
		double data[16];
		
	public:
		Matrice4x4();
		Matrice4x4(double a00, double a01, double a02, double a03,
			   double a10, double a11, double a12, double a13,
			   double a20, double a21, double a22, double a23,
			   double a30, double a31, double a32, double a33);
		Matrice4x4(double* d);
		Matrice4x4(const vector<double>& d);
		Matrice4x4(const Matrice4x4& mat);
		
		Matrice4x4& operator=(const Matrice4x4& mat);
		
		double& operator() (unsigned row, unsigned col);
		double  operator() (unsigned row, unsigned col) const;
		
		double det();
		Matrice4x4 inv();
		//void inv(const Matrice4x4& mat);
		
		Matrice4x4 operator*(double cst);
		Vector3 operator*(const Vector3& vect) const;
		
		bool operator==(const Matrice4x4& vect) const;
		
		friend ostream& operator<<(ostream& stream, const Matrice4x4& mat);
};

Matrice3x3 outerProduct(const Vector3& vect1, const Vector3& vect2);

class Matrice3x3_complex
{
	private:
		complex<double> data[9];
		
	public:
		Matrice3x3_complex();
		Matrice3x3_complex(complex<double> a00, complex<double> a01, complex<double> a02,
				   complex<double> a10, complex<double> a11, complex<double> a12,
				   complex<double> a20, complex<double> a21, complex<double> a22);
		Matrice3x3_complex(const vector< complex<double> >& d);
		
		complex<double> operator() (unsigned row, unsigned col) const;
		
		Matrice3x3_complex operator*(complex<double> cst);
		
		complex<double> det();
		Matrice3x3_complex inv();
		
		friend ostream& operator<<(ostream& stream, const Matrice3x3_complex& mat);
};

class Matrice4x4_complex
{
	private:
		complex<double> data[16];
		
	public:
		Matrice4x4_complex();
		Matrice4x4_complex(complex<double> a00, complex<double> a01, complex<double> a02, complex<double> a03,
				   complex<double> a10, complex<double> a11, complex<double> a12, complex<double> a13,
				   complex<double> a20, complex<double> a21, complex<double> a22, complex<double> a23,
				   complex<double> a30, complex<double> a31, complex<double> a32, complex<double> a33);
		Matrice4x4_complex(const vector< complex<double> >& d);
		
		complex<double> operator() (unsigned row, unsigned col) const;
		
		Matrice4x4_complex operator*(complex<double> cst);
		
		complex<double> det();
		Matrice4x4_complex inv();
		
		friend ostream& operator<<(ostream& stream, const Matrice4x4_complex& mat);
};

template <typename T> class Point3
{
	private:
		T X;
		T Y;
		T Z;
		
	public:
		Point3()
		{
			X = 0;
			Y = 0;
			Z = 0;
		};
		
		Point3(T x, T y, T z)
		{
			X = x;
			Y = y;
			Z = z;
		};
		
		Point3(const Vector3& v)
		{
			X = (T)v[0];
			Y = (T)v[1];
			Z = (T)v[2];
		};
		
		Point3(const Point3<T>& point)
		{
			X = point.x();
			Y = point.y();
			Z = point.z();
		};
		
		T x() const { return X; };
		T y() const { return Y; };
		T z() const { return Z; };
		
		void round_coords()
		{
			X = round(X);
			Y = round(Y);
			Z = round(Z);
		}
		
		double normDistance2(const Point3<T>& point)
		{
			return (X-point.x()) * (X-point.x()) + (Y-point.y()) * (Y-point.y()) + (Z-point.z()) * (Z-point.z());
		};
		
		double distanceTo(const Point3<T>& point)
		{
			return sqrt( pow(X - point.x(), 2) + pow(Y - point.y(), 2) + pow(Z - point.z(), 2) );
		};
		
		Vector3 vectTo(const Point3<T>& point) const
		{
			double dist_x = (double)(point.x() - X);
			double dist_y = (double)(point.y() - Y);
			double dist_z = (double)(point.z() - Z);
			
			return Vector3(dist_x, dist_y, dist_z);
		};
		
		Point3<T> operator+(const Point3<T>& point) const
		{
			return Point3<T>(X + point.x(),
					 Y + point.y(),
					 Z + point.z());
		}
		
		Vector3 operator-(const Point3<T>& point) const
		{
			double dist_x = (double)(X - point.x());
			double dist_y = (double)(Y - point.y());
			double dist_z = (double)(Z - point.z());
			
			return Vector3(dist_x, dist_y, dist_z);
		}
		
		Point3<T> operator/(float coef) const
		{
			T dist_x = (T)(X / coef);
			T dist_y = (T)(Y / coef);
			T dist_z = (T)(Z / coef);
			
			return Point3<T>(dist_x, dist_y, dist_z);
		}
		
		Point3<T> operator*(float coef) const
		{
			T dist_x = (T)(X * coef);
			T dist_y = (T)(Y * coef);
			T dist_z = (T)(Z * coef);
			
			return Point3<T>(dist_x, dist_y, dist_z);
		}
		
		bool operator==(const Point3<T>& point) const
		{
			if (point.x() == X && point.y() == Y && point.z() == Z) return true;
			else return false;
		};
		
		bool operator!=(const Point3<T>& point) const
		{
			if (point.x() == X && point.y() == Y && point.z() == Z) return false;
			else return true;
		};
		
		Point3<T>& operator=(const Point3<T>& point)
		{
			X = point.x();
			Y = point.y();
			Z = point.z();
			
			return *this;
		}
		
		Point3<T>& operator+=(const Point3<T>& point)
		{
			X += point.x();
			Y += point.y();
			Z += point.z();
			
			return *this;
		}
		
		Point3<T>& operator-=(const Point3<T>& point)
		{
			X -= point.x();
			Y -= point.y();
			Z -= point.z();
			
			return *this;
		}
		
		Point3<T>& operator/=(float coef)
		{
			X /= coef;
			Y /= coef;
			Z /= coef;
			
			return *this;
		}
		
		bool operator<(const Point3<T>& point) const
		{
			/*if (Z < point.z()) return true;
			else
				if (Z > point.z()) return false;
				else
				{
					if (Y < point.y()) return true;
					else
						if (Y > point.y()) return false;
						else
						{
							if (X < point.x()) return true;
							else return false;
						}
				}
			*/
			
			if (X < point.x()) return true;
			else if (X > point.x()) return false;
			else
			{
				if (Y < point.y()) return true;
				else if (Y > point.y()) return false;
				else
				{
					if (Z < point.z()) return true;
					else return false;
				}
			}
		}
};

template <typename T> ostream& operator<<(ostream& stream, const Point3<T>& point)
{
	return stream << "(" << point.x() << ", " << point.y() << ", " << point.z() << ")";
}

template <typename T> class Point4
{
	private:
		T X;
		T Y;
		T Z;
		T M;
		
	public:
		Point4()
		{
			X = 0;
			Y = 0;
			Z = 0;
			M = 0;
		};
		
		Point4(T x, T y, T z, T m)
		{
			X = x;
			Y = y;
			Z = z;
			M = m;
		};

		Point4(const Vector4& v)
		{
			X = (T)v[0];
			Y = (T)v[1];
			Z = (T)v[2];
			M = (T)v[3];
		};
		
		Point4(const Point4<T>& point)
		{
			X = point.x();
			Y = point.y();
			Z = point.z();
			M = point.t();
		};
		
		T x() const { return X; };
		T y() const { return Y; };
		T z() const { return Z; };
		T t() const { return M; };
		
		double normDistance2(const Point4<T>& point)
		{
			return	(X-point.x()) * (X-point.x()) +
				(Y-point.y()) * (Y-point.y()) +
				(Z-point.z()) * (Z-point.z()) +
				(M-point.m()) * (M-point.m());
		};

		double distanceTo(const Point4<T>& point)
		{
			return sqrt( (X-point.x()) * (X-point.x()) +
				(Y-point.y()) * (Y-point.y()) +
				(Z-point.z()) * (Z-point.z()) +
				(M-point.m()) * (M-point.m()) );
		};
		
		Vector4 operator-(const Point4<T>& point) const
		{
			double dist_x = (double)(X - point.x());
			double dist_y = (double)(Y - point.y());
			double dist_z = (double)(Z - point.z());
			double dist_m = (double)(M - point.m());
			
			return Vector4(dist_x, dist_y, dist_z, dist_m);
		}
		
		Point4<T> operator/(float coef) const
		{
			T dist_x = (T)(X / coef);
			T dist_y = (T)(Y / coef);
			T dist_z = (T)(Z / coef);
			T dist_m = (T)(M / coef);
			
			return Point4<T>(dist_x, dist_y, dist_z, dist_m);
		}
		
		bool operator==(const Point4<T>& point) const
		{
			if (point.x() == X && point.y() == Y && point.z() == Z && point.m() == M) return true;
			else return false;
		};
		
		bool operator!=(const Point4<T>& point) const
		{
			if (point.x() == X && point.y() == Y && point.z() == Z && point.m() == M) return false;
			else return true;
		};
		
		Point4<T>& operator+=(const Point4<T>& point)
		{
			X += point.x();
			Y += point.y();
			Z += point.z();
			M += point.m();
			
			return *this;
		}
		
		Point4<T>& operator-=(const Point4<T>& point)
		{
			X -= point.x();
			Y -= point.y();
			Z -= point.z();
			M -= point.m();
			
			return *this;
		}
		
		bool operator<(const Point4<T>& point) const
		{
			if (M < point.t()) return true;
			else
				if (M > point.t()) return false;
				else
				{
					if (Z < point.z()) return true;
					else
						if (Z > point.z()) return false;
						else
						{
							if (Y < point.y()) return true;
							else
								if (Y > point.y()) return false;
								else
								{
									if (X < point.x()) return true;
									else return false;
								}
						}
				}
		}
};

template <typename T> ostream& operator<<(ostream& stream, const Point4<T>& point)
{
	return stream << "(" << point.x() << ", " << point.y() << ", " << point.z() << ", " << point.t() << ")";
}

inline float dirac(float x, float alpha);

float dirac(float x, float alpha)
{
	return 1. / ( (M_PI * alpha) * (1 + x*x / (alpha*alpha) ) );
}

inline int computeIndex(const pair<int, int>& size, int r, int c);
inline int computeIndex(const int* size, int x, int y, int z);
inline int computeIndex(const int* size, int x, int y, int z, int m);

inline void computeCoords(int indice, const pair<int, int>& size, int& r, int& c);
inline void computeCoords(int indice, const int* size, int& x, int& y, int& z);
inline void computeCoords(int indice, const int* size, int& x, int& y, int& z, int& m);

inline int computeIndex(const pair<int, int>& size, int r, int c)
{
	if (r < 0 || r >= size.first || c < 0 || c >= size.second)
	{
		cout << "point (" << r << ", " << c << ") et taille volume = " << size.first << " x " << size.second << endl;
		vector<int> test;
		cout << test[10] << endl;
		throw( string("Erreur dans computeIndex : coordonnees hors du volume !") );
	}
	
	return c + size.second * r;
}

inline int computeIndex(const int* size, int x, int y, int z)
{
	if (x < 0 || x >= size[0] || y < 0 || y >= size[1] || z < 0 || z >= size[2])
	{
		cout << "point " << Point3<int>(x, y, z) << " et taille volume = " << size[0] << " " << size[1] << " " << size[2] << endl;
		vector<int> test;
		cout << test[10] << endl;
		throw( string("Erreur dans computeIndex : coordonnees hors du volume !") );
	}
	
	return x + size[0] * (y + size[1] * z);
}

inline int computeIndex(const int* size, int x, int y, int z, int m)
{
	if (x < 0 || x >= size[0] || y < 0 || y >= size[1] || z < 0 || z >= size[2] || m < 0 || m >= size[3])
	{
		cout << "point " << Point4<int>(x, y, z, m) << " et taille volume = " << size[0] << " " << size[1] << " " << size[2] << " " << size[3] << endl;
		vector<int> test;
		cout << test[10] << endl;
		throw( string("Erreur dans computeIndex : coordonnees hors du volume !") );
	}
	
	return x + size[0] * (y + size[1] * (z + size[2] * m));
}


inline void computeCoords(int indice, const pair<int, int>& size, int& r, int& c)
{
	if (indice < 0 || indice >= size.first * size.second)
		throw( string("Erreur dans computeCoords : coordonnees hors du volume !") );
	
	r = (int)floor( (double)indice / (double)size.second );
	c = indice - r * size.second;
}

inline void computeCoords(int indice, const int* size, int& x, int& y, int& z)
{
	if (indice < 0 || indice >= size[0] * size[1] * size[2])
		throw( string("Erreur dans computeCoords : coordonnees hors du volume !") );
	
	z = (int)floor((double)indice / (double)(size[0] * size[1]));
	y = (int)floor((double)(indice - z * size[0] * size[1]) / (double)size[0]);
	x = indice - z * size[0] * size[1] - y * size[0];
}

inline void computeCoords(int indice, const int* size, int& x, int& y, int& z, int& m)
{
	if (indice < 0 || indice >= size[0] * size[1] * size[2] * size[3])
		throw( string("Erreur dans computeCoords : coordonnees hors du volume !") );
	
	m = (int)floor((double)indice / (double)(size[0] * size[1] * size[2]));
	z = (int)floor((double)(indice - m * size[0] * size[1] * size[2]) / (double)(size[0] * size[1]));
	y = (int)floor((double)(indice - m * size[0] * size[1] * size[2] - z * size[0] * size[1]) / (double)size[0]);
	x = indice - m * size[0] * size[1] * size[2] - z * size[0] * size[1] - y * size[0];
}

void calculNeighboursCoords(int x, int y, int z, int m, int* xN, int* yN, int* zN, int* mN);

#endif
