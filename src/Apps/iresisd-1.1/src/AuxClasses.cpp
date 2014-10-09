#include "AuxClasses.h"

using namespace std;

void calculNeighboursCoords(int x, int y, int z, int m, int* xN, int* yN, int* zN, int* mN)
{
	int xNtmp[80] = {
		x-1, x, x+1, x-1, x, x+1, x-1, x, x+1, x-1, x, x+1, x-1, x, x+1, x-1, x, x+1, x-1, x, x+1, x-1, x, x+1, x-1, x, x+1,
		x-1, x, x+1, x-1, x, x+1, x-1, x, x+1, x-1, x, x+1, x-1, x+1, x-1, x, x+1, x-1, x, x+1, x-1, x, x+1, x-1, x, x+1,
		x-1, x, x+1, x-1, x, x+1, x-1, x, x+1, x-1, x, x+1, x-1, x, x+1, x-1, x, x+1, x-1, x, x+1, x-1, x, x+1, x-1, x, x+1};
	
	int yNtmp[80] = {
		y-1, y-1, y-1, y, y, y, y+1, y+1, y+1, y-1, y-1, y-1, y, y, y, y+1, y+1, y+1, y-1, y-1, y-1, y, y, y, y+1, y+1, y+1,
		y-1, y-1, y-1, y, y, y, y+1, y+1, y+1, y-1, y-1, y-1, y, y, y+1, y+1, y+1, y-1, y-1, y-1, y, y, y, y+1, y+1, y+1,
		y-1, y-1, y-1, y, y, y, y+1, y+1, y+1, y-1, y-1, y-1, y, y, y, y+1, y+1, y+1, y-1, y-1, y-1, y, y, y, y+1, y+1, y+1};
	
	int zNtmp[80] = {
		z-1, z-1, z-1, z-1, z-1, z-1, z-1, z-1, z-1, z, z, z, z, z, z, z, z, z, z+1, z+1, z+1, z+1, z+1, z+1, z+1, z+1, z+1,
		z-1, z-1, z-1, z-1, z-1, z-1, z-1, z-1, z-1, z, z, z, z, z, z, z, z, z+1, z+1, z+1, z+1, z+1, z+1, z+1, z+1, z+1,
		z-1, z-1, z-1, z-1, z-1, z-1, z-1, z-1, z-1, z, z, z, z, z, z, z, z, z, z+1, z+1, z+1, z+1, z+1, z+1, z+1, z+1, z+1};
	
	int mNtmp[80] = {
		m-1, m-1, m-1, m-1, m-1, m-1, m-1, m-1, m-1, m-1, m-1, m-1, m-1, m-1, m-1, m-1, m-1, m-1, m-1, m-1, m-1, m-1, m-1, m-1, m-1, m-1, m-1,
		m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m,
		m+1, m+1, m+1, m+1, m+1, m+1, m+1, m+1, m+1, m+1, m+1, m+1, m+1, m+1, m+1, m+1, m+1, m+1, m+1, m+1, m+1, m+1, m+1, m+1, m+1, m+1, m+1};
	
	for (int i=0; i<80; i++)
	{
		xN[i] = xNtmp[i];
		yN[i] = yNtmp[i];
		zN[i] = zNtmp[i];
		mN[i] = mNtmp[i];
	}
}

/** Vector3 **/

Vector3::Vector3()
{
	a = 0;
	b = 0;
	c = 0;
	//isnull = true;
}

Vector3::Vector3(double A, double B, double C)
{
	a = A;
	b = B;
	c = C;
	//isnull = false;
}

Vector3::Vector3(double* v)
{
	a = v[0];
	b = v[1];
	c = v[2];
	//isnull = false;
}

Vector3::Vector3(const Vector3& vect)
{
	a = vect[0];
	b = vect[1];
	c = vect[2];
	//isnull = vect.isNull();
}

Vector3& Vector3::operator=(const Vector3& vect)
{
	a = vect[0];
	b = vect[1];
	c = vect[2];
	//isnull = vect.isNull();
	
	return *this;
}

bool Vector3::operator!=(const Vector3& vect) const
{
	if (a != vect[0] || b != vect[1] || c != vect[2]) return true;
	else return false;
}

bool Vector3::operator==(const Vector3& vect) const
{
	if (a == vect[0] && b == vect[1] && c == vect[2]) return true;
	else return false;
}

/*bool Vector3::isNull() const
{
	//return isnull;
	if (!a && !b && !c) return true;
	else return false;
}*/

double& Vector3::operator[](int i)
{
	switch(i)
	{
		case 0:
			//isnull = false;
			return a;
			
		case 1:
			//isnull = false;
			return b;
			
		case 2:
			//isnull = false;
			return c;
			
		default:
			throw(string("Error: indice out of bound"));
	}
}

double Vector3::operator[](int i) const
{
	switch(i)
	{
		case 0:
			return a;
			
		case 1:
			return b;
			
		case 2:
			return c;
			
		default:
			throw(string("Error: indice out of bound"));
	}
}

double Vector3::operator*(const Vector3& vect) const
{
	return	a * vect[0] +
		b * vect[1] +
		c * vect[2];
}

Vector3& Vector3::operator+=(const Vector3& vect)
{
	a += vect[0];
	b += vect[1];
	c += vect[2];
	
	//isnull = false;
	
	return *this;
}

Vector3& Vector3::operator-=(const Vector3& vect)
{
	a -= vect[0];
	b -= vect[1];
	c -= vect[2];
	
	//isnull = false;
	
	return *this;
}

Vector3& Vector3::operator*=(double fac)
{
	a *= fac;
	b *= fac;
	c *= fac;
	
	return *this;
}

Vector3& Vector3::operator/=(double div)
{
	if (fabs(div) > 0.0001)
	{
		a /= div;
		b /= div;
		c /= div;
	
		//isnull = false;
	}
	else
	{
		cout << div << " Vector3::operator/= : Warning: dividing by zero! -> Skiping" << endl;
		//cv::waitKey();
	}
	
	return *this;
}

Vector3 Vector3::operator-() const
{
    //cout << "operator-()" << endl;

	return Vector3(-a, -b, -c);
}

Vector3 Vector3::operator-(const Vector3& vect) const
{
    //cout << "operator-(const Vector3& vect)" << endl;

	return Vector3(a - vect[0], b - vect[1], c - vect[2]);
}

Vector3 Vector3::operator+(const Vector3& vect) const
{
	return Vector3(a + vect[0], b + vect[1], c + vect[2]);
}

Vector3 Vector3::operator*(double cst) const
{
	return Vector3(a*cst, b*cst, c*cst);
}

Vector3 Vector3::operator/(double div) const
{
	double aTmp, bTmp, cTmp;
	
	if (div)
	{
		aTmp = a / div;
		bTmp = b / div;
		cTmp = c / div;
		
		return Vector3(aTmp, bTmp, cTmp);
	}
	else
		cout << "Warning: dividing by zero! -> Skiping" << endl;
	
	return *this;
}

Vector3 Vector3::operator/(Vector3 div) const
{
	double aTmp, bTmp, cTmp;
	
	if (div!=VECT_NULL)
	{
		aTmp = a / div[0];
		bTmp = b / div[1];
		cTmp = c / div[2];
		
		return Vector3(aTmp, bTmp, cTmp);
	}
	else
		cout << "Warning: dividing by zero! -> Skiping" << endl;
	
	return *this;
}

double Vector3::norm2() const
{
	return a*a + b*b + c*c;
}

double Vector3::norm() const
{
	double _norm2 = norm2();
	
	if (_norm2 == 0) return 0;
	else return sqrt(_norm2);
}

ostream& operator<<(ostream& stream, const Vector3& vect)
{
	return stream << "[" << vect[0] << ", " << vect[1] << ", " << vect[2] << "]";
}

/** Vector4 **/

Vector4::Vector4()
{
	a = 0;
	b = 0;
	c = 0;
	d = 0;
}

Vector4::Vector4(double A, double B, double C, double D)
{
	a = A;
	b = B;
	c = C;
	d = D;
}

Vector4::Vector4(double* v)
{
	a = v[0];
	b = v[1];
	c = v[2];
	d = v[3];
}

Vector4::Vector4(const Vector4& vect)
{
	a = vect[0];
	b = vect[1];
	c = vect[2];
	d = vect[3];
}

Vector4& Vector4::operator=(const Vector4& vect)
{
	a = vect[0];
	b = vect[1];
	c = vect[2];
	d = vect[3];
	
	return *this;
}

bool Vector4::operator!=(const Vector4& vect) const
{
	if (a != vect[0] || b != vect[1] || c != vect[2] || d != vect[3]) return true;
	else return false;
}

bool Vector4::operator==(const Vector4& vect) const
{
	if (a == vect[0] && b == vect[1] && c == vect[2] && d == vect[3]) return true;
	else return false;
}

double& Vector4::operator[](int i)
{
	switch(i)
	{
		case 0:
			return a;
			
		case 1:
			return b;
			
		case 2:
			return c;
			
		case 3:
			return d;
			
		default:
			throw(string("Error: indice out of bound"));
	}
}

double Vector4::operator[](int i) const
{
	switch(i)
	{
		case 0:
			return a;
			
		case 1:
			return b;
			
		case 2:
			return c;
			
		case 3:
			return d;
			
		default:
			throw(string("Error: indice out of bound"));
	}
}

double Vector4::operator*(const Vector4& vect) const
{
	return	a * vect[0] +
		b * vect[1] +
		c * vect[2] +
		d * vect[3];
}

Vector4& Vector4::operator+=(const Vector4& vect)
{
	a += vect[0];
	b += vect[1];
	c += vect[2];
	d += vect[3];
	
	return *this;
}

Vector4& Vector4::operator-=(const Vector4& vect)
{
	a -= vect[0];
	b -= vect[1];
	c -= vect[2];
	d -= vect[3];
	
	return *this;
}

Vector4& Vector4::operator*=(double fac)
{
	a *= fac;
	b *= fac;
	c *= fac;
	d *= fac;
	
	return *this;
}

Vector4& Vector4::operator/=(double div)
{
	if (fabs(div) > 0.0001)
	{
		a /= div;
		b /= div;
		c /= div;
		d /= div;
	}
	else
	{
		cout << div << " Vector3::operator/= : Warning: dividing by zero! -> Skiping" << endl;
		//cv::waitKey();
	}
	
	return *this;
}

Vector4 Vector4::operator-() const
{
	return Vector4(-a, -b, -c, -d);
}

Vector4 Vector4::operator-(const Vector4& vect) const
{
    //cout << "operator-(const Vector3& vect)" << endl;

	return Vector4(a - vect[0], b - vect[1], c - vect[2], d - vect[3]);
}

Vector4 Vector4::operator+(const Vector4& vect) const
{
	return Vector4(a + vect[0], b + vect[1], c + vect[2], d + vect[3]);
}

Vector4 Vector4::operator*(double cst) const
{
	return Vector4(a*cst, b*cst, c*cst, d*cst);
}

Vector4 Vector4::operator/(double div) const
{
	double aTmp, bTmp, cTmp, dTmp;
	
	if (div)
	{
		aTmp = a / div;
		bTmp = b / div;
		cTmp = c / div;
		dTmp = d / div;
		
		return Vector4(aTmp, bTmp, cTmp, dTmp);
	}
	else
		cout << "Warning: dividing by zero! -> Skiping" << endl;
	
	return *this;
}

double Vector4::norm2() const
{
	return a*a + b*b + c*c + d*d;
}

double Vector4::norm() const
{
	double _norm2 = norm2();
	
	if (_norm2 == 0) return 0;
	else return sqrt(_norm2);
}

bool Vector4::operator<(const Vector4& vect) const
{
	if (a < vect[0]) return true;
	if (b < vect[1]) return true;
	if (c < vect[2]) return true;
	if (d < vect[3]) return true;
	return false;
}

ostream& operator<<(ostream& stream, const Vector4& vect)
{
	return stream << "[" << vect[0] << ", " << vect[1] << ", " << vect[2] << ", " << vect[3] << "]";
}


/** VectorN **/

VectorN::VectorN()
{
	dimension = 0;
	val = NULL;
	//isnull = true;
}

VectorN::VectorN(int d, double* v)
{
	dimension = d;
	
	val = new double[dimension];
	for (int i=0; i<dimension; i++) val[i] = v[i];
	
	//isnull = false;
}

VectorN::VectorN(const VectorN& vect)
{
	dimension = vect.dim();
	
	val = new double[dimension];
	for (int i=0; i<dimension; i++) val[i] = vect[i];
	
	//isnull = vect.isNull();
}

VectorN::~VectorN()
{
	if (val != NULL) delete val;
}

int VectorN::dim() const
{
	return dimension;
}

VectorN& VectorN::operator=(const VectorN& vect)
{
	if (dimension != vect.dim())
		if (val != NULL) delete val;
	
	dimension = vect.dim();
	for (int i=0; i<dimension; i++) val[i] = vect[i];
	
	//isnull = vect.isNull();
	
	return *this;
}

bool VectorN::operator!=(const VectorN& vect) const
{
	if (dimension != vect.dim()) throw string("VectorN::operator!= Error : operation sur des vecteurs de dimensions differentes");
	
	for (int i=0; i<dimension; i++)
		if (val[i] != vect[i]) return true;
	
	return false;
}

bool VectorN::operator==(const VectorN& vect) const
{
	if (dimension != vect.dim()) throw string("VectorN::operator== Error : operation sur des vecteurs de dimensions differentes");
	
	for (int i=0; i<dimension; i++)
		if (val[i] != vect[i]) return false;
	
	return true;
}

/*bool Vector3::isNull() const
{
	//return isnull;
	if (!a && !b && !c) return true;
	else return false;
}*/

double& VectorN::operator[](unsigned int i)
{
	if (i >= dimension) throw string("VectorN::operator[] Error : indice depasse la taille du vecteur");
	
	return val[i];
}

double VectorN::operator[](unsigned int i) const
{
	if (i >= dimension) throw string("VectorN::operator[] Error : indice depasse la taille du vecteur");
	
	return val[i];
}

double VectorN::operator*(const VectorN& vect) const
{
	if (dimension != vect.dim()) throw string("VectorN::operator* Error : operation sur des vecteurs de dimensions differentes");
	
	double res = 0;
	for (int i=0; i<dimension; i++) res += val[i] * vect[i];
	
	return res;
}

VectorN& VectorN::operator+=(const VectorN& vect)
{
	if (dimension != vect.dim()) throw string("VectorN::operator+= Error : operation sur des vecteurs de dimensions differentes");
	
	for (int i=0; i<dimension; i++) val[i] += vect[i];
	
	//isnull = false;
	
	return *this;
}

VectorN& VectorN::operator-=(const VectorN& vect)
{
	if (dimension != vect.dim()) throw string("VectorN::operator-= Error : operation sur des vecteurs de dimensions differentes");
	
	for (int i=0; i<dimension; i++) val[i] -= vect[i];
	
	//isnull = false;
	
	return *this;
}

VectorN& VectorN::operator*=(double fac)
{
	for (int i=0; i<dimension; i++) val[i] *= fac;
	
	return *this;
}

VectorN& VectorN::operator/=(double div)
{
	if (fabs(div) > 0.0001)
	{
		for (int i=0; i<dimension; i++) val[i] /= div;
	
		//isnull = false;
	}
	else
		cout << "Warning: dividing by zero! -> Skiping" << endl;
	
	return *this;
}

VectorN VectorN::operator-() const
{
	VectorN res(*this);
	res *= -1;
	
	return res;
}

VectorN VectorN::operator-(const VectorN& vect) const
{
	if (dimension != vect.dim()) throw string("VectorN::operator- Error : operation sur des vecteurs de dimensions differentes");
	
	VectorN res(*this);
	for (int i=0; i<dimension; i++) res[i] -= vect[i];
	
	return res;
}

VectorN VectorN::operator+(const VectorN& vect) const
{
	if (dimension != vect.dim()) throw string("VectorN::operator+ Error : operation sur des vecteurs de dimensions differentes");
	
	VectorN res(*this);
	for (int i=0; i<dimension; i++) res[i] += vect[i];
	
	return res;
}

VectorN VectorN::operator*(double cst) const
{
	VectorN res(*this);
	res *= cst;
	
	return res;
}

VectorN VectorN::operator/(double div) const
{
	if (div)
	{
		VectorN res(*this);
		res /= div;
		
		return res;
	}
	else
		cout << "Warning: dividing by zero! -> Skiping" << endl;
	
	return *this;
}

double VectorN::norm2() const
{
	double tmp = 0;
	for (int i=0; i<dimension; i++) tmp += val[i] * val[i];
	
	return tmp;
}

double VectorN::norm() const
{
	double tmp = 0;
	for (int i=0; i<dimension; i++) tmp += val[i] * val[i];
	
	return sqrt(tmp);
}

ostream& operator<<(ostream& stream, const VectorN& vect)
{
	stream << "[";
	
	for (int i=0; i<vect.dim()-1; i++) stream << vect[i] << " ";
	
	return stream << vect[vect.dim()-1] << "]";
}

/** Matrice3x3 **/

Matrice3x3::Matrice3x3()
{
	for (int i=0; i<9; i++)
		data[i] = 0;
}

Matrice3x3::Matrice3x3(double a00, double a01, double a02,
		       double a10, double a11, double a12,
		       double a20, double a21, double a22)
{
	data[0] = a00;
	data[1] = a01;
	data[2] = a02;
	data[3] = a10;
	data[4] = a11;
	data[5] = a12;
	data[6] = a20;
	data[7] = a21;
	data[8] = a22;
}

Matrice3x3::Matrice3x3(Vector3 col0, Vector3 col1, Vector3 col2)
{
	data[0] = col0[0];
	data[3] = col0[1];
	data[6] = col0[2];
	data[1] = col1[0];
	data[4] = col1[1];
	data[7] = col1[2];
	data[2] = col2[0];
	data[5] = col2[1];
	data[8] = col2[2];
}

Matrice3x3::Matrice3x3(double* d)
{
	try
	{
		for (int i=0; i<9; i++)
			data[i] = d[i];
	}
	catch(...)
	{
		throw("Error: Please, check array size (must be 9)");
	}
}

Matrice3x3::Matrice3x3(const Matrice3x3& mat)
{
	for (int r=0; r<3; r++)
	for (int c=0; c<3; c++)
		data[c + 3 * r] = mat(r, c);
}

Matrice3x3& Matrice3x3::operator=(const Matrice3x3& mat)
{
	for (int r=0; r<3; r++)
	for (int c=0; c<3; c++)
		data[c + 3 * r] = mat(r, c);
	
	return *this;
}

double& Matrice3x3::operator() (unsigned row, unsigned col)
{
	if (row >= 3 || col >= 3)
		throw string("Error: matrix subscript out of bounds");
	
	return data[col + 3 * row];
}

double Matrice3x3::operator() (unsigned row, unsigned col) const
{
	if (row >= 3 || col >= 3)
		throw string("Error: matrix subscript out of bounds");
	
	return data[col + 3 * row];
}

double Matrice3x3::det()
{
	return 	data[0] * (data[4] * data[8] - data[5] * data[7]) -
		data[1] * (data[3] * data[8] - data[5] * data[6]) +
		data[2] * (data[3] * data[7] - data[4] * data[6]);
}

Matrice3x3 Matrice3x3::inv()
{
	return Matrice3x3(data[4] * data[8] - data[5] * data[7], data[2] * data[7] - data[1] * data[8], data[1] * data[5] - data[2] * data[4],
			  data[5] * data[6] - data[3] * data[8], data[0] * data[8] - data[2] * data[6], data[2] * data[3] - data[0] * data[5],
			  data[3] * data[7] - data[4] * data[6], data[1] * data[6] - data[0] * data[7], data[0] * data[4] - data[1] * data[3])
			* (1. / det());
	
	/*//Gauss-Jordan inversion algorithm
	
	int r, c;
	double ligne[6];
	double diviseur, multiplieur;
	
	double matTmp[18];
	for (r=0; r<3; r++)
	{
		for (c=0; c<3; c++) matTmp[c + 6 * r] = data[c + 3 * r];
		for (c=3; c<6; c++) matTmp[c + 6 * r] = 0;
	}
	
	matTmp[3] = 1;
	matTmp[10] = 1;
	matTmp[17] = 1;
	
	for (int k=0; k<3; k++)
	{
		for (r=k; r<3; r++)
		{
			if (matTmp[k + 6 * r]) break;
			
			if (r == 2) throw string("Error: matrice not invertible!");
		}
		
		if (r != k)
		{
			for (c=0; c<6; c++) ligne[c] = matTmp[c + 6 * k];
			for (c=0; c<6; c++) matTmp[c + 6 * k] = matTmp[c + 6 * r];
			for (c=0; c<6; c++) matTmp[c + 6 * r] = ligne[c];
		}
		
		diviseur = matTmp[k + 6 * k];
		for (c=0; c<6; c++) matTmp[c + 6 * k] /= diviseur;
		
		for (int i=0; i<3; i++)
		{
			if (i == k) continue;
			
			multiplieur = matTmp[k + 6 * i];
			for (c=0; c<6; c++) matTmp[c + 6 * i] -= multiplieur * matTmp[c + 6 * k];
		}
	}
	
	double matTmp2[9];
	for (r=0; r<3; r++)
	for (c=0; c<3; c++)
		matTmp2[c + 3 * r] = matTmp[c+3 + 6 * r];
	
	return Matrice3x3(matTmp2);*/
}

Matrice3x3 Matrice3x3::operator*(double cst)
{
	double tmp[9];
	for (int i=0; i<9; i++) tmp[i] = cst * data[i];
	
	return Matrice3x3(tmp);
}

Vector3 Matrice3x3::operator*(const Vector3& vect)
{
	return Vector3(vect[0]*data[0] + vect[1]*data[1] + vect[2]*data[2],
		       vect[0]*data[3] + vect[1]*data[4] + vect[2]*data[5],
		       vect[0]*data[6] + vect[1]*data[7] + vect[2]*data[8]);
}

Matrice3x3& Matrice3x3::operator+=(const Matrice3x3& mat)
{
	for (int col=0; col<3; col++)
	for (int row=0; row<3; row++)
		data[col + 3 * row] += mat(row, col);
	
	return *this;
}

ostream& operator<<(ostream& stream, const Matrice3x3& mat)
{
	for (int r=0; r<3; r++)
	{
		for (int c=0; c<3; c++)
			stream << mat(r, c) << " ";
	
		stream << endl;
	}
	
	return stream;
}

Matrice3x3 product_vectvectT(Vector3& v1, Vector3& v2)
{
	return Matrice3x3(v1[0]*v2[0], v1[0]*v2[1], v1[0]*v2[2],
			  v1[1]*v2[0], v1[1]*v2[1], v1[1]*v2[2],
			  v1[2]*v2[0], v1[2]*v2[1], v1[2]*v2[2]);
}

/** Matrice4x4 **/

Matrice4x4::Matrice4x4()
{
	for (int i=0; i<16; i++)
		data[i] = 0;
}

Matrice4x4::Matrice4x4(double a00, double a01, double a02, double a03,
		       double a10, double a11, double a12, double a13,
		       double a20, double a21, double a22, double a23,
		       double a30, double a31, double a32, double a33)
{
	data[0] = a00;
	data[1] = a01;
	data[2] = a02;
	data[3] = a03;
	data[4] = a10;
	data[5] = a11;
	data[6] = a12;
	data[7] = a13;
	data[8] = a20;
	data[9] = a21;
	data[10] = a22;
	data[11] = a23;
	data[12] = a30;
	data[13] = a31;
	data[14] = a32;
	data[15] = a33;
}

Matrice4x4::Matrice4x4(double* d)
{
	try
	{
		for (int i=0; i<16; i++) data[i] = d[i];
	}
	catch(...)
	{
		throw("Error: Please, check array size (must be 16)");
	}
}

Matrice4x4::Matrice4x4(const vector<double>& d)
{
	if (d.size() != 16) throw("Error: Please, check array size (must be 16)");
	
	for (int i=0; i<16; i++) data[i] = d[i];
}

Matrice4x4::Matrice4x4(const Matrice4x4& mat)
{
	for (int r=0; r<4; r++)
	for (int c=0; c<4; c++)
		data[c + 4 * r] = mat(r, c);
}

Matrice4x4& Matrice4x4::operator=(const Matrice4x4& mat)
{
	for (int r=0; r<4; r++)
	for (int c=0; c<4; c++)
		data[c + 4 * r] = mat(r, c);
	
	return *this;
}

double& Matrice4x4::operator() (unsigned row, unsigned col)
{
	if (row >= 4 || col >= 4)
		throw string("Error: matrix subscript out of bounds");
	
	return data[col + 4 * row];
}

double Matrice4x4::operator() (unsigned row, unsigned col) const
{
	if (row >= 4 || col >= 4)
		throw string("Error: matrix subscript out of bounds");
	
	return data[col + 4 * row];
}

double Matrice4x4::det()
{
	Matrice3x3 minor0(data[5], data[6], data[7],
			  data[9], data[10], data[11],
			  data[13], data[14], data[15]);
	
	Matrice3x3 minor1(data[4], data[6], data[7],
			  data[8], data[10], data[11],
			  data[12], data[14], data[15]);
	
	Matrice3x3 minor2(data[4], data[5], data[7],
			  data[8], data[9], data[11],
			  data[12], data[13], data[15]);
	
	Matrice3x3 minor3(data[4], data[5], data[6],
			  data[8], data[9], data[10],
			  data[12], data[13], data[14]);
	
	return 	data[0] * minor0.det() -
		data[1] * minor1.det() +
		data[2] * minor2.det() -
		data[3] * minor3.det();
}

Matrice4x4 Matrice4x4::inv()
{
	Matrice3x3 minor0(data[5], data[6], data[7],
			  data[9], data[10], data[11],
			  data[13], data[14], data[15]);
	
	Matrice3x3 minor1(data[4], data[6], data[7],
			  data[8], data[10], data[11],
			  data[12], data[14], data[15]);
	
	Matrice3x3 minor2(data[4], data[5], data[7],
			  data[8], data[9], data[11],
			  data[12], data[13], data[15]);
	
	Matrice3x3 minor3(data[4], data[5], data[6],
			  data[8], data[9], data[10],
			  data[12], data[13], data[14]);
	
	Matrice3x3 minor4(data[1], data[2], data[3],
			  data[9], data[10], data[11],
			  data[13], data[14], data[15]);
	
	Matrice3x3 minor5(data[0], data[2], data[3],
			  data[8], data[10], data[11],
			  data[12], data[14], data[15]);
	
	Matrice3x3 minor6(data[0], data[1], data[3],
			  data[8], data[9], data[11],
			  data[12], data[13], data[15]);
	
	Matrice3x3 minor7(data[0], data[1], data[2],
			  data[8], data[9], data[10],
			  data[12], data[13], data[14]);
	
	Matrice3x3 minor8(data[1], data[2], data[3],
			  data[5], data[6], data[7],
			  data[13], data[14], data[15]);
	
	Matrice3x3 minor9(data[0], data[2], data[3],
			  data[4], data[6], data[7],
			  data[12], data[14], data[15]);
	
	Matrice3x3 minor10(data[0], data[1], data[3],
			  data[4], data[5], data[7],
			  data[12], data[13], data[15]);
	
	Matrice3x3 minor11(data[0], data[1], data[2],
			  data[4], data[5], data[6],
			  data[12], data[13], data[14]);
	
	Matrice3x3 minor12(data[1], data[2], data[3],
			  data[5], data[6], data[7],
			  data[9], data[10], data[11]);
	
	Matrice3x3 minor13(data[0], data[2], data[3],
			  data[4], data[6], data[7],
			  data[8], data[10], data[11]);
	
	Matrice3x3 minor14(data[0], data[1], data[3],
			   data[4], data[5], data[7],
			   data[8], data[9], data[11]);
	
	Matrice3x3 minor15(data[0], data[1], data[2],
			   data[4], data[5], data[6],
			   data[8], data[9], data[10]);
	
	return Matrice4x4( minor0.det(), - minor4.det(), minor8.det(), - minor12.det(),
			   - minor1.det(), minor5.det(), - minor9.det(), minor13.det(),
			   minor2.det(), - minor6.det(), minor10.det(), - minor14.det(),
			   - minor3.det(), minor7.det(), - minor11.det(), minor15.det() )
			* (1. / det());
	
	/*//Gauss-Jordan inversion algorithm
	
	int r, c;
	double ligne[8];
	double diviseur, multiplieur;
	
	double matTmp[32];
	for (r=0; r<4; r++)
	{
		for (c=0; c<4; c++) matTmp[c + 8 * r] = data[c + 4 * r];
		for (c=4; c<8; c++) matTmp[c + 8 * r] = 0;
	}
	
	matTmp[4] = 1;
	matTmp[13] = 1;
	matTmp[22] = 1;
	matTmp[31] = 1;
	
	for (int k=0; k<4; k++)
	{
		for (r=k; r<4; r++)
		{
			if (matTmp[k + 8 * r]) break;
			
			if (r == 3) throw string("Error: matrice not invertible!");
		}
		
		if (r != k)
		{
			for (c=0; c<8; c++) ligne[c] = matTmp[c + 8 * k];
			for (c=0; c<8; c++) matTmp[c + 8 * k] = matTmp[c + 8 * r];
			for (c=0; c<8; c++) matTmp[c + 8 * r] = ligne[c];
		}
		
		diviseur = matTmp[k + 8 * k];
		for (c=0; c<8; c++) matTmp[c + 8 * k] /= diviseur;
		
		for (int i=0; i<4; i++)
		{
			if (i == k) continue;
			
			multiplieur = matTmp[k + 8 * i];
			for (c=0; c<8; c++) matTmp[c + 8 * i] -= multiplieur * matTmp[c + 8 * k];
		}
	}
	
	double datatmp[16];
	
	for (r=0; r<4; r++)
	for (c=0; c<4; c++)
		datatmp[c + 4 * r] = matTmp[c+4 + 8 * r];
	
	return Matrice4x4(datatmp);*/
}

Matrice4x4 Matrice4x4::operator*(double cst)
{
	double tmp[16];
	for (int i=0; i<16; i++) tmp[i] = cst * data[i];
	
	return Matrice4x4(tmp);
}

Vector3 Matrice4x4::operator*(const Vector3& V) const
{
	double a = V[0] * data[0] + V[1] * data[1] + V[2] * data[2] + data[3];
	double b = V[0] * data[4] + V[1] * data[5] + V[2] * data[6] + data[7];
	double c = V[0] * data[8] + V[1] * data[9] + V[2] * data[10] + data[11];
	
	return Vector3(a, b, c);
}

bool Matrice4x4::operator==(const Matrice4x4& mat) const
{
	bool same = true;
	
	for (int i=0; i<16; i++)
	{
		if (data[i] != mat.data[i])
		{
			same = false;
			break;
		}
	}
	
	return same;
}

ostream& operator<<(ostream& stream, const Matrice4x4& mat)
{
	for (int r=0; r<4; r++)
	{
		for (int c=0; c<4; c++)
			stream << mat(r, c) << " ";
	
		stream << endl;
	}
	
	return stream;
}

Vector3 vectProduct(const Vector3& V1, const Vector3& V2)
{
	double a = V1[1] * V2[2] - V1[2] * V2[1];
	double b = V1[2] * V2[0] - V1[0] * V2[2];
	double c = V1[0] * V2[1] - V1[1] * V2[0];
	
	return Vector3(a, b, c);
}

Vector3 projectionProduct(const Vector3& V1, const Vector3& V2)
{
	return Vector3(V1[0] * V2[0], V1[1] * V2[1], V1[2] * V2[2]);
}

Matrice3x3 outerProduct(const Vector3& A, const Vector3& B)
{
	return Matrice3x3(A[0]*B[0], A[0]*B[1], A[0]*B[2],
			  A[1]*B[0], A[1]*B[1], A[1]*B[2],
			  A[2]*B[0], A[2]*B[1], A[2]*B[2]);
}

/** Matrice3x3_complex **/

Matrice3x3_complex::Matrice3x3_complex()
{
	for (int i=0; i<9; i++)
		data[i] = complex<double>(0, 0);
}

Matrice3x3_complex::Matrice3x3_complex(complex<double> a00, complex<double> a01, complex<double> a02,
				       complex<double> a10, complex<double> a11, complex<double> a12,
				       complex<double> a20, complex<double> a21, complex<double> a22)
{
	data[0] = a00;
	data[1] = a01;
	data[2] = a02;
	data[3] = a10;
	data[4] = a11;
	data[5] = a12;
	data[6] = a20;
	data[7] = a21;
	data[8] = a22;
}

Matrice3x3_complex::Matrice3x3_complex(const vector< complex<double> >& d)
{
	if (d.size() != 9) throw("Error: Please, check array size (must be 9)");
	
	for (int i=0; i<9; i++) data[i] = d[i];
}

complex<double> Matrice3x3_complex::operator() (unsigned row, unsigned col) const
{
	if (row >= 3 || col >= 3)
		throw string("Error: matrix subscript out of bounds");
	
	return data[col + 3 * row];
}

Matrice3x3_complex Matrice3x3_complex::operator*(complex<double> cst)
{
	vector< complex<double> > tmp(9);
	for (int i=0; i<9; i++) tmp[i] = cst * data[i];
	
	return Matrice3x3_complex(tmp);
}

complex<double> Matrice3x3_complex::det()
{
	return 	data[0] * (data[4] * data[8] - data[5] * data[7]) -
		data[1] * (data[3] * data[8] - data[5] * data[6]) +
		data[2] * (data[3] * data[7] - data[4] * data[6]);
}

Matrice3x3_complex Matrice3x3_complex::inv()
{
	return Matrice3x3_complex(
			data[4] * data[8] - data[5] * data[7], data[2] * data[7] - data[1] * data[8], data[1] * data[5] - data[2] * data[4],
			data[5] * data[6] - data[3] * data[8], data[0] * data[8] - data[2] * data[6], data[2] * data[3] - data[0] * data[5],
			data[3] * data[7] - data[4] * data[6], data[1] * data[6] - data[0] * data[7], data[0] * data[4] - data[1] * data[3]
			) * (1. / det());
}

ostream& operator<<(ostream& stream, const Matrice3x3_complex& mat)
{
	for (int r=0; r<3; r++)
	{
		for (int c=0; c<3; c++)
			stream << mat(r, c) << " ";
	
		stream << endl;
	}
	
	return stream;
}


/** Matrice4x4_complex **/

Matrice4x4_complex::Matrice4x4_complex()
{
	for (int i=0; i<16; i++)
		data[i] = 0;
}

Matrice4x4_complex::Matrice4x4_complex(complex<double> a00, complex<double> a01, complex<double> a02, complex<double> a03,
				       complex<double> a10, complex<double> a11, complex<double> a12, complex<double> a13,
				       complex<double> a20, complex<double> a21, complex<double> a22, complex<double> a23,
				       complex<double> a30, complex<double> a31, complex<double> a32, complex<double> a33)
{
	data[0] = a00;
	data[1] = a01;
	data[2] = a02;
	data[3] = a03;
	data[4] = a10;
	data[5] = a11;
	data[6] = a12;
	data[7] = a13;
	data[8] = a20;
	data[9] = a21;
	data[10] = a22;
	data[11] = a23;
	data[12] = a30;
	data[13] = a31;
	data[14] = a32;
	data[15] = a33;
}

Matrice4x4_complex::Matrice4x4_complex(const vector< complex<double> >& d)
{
	if (d.size() != 16) throw("Error: Please, check array size (must be 16)");
	
	for (int i=0; i<16; i++) data[i] = d[i];
}

complex<double> Matrice4x4_complex::operator() (unsigned row, unsigned col) const
{
	if (row >= 4 || col >= 4)
		throw string("Error: matrix subscript out of bounds");
	
	return data[col + 4 * row];
}

Matrice4x4_complex Matrice4x4_complex::operator*(complex<double> cst)
{
	vector< complex<double> > tmp(16);
	for (int i=0; i<16; i++) tmp[i] = cst * data[i];
	
	return Matrice4x4_complex(tmp);
}

complex<double> Matrice4x4_complex::det()
{
	Matrice3x3_complex minor0(data[5], data[6], data[7],
			  data[9], data[10], data[11],
			  data[13], data[14], data[15]);
	
	Matrice3x3_complex minor1(data[4], data[6], data[7],
			  data[8], data[10], data[11],
			  data[12], data[14], data[15]);
	
	Matrice3x3_complex minor2(data[4], data[5], data[7],
			  data[8], data[9], data[11],
			  data[12], data[13], data[15]);
	
	Matrice3x3_complex minor3(data[4], data[5], data[6],
			  data[8], data[9], data[10],
			  data[12], data[13], data[14]);
	
	return 	data[0] * minor0.det() -
			data[1] * minor1.det() +
			data[2] * minor2.det() -
			data[3] * minor3.det();
}

Matrice4x4_complex Matrice4x4_complex::inv()
{
	Matrice3x3_complex minor0(data[5], data[6], data[7],
			  data[9], data[10], data[11],
			  data[13], data[14], data[15]);
	
	Matrice3x3_complex minor1(data[4], data[6], data[7],
			  data[8], data[10], data[11],
			  data[12], data[14], data[15]);
	
	Matrice3x3_complex minor2(data[4], data[5], data[7],
			  data[8], data[9], data[11],
			  data[12], data[13], data[15]);
	
	Matrice3x3_complex minor3(data[4], data[5], data[6],
			  data[8], data[9], data[10],
			  data[12], data[13], data[14]);
	
	Matrice3x3_complex minor4(data[1], data[2], data[3],
			  data[9], data[10], data[11],
			  data[13], data[14], data[15]);
	
	Matrice3x3_complex minor5(data[0], data[2], data[3],
			  data[8], data[10], data[11],
			  data[12], data[14], data[15]);
	
	Matrice3x3_complex minor6(data[0], data[1], data[3],
			  data[8], data[9], data[11],
			  data[12], data[13], data[15]);
	
	Matrice3x3_complex minor7(data[0], data[1], data[2],
			  data[8], data[9], data[10],
			  data[12], data[13], data[14]);
	
	Matrice3x3_complex minor8(data[1], data[2], data[3],
			  data[5], data[6], data[7],
			  data[13], data[14], data[15]);
	
	Matrice3x3_complex minor9(data[0], data[2], data[3],
			  data[4], data[6], data[7],
			  data[12], data[14], data[15]);
	
	Matrice3x3_complex minor10(data[0], data[1], data[3],
			   data[4], data[5], data[7],
			   data[12], data[13], data[15]);
	
	Matrice3x3_complex minor11(data[0], data[1], data[2],
			   data[4], data[5], data[6],
			   data[12], data[13], data[14]);
	
	Matrice3x3_complex minor12(data[1], data[2], data[3],
			   data[5], data[6], data[7],
			   data[9], data[10], data[11]);
	
	Matrice3x3_complex minor13(data[0], data[2], data[3],
			   data[4], data[6], data[7],
			   data[8], data[10], data[11]);
	
	Matrice3x3_complex minor14(data[0], data[1], data[3],
			   data[4], data[5], data[7],
			   data[8], data[9], data[11]);
	
	Matrice3x3_complex minor15(data[0], data[1], data[2],
			   data[4], data[5], data[6],
			   data[8], data[9], data[10]);
	
	return Matrice4x4_complex( minor0.det(), - minor4.det(), minor8.det(), - minor12.det(),
				   - minor1.det(), minor5.det(), - minor9.det(), minor13.det(),
				   minor2.det(), - minor6.det(), minor10.det(), - minor14.det(),
				   - minor3.det(), minor7.det(), - minor11.det(), minor15.det() )
			* (1. / det());
}

ostream& operator<<(ostream& stream, const Matrice4x4_complex& mat)
{
	for (int r=0; r<4; r++)
	{
		for (int c=0; c<4; c++)
			stream << mat(r, c) << " ";
	
		stream << endl;
	}
	
	return stream;
}
