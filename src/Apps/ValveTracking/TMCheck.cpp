#include <iostream>

#include <SVMClassifier.h>
#include <MatrixCommon.h>

int main(int argc, char ** argv)
{
	typedef utils::DoubleMatrixType MatrixType;
	typedef utils::IntMatrixType IntMatrixType;

	MatrixType X = MatrixType::Random(1000,100);
	IntMatrixType y = IntMatrixType::Ones(1000,1);

	for(unsigned int i = 0; i < 500; i++)
	{
		y(i,0) = 0;		
	}


	vt::SVMClassifier::Pointer cls = vt::SVMClassifier::New();
	cls->Train(X,y);

	return 0;



}
