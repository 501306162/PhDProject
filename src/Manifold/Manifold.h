#ifndef MANIFOLD_H
#define MANIFOLD_H

#include <itkObject.h>
#include <itkObjectFactory.h>
#include <Eigen/Dense>


#include "MatrixCommon.h"

namespace manifold
{
class Manifold : public itk::Object
{
public:
	typedef Manifold 						Self;
	typedef itk::Object 					Superclass;
	typedef itk::SmartPointer<Self> 		Pointer;
	typedef itk::SmartPointer<const Self> 	ConstPointer;

	itkTypeMacro(Manifold, itk::Object);

	typedef utils::DoubleMatrixType MatrixType;

	virtual void Fit(const MatrixType &X) = 0;
	virtual void Transform(const MatrixType &X, MatrixType &Y) = 0;
	virtual void GetEmbedding(MatrixType &X) = 0;

protected:

};
}





#endif
