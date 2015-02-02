#ifndef INITIAL_PARTICLE_GENERATOR_H
#define INITIAL_PARTICLE_GENERATOR_H

#include <itkImage.h>
#include <MatrixCommon.h>
#include <itkSimilarity3DTransform.h>

namespace vt
{
class InitialParticleGenerator : public itk::Object
{
public:
	typedef InitialParticleGenerator Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkTypeMacro(InitialParticleGenerator, Object);
	itkNewMacro(Self);

	typedef itk::Image<unsigned short, 3> ImageType;
	typedef ImageType::PointType PointType;
	typedef itk::Vector<double, 3> VectorType;
	typedef itk::Similarity3DTransform<double> TransformType;

	typedef struct particle_
	{
		PointType p1;
		PointType p2;
		VectorType direction;
		double length;
	} Particle;


	typedef utils::DoubleMatrixType MatrixType;
	typedef std::pair<MatrixType, MatrixType> PointPair;
	typedef std::pair<ImageType::Pointer, PointPair> ImagePointPair;

	typedef std::vector<Particle> ParticleList;
	typedef std::vector<PointType> PointList;
	typedef std::map<std::string, ParticleList> ParticleMap;


	void SetImage(const ImageType::Pointer &image) { m_Image = image; }
	void SetPoints(const MatrixType &points1, const MatrixType &points2) { m_Points1 = points1; m_Points2 = points2; }
	void SetTransform(const TransformType::Pointer &transform) { m_Transform = transform; }
	void SetNumParticles(const unsigned int num) { m_NumParticles = num; }
	void AddInput(const std::string &type, const ImagePointPair &data) { m_Data[type] = data; }

	void Generate();
	ParticleMap GetMapOutput() const { return m_MapOut; }

	ParticleList GetOutput() const { return m_Output; }

protected:
	InitialParticleGenerator() {}
	virtual ~InitialParticleGenerator() {}

private:

	void GetAlignedPoints(const MatrixType &input, MatrixType &aligned);
	void GetAlignedPoints3(const MatrixType &input, ImageType::Pointer &image, MatrixType &aligned);
	void GetAlignedPoints2(const MatrixType &input, MatrixType &aligned);
	PointType ProjectPoint(const PointType &p);
	PointType ProjectPoint2(const PointType &p, const ImageType::Pointer &image);
	void GetSamples(const MatrixType &aligned, PointList &samples);
	void GetSamples2(const MatrixType &aligned1, const MatrixType &aligned2, PointList &samples1, PointList &samples2);

	TransformType::Pointer m_Transform;
	MatrixType m_Points1;
	MatrixType m_Points2;
	unsigned int m_NumParticles;
	ImageType::Pointer m_Image;
	ParticleList m_Output;


	ParticleMap m_MapOut;
	std::map<std::string, ImagePointPair> m_Data;

	InitialParticleGenerator(const Self&);
	void operator=(const Self&);
};



}

#endif
