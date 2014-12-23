#include "ValveIO.h"

#include <CommonDefinitions.h>

#include <vtkSmartPointer.h>
#include <vtkLineSource.h>
#include <vtkPolyDataWriter.h>


#include <QTextStream>
#include <QFile>
#include <QVariant>
#include <qjson/serializer.h>
#include <qjson/parser.h>

#include <vtkPolyDataReader.h>
#include <vtkPoints.h>


namespace vt 
{
// ------------------------------------------------------------------------
template<int VDimensions>
void ValveWriter<VDimensions>::Write()
{
	// construct the names of the parts 
	std::string linePart = m_Filename + "_line.vtk";
	std::string imagePart = m_Filename + "_image.nrrd";
	std::string dataPart = m_Filename + ".txt";


	typename ValveLineType::ImageType::Pointer image = m_Valve->GetImage();
	utils::ImageIO<typename ValveLineType::ImageType>::Write(imagePart, image);

	vtkSmartPointer<vtkPolyData> line = m_Valve->GetPolyData();
	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetFileName(linePart.c_str());
	writer->SetInputData(line);
	writer->Write();
	writer->Update();


}

// ------------------------------------------------------------------------
template<int VDimensions>
typename ValveReader<VDimensions>::ValveLineType::Pointer 
ValveReader<VDimensions>::GetOutput() const 
{
	std::string imagePart = m_Filename + "_image.nrrd";
	std::string linePart = m_Filename + "_line.vtk";


	// load the polydata 
	vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
	reader->SetFileName(linePart.c_str());
	reader->Update();

	vtkSmartPointer<vtkPolyData> poly = reader->GetOutput();
	typename ValveLineType::Pointer valve = ValveLineType::New();


	vtkPoints * points = poly->GetPoints();
	typename ValveLineType::PointType p1,p2;
	for(unsigned int i = 0; i < VDimensions; i++)
	{
		p1[i] = points->GetPoint(0)[i];
		p2[i] = points->GetPoint(1)[i];
	}
	valve->SetP1(p1);
	valve->SetP2(p2);


	typename ValveLineType::ImageType::Pointer image = 
		utils::ImageIO<typename ValveLineType::ImageType>::Read(imagePart);
	valve->SetImage(image);


	typename ValveLineType::ContIndexType ind1, ind2;
	image->TransformPhysicalPointToContinuousIndex(p1, ind1);
	image->TransformPhysicalPointToContinuousIndex(p2, ind2);
	valve->SetInd1(ind1);
	valve->SetInd2(ind2);

	return valve;

}


// ------------------------------------------------------------------------
template<int VDimensions>
typename ValveSequenceReader<VDimensions>::ValveSequenceType::Pointer 
ValveSequenceReader<VDimensions>::GetOutput() const 
{
	// read the file 
	QFile * file = new QFile(QString::fromStdString(m_Filename));
	file->open(QIODevice::ReadOnly | QIODevice::Text);


	QJson::Parser parser;
	bool ok;
	QVariantList list = parser.parse(file, &ok).toList();
	file->close();
	delete file;

	typename ValveSequenceType::Pointer valve = ValveSequenceType::New();

	for(int i = 0; i < list.size(); i++)
	{
		QString fname = list[i].toString();
		std::string name = fname.toStdString();

		typename ValveReader<VDimensions>::Pointer reader = ValveReader<VDimensions>::New();
		reader->SetFileName(name);
		valve->AddValveLine(reader->GetOutput());
	}

	return valve;

}


// ------------------------------------------------------------------------
template<int VDimensions>
void ValveSequenceWriter<VDimensions>::Write()
{
	std::string dataPart = m_Filename + ".txt";
	
	unsigned int timeSteps = m_Sequence->GetNumberOfLines();

	QVariantList fileList;

	for(unsigned int i = 0; i < timeSteps; i++)
	{
		std::stringstream ss;
		ss << m_Filename << "_time_" << i;
		std::string partFilename = ss.str();

		typename ValveWriter<VDimensions>::Pointer writer = ValveWriter<VDimensions>::New();
		writer->SetFileName(partFilename);
		writer->SetInput(m_Sequence->GetValveLine(i));
		writer->Write();

		QString fname = QString::fromStdString(partFilename);

		fileList << fname;
	}

	QJson::Serializer serialiser;
	bool ok;
	QString json = serialiser.serialize(fileList, &ok);
	
	QFile file(QString::fromStdString(dataPart));
	file.open(QIODevice::WriteOnly | QIODevice::Text);
	QTextStream out(&file);

	out << json;
	file.close();
	
}

}


template class vt::ValveReader<2>;
template class vt::ValveReader<3>;
template class vt::ValveWriter<2>;
template class vt::ValveWriter<3>;
template class vt::ValveSequenceReader<2>;
template class vt::ValveSequenceReader<3>;
template class vt::ValveSequenceWriter<2>;
template class vt::ValveSequenceWriter<3>;



