#include <iostream>
#include <MatrixReader.h>
#include <DiffusionMap.h>
#include <MatrixWriter.h>

int main(int argc, char *argv[])
{
	// load the data
	typedef utils::MatrixReader MatReader;
	MatReader::Pointer reader = MatReader::New();
	reader->SetFilename(argv[1]);
	reader->Read();

	utils::MatrixDataSet::Pointer data = reader->GetOutput();

	typedef manifold::DiffusionMap DiffMap;
	DiffMap::Pointer diffMap = DiffMap::New();
	diffMap->Fit(data->doubleData["swiss"].transpose());

	DiffMap::MatrixType E;
	diffMap->GetEmbedding(E);


	utils::MatrixDataSet::Pointer output = utils::MatrixDataSet::New();
	output->AddData("E",E);



	typedef utils::MatrixWriter Writer;
	Writer::Pointer writer = Writer::New();
	writer->SetInput(output);
	writer->SetFilename("E.mat");
	writer->Write();

	

	return 0;
}
