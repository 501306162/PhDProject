#include <iostream>

#include <Directory.h>
#include <ValveIO.h>
#include <ValveNormaliser.h>
#include <QString>
#include <FlipChecker.h>
#include <MatrixCommon.h>
#include <MatrixWriter.h>

using namespace vt;
int main(int argc, char ** argv)
{
	utils::Directory::FilenamesType folderNames = utils::Directory::GetDirectories(argv[1]);
	const std::string outputDirectory = argv[2];
	const unsigned int numTimeSteps = 25;

	typedef utils::DoubleMatrixType MatrixType;
	typedef utils::IntMatrixType IntMatrixType;

	

	for(unsigned int i = 0; i < folderNames.size(); i++)
	{
		// get the valve type 
		const std::string folder = folderNames[i];
		std::string type = utils::Directory::GetFileName(folder);

		// get the line files in the directory 
		utils::Directory::FilenamesType lineFiles = utils::Directory::GetFiles(folder, ".txt");
		MatrixType lengths = MatrixType(lineFiles.size(), numTimeSteps);
		IntMatrixType owners = IntMatrixType(lineFiles.size(), 1);

		for(unsigned int lineNum = 0; lineNum < lineFiles.size(); lineNum++)
		{
			const std::string lineFileName = lineFiles[lineNum];
			QString fname = QString::fromStdString(utils::Directory::GetFileName(lineFileName));
			const unsigned int id = fname.replace("d","").replace(".txt","").toInt();

			// load the valve sequence
			ValveSequenceReader<3>::Pointer reader = ValveSequenceReader<3>::New();
			reader->SetFileName(lineFileName);
			ValveSequence<3>::Pointer sequence = reader->GetOutput();

			owners(lineNum, 0) = id;

			// loop through the time steps
			for(unsigned int timeStep = 0; timeStep < numTimeSteps; timeStep++)
			{
				ValveLine<3>::Pointer line = sequence->GetValveLine(timeStep);
				double dist = line->GetP1().EuclideanDistanceTo(line->GetP2());
				lengths(lineNum, timeStep) = dist;
			}
		}
	
		// save the data
		utils::MatrixDataSet::Pointer dset = utils::MatrixDataSet::New();
		std::string filename = utils::Directory::GetPath(outputDirectory, type+".hdf5");
		std::cout << filename << std::endl;
		dset->AddData("lengths", lengths);
		dset->AddData("owners", owners);
		utils::MatrixWriter::Pointer writer = utils::MatrixWriter::New();
		writer->SetInput(dset);


		writer->SetFilename(filename);
		writer->Write();
			

	}



	return 0;
}
