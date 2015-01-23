#include <iostream>
#include <Directory.h>

#include <ValveIO.h>
#include <FlipChecker.h>

using namespace vt;
int main(int argc, char ** argv)
{
	utils::Directory::FilenamesType filenames = utils::Directory::GetFiles(argv[1], ".txt");
	for(unsigned int i = 0; i < filenames.size(); i++)
	{
		// get the id 
		const std::string filename = filenames[i];
		std::string fname = utils::Directory::GetFileName(filename);
		unsigned int id = QString::fromStdString(fname).replace("d","").replace(".txt","").toInt();


		FlipChecker::Pointer checker = FlipChecker::New();
		bool flipPoints = checker->FlipPoints(argv[2], id);

		if(!flipPoints) continue;

		ValveSequenceReader<3>::Pointer reader = ValveSequenceReader<3>::New();
		reader->SetFileName(filename);
		ValveSequence<3>::Pointer sequence = reader->GetOutput();
		ValveSequence<3>::Pointer output = ValveSequence<3>::New();
		for(unsigned int j = 0; j < sequence->GetNumberOfLines(); j++)
		{
			ValveLine<3>::Pointer line = sequence->GetValveLine(j);			
			ValveLine<3>::PointType p1 = line->GetP1();
			ValveLine<3>::PointType p2 = line->GetP2();
			line->SetP1(p2);
			line->SetP2(p1);
			line->UpdateIndexs();

			output->AddValveLine(line);
		}

		std::stringstream ss;
		ss << utils::Directory::GetPath(argv[1], "d") << id;

		std::cout << ss.str() << std::endl;
		ValveSequenceWriter<3>::Pointer writer = ValveSequenceWriter<3>::New();
		writer->SetInput(output);
		writer->SetFileName(ss.str());
		writer->Write();



	}



	return 0;
}
