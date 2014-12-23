#include <iostream>

#include <Directory.h>
#include <ValveIO.h>

#include <ValveHelpers.h>
#include <ValveAligner.h>


using namespace vt;

typedef ValveLine<2> FlatValve;
typedef ValveSequence<2> FlatSequence;


int main(int argc, char ** argv)
{
	utils::Directory::Pointer dir = utils::Directory::New();
	dir->SetDirectory(argv[1]);
	dir->SetExtension(".txt");
	utils::Directory::FilenamesType filenames = dir->GetOutput();

	std::vector<FlatSequence::Pointer> valves;
	for(unsigned int i = 0; i < filenames.size(); i++)
	{
		ValveSequenceReader<3>::Pointer reader = ValveSequenceReader<3>::New();
		reader->SetFileName(filenames[i]);

		// flatten the valves
		FlatSequence::Pointer valveSequence = FlatSequence::New();
		FlattenValveSequence(reader->GetOutput(), valveSequence);

		valves.push_back(valveSequence);

		ValveSequenceAligner::Pointer aligner = ValveSequenceAligner::New();
		aligner->SetFixedSequence(valves[0]);
		aligner->SetMovingSequence(valves[i]);
		aligner->Align();

		std::stringstream ss;
		ss << argv[2] << "/d" << i;

		ValveSequenceWriter<2>::Pointer writer = ValveSequenceWriter<2>::New();
		writer->SetInput(aligner->GetOutput());
		writer->SetFileName(ss.str());
		writer->Write();


	}


	return 0;

}
