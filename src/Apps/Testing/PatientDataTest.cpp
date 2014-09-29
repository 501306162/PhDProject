#include <PatientData.h>
#include <iostream>

int main(int argc, char *argv[])
{
	
	utils::PatientData::Pointer data = utils::PatientData::Load(argv[1]);


	return 0;
}
