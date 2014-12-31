#include <iostream>

#include <BoundingBox.h>


int main(int argc, char ** argv)
{
	// load the image and the bounding box
	vt::BoundingBox::Pointer boundingBox = vt::BoundingBox::New();
	boundingBox->Load(argv[1]);
	


	return 0;
}
