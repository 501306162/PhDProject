#ifndef HELPERS_H
#define HELPERS_H

namespace utils
{
template<typename TInputType, typename TOutputType>
void
ExtractSlice(const TInputType * input, const unsigned int slice,
	   typename TOutputType::Pointer &output);	
} /* utils */ 


#include "Helpers.hpp"

#endif
