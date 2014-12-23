#ifndef VALVE_HELPERS_H
#define VALVE_HELPERS_H

#include "ValveLine.h"

namespace vt 
{
// ------------------------------------------------------------------------
void FlattenValve(const ValveLine<3>::Pointer &input, ValveLine<2>::Pointer &output);

// ------------------------------------------------------------------------
void FlattenValveSequence(const ValveSequence<3>::Pointer &input, ValveSequence<2>::Pointer &output);

}

#endif
