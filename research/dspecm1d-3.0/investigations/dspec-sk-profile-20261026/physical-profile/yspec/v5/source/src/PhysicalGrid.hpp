#pragma once

namespace YSpec {
// Set only by the isolated physical-profile driver. Production code and the
// normal adapter leave the solver's historical DC skip in place.
extern bool PhysicalIncludeDC;
extern unsigned PhysicalFftPadding;
}
