#include <float.h>
#ifdef WIN32
void BMDS_Win_fpset(void) {
  _fpreset();
  _controlfp(_MCW_EM, _MCW_EM);
  _controlfp(_PC_64, _MCW_PC);
}

#endif
