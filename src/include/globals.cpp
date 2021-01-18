
// for declaration, see globals.h
#ifdef _debugCompile_
int    globals::debugLevel_ = 1;
#else
int    globals::debugLevel_ = 0;
#endif
size_t globals::typeCounter_ = 0;