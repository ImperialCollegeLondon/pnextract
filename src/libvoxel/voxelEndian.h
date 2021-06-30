//stackoverflow : Kevin: https://stackoverflow.com/a/105342/2718352
//How do I convert between big-endian and little-endian values in C++?

template<typename T> inline T be2le(T vv) {
	std::cout<<"err"<<std::endl; alert("not implemented",-1);	return vv; }
template<> inline unsigned short  be2le(unsigned short vv) {
	return (vv >> 8) | (vv << 8);  }
template<> inline unsigned int    be2le(unsigned int ui) {
    return (ui >> 24) | ((ui<<8) & 0x00FF0000) | ((ui>>8) & 0x0000FF00) | (ui << 24);  }
template<> inline unsigned long long  be2le(unsigned long long ull)  {
    return (ull >> 56) |
          ((ull<<40) & 0x00FF000000000000) | ((ull<<24) & 0x0000FF0000000000) |
          ((ull<<8)  & 0x000000FF00000000) |  ((ull>>8) & 0x00000000FF000000) |
          ((ull>>24) & 0x0000000000FF0000) | ((ull>>40) & 0x000000000000FF00) |
          (ull << 56);  }


template<typename T, template<typename ...> class C>
void flipEndian(C<T>& vImage)
{
	if(sizeof(T)>1)
		std::transform(vImage.begin(), vImage.end(), vImage.begin(), be2le<T>);
}
