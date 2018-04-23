/*
 * Test program for gzifstream and gzofstream
 *
 * by Ludwig Schwardt <schwardt@sun.ac.za>
 * original version by Kevin Ruland <kevin@rodin.wustl.edu>
 */

#include "zfstream.h"
#include <iostream>      // for cout

int main() {

  gzofstream outf;
  gzifstream inf;
  char buf[80];

  outf.open("test1.txt.gz");
  outf 
  << "baaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaae\n" 
  << "baaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaae\n" 
  << "baaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaae\n" 
  << "baaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaae\n" 
  << "baaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaae\n" 
  << 1.3 << "\nPlan " << 9 << std::endl;
  outf.close();

  std::cout << "\nReading 'test1.txt.gz' (buffered) produces:\n";
  inf.open("test1.txt.gz");
  (inf.read(buf,50)) ;
  {
    std::cout << buf << "\t(" << inf.rdbuf()->in_avail() << " chars left in buffer)\n";
  }
  inf.close();

  outf.rdbuf()->pubsetbuf(0,0);
  outf.open("test11.txt.gz");
  outf << setcompression(Z_DEFAULT_COMPRESSION)
       << buf << 1.3 << "aa\nPlan " << 9 << std::endl;
  outf.close();

  std::cout << "\nReading 'test2.txt.gz' (unbuffered) produces:\n";
  inf.rdbuf()->pubsetbuf(0,0);
  inf.open("test11.txt.gz");
  while (inf.getline(buf,80)) { 
    std::cout << buf <<a<< "\t(" << inf.rdbuf()->in_avail() << " chars left in buffer)\n";
  }
  inf.close();

  return 0;

}
