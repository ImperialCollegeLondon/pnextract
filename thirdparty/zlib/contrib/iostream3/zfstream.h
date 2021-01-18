/*
 * A C++ I/O streams interface to the zlib gz* functions
 *
 * by Ludwig Schwardt <schwardt@sun.ac.za>
 * original version by Kevin Ruland <kevin@rodin.wustl.edu>
 *
 * This version is standard-compliant and compatible with gcc 3.x.
 */


// Ali Q Raeini (2019), merged zfstream.cc with this .h file, to make it header only. 
//The .cc file contents are added in the end to be able to recover the original structure


#ifndef ZFSTREAM_H
#define ZFSTREAM_H

#include <istream>  // not iostream, since we don't need cin/cout
#include <ostream>
#include "zlib.h"

/*****************************************************************************/

/**
 *  @brief  Gzipped file stream buffer class.
 *
 *  This class implements basic_filebuf for gzipped files. It doesn't yet support
 *  seeking (allowed by zlib but slow/limited), putback and read/write access
 *  (tricky). Otherwise, it attempts to be a drop-in replacement for the standard
 *  file streambuf.
*/
class gzfilebuf : public std::streambuf
{
public:
  //  Default constructor.
  gzfilebuf();

  //  Destructor.
  virtual
  ~gzfilebuf();

  /**
   *  @brief  Set compression level and strategy on the fly.
   *  @param  comp_level  Compression level (see zlib.h for allowed values)
   *  @param  comp_strategy  Compression strategy (see zlib.h for allowed values)
   *  @return  Z_OK on success, Z_STREAM_ERROR otherwise.
   *
   *  Unfortunately, these parameters cannot be modified separately, as the
   *  previous zfstream version assumed. Since the strategy is seldom changed,
   *  it can default and setcompression(level) then becomes like the old
   *  setcompressionlevel(level).
  */
  int
  setcompression(int comp_level,
                 int comp_strategy = Z_DEFAULT_STRATEGY);

  /**
   *  @brief  Check if file is open.
   *  @return  True if file is open.
  */
  bool
  is_open() const { return (file != NULL); }

  /**
   *  @brief  Open gzipped file.
   *  @param  name  File name.
   *  @param  mode  Open mode flags.
   *  @return  @c this on success, NULL on failure.
  */
  gzfilebuf*
  open(const char* name,
       std::ios_base::openmode mode);

  /**
   *  @brief  Attach to already open gzipped file.
   *  @param  fd  File descriptor.
   *  @param  mode  Open mode flags.
   *  @return  @c this on success, NULL on failure.
  */
  gzfilebuf*
  attach(int fd,
         std::ios_base::openmode mode);

  /**
   *  @brief  Close gzipped file.
   *  @return  @c this on success, NULL on failure.
  */
  gzfilebuf*
  close();

protected:
  /**
   *  @brief  Convert ios open mode int to mode string used by zlib.
   *  @return  True if valid mode flag combination.
  */
  bool
  open_mode(std::ios_base::openmode mode,
            char* c_mode) const;

  /**
   *  @brief  Number of characters available in stream buffer.
   *  @return  Number of characters.
   *
   *  This indicates number of characters in get area of stream buffer.
   *  These characters can be read without accessing the gzipped file.
  */
  virtual std::streamsize
  showmanyc();

  /**
   *  @brief  Fill get area from gzipped file.
   *  @return  First character in get area on success, EOF on error.
   *
   *  This actually reads characters from gzipped file to stream
   *  buffer. Always buffered.
  */
  virtual int_type
  underflow();

  /**
   *  @brief  Write put area to gzipped file.
   *  @param  c  Extra character to add to buffer contents.
   *  @return  Non-EOF on success, EOF on error.
   *
   *  This actually writes characters in stream buffer to
   *  gzipped file. With unbuffered output this is done one
   *  character at a time.
  */
  virtual int_type
  overflow(int_type c = traits_type::eof());

  /**
   *  @brief  Installs external stream buffer.
   *  @param  p  Pointer to char buffer.
   *  @param  n  Size of external buffer.
   *  @return  @c this on success, NULL on failure.
   *
   *  Call setbuf(0,0) to enable unbuffered output.
  */
  virtual std::streambuf*
  setbuf(char_type* p,
         std::streamsize n);

  /**
   *  @brief  Flush stream buffer to file.
   *  @return  0 on success, -1 on error.
   *
   *  This calls underflow(EOF) to do the job.
  */
  virtual int
  sync();

//
// Some future enhancements
//
//  virtual int_type uflow();
//  virtual int_type pbackfail(int_type c = traits_type::eof());
//  virtual pos_type
//  seekoff(off_type off,
//          std::ios_base::seekdir way,
//          std::ios_base::openmode mode = std::ios_base::in|std::ios_base::out);
//  virtual pos_type
//  seekpos(pos_type sp,
//          std::ios_base::openmode mode = std::ios_base::in|std::ios_base::out);

private:
  /**
   *  @brief  Allocate internal buffer.
   *
   *  This function is safe to call multiple times. It will ensure
   *  that a proper internal buffer exists if it is required. If the
   *  buffer already exists or is external, the buffer pointers will be
   *  reset to their original state.
  */
  void
  enable_buffer();

  /**
   *  @brief  Destroy internal buffer.
   *
   *  This function is safe to call multiple times. It will ensure
   *  that the internal buffer is deallocated if it exists. In any
   *  case, it will also reset the buffer pointers.
  */
  void
  disable_buffer();

  /**
   *  Underlying file pointer.
  */
  gzFile file;

  /**
   *  Mode in which file was opened.
  */
  std::ios_base::openmode io_mode;

  /**
   *  @brief  True if this object owns file descriptor.
   *
   *  This makes the class responsible for closing the file
   *  upon destruction.
  */
  bool own_fd;

  /**
   *  @brief  Stream buffer.
   *
   *  For simplicity this remains allocated on the free store for the
   *  entire life span of the gzfilebuf object, unless replaced by setbuf.
  */
  char_type* buffer;

  /**
   *  @brief  Stream buffer size.
   *
   *  Defaults to system default buffer size (typically 8192 bytes).
   *  Modified by setbuf.
  */
  std::streamsize buffer_size;

  /**
   *  @brief  True if this object owns stream buffer.
   *
   *  This makes the class responsible for deleting the buffer
   *  upon destruction.
  */
  bool own_buffer;
};

/*****************************************************************************/

/**
 *  @brief  Gzipped file input stream class.
 *
 *  This class implements ifstream for gzipped files. Seeking and putback
 *  is not supported yet.
*/
class gzifstream : public std::istream
{
public:
  //  Default constructor
  gzifstream();

  /**
   *  @brief  Construct stream on gzipped file to be opened.
   *  @param  name  File name.
   *  @param  mode  Open mode flags (forced to contain ios::in).
  */
  explicit
  gzifstream(const char* name,
             std::ios_base::openmode mode = std::ios_base::in);

  /**
   *  @brief  Construct stream on already open gzipped file.
   *  @param  fd    File descriptor.
   *  @param  mode  Open mode flags (forced to contain ios::in).
  */
  explicit
  gzifstream(int fd,
             std::ios_base::openmode mode = std::ios_base::in);

  /**
   *  Obtain underlying stream buffer.
  */
  gzfilebuf*
  rdbuf() const
  { return const_cast<gzfilebuf*>(&sb); }

  /**
   *  @brief  Check if file is open.
   *  @return  True if file is open.
  */
  bool
  is_open() { return sb.is_open(); }

  /**
   *  @brief  Open gzipped file.
   *  @param  name  File name.
   *  @param  mode  Open mode flags (forced to contain ios::in).
   *
   *  Stream will be in state good() if file opens successfully;
   *  otherwise in state fail(). This differs from the behavior of
   *  ifstream, which never sets the state to good() and therefore
   *  won't allow you to reuse the stream for a second file unless
   *  you manually clear() the state. The choice is a matter of
   *  convenience.
  */
  void
  open(const char* name,
       std::ios_base::openmode mode = std::ios_base::in);

  /**
   *  @brief  Attach to already open gzipped file.
   *  @param  fd  File descriptor.
   *  @param  mode  Open mode flags (forced to contain ios::in).
   *
   *  Stream will be in state good() if attach succeeded; otherwise
   *  in state fail().
  */
  void
  attach(int fd,
         std::ios_base::openmode mode = std::ios_base::in);

  /**
   *  @brief  Close gzipped file.
   *
   *  Stream will be in state fail() if close failed.
  */
  void
  close();

private:
  /**
   *  Underlying stream buffer.
  */
  gzfilebuf sb;
};

/*****************************************************************************/

/**
 *  @brief  Gzipped file output stream class.
 *
 *  This class implements ofstream for gzipped files. Seeking and putback
 *  is not supported yet.
*/
class gzofstream : public std::ostream
{
public:
  //  Default constructor
  gzofstream();

  /**
   *  @brief  Construct stream on gzipped file to be opened.
   *  @param  name  File name.
   *  @param  mode  Open mode flags (forced to contain ios::out).
  */
  explicit
  gzofstream(const char* name,
             std::ios_base::openmode mode = std::ios_base::out);

  /**
   *  @brief  Construct stream on already open gzipped file.
   *  @param  fd    File descriptor.
   *  @param  mode  Open mode flags (forced to contain ios::out).
  */
  explicit
  gzofstream(int fd,
             std::ios_base::openmode mode = std::ios_base::out);

  /**
   *  Obtain underlying stream buffer.
  */
  gzfilebuf*
  rdbuf() const
  { return const_cast<gzfilebuf*>(&sb); }

  /**
   *  @brief  Check if file is open.
   *  @return  True if file is open.
  */
  bool
  is_open() { return sb.is_open(); }

  /**
   *  @brief  Open gzipped file.
   *  @param  name  File name.
   *  @param  mode  Open mode flags (forced to contain ios::out).
   *
   *  Stream will be in state good() if file opens successfully;
   *  otherwise in state fail(). This differs from the behavior of
   *  ofstream, which never sets the state to good() and therefore
   *  won't allow you to reuse the stream for a second file unless
   *  you manually clear() the state. The choice is a matter of
   *  convenience.
  */
  void
  open(const char* name,
       std::ios_base::openmode mode = std::ios_base::out);

  /**
   *  @brief  Attach to already open gzipped file.
   *  @param  fd  File descriptor.
   *  @param  mode  Open mode flags (forced to contain ios::out).
   *
   *  Stream will be in state good() if attach succeeded; otherwise
   *  in state fail().
  */
  void
  attach(int fd,
         std::ios_base::openmode mode = std::ios_base::out);

  /**
   *  @brief  Close gzipped file.
   *
   *  Stream will be in state fail() if close failed.
  */
  void
  close();

private:
  /**
   *  Underlying stream buffer.
  */
  gzfilebuf sb;
};

/*****************************************************************************/

/**
 *  @brief  Gzipped file output stream manipulator class.
 *
 *  This class defines a two-argument manipulator for gzofstream. It is used
 *  as base for the setcompression(int,int) manipulator.
*/
template<typename T1, typename T2>
  class gzomanip2
  {
  public:
    // Allows insertor to peek at internals
    template <typename Ta, typename Tb>
      friend gzofstream&
      operator<<(gzofstream&,
                 const gzomanip2<Ta,Tb>&);

    // Constructor
    gzomanip2(gzofstream& (*f)(gzofstream&, T1, T2),
              T1 v1,
              T2 v2);
  private:
    // Underlying manipulator function
    gzofstream&
    (*func)(gzofstream&, T1, T2);

    // Arguments for manipulator function
    T1 val1;
    T2 val2;
  };

/*****************************************************************************/

// Manipulator function thunks through to stream buffer
inline gzofstream&
setcompression(gzofstream &gzs, int l, int s = Z_DEFAULT_STRATEGY)
{
  (gzs.rdbuf())->setcompression(l, s);
  return gzs;
}

// Manipulator constructor stores arguments
template<typename T1, typename T2>
  inline
  gzomanip2<T1,T2>::gzomanip2(gzofstream &(*f)(gzofstream &, T1, T2),
                              T1 v1,
                              T2 v2)
  : func(f), val1(v1), val2(v2)
  { }

// Insertor applies underlying manipulator function to stream
template<typename T1, typename T2>
  inline gzofstream&
  operator<<(gzofstream& s, const gzomanip2<T1,T2>& m)
  { return (*m.func)(s, m.val1, m.val2); }

// Insert this onto stream to simplify setting of compression level
inline gzomanip2<int,int>
setcompression(int l, int s = Z_DEFAULT_STRATEGY)
{ return gzomanip2<int,int>(&setcompression, l, s); }



// ***************************************************************************
// contents of zfstream.cc ***************************************************


// Internal buffer sizes (default and "unbuffered" versions)
#define BIGBUFSIZE BUFSIZ
#define SMALLBUFSIZE 1

/*****************************************************************************/

// Default constructor
inline gzfilebuf::gzfilebuf()
: file(NULL), io_mode(std::ios_base::openmode(0)), own_fd(false),
  buffer(NULL), buffer_size(BIGBUFSIZE), own_buffer(true)
{
  // No buffers to start with
  this->disable_buffer();
}

// Destructor
inline gzfilebuf::~gzfilebuf()
{
  // Sync output buffer and close only if responsible for file
  // (i.e. attached streams should be left open at this stage)
  this->sync();
  if (own_fd)
    this->close();
  // Make sure internal buffer is deallocated
  this->disable_buffer();
}

// Set compression level and strategy
inline int
gzfilebuf::setcompression(int comp_level,
                          int comp_strategy)
{
  return gzsetparams(file, comp_level, comp_strategy);
}

// Open gzipped file
inline gzfilebuf*
gzfilebuf::open(const char *name,
                std::ios_base::openmode mode)
{
  // Fail if file already open
  if (this->is_open())
    return NULL;
  // Don't support simultaneous read/write access (yet)
  if ((mode & std::ios_base::in) && (mode & std::ios_base::out))
    return NULL;

  // Build mode string for gzopen and check it [27.8.1.3.2]
  char char_mode[6] = "\0\0\0\0\0";
  if (!this->open_mode(mode, char_mode))
    return NULL;

  // Attempt to open file
  if ((file = gzopen(name, char_mode)) == NULL)
    return NULL;

  // On success, allocate internal buffer and set flags
  this->enable_buffer();
  io_mode = mode;
  own_fd = true;
  return this;
}

// Attach to gzipped file
inline gzfilebuf*
gzfilebuf::attach(int fd,
                  std::ios_base::openmode mode)
{
  // Fail if file already open
  if (this->is_open())
    return NULL;
  // Don't support simultaneous read/write access (yet)
  if ((mode & std::ios_base::in) && (mode & std::ios_base::out))
    return NULL;

  // Build mode string for gzdopen and check it [27.8.1.3.2]
  char char_mode[6] = "\0\0\0\0\0";
  if (!this->open_mode(mode, char_mode))
    return NULL;

  // Attempt to attach to file
  if ((file = gzdopen(fd, char_mode)) == NULL)
    return NULL;

  // On success, allocate internal buffer and set flags
  this->enable_buffer();
  io_mode = mode;
  own_fd = false;
  return this;
}

// Close gzipped file
inline gzfilebuf*
gzfilebuf::close()
{
  // Fail immediately if no file is open
  if (!this->is_open())
    return NULL;
  // Assume success
  gzfilebuf* retval = this;
  // Attempt to sync and close gzipped file
  if (this->sync() == -1)
    retval = NULL;
  if (gzclose(file) < 0)
    retval = NULL;
  // File is now gone anyway (postcondition [27.8.1.3.8])
  file = NULL;
  own_fd = false;
  // Destroy internal buffer if it exists
  this->disable_buffer();
  return retval;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Convert int open mode to mode string
inline bool
gzfilebuf::open_mode(std::ios_base::openmode mode,
                     char* c_mode) const
{
  bool testb = mode & std::ios_base::binary;
  bool testi = mode & std::ios_base::in;
  bool testo = mode & std::ios_base::out;
  bool testt = mode & std::ios_base::trunc;
  bool testa = mode & std::ios_base::app;

  // Check for valid flag combinations - see [27.8.1.3.2] (Table 92)
  // Original zfstream hardcoded the compression level to maximum here...
  // Double the time for less than 1% size improvement seems
  // excessive though - keeping it at the default level
  // To change back, just append "9" to the next three mode strings
  if (!testi && testo && !testt && !testa)
    strcpy(c_mode, "w");
  if (!testi && testo && !testt && testa)
    strcpy(c_mode, "a");
  if (!testi && testo && testt && !testa)
    strcpy(c_mode, "w");
  if (testi && !testo && !testt && !testa)
    strcpy(c_mode, "r");
  // No read/write mode yet
//  if (testi && testo && !testt && !testa)
//    strcpy(c_mode, "r+");
//  if (testi && testo && testt && !testa)
//    strcpy(c_mode, "w+");

  // Mode string should be empty for invalid combination of flags
  if (strlen(c_mode) == 0)
    return false;
  if (testb)
    strcat(c_mode, "b");
  return true;
}

// Determine number of characters in internal get buffer
inline std::streamsize
gzfilebuf::showmanyc()
{
  // Calls to underflow will fail if file not opened for reading
  if (!this->is_open() || !(io_mode & std::ios_base::in))
    return -1;
  // Make sure get area is in use
  if (this->gptr() && (this->gptr() < this->egptr()))
    return std::streamsize(this->egptr() - this->gptr());
  else
    return 0;
}

// Fill get area from gzipped file
inline gzfilebuf::int_type
gzfilebuf::underflow()
{
  // If something is left in the get area by chance, return it
  // (this shouldn't normally happen, as underflow is only supposed
  // to be called when gptr >= egptr, but it serves as error check)
  if (this->gptr() && (this->gptr() < this->egptr()))
    return traits_type::to_int_type(*(this->gptr()));

  // If the file hasn't been opened for reading, produce error
  if (!this->is_open() || !(io_mode & std::ios_base::in))
    return traits_type::eof();

  // Attempt to fill internal buffer from gzipped file
  // (buffer must be guaranteed to exist...)
  int bytes_read = gzread(file, buffer, buffer_size);
  // Indicates error or EOF
  if (bytes_read <= 0)
  {
    // Reset get area
    this->setg(buffer, buffer, buffer);
    return traits_type::eof();
  }
  // Make all bytes read from file available as get area
  this->setg(buffer, buffer, buffer + bytes_read);

  // Return next character in get area
  return traits_type::to_int_type(*(this->gptr()));
}

// Write put area to gzipped file
inline gzfilebuf::int_type
gzfilebuf::overflow(int_type c)
{
  // Determine whether put area is in use
  if (this->pbase())
  {
    // Double-check pointer range
    if (this->pptr() > this->epptr() || this->pptr() < this->pbase())
      return traits_type::eof();
    // Add extra character to buffer if not EOF
    if (!traits_type::eq_int_type(c, traits_type::eof()))
    {
      *(this->pptr()) = traits_type::to_char_type(c);
      this->pbump(1);
    }
    // Number of characters to write to file
    int bytes_to_write = this->pptr() - this->pbase();
    // Overflow doesn't fail if nothing is to be written
    if (bytes_to_write > 0)
    {
      // If the file hasn't been opened for writing, produce error
      if (!this->is_open() || !(io_mode & std::ios_base::out))
        return traits_type::eof();
      // If gzipped file won't accept all bytes written to it, fail
      if (gzwrite(file, this->pbase(), bytes_to_write) != bytes_to_write)
        return traits_type::eof();
      // Reset next pointer to point to pbase on success
      this->pbump(-bytes_to_write);
    }
  }
  // Write extra character to file if not EOF
  else if (!traits_type::eq_int_type(c, traits_type::eof()))
  {
    // If the file hasn't been opened for writing, produce error
    if (!this->is_open() || !(io_mode & std::ios_base::out))
      return traits_type::eof();
    // Impromptu char buffer (allows "unbuffered" output)
    char_type last_char = traits_type::to_char_type(c);
    // If gzipped file won't accept this character, fail
    if (gzwrite(file, &last_char, 1) != 1)
      return traits_type::eof();
  }

  // If you got here, you have succeeded (even if c was EOF)
  // The return value should therefore be non-EOF
  if (traits_type::eq_int_type(c, traits_type::eof()))
    return traits_type::not_eof(c);
  else
    return c;
}

// Assign new buffer
inline std::streambuf*
gzfilebuf::setbuf(char_type* p,
                  std::streamsize n)
{
  // First make sure stuff is sync'ed, for safety
  if (this->sync() == -1)
    return NULL;
  // If buffering is turned off on purpose via setbuf(0,0), still allocate one...
  // "Unbuffered" only really refers to put [27.8.1.4.10], while get needs at
  // least a buffer of size 1 (very inefficient though, therefore make it bigger?)
  // This follows from [27.5.2.4.3]/12 (gptr needs to point at something, it seems)
  if (!p || !n)
  {
    // Replace existing buffer (if any) with small internal buffer
    this->disable_buffer();
    buffer = NULL;
    buffer_size = 0;
    own_buffer = true;
    this->enable_buffer();
  }
  else
  {
    // Replace existing buffer (if any) with external buffer
    this->disable_buffer();
    buffer = p;
    buffer_size = n;
    own_buffer = false;
    this->enable_buffer();
  }
  return this;
}

// Write put area to gzipped file (i.e. ensures that put area is empty)
inline int
gzfilebuf::sync()
{
  return traits_type::eq_int_type(this->overflow(), traits_type::eof()) ? -1 : 0;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Allocate internal buffer
inline void
gzfilebuf::enable_buffer()
{
  // If internal buffer required, allocate one
  if (own_buffer && !buffer)
  {
    // Check for buffered vs. "unbuffered"
    if (buffer_size > 0)
    {
      // Allocate internal buffer
      buffer = new char_type[buffer_size];
      // Get area starts empty and will be expanded by underflow as need arises
      this->setg(buffer, buffer, buffer);
      // Setup entire internal buffer as put area.
      // The one-past-end pointer actually points to the last element of the buffer,
      // so that overflow(c) can safely add the extra character c to the sequence.
      // These pointers remain in place for the duration of the buffer
      this->setp(buffer, buffer + buffer_size - 1);
    }
    else
    {
      // Even in "unbuffered" case, (small?) get buffer is still required
      buffer_size = SMALLBUFSIZE;
      buffer = new char_type[buffer_size];
      this->setg(buffer, buffer, buffer);
      // "Unbuffered" means no put buffer
      this->setp(0, 0);
    }
  }
  else
  {
    // If buffer already allocated, reset buffer pointers just to make sure no
    // stale chars are lying around
    this->setg(buffer, buffer, buffer);
    this->setp(buffer, buffer + buffer_size - 1);
  }
}

// Destroy internal buffer
inline void
gzfilebuf::disable_buffer()
{
  // If internal buffer exists, deallocate it
  if (own_buffer && buffer)
  {
    // Preserve unbuffered status by zeroing size
    if (!this->pbase())
      buffer_size = 0;
    delete[] buffer;
    buffer = NULL;
    this->setg(0, 0, 0);
    this->setp(0, 0);
  }
  else
  {
    // Reset buffer pointers to initial state if external buffer exists
    this->setg(buffer, buffer, buffer);
    if (buffer)
      this->setp(buffer, buffer + buffer_size - 1);
    else
      this->setp(0, 0);
  }
}

/*****************************************************************************/

// Default constructor initializes stream buffer
inline gzifstream::gzifstream()
: std::istream(NULL), sb()
{ this->init(&sb); }

// Initialize stream buffer and open file
inline gzifstream::gzifstream(const char* name,
                       std::ios_base::openmode mode)
: std::istream(NULL), sb()
{
  this->init(&sb);
  this->open(name, mode);
}

// Initialize stream buffer and attach to file
inline gzifstream::gzifstream(int fd,
                       std::ios_base::openmode mode)
: std::istream(NULL), sb()
{
  this->init(&sb);
  this->attach(fd, mode);
}

// Open file and go into fail() state if unsuccessful
inline void
gzifstream::open(const char* name,
                 std::ios_base::openmode mode)
{
  if (!sb.open(name, mode | std::ios_base::in))
    this->setstate(std::ios_base::failbit);
  else
    this->clear();
}

// Attach to file and go into fail() state if unsuccessful
inline void
gzifstream::attach(int fd,
                   std::ios_base::openmode mode)
{
  if (!sb.attach(fd, mode | std::ios_base::in))
    this->setstate(std::ios_base::failbit);
  else
    this->clear();
}

// Close file
inline void
gzifstream::close()
{
  if (!sb.close())
    this->setstate(std::ios_base::failbit);
}

/*****************************************************************************/

// Default constructor initializes stream buffer
inline gzofstream::gzofstream()
: std::ostream(NULL), sb()
{ this->init(&sb); }

// Initialize stream buffer and open file
inline gzofstream::gzofstream(const char* name,
                       std::ios_base::openmode mode)
: std::ostream(NULL), sb()
{
  this->init(&sb);
  this->open(name, mode);
}

// Initialize stream buffer and attach to file
inline gzofstream::gzofstream(int fd,
                       std::ios_base::openmode mode)
: std::ostream(NULL), sb()
{
  this->init(&sb);
  this->attach(fd, mode);
}

// Open file and go into fail() state if unsuccessful
inline void
gzofstream::open(const char* name,
                 std::ios_base::openmode mode)
{
  if (!sb.open(name, mode | std::ios_base::out))
    this->setstate(std::ios_base::failbit);
  else
    this->clear();
}

// Attach to file and go into fail() state if unsuccessful
inline void
gzofstream::attach(int fd,
                   std::ios_base::openmode mode)
{
  if (!sb.attach(fd, mode | std::ios_base::out))
    this->setstate(std::ios_base::failbit);
  else
    this->clear();
}

// Close file
inline void
gzofstream::close()
{
  if (!sb.close())
    this->setstate(std::ios_base::failbit);
}

#endif // ZFSTREAM_H
