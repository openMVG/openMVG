// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image.hpp"

#include <cstring>
#include <iostream>
#include <cmath>

extern "C" {
  #include "jpeglib.h"
  #include "png.h"
}

using namespace std;

namespace openMVG {

static bool CmpFormatExt(const char *a, const char *b) {
  size_t len_a = strlen(a);
  size_t len_b = strlen(b);
  if (len_a != len_b) return false;
  for (size_t i = 0; i < len_a; ++i)
    if (tolower(a[i]) != tolower(b[i]))
      return false;
  return true;
}

Format GetFormat(const char *c) {
  const char *p = strrchr (c, '.');

  if (p == NULL)
    return Unknown;

  if (CmpFormatExt(p, ".png")) return Png;
  if (CmpFormatExt(p, ".ppm")) return Pnm;
  if (CmpFormatExt(p, ".pgm")) return Pnm;
  if (CmpFormatExt(p, ".pbm")) return Pnm;
  if (CmpFormatExt(p, ".pnm")) return Pnm;
  if (CmpFormatExt(p, ".jpg")) return Jpg;
  if (CmpFormatExt(p, ".jpeg")) return Jpg;

  cerr << "Error: Couldn't open " << c << " Unknown file format" << std::endl;
  return Unknown;
}

int ReadImage(const char *filename,
              vector<unsigned char> * ptr,
              int * w,
              int * h,
              int * depth){
  Format f = GetFormat(filename);

  switch (f) {
    case Pnm:
      return ReadPnm(filename, ptr, w, h, depth);
    case Png:
      return ReadPng(filename, ptr, w, h, depth);
    case Jpg:
      return ReadJpg(filename, ptr, w, h, depth);
    default:
      return 0;
  };
}

int WriteImage(const char * filename,
              const vector<unsigned char> & ptr,
              int w,
              int h,
              int depth){
  Format f = GetFormat(filename);

  switch (f) {
    case Pnm:
      return WritePnm(filename, ptr, w, h, depth);
    case Png:
      return WritePng(filename, ptr, w, h, depth);
    case Jpg:
      return WriteJpg(filename, ptr, w, h, depth);
    default:
      return 0;
  };
}

int ReadJpg(const char * filename,
            vector<unsigned char> * ptr,
            int * w,
            int * h,
            int * depth) {

  FILE *file = fopen(filename, "rb");
  if (!file) {
    cerr << "Error: Couldn't open " << filename << " fopen returned 0";
    return 0;
  }
  int res = ReadJpgStream(file, ptr, w, h, depth);
  fclose(file);
  return res;
}

struct my_error_mgr {
  struct jpeg_error_mgr pub;
  jmp_buf setjmp_buffer;
};

METHODDEF(void)
jpeg_error (j_common_ptr cinfo)
{
  my_error_mgr *myerr = (my_error_mgr*) (cinfo->err);
  (*cinfo->err->output_message) (cinfo);
  longjmp(myerr->setjmp_buffer, 1);
}

int ReadJpgStream(FILE * file,
                  vector<unsigned char> * ptr,
                  int * w,
                  int * h,
                  int * depth) {
  jpeg_decompress_struct cinfo;
  struct my_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr.pub);
  jerr.pub.error_exit = &jpeg_error;

  if (setjmp(jerr.setjmp_buffer)) {
    jpeg_destroy_decompress(&cinfo);
    return 0;
  }

  jpeg_create_decompress(&cinfo);
  jpeg_stdio_src(&cinfo, file);
  jpeg_read_header(&cinfo, TRUE);
  jpeg_start_decompress(&cinfo);

  int row_stride = cinfo.output_width * cinfo.output_components;

  *h = cinfo.output_height;
  *w = cinfo.output_width;
  *depth = cinfo.output_components;
  ptr->resize((*h)*(*w)*(*depth));

  unsigned char *ptrCpy = &(*ptr)[0];

  while (cinfo.output_scanline < cinfo.output_height) {
    JSAMPROW scanline[1] = { ptrCpy };
    jpeg_read_scanlines(&cinfo, scanline, 1);
    ptrCpy += row_stride;
  }

  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);
  return 1;
}


int WriteJpg(const char * filename,
             const vector<unsigned char> & array,
             int w,
             int h,
             int depth,
             int quality) {
  FILE *file = fopen(filename, "wb");
  if (!file) {
    cerr << "Error: Couldn't open " << filename << " fopen returned 0";
    return 0;
  }
  int res = WriteJpgStream(file, array, w, h, depth, quality);
  fclose(file);
  return res;
}

int WriteJpgStream(FILE *file,
                   const vector<unsigned char> & array,
                   int w,
                   int h,
                   int depth,
                   int quality) {
  if (quality < 0 || quality > 100)
    cerr << "Error: The quality parameter should be between 0 and 100";

  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;

  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo, file);

  cinfo.image_width = w;
  cinfo.image_height = h;
  cinfo.input_components = depth;

  if (cinfo.input_components==3) {
    cinfo.in_color_space = JCS_RGB;
  } else if (cinfo.input_components==1) {
    cinfo.in_color_space = JCS_GRAYSCALE;
  } else {
    cerr << "Error: Unsupported number of channels in file";
    jpeg_destroy_compress(&cinfo);
    return 0;
  }

  jpeg_set_defaults(&cinfo);
  jpeg_set_quality(&cinfo, quality, TRUE);
  jpeg_start_compress(&cinfo, TRUE);

  const unsigned char *ptr = &array[0];
  int row_bytes = cinfo.image_width*cinfo.input_components;

  JSAMPLE *row = new JSAMPLE[row_bytes];

  while (cinfo.next_scanline < cinfo.image_height) {
    memcpy(&row[0], &ptr[0], row_bytes*sizeof(unsigned char));
    //int i;
    //for (i = 0; i < row_bytes; ++i)
    //	row[i] = ptr[i];
    jpeg_write_scanlines(&cinfo, &row, 1);
    ptr += row_bytes;
  }

  delete [] row;

  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);
  return 1;
}

int ReadPng(const char *filename,
            vector<unsigned char> * ptr,
            int * w,
            int * h,
            int * depth) {
  FILE *file = fopen(filename, "rb");
  if (!file) {
    cerr << "Error: Couldn't open " << filename << " fopen returned 0";
    return 0;
  }
  int res = ReadPngStream(file, ptr, w, h, depth);
  fclose(file);
  return res;
}

int ReadPngStream(FILE *file,
                  vector<unsigned char> * ptr,
                  int * w,
                  int * h,
                  int * depth)  {

  // first check the eight byte PNG signature
  png_byte  pbSig[8];
  size_t readcnt = fread(pbSig, 1, 8, file);
  (void) readcnt;
  if (png_sig_cmp(pbSig, 0, 8))
  {
    return 0;
  }

  // create the two png(-info) structures
  png_structp png_ptr = NULL;
  png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL,
    (png_error_ptr)NULL, (png_error_ptr)NULL);
  if (!png_ptr)
  {
    return 0;
  }
  png_infop info_ptr = NULL;
  info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr)
  {
    png_destroy_read_struct(&png_ptr, NULL, NULL);
    return 0;
  }

  // initialize the png structure
  png_init_io(png_ptr, file);
  png_set_sig_bytes(png_ptr, 8);

  // read all PNG info up to image data

  png_read_info(png_ptr, info_ptr);

  // get width, height, bit-depth and color-type
  png_uint_32 wPNG, hPNG;
  int                 iBitDepth;
  int                 iColorType;
  png_get_IHDR(png_ptr, info_ptr, &wPNG, &hPNG, &iBitDepth,
    &iColorType, NULL, NULL, NULL);

  // expand images of all color-type to 8-bit

  if (iColorType == PNG_COLOR_TYPE_PALETTE)
    png_set_expand(png_ptr);
  if (iBitDepth < 8)
    png_set_expand(png_ptr);
  if (png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS))
    png_set_expand(png_ptr);
  if (iBitDepth == 16) // convert 16-bit to 8-bit on the fly
    png_set_strip_16(png_ptr);

  double dGamma;
  // if required set gamma conversion
  if (png_get_gAMA(png_ptr, info_ptr, &dGamma))
    png_set_gamma(png_ptr, (double) 2.2, dGamma);

  // after the transformations are registered, update info_ptr data

  png_read_update_info(png_ptr, info_ptr);

  // get again width, height and the new bit-depth and color-type

  png_get_IHDR(png_ptr, info_ptr, &wPNG, &hPNG, &iBitDepth,
    &iColorType, NULL, NULL, NULL);

  // Get number of byte along a tow
  png_uint_32         ulRowBytes;
  ulRowBytes = png_get_rowbytes(png_ptr, info_ptr);

  // and allocate memory for an array of row-pointers
  png_byte   **ppbRowPointers = NULL;
  if ((ppbRowPointers = (png_bytepp) malloc(hPNG
    * sizeof(png_bytep))) == NULL)
  {
    std::cerr << "PNG: out of memory" << std::endl;
    return 0;
  }

  *w = wPNG;
  *h = hPNG;
  *depth = png_get_channels(png_ptr, info_ptr);

  // now we can allocate memory to store the image
  ptr->resize((*h)*(*w)*(*depth));

  // set the individual row-pointers to point at the correct offsets
  for (int i = 0; i < hPNG; i++)
    ppbRowPointers[i] = &((*ptr)[0]) + i * ulRowBytes;

  // now we can go ahead and just read the whole image
  png_read_image(png_ptr, ppbRowPointers);

  // read the additional chunks in the PNG file (not really needed)
  png_read_end(png_ptr, NULL);

  free (ppbRowPointers);

  png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
  return 1;
}

int WritePng(const char * filename,
             const vector<unsigned char> & ptr,
             int w,
             int h,
             int depth) {
  FILE *file = fopen(filename, "wb");
  if (!file) {
    cerr << "Error: Couldn't open " << filename << " fopen returned 0";
    return 0;
  }
  int res = WritePngStream(file, ptr, w, h, depth);
  fclose(file);
  return res;
}

int WritePngStream(FILE * file,
                   const vector<unsigned char> & ptr,
                   int w,
                   int h,
                   int depth) {
  png_structp png_ptr =
      png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

  if (!png_ptr)
    return 0;

  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr)
    return 0;

  png_init_io(png_ptr, file);

  // color types are defined at png.h:841+.
  char colour;
  switch(depth)
  {
    case 4: colour = PNG_COLOR_TYPE_RGBA;
      break;
    case 3: colour = PNG_COLOR_TYPE_RGB;
      break;
    case 1: colour = PNG_COLOR_TYPE_GRAY;
      break;
    default:
      return 0;
  }

  png_set_IHDR(png_ptr, info_ptr, w, h,
      8, colour, PNG_INTERLACE_NONE,
      PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

  png_write_info(png_ptr, info_ptr);

  png_bytep *row_pointers =
      (png_bytep*) malloc(sizeof(png_bytep) * depth * h);

  for (int y = 0; y < h; ++y)
    row_pointers[y] = (png_byte*) (&ptr[0]) + w * depth * y;

  png_write_image(png_ptr, row_pointers);
  free(row_pointers);
  png_write_end(png_ptr, NULL);
  png_destroy_write_struct(&png_ptr, &info_ptr);
  return 1;
}

int ReadPnm(const char * filename,
            vector<unsigned char> * array,
            int * w,
            int * h,
            int * depth)  {
  FILE *file = fopen(filename, "rb");
  if (!file) {
    cerr << "Error: Couldn't open " << filename << " fopen returned 0";
    return 0;
  }
  int res = ReadPnmStream(file, array, w, h, depth);
  fclose(file);
  return res;
}


// Comment handling as per the description provided at
//   http://netpbm.sourceforge.net/doc/pgm.html
// and http://netpbm.sourceforge.net/doc/pbm.html
int ReadPnmStream(FILE *file,
                  vector<unsigned char> * array,
                  int * w,
                  int * h,
                  int * depth) {

  const int NUM_VALUES = 3;
  const int INT_BUFFER_SIZE = 256;

  int magicnumber;
  char intBuffer[INT_BUFFER_SIZE];
  int values[NUM_VALUES], valuesIndex = 0, intIndex = 0, inToken = 0;
  size_t res;

  // Check magic number.
  res = size_t(fscanf(file, "P%d", &magicnumber));
  if (res != 1) {
    return 0;
  }
  if (magicnumber == 5) {
    *depth = 1;
  } else if (magicnumber == 6) {
    *depth = 3;
  } else {
    return 0;
  }

  // the following loop parses the PNM header one character at a time, looking
  // for the int tokens width, height and maxValues (in that order), and
  // discarding all comment (everything from '#' to '\n' inclusive), where
  // comments *may occur inside tokens*. Each token must be terminate with a
  // whitespace character, and only one whitespace char is eaten after the
  // third int token is parsed.
  while (valuesIndex < NUM_VALUES) {
    char nextChar ;
    res = fread(&nextChar,1,1,file);
    if (res == 0) return 0; // read failed, EOF?

    if (isspace(nextChar)) {
      if (inToken) { // we were reading a token, so this white space delimits it
        inToken = 0;
        intBuffer[intIndex] = 0 ; // NULL-terminate the string
        values[valuesIndex++] = atoi(intBuffer);
        intIndex = 0; // reset for next int token
        // to conform with current image class
        if (valuesIndex == 3 && values[2] > 255) return 0;
      }
    }
    else if (isdigit(nextChar)) {
      inToken = 1 ; // in case it's not already set
      intBuffer[intIndex++] = nextChar ;
      if (intIndex == INT_BUFFER_SIZE) // tokens should never be this long
        return 0;
    }
    else if (nextChar == '#') {
      do { // eat all characters from input stream until newline
        res = fread(&nextChar,1,1,file);
      } while (res == 1 && nextChar != '\n');
      if (res == 0) return 0; // read failed, EOF?
    }
    else {
      // Encountered a non-whitespace, non-digit outside a comment - bail out.
      return 0;
    }
  }

  // Read pixels.
  (*array).resize( values[1] * values[0] * (*depth));
  *w = values[0];
  *h = values[1];
  res = fread( &(*array)[0], 1, array->size(), file);
  if (res != array->size()) {
    return 0;
  }
  return 1;
}

int WritePnm(const char * filename,
              const vector<unsigned char> & array,
              int w,
              int h,
              int depth) {
  FILE *file = fopen(filename, "wb");
  if (!file) {
    cerr << "Error: Couldn't open " << filename << " fopen returned 0";
    return 0;
  }
  int res = WritePnmStream(file, array, w, h, depth);
  fclose(file);
  return res;
}


int WritePnmStream(FILE * file,
                   const vector<unsigned char> & array,
                   int w,
                   int h,
                   int depth) {
  int res;

  // Write magic number.
  if (depth == 1) {
    fprintf(file, "P5\n");
  } else if (depth == 3) {
    fprintf(file, "P6\n");
  } else {
    return 0;
  }

  // Write sizes.
  fprintf(file, "%d %d %d\n", w, h, 255);

  // Write pixels.
  res = fwrite( &array[0], 1, static_cast<int>(array.size()), file);
  if (res != array.size()) {
    return 0;
  }
  return 1;
}

}  // namespace openMVG
