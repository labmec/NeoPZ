/**
 * @file pzstring.h
 * @brief String implementation.
 */
// $Id: pzstring.h,v 1.1.1.1 2003-02-04 16:45:27 cantao Exp $

#ifndef PZSTRING_H
#define PZSTRING_H

#include "pzstack.h" 
#include "pzmanvector.h"

#include <string.h>

class TPZString: public TPZStack< char >
{
   public:

      /** Default Constructor. */
      inline TPZString();

      /** Constructs a TPZString object based on a null char ended string. */
      inline TPZString(const char * source);

      /** Initializes the TPZString object with a predefined size. */
      inline TPZString(const int size);

      /** Single char constructor. */
      inline TPZString(const char chr);
	
      /** Destructor. */
      ~TPZString(){}

      /**
       * Concatenates the 'this' string with another character string.
       * Original strings are left unchanged.
       */
      TPZString operator+(const char * increment) const;

      /**
       * Concatenates the 'this' string with a character.
       * Original string and character are left unchanged.
       */
      TPZString operator+(const TPZString & increment) const;

      /**
       * Appends a string at the tail. Resizes the TPZString if necessary.
       */
      void operator+=(const char * increment);

      /**
       * Appends a character at the end. Resizes if necessary.
       */
      void operator+=(const char increment);

      /**
       *  operator =. Resizes if necessary.
       */
      inline void operator=(const char * source)
      {
	 int len = strlen(source);

	 if(NElements()<(len+1))
	 {
	    Resize(len+1);
	 }

	 strcpy(fStore, source);
      }

      /**
       * Appends a character at the end. Resizes if necessary.
       */
      void Append(const char TailIncrement);

      /**
       * Appends a string at the tail. Resizes the TPZString if necessary.
       */
      void Append(const char * TailIncrement);

      /**
       * Explicitly convertes a TPZString into a const null ended
       * char string .
       */
      inline const char * Str() const;

      /**
       * Implicit conversion
       */
      inline operator const char *() const;

      /**
       * Similar to strlen(string). Also returns the number of
       * non-null characters. 
       */
      inline int Length() const;

      /**
       * Returns a subset string.
       *
       * If start and end are the same, a substring containing one
       * character is returned.
       * <li>
       * <ul> if start <=0, start is assumed to be 0;
       * <ul> if end < start, null pointer is returned;
       * <ul> if end > length, end is assumed to be equal to length.
       * </li>
       *
       * @param start start of the string (zero-based)
       * @param end end of the string (zero-based)
       */
      TPZString SubStr(const int start,const int end) const;

      /** Empties the string. */
      void Empty();

      /**
       * Internally allocates the exact string size (length + 1)
       * to store it.
       */
      void Optimize();
};

typedef TPZVec< TPZString > TPZText;

#endif

//--| PZ |----------------------------------------------------------------------
