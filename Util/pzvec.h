// Generated by Together

/**
 * @file pzvec.h
 * @brief Templated vector implementation.
 */
// $Id: pzvec.h,v 1.17 2010-08-12 12:55:32 phil Exp $

#ifndef TVEC_H
#define TVEC_H

//#pragma interface

#include "pzreal.h"
#include "pzerror.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include <stdlib.h>


/**
 * @ingroup util
 *
 *  This class implements a simple vector storage scheme for a
 *  templated class T. The copy constructor and operator= requires the
 *  operator= to be implemented on the class T.
 */
template< class T >
class TPZVec {
   public:
      /** Creates a vector with size 0. */
      TPZVec();

      /**
       * Creates a vector of a given size.
       *
       * It will call the empty constructor on all objects of type T
       * created.
       *
       * @param size Size of the new vector.
       */
      TPZVec(const int size);

      /**
       * Creates a vector of a given size.
       *
       * It will call the empty constructor on all objects of type T
       * created copies the object copy to all elements.
       *
       * @param size Size of the new vector.
       * @param copy Model object to initialize the other objects.
       */
      TPZVec(const int size, const T& copy);

      /**
       * Creates a vector with copy constructor
       * will call the empty constructor on all objects of type T created
       * @param copy : original vector
       */
      TPZVec(const TPZVec<T> &copy);

      /**
       * destructor, will delete the storage allocated
       */
      virtual ~TPZVec();

      /**
       * will copy the vector into the current vector. Will first
       * delete the allocated storage before allocating storage for
       * the copy
       *
       * @param copy vector which will be copied
       * @return reference to the current object
       */
      TPZVec<T> &operator=(const TPZVec<T> &copy);

      /**
       * Fills the vector with a value of type T.
       *
       * @param a Element to fill the vector with.
       * @return Reference to the current object.
       */
      TPZVec<T>& operator=(const T& a);

      /**
       * access operator, will perform bounds checking unless the
       * variable NODEBUG is defined
       *
       * @param index element in the vector which will be acessed
       * @return a reference to the element specified by index\
       */
      T& operator[]( const int index ) const
      {
#ifdef DEBUG
         if( index < 0 || index >= fNElements )
         {
            PZError << __PRETTY_FUNCTION__ << " acessing element out of range.";
            PZError << "|" << std::endl;
            PZError << "+-> NElements = " << NElements() << std::endl;
            PZError << "|" << std::endl;
            PZError << "+-> Index = " << index << std::endl;
            DebugStop();
            exit( -1 );
         }
#endif
         return fStore[ index ];
      }

      /// Extraction operator.

			friend std::ostream& operator<<( std::ostream& Out, const TPZVec< T >& v )
      {
         std::streamsize width = Out.width();

	 const char* sep = ( width == 0 ? ", " : "" );

         int size = v.NElements();
		  
		  if(size) Out << std::setw(width) << v.fStore[0];

         for( int ii = 1; ii < size; ii++ )
         {
            Out << std::setw( width ) << sep << v.fStore[ ii ];
         }

         return Out;
      }

      /** Casting operator.
       *
       *  @return The fStore pointer.
       */
      operator T*() const { return fStore; }

      /**
       * Will fill the elements of the vector with a copy object.
       * @param copy object which will be copied
       * @param from first index which will be overwritten
       * @param numelem number of elements which will be overwritten
       */
      void Fill(const T& copy, const int from=0, const int numelem=-1);

      /**
       * Returns the number of elements of the vector
       * @return number of elements used by the vector
       */
      inline int NElements() const { return fNElements; }

	/**
	 * Returns the number of elements of the vector
	 * @return number of elements used by the vector
	 */
	inline int size() const { return fNElements; }
	
	/**
       * Resizes the vector object reallocating the necessary storage,
       * copying the existing objects to the new storage.
       *
       * @param newsize size of the vector
       * @param object object used to initialize the new members
       */
      virtual void Resize(const int newsize,const T& object);

      /**
       * Resizes the vector object reallocating the necessary storage,
       * copying the existing objects to the new storage. The new
       * members are not initialized.
       *
       * @param newsize size of the vector
       */
      virtual void Resize(const int newsize);

      /**
       * Print the structural information of the vector object to the
       * output stream. This method will not print the objects
       * themselves!
       */
	void Print(std::ostream &out = std::cout);

   protected:
      /**
       * Allocated storage for the vector object
       */
      T* fStore;

      /**
       * Number of elements of the vector object
       */
      int fNElements;
};

//--| IMPLEMENTATION |----------------------------------------------------------

/* template <class T> */
/* T &TPZVec<T>::operator[](const int index) const{ */
/* #ifndef NODEBUG */
/*        if(index <0 || index >= fNElements) { */
/*          cout << "TPZVec acessing element out of range\n"; */
/*          exit(-1); */
/*        } */
/* #endif */
/*              return fStore[index]; */
/* } */

template< class T >
inline TPZVec<T>::TPZVec() : fStore( 0 ), fNElements( 0 )
{
   // NOTHING TO DO HERE!
}

template< class T >
TPZVec<T>::TPZVec( const int size ) : fStore( 0 )
{
#ifndef NOTDEBUG
   if( size < 0 )
   {
      PZError << "TPZVec constructor. Bad parameter size, then size = 0."
              << std::endl;
      PZError.flush();
      fNElements = 0;
      return;
   }
#endif

   // If a positive value was requested, allocate it.
   if( size > 0 )
   {
      fStore = new T[ size ];
   }

   // Note that even 0 sized vectors are allowed.
   fNElements = size;
}

template< class T >
TPZVec<T>::TPZVec( const int size, const T& copy ) : fStore( 0 )
{
#ifndef NOTDEBUG
   if( size < 0 )
   {
      PZError << "TPZVec constructor. Bad parameter size, then size = 0."
              << std::endl;
      PZError.flush();
      fNElements = 0;
      return;
   }
#endif

   if( size )
   {
      fStore = new T[size];
   }

   fNElements = size;

   for( int i = 0; i < size; i++ )
   {
      fStore[i] = copy;
   }
}

template< class T >
TPZVec<T>::TPZVec(const TPZVec<T> &copy){
   fStore = 0;

   if( copy.fNElements > 0 )
      fStore = new T[copy.fNElements];
   else
      fStore = 0;

   for(int i=0; i<copy.fNElements; i++)
      fStore[i]=copy.fStore[i];

   fNElements = copy.fNElements;
}


template<class T>
inline TPZVec<T>::~TPZVec() {
   if( fStore )
   {
      delete [] fStore;
   }
}

template< class T >
TPZVec<T> &TPZVec<T>::operator=(const TPZVec<T> &copy){
   if(this == &copy) return *this;

   Resize(copy.NElements());

   for(int i=0; i<copy.fNElements; i++)
      fStore[i]=copy.fStore[i];

   fNElements = copy.fNElements;

   return *this;
}

// OPENED QUESTION: what to do with 0 size vectors??? Cantao (2002.01.09)
template< class T >
TPZVec<T>& TPZVec<T>::operator=( const T& a )
{
   T* end = fStore + fNElements;

   for( T* walk = fStore; walk < end; *walk++ = a );

   return *this;
}


template< class T >
void TPZVec<T>::Resize(const int newsize,const T& object){
#ifndef NOTDEBUG
   if(newsize<0) {
      PZError << "TPZVec::Resize. Bad parameter newsize." << std::endl;
      PZError.flush();
      return;
   }
#endif
   if(newsize == fNElements) return;
   T* newstore;
   if(newsize) newstore = new T[newsize];
   else newstore = 0;
   int large = (fNElements < newsize) ? fNElements : newsize;
   int i;
   for(i=0; i<large; i++) {
      newstore[i] = fStore[i];
   }
   for(;i<newsize;i++) {   // then only to case : large=fNElement < newsize
      newstore[i] = object;
   }
   delete[] fStore;
   fStore = newstore;
   fNElements = newsize;//cedric 20/11/99 e 29/04/00
}

template< class T >
void TPZVec<T>::Resize(const int newsize){
#ifndef NOTDEBUG
   if(newsize<0) {
      PZError << "TPZVec::Resize. Bad parameter newsize." << newsize <<  std::endl;
      PZError.flush();
      DebugStop();
      return;
   }
#endif

   if(newsize == fNElements) return;
   if (newsize == 0) {
      fNElements = 0;
      delete[] fStore;
      fStore = 0;
      return;
   }
   T *newstore = new T[newsize];
   int large = (fNElements < newsize) ? fNElements : newsize;
   int i;
   for(i=0; i<large; i++) {
      newstore[i] = fStore[i];
   }
   if (fStore) delete[] fStore;
   fStore = newstore;
   fNElements = newsize;
}


template< class T >
void TPZVec<T>::Fill(const T& copy, const int from, const int numelem){
#ifndef NOTDEBUG
   if(numelem<0 && numelem != -1) {
      PZError << "TPZVec::Fill" << std::endl
              << "It's negative parameter numelem, then numelem = "
              << fNElements << std::endl;

      PZError.flush();
   }
#endif

   int first = (from < 0) ? 0 : from;
   int nel = numelem;
   first = (first > fNElements) ? fNElements : first;
   if (nel < 0) nel = fNElements;
   int last = (from+nel > fNElements) ? fNElements : from+nel;

   for(int i=first; i<last; i++)
      fStore[i] = copy;
}

template< class T >
inline void TPZVec<T>::Print(std::ostream &out)
{
   out << std::endl << "Number of elements = " << fNElements;

   //   for(int i=0;i<fNElements;i++)
   //   out << setw( 14 ) << setprecision( 6 ) << fStore[i];
}

inline std::ostream &operator<<(std::ostream &out, std::pair<int,int> &element)
{
	out << element.first << "|" << element.second;
	return out;
}

#endif

//--| PZ |----------------------------------------------------------------------
