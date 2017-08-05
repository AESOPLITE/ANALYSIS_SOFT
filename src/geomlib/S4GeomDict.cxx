// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME S4GeomDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "THelicalTrack.h"
#include "TCircle.h"
#include "TCutCone.h"
#include "TVCurve.h"
#include "TCylinder.h"
#include "TPlane.h"
#include "THype.h"
#include "TStraightTrack.h"
#include "TTube.h"
#include "TVSolid.h"
#include "TVSurface.h"
#include "TVTrack.h"
#include "TBField.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TVCurve(void *p = 0);
   static void *newArray_TVCurve(Long_t size, void *p);
   static void delete_TVCurve(void *p);
   static void deleteArray_TVCurve(void *p);
   static void destruct_TVCurve(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TVCurve*)
   {
      ::TVCurve *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TVCurve >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TVCurve", ::TVCurve::Class_Version(), "TVCurve.h", 25,
                  typeid(::TVCurve), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TVCurve::Dictionary, isa_proxy, 4,
                  sizeof(::TVCurve) );
      instance.SetNew(&new_TVCurve);
      instance.SetNewArray(&newArray_TVCurve);
      instance.SetDelete(&delete_TVCurve);
      instance.SetDeleteArray(&deleteArray_TVCurve);
      instance.SetDestructor(&destruct_TVCurve);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TVCurve*)
   {
      return GenerateInitInstanceLocal((::TVCurve*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TVCurve*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_TVTrack(void *p);
   static void deleteArray_TVTrack(void *p);
   static void destruct_TVTrack(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TVTrack*)
   {
      ::TVTrack *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TVTrack >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TVTrack", ::TVTrack::Class_Version(), "TVTrack.h", 40,
                  typeid(::TVTrack), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TVTrack::Dictionary, isa_proxy, 4,
                  sizeof(::TVTrack) );
      instance.SetDelete(&delete_TVTrack);
      instance.SetDeleteArray(&deleteArray_TVTrack);
      instance.SetDestructor(&destruct_TVTrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TVTrack*)
   {
      return GenerateInitInstanceLocal((::TVTrack*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TVTrack*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_THelicalTrack(void *p = 0);
   static void *newArray_THelicalTrack(Long_t size, void *p);
   static void delete_THelicalTrack(void *p);
   static void deleteArray_THelicalTrack(void *p);
   static void destruct_THelicalTrack(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::THelicalTrack*)
   {
      ::THelicalTrack *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::THelicalTrack >(0);
      static ::ROOT::TGenericClassInfo 
         instance("THelicalTrack", ::THelicalTrack::Class_Version(), "THelicalTrack.h", 39,
                  typeid(::THelicalTrack), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::THelicalTrack::Dictionary, isa_proxy, 4,
                  sizeof(::THelicalTrack) );
      instance.SetNew(&new_THelicalTrack);
      instance.SetNewArray(&newArray_THelicalTrack);
      instance.SetDelete(&delete_THelicalTrack);
      instance.SetDeleteArray(&deleteArray_THelicalTrack);
      instance.SetDestructor(&destruct_THelicalTrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::THelicalTrack*)
   {
      return GenerateInitInstanceLocal((::THelicalTrack*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::THelicalTrack*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TCircle(void *p = 0);
   static void *newArray_TCircle(Long_t size, void *p);
   static void delete_TCircle(void *p);
   static void deleteArray_TCircle(void *p);
   static void destruct_TCircle(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TCircle*)
   {
      ::TCircle *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TCircle >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TCircle", ::TCircle::Class_Version(), "TCircle.h", 29,
                  typeid(::TCircle), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TCircle::Dictionary, isa_proxy, 4,
                  sizeof(::TCircle) );
      instance.SetNew(&new_TCircle);
      instance.SetNewArray(&newArray_TCircle);
      instance.SetDelete(&delete_TCircle);
      instance.SetDeleteArray(&deleteArray_TCircle);
      instance.SetDestructor(&destruct_TCircle);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TCircle*)
   {
      return GenerateInitInstanceLocal((::TCircle*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TCircle*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_TVSurface(void *p);
   static void deleteArray_TVSurface(void *p);
   static void destruct_TVSurface(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TVSurface*)
   {
      ::TVSurface *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TVSurface >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TVSurface", ::TVSurface::Class_Version(), "TVSurface.h", 33,
                  typeid(::TVSurface), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TVSurface::Dictionary, isa_proxy, 4,
                  sizeof(::TVSurface) );
      instance.SetDelete(&delete_TVSurface);
      instance.SetDeleteArray(&deleteArray_TVSurface);
      instance.SetDestructor(&destruct_TVSurface);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TVSurface*)
   {
      return GenerateInitInstanceLocal((::TVSurface*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TVSurface*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TCutCone(void *p = 0);
   static void *newArray_TCutCone(Long_t size, void *p);
   static void delete_TCutCone(void *p);
   static void deleteArray_TCutCone(void *p);
   static void destruct_TCutCone(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TCutCone*)
   {
      ::TCutCone *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TCutCone >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TCutCone", ::TCutCone::Class_Version(), "TCutCone.h", 34,
                  typeid(::TCutCone), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TCutCone::Dictionary, isa_proxy, 4,
                  sizeof(::TCutCone) );
      instance.SetNew(&new_TCutCone);
      instance.SetNewArray(&newArray_TCutCone);
      instance.SetDelete(&delete_TCutCone);
      instance.SetDeleteArray(&deleteArray_TCutCone);
      instance.SetDestructor(&destruct_TCutCone);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TCutCone*)
   {
      return GenerateInitInstanceLocal((::TCutCone*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TCutCone*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TCylinder(void *p = 0);
   static void *newArray_TCylinder(Long_t size, void *p);
   static void delete_TCylinder(void *p);
   static void deleteArray_TCylinder(void *p);
   static void destruct_TCylinder(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TCylinder*)
   {
      ::TCylinder *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TCylinder >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TCylinder", ::TCylinder::Class_Version(), "TCylinder.h", 33,
                  typeid(::TCylinder), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TCylinder::Dictionary, isa_proxy, 4,
                  sizeof(::TCylinder) );
      instance.SetNew(&new_TCylinder);
      instance.SetNewArray(&newArray_TCylinder);
      instance.SetDelete(&delete_TCylinder);
      instance.SetDeleteArray(&deleteArray_TCylinder);
      instance.SetDestructor(&destruct_TCylinder);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TCylinder*)
   {
      return GenerateInitInstanceLocal((::TCylinder*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TCylinder*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TPlane(void *p = 0);
   static void *newArray_TPlane(Long_t size, void *p);
   static void delete_TPlane(void *p);
   static void deleteArray_TPlane(void *p);
   static void destruct_TPlane(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TPlane*)
   {
      ::TPlane *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TPlane >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TPlane", ::TPlane::Class_Version(), "TPlane.h", 31,
                  typeid(::TPlane), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TPlane::Dictionary, isa_proxy, 4,
                  sizeof(::TPlane) );
      instance.SetNew(&new_TPlane);
      instance.SetNewArray(&newArray_TPlane);
      instance.SetDelete(&delete_TPlane);
      instance.SetDeleteArray(&deleteArray_TPlane);
      instance.SetDestructor(&destruct_TPlane);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TPlane*)
   {
      return GenerateInitInstanceLocal((::TPlane*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TPlane*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_THype(void *p = 0);
   static void *newArray_THype(Long_t size, void *p);
   static void delete_THype(void *p);
   static void deleteArray_THype(void *p);
   static void destruct_THype(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::THype*)
   {
      ::THype *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::THype >(0);
      static ::ROOT::TGenericClassInfo 
         instance("THype", ::THype::Class_Version(), "THype.h", 31,
                  typeid(::THype), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::THype::Dictionary, isa_proxy, 4,
                  sizeof(::THype) );
      instance.SetNew(&new_THype);
      instance.SetNewArray(&newArray_THype);
      instance.SetDelete(&delete_THype);
      instance.SetDeleteArray(&deleteArray_THype);
      instance.SetDestructor(&destruct_THype);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::THype*)
   {
      return GenerateInitInstanceLocal((::THype*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::THype*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TStraightTrack(void *p = 0);
   static void *newArray_TStraightTrack(Long_t size, void *p);
   static void delete_TStraightTrack(void *p);
   static void deleteArray_TStraightTrack(void *p);
   static void destruct_TStraightTrack(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TStraightTrack*)
   {
      ::TStraightTrack *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TStraightTrack >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TStraightTrack", ::TStraightTrack::Class_Version(), "TStraightTrack.h", 35,
                  typeid(::TStraightTrack), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TStraightTrack::Dictionary, isa_proxy, 4,
                  sizeof(::TStraightTrack) );
      instance.SetNew(&new_TStraightTrack);
      instance.SetNewArray(&newArray_TStraightTrack);
      instance.SetDelete(&delete_TStraightTrack);
      instance.SetDeleteArray(&deleteArray_TStraightTrack);
      instance.SetDestructor(&destruct_TStraightTrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TStraightTrack*)
   {
      return GenerateInitInstanceLocal((::TStraightTrack*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TStraightTrack*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TVSolid(void *p = 0);
   static void *newArray_TVSolid(Long_t size, void *p);
   static void delete_TVSolid(void *p);
   static void deleteArray_TVSolid(void *p);
   static void destruct_TVSolid(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TVSolid*)
   {
      ::TVSolid *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TVSolid >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TVSolid", ::TVSolid::Class_Version(), "TVSolid.h", 25,
                  typeid(::TVSolid), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TVSolid::Dictionary, isa_proxy, 4,
                  sizeof(::TVSolid) );
      instance.SetNew(&new_TVSolid);
      instance.SetNewArray(&newArray_TVSolid);
      instance.SetDelete(&delete_TVSolid);
      instance.SetDeleteArray(&deleteArray_TVSolid);
      instance.SetDestructor(&destruct_TVSolid);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TVSolid*)
   {
      return GenerateInitInstanceLocal((::TVSolid*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TVSolid*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TTube(void *p = 0);
   static void *newArray_TTube(Long_t size, void *p);
   static void delete_TTube(void *p);
   static void deleteArray_TTube(void *p);
   static void destruct_TTube(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TTube*)
   {
      ::TTube *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TTube >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TTube", ::TTube::Class_Version(), "TTube.h", 29,
                  typeid(::TTube), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TTube::Dictionary, isa_proxy, 4,
                  sizeof(::TTube) );
      instance.SetNew(&new_TTube);
      instance.SetNewArray(&newArray_TTube);
      instance.SetDelete(&delete_TTube);
      instance.SetDeleteArray(&deleteArray_TTube);
      instance.SetDestructor(&destruct_TTube);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TTube*)
   {
      return GenerateInitInstanceLocal((::TTube*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TTube*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TBField(void *p = 0);
   static void *newArray_TBField(Long_t size, void *p);
   static void delete_TBField(void *p);
   static void deleteArray_TBField(void *p);
   static void destruct_TBField(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TBField*)
   {
      ::TBField *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TBField >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TBField", ::TBField::Class_Version(), "TBField.h", 34,
                  typeid(::TBField), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TBField::Dictionary, isa_proxy, 4,
                  sizeof(::TBField) );
      instance.SetNew(&new_TBField);
      instance.SetNewArray(&newArray_TBField);
      instance.SetDelete(&delete_TBField);
      instance.SetDeleteArray(&deleteArray_TBField);
      instance.SetDestructor(&destruct_TBField);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TBField*)
   {
      return GenerateInitInstanceLocal((::TBField*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TBField*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TVCurve::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TVCurve::Class_Name()
{
   return "TVCurve";
}

//______________________________________________________________________________
const char *TVCurve::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TVCurve*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TVCurve::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TVCurve*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TVCurve::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TVCurve*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TVCurve::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TVCurve*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TVTrack::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TVTrack::Class_Name()
{
   return "TVTrack";
}

//______________________________________________________________________________
const char *TVTrack::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TVTrack*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TVTrack::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TVTrack*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TVTrack::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TVTrack*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TVTrack::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TVTrack*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr THelicalTrack::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *THelicalTrack::Class_Name()
{
   return "THelicalTrack";
}

//______________________________________________________________________________
const char *THelicalTrack::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::THelicalTrack*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int THelicalTrack::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::THelicalTrack*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *THelicalTrack::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::THelicalTrack*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *THelicalTrack::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::THelicalTrack*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TCircle::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TCircle::Class_Name()
{
   return "TCircle";
}

//______________________________________________________________________________
const char *TCircle::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TCircle*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TCircle::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TCircle*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TCircle::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TCircle*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TCircle::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TCircle*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TVSurface::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TVSurface::Class_Name()
{
   return "TVSurface";
}

//______________________________________________________________________________
const char *TVSurface::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TVSurface*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TVSurface::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TVSurface*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TVSurface::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TVSurface*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TVSurface::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TVSurface*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TCutCone::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TCutCone::Class_Name()
{
   return "TCutCone";
}

//______________________________________________________________________________
const char *TCutCone::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TCutCone*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TCutCone::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TCutCone*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TCutCone::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TCutCone*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TCutCone::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TCutCone*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TCylinder::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TCylinder::Class_Name()
{
   return "TCylinder";
}

//______________________________________________________________________________
const char *TCylinder::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TCylinder*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TCylinder::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TCylinder*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TCylinder::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TCylinder*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TCylinder::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TCylinder*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TPlane::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TPlane::Class_Name()
{
   return "TPlane";
}

//______________________________________________________________________________
const char *TPlane::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TPlane*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TPlane::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TPlane*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TPlane::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TPlane*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TPlane::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TPlane*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr THype::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *THype::Class_Name()
{
   return "THype";
}

//______________________________________________________________________________
const char *THype::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::THype*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int THype::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::THype*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *THype::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::THype*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *THype::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::THype*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TStraightTrack::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TStraightTrack::Class_Name()
{
   return "TStraightTrack";
}

//______________________________________________________________________________
const char *TStraightTrack::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TStraightTrack*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TStraightTrack::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TStraightTrack*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TStraightTrack::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TStraightTrack*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TStraightTrack::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TStraightTrack*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TVSolid::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TVSolid::Class_Name()
{
   return "TVSolid";
}

//______________________________________________________________________________
const char *TVSolid::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TVSolid*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TVSolid::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TVSolid*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TVSolid::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TVSolid*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TVSolid::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TVSolid*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TTube::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TTube::Class_Name()
{
   return "TTube";
}

//______________________________________________________________________________
const char *TTube::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TTube*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TTube::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TTube*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TTube::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TTube*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TTube::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TTube*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TBField::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TBField::Class_Name()
{
   return "TBField";
}

//______________________________________________________________________________
const char *TBField::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TBField*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TBField::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TBField*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TBField::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TBField*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TBField::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TBField*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TVCurve::Streamer(TBuffer &R__b)
{
   // Stream an object of class TVCurve.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TVCurve::Class(),this);
   } else {
      R__b.WriteClassBuffer(TVCurve::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TVCurve(void *p) {
      return  p ? new(p) ::TVCurve : new ::TVCurve;
   }
   static void *newArray_TVCurve(Long_t nElements, void *p) {
      return p ? new(p) ::TVCurve[nElements] : new ::TVCurve[nElements];
   }
   // Wrapper around operator delete
   static void delete_TVCurve(void *p) {
      delete ((::TVCurve*)p);
   }
   static void deleteArray_TVCurve(void *p) {
      delete [] ((::TVCurve*)p);
   }
   static void destruct_TVCurve(void *p) {
      typedef ::TVCurve current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TVCurve

//______________________________________________________________________________
void TVTrack::Streamer(TBuffer &R__b)
{
   // Stream an object of class TVTrack.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TVTrack::Class(),this);
   } else {
      R__b.WriteClassBuffer(TVTrack::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_TVTrack(void *p) {
      delete ((::TVTrack*)p);
   }
   static void deleteArray_TVTrack(void *p) {
      delete [] ((::TVTrack*)p);
   }
   static void destruct_TVTrack(void *p) {
      typedef ::TVTrack current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TVTrack

//______________________________________________________________________________
void THelicalTrack::Streamer(TBuffer &R__b)
{
   // Stream an object of class THelicalTrack.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(THelicalTrack::Class(),this);
   } else {
      R__b.WriteClassBuffer(THelicalTrack::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_THelicalTrack(void *p) {
      return  p ? new(p) ::THelicalTrack : new ::THelicalTrack;
   }
   static void *newArray_THelicalTrack(Long_t nElements, void *p) {
      return p ? new(p) ::THelicalTrack[nElements] : new ::THelicalTrack[nElements];
   }
   // Wrapper around operator delete
   static void delete_THelicalTrack(void *p) {
      delete ((::THelicalTrack*)p);
   }
   static void deleteArray_THelicalTrack(void *p) {
      delete [] ((::THelicalTrack*)p);
   }
   static void destruct_THelicalTrack(void *p) {
      typedef ::THelicalTrack current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::THelicalTrack

//______________________________________________________________________________
void TCircle::Streamer(TBuffer &R__b)
{
   // Stream an object of class TCircle.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TCircle::Class(),this);
   } else {
      R__b.WriteClassBuffer(TCircle::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TCircle(void *p) {
      return  p ? new(p) ::TCircle : new ::TCircle;
   }
   static void *newArray_TCircle(Long_t nElements, void *p) {
      return p ? new(p) ::TCircle[nElements] : new ::TCircle[nElements];
   }
   // Wrapper around operator delete
   static void delete_TCircle(void *p) {
      delete ((::TCircle*)p);
   }
   static void deleteArray_TCircle(void *p) {
      delete [] ((::TCircle*)p);
   }
   static void destruct_TCircle(void *p) {
      typedef ::TCircle current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TCircle

//______________________________________________________________________________
void TVSurface::Streamer(TBuffer &R__b)
{
   // Stream an object of class TVSurface.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TVSurface::Class(),this);
   } else {
      R__b.WriteClassBuffer(TVSurface::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_TVSurface(void *p) {
      delete ((::TVSurface*)p);
   }
   static void deleteArray_TVSurface(void *p) {
      delete [] ((::TVSurface*)p);
   }
   static void destruct_TVSurface(void *p) {
      typedef ::TVSurface current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TVSurface

//______________________________________________________________________________
void TCutCone::Streamer(TBuffer &R__b)
{
   // Stream an object of class TCutCone.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TCutCone::Class(),this);
   } else {
      R__b.WriteClassBuffer(TCutCone::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TCutCone(void *p) {
      return  p ? new(p) ::TCutCone : new ::TCutCone;
   }
   static void *newArray_TCutCone(Long_t nElements, void *p) {
      return p ? new(p) ::TCutCone[nElements] : new ::TCutCone[nElements];
   }
   // Wrapper around operator delete
   static void delete_TCutCone(void *p) {
      delete ((::TCutCone*)p);
   }
   static void deleteArray_TCutCone(void *p) {
      delete [] ((::TCutCone*)p);
   }
   static void destruct_TCutCone(void *p) {
      typedef ::TCutCone current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TCutCone

//______________________________________________________________________________
void TCylinder::Streamer(TBuffer &R__b)
{
   // Stream an object of class TCylinder.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TCylinder::Class(),this);
   } else {
      R__b.WriteClassBuffer(TCylinder::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TCylinder(void *p) {
      return  p ? new(p) ::TCylinder : new ::TCylinder;
   }
   static void *newArray_TCylinder(Long_t nElements, void *p) {
      return p ? new(p) ::TCylinder[nElements] : new ::TCylinder[nElements];
   }
   // Wrapper around operator delete
   static void delete_TCylinder(void *p) {
      delete ((::TCylinder*)p);
   }
   static void deleteArray_TCylinder(void *p) {
      delete [] ((::TCylinder*)p);
   }
   static void destruct_TCylinder(void *p) {
      typedef ::TCylinder current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TCylinder

//______________________________________________________________________________
void TPlane::Streamer(TBuffer &R__b)
{
   // Stream an object of class TPlane.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TPlane::Class(),this);
   } else {
      R__b.WriteClassBuffer(TPlane::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TPlane(void *p) {
      return  p ? new(p) ::TPlane : new ::TPlane;
   }
   static void *newArray_TPlane(Long_t nElements, void *p) {
      return p ? new(p) ::TPlane[nElements] : new ::TPlane[nElements];
   }
   // Wrapper around operator delete
   static void delete_TPlane(void *p) {
      delete ((::TPlane*)p);
   }
   static void deleteArray_TPlane(void *p) {
      delete [] ((::TPlane*)p);
   }
   static void destruct_TPlane(void *p) {
      typedef ::TPlane current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TPlane

//______________________________________________________________________________
void THype::Streamer(TBuffer &R__b)
{
   // Stream an object of class THype.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(THype::Class(),this);
   } else {
      R__b.WriteClassBuffer(THype::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_THype(void *p) {
      return  p ? new(p) ::THype : new ::THype;
   }
   static void *newArray_THype(Long_t nElements, void *p) {
      return p ? new(p) ::THype[nElements] : new ::THype[nElements];
   }
   // Wrapper around operator delete
   static void delete_THype(void *p) {
      delete ((::THype*)p);
   }
   static void deleteArray_THype(void *p) {
      delete [] ((::THype*)p);
   }
   static void destruct_THype(void *p) {
      typedef ::THype current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::THype

//______________________________________________________________________________
void TStraightTrack::Streamer(TBuffer &R__b)
{
   // Stream an object of class TStraightTrack.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TStraightTrack::Class(),this);
   } else {
      R__b.WriteClassBuffer(TStraightTrack::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TStraightTrack(void *p) {
      return  p ? new(p) ::TStraightTrack : new ::TStraightTrack;
   }
   static void *newArray_TStraightTrack(Long_t nElements, void *p) {
      return p ? new(p) ::TStraightTrack[nElements] : new ::TStraightTrack[nElements];
   }
   // Wrapper around operator delete
   static void delete_TStraightTrack(void *p) {
      delete ((::TStraightTrack*)p);
   }
   static void deleteArray_TStraightTrack(void *p) {
      delete [] ((::TStraightTrack*)p);
   }
   static void destruct_TStraightTrack(void *p) {
      typedef ::TStraightTrack current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TStraightTrack

//______________________________________________________________________________
void TVSolid::Streamer(TBuffer &R__b)
{
   // Stream an object of class TVSolid.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TVSolid::Class(),this);
   } else {
      R__b.WriteClassBuffer(TVSolid::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TVSolid(void *p) {
      return  p ? new(p) ::TVSolid : new ::TVSolid;
   }
   static void *newArray_TVSolid(Long_t nElements, void *p) {
      return p ? new(p) ::TVSolid[nElements] : new ::TVSolid[nElements];
   }
   // Wrapper around operator delete
   static void delete_TVSolid(void *p) {
      delete ((::TVSolid*)p);
   }
   static void deleteArray_TVSolid(void *p) {
      delete [] ((::TVSolid*)p);
   }
   static void destruct_TVSolid(void *p) {
      typedef ::TVSolid current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TVSolid

//______________________________________________________________________________
void TTube::Streamer(TBuffer &R__b)
{
   // Stream an object of class TTube.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TTube::Class(),this);
   } else {
      R__b.WriteClassBuffer(TTube::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TTube(void *p) {
      return  p ? new(p) ::TTube : new ::TTube;
   }
   static void *newArray_TTube(Long_t nElements, void *p) {
      return p ? new(p) ::TTube[nElements] : new ::TTube[nElements];
   }
   // Wrapper around operator delete
   static void delete_TTube(void *p) {
      delete ((::TTube*)p);
   }
   static void deleteArray_TTube(void *p) {
      delete [] ((::TTube*)p);
   }
   static void destruct_TTube(void *p) {
      typedef ::TTube current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TTube

//______________________________________________________________________________
void TBField::Streamer(TBuffer &R__b)
{
   // Stream an object of class TBField.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TBField::Class(),this);
   } else {
      R__b.WriteClassBuffer(TBField::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TBField(void *p) {
      return  p ? new(p) ::TBField : new ::TBField;
   }
   static void *newArray_TBField(Long_t nElements, void *p) {
      return p ? new(p) ::TBField[nElements] : new ::TBField[nElements];
   }
   // Wrapper around operator delete
   static void delete_TBField(void *p) {
      delete ((::TBField*)p);
   }
   static void deleteArray_TBField(void *p) {
      delete [] ((::TBField*)p);
   }
   static void destruct_TBField(void *p) {
      typedef ::TBField current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TBField

namespace {
  void TriggerDictionaryInitialization_S4GeomDict_Impl() {
    static const char* headers[] = {
"THelicalTrack.h",
"TCircle.h",
"TCutCone.h",
"TVCurve.h",
"TCylinder.h",
"TPlane.h",
"THype.h",
"TStraightTrack.h",
"TTube.h",
"TVSolid.h",
"TVSurface.h",
"TVTrack.h",
"TBField.h",
0
    };
    static const char* includePaths[] = {
"../../include",
"/opt/root/include",
"/home/smechbal/ANALYSIS_SOFT/src/geomlib/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "S4GeomDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(Base class for any curve)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TVCurve.h")))  __attribute__((annotate("$clingAutoload$THelicalTrack.h")))  TVCurve;
class __attribute__((annotate(R"ATTRDUMP(Base class for any track)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TVTrack.h")))  __attribute__((annotate("$clingAutoload$THelicalTrack.h")))  TVTrack;
class __attribute__((annotate(R"ATTRDUMP(circle class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$THelicalTrack.h")))  THelicalTrack;
class __attribute__((annotate(R"ATTRDUMP(circle class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TCircle.h")))  TCircle;
class __attribute__((annotate(R"ATTRDUMP(Base class for any surface)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TVSurface.h")))  __attribute__((annotate("$clingAutoload$TCutCone.h")))  TVSurface;
class __attribute__((annotate(R"ATTRDUMP(hype class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TCutCone.h")))  TCutCone;
class __attribute__((annotate(R"ATTRDUMP(cylinder class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TCylinder.h")))  TCylinder;
class __attribute__((annotate(R"ATTRDUMP(plane class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TPlane.h")))  TPlane;
class __attribute__((annotate(R"ATTRDUMP(hype class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$THype.h")))  THype;
class __attribute__((annotate(R"ATTRDUMP(circle class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TStraightTrack.h")))  TStraightTrack;
class __attribute__((annotate(R"ATTRDUMP(Base class for any solid)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TVSolid.h")))  __attribute__((annotate("$clingAutoload$TTube.h")))  TVSolid;
class __attribute__((annotate(R"ATTRDUMP(TTube class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TTube.h")))  TTube;
class __attribute__((annotate(R"ATTRDUMP(Base class for detector system)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TBField.h")))  TBField;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "S4GeomDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "THelicalTrack.h"
#include "TCircle.h"
#include "TCutCone.h"
#include "TVCurve.h"
#include "TCylinder.h"
#include "TPlane.h"
#include "THype.h"
#include "TStraightTrack.h"
#include "TTube.h"
#include "TVSolid.h"
#include "TVSurface.h"
#include "TVTrack.h"
#include "TBField.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TBField", payloadCode, "@",
"TCircle", payloadCode, "@",
"TCutCone", payloadCode, "@",
"TCylinder", payloadCode, "@",
"THelicalTrack", payloadCode, "@",
"THype", payloadCode, "@",
"TPlane", payloadCode, "@",
"TStraightTrack", payloadCode, "@",
"TTube", payloadCode, "@",
"TVCurve", payloadCode, "@",
"TVSolid", payloadCode, "@",
"TVSurface", payloadCode, "@",
"TVTrack", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("S4GeomDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_S4GeomDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_S4GeomDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_S4GeomDict() {
  TriggerDictionaryInitialization_S4GeomDict_Impl();
}
