// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME S4KalTrackDict

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
#include "TVTrackHit.h"
#include "TVMeasLayer.h"
#include "TKalTrackSite.h"
#include "TKalTrackState.h"
#include "TKalTrack.h"
#include "TTrackFrame.h"
#include "TKalDetCradle.h"
#include "TVKalDetector.h"
#include "TKalFilterCond.h"
#include "KalTrackDim.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void delete_TVMeasLayer(void *p);
   static void deleteArray_TVMeasLayer(void *p);
   static void destruct_TVMeasLayer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TVMeasLayer*)
   {
      ::TVMeasLayer *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TVMeasLayer >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TVMeasLayer", ::TVMeasLayer::Class_Version(), "TVMeasLayer.h", 40,
                  typeid(::TVMeasLayer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TVMeasLayer::Dictionary, isa_proxy, 4,
                  sizeof(::TVMeasLayer) );
      instance.SetDelete(&delete_TVMeasLayer);
      instance.SetDeleteArray(&deleteArray_TVMeasLayer);
      instance.SetDestructor(&destruct_TVMeasLayer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TVMeasLayer*)
   {
      return GenerateInitInstanceLocal((::TVMeasLayer*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TVMeasLayer*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_TVTrackHit(void *p);
   static void deleteArray_TVTrackHit(void *p);
   static void destruct_TVTrackHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TVTrackHit*)
   {
      ::TVTrackHit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TVTrackHit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TVTrackHit", ::TVTrackHit::Class_Version(), "TVTrackHit.h", 25,
                  typeid(::TVTrackHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TVTrackHit::Dictionary, isa_proxy, 4,
                  sizeof(::TVTrackHit) );
      instance.SetDelete(&delete_TVTrackHit);
      instance.SetDeleteArray(&deleteArray_TVTrackHit);
      instance.SetDestructor(&destruct_TVTrackHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TVTrackHit*)
   {
      return GenerateInitInstanceLocal((::TVTrackHit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TVTrackHit*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TTrackFrame(void *p = 0);
   static void *newArray_TTrackFrame(Long_t size, void *p);
   static void delete_TTrackFrame(void *p);
   static void deleteArray_TTrackFrame(void *p);
   static void destruct_TTrackFrame(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TTrackFrame*)
   {
      ::TTrackFrame *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TTrackFrame >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TTrackFrame", ::TTrackFrame::Class_Version(), "TTrackFrame.h", 24,
                  typeid(::TTrackFrame), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TTrackFrame::Dictionary, isa_proxy, 4,
                  sizeof(::TTrackFrame) );
      instance.SetNew(&new_TTrackFrame);
      instance.SetNewArray(&newArray_TTrackFrame);
      instance.SetDelete(&delete_TTrackFrame);
      instance.SetDeleteArray(&deleteArray_TTrackFrame);
      instance.SetDestructor(&destruct_TTrackFrame);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TTrackFrame*)
   {
      return GenerateInitInstanceLocal((::TTrackFrame*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TTrackFrame*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TKalTrackSite(void *p = 0);
   static void *newArray_TKalTrackSite(Long_t size, void *p);
   static void delete_TKalTrackSite(void *p);
   static void deleteArray_TKalTrackSite(void *p);
   static void destruct_TKalTrackSite(void *p);
   static Long64_t merge_TKalTrackSite(void *obj, TCollection *coll,TFileMergeInfo *info);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TKalTrackSite*)
   {
      ::TKalTrackSite *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TKalTrackSite >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TKalTrackSite", ::TKalTrackSite::Class_Version(), "TKalTrackSite.h", 38,
                  typeid(::TKalTrackSite), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TKalTrackSite::Dictionary, isa_proxy, 4,
                  sizeof(::TKalTrackSite) );
      instance.SetNew(&new_TKalTrackSite);
      instance.SetNewArray(&newArray_TKalTrackSite);
      instance.SetDelete(&delete_TKalTrackSite);
      instance.SetDeleteArray(&deleteArray_TKalTrackSite);
      instance.SetDestructor(&destruct_TKalTrackSite);
      instance.SetMerge(&merge_TKalTrackSite);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TKalTrackSite*)
   {
      return GenerateInitInstanceLocal((::TKalTrackSite*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TKalTrackSite*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TKalTrackState(void *p = 0);
   static void *newArray_TKalTrackState(Long_t size, void *p);
   static void delete_TKalTrackState(void *p);
   static void deleteArray_TKalTrackState(void *p);
   static void destruct_TKalTrackState(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TKalTrackState*)
   {
      ::TKalTrackState *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TKalTrackState >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TKalTrackState", ::TKalTrackState::Class_Version(), "TKalTrackState.h", 38,
                  typeid(::TKalTrackState), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TKalTrackState::Dictionary, isa_proxy, 4,
                  sizeof(::TKalTrackState) );
      instance.SetNew(&new_TKalTrackState);
      instance.SetNewArray(&newArray_TKalTrackState);
      instance.SetDelete(&delete_TKalTrackState);
      instance.SetDeleteArray(&deleteArray_TKalTrackState);
      instance.SetDestructor(&destruct_TKalTrackState);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TKalTrackState*)
   {
      return GenerateInitInstanceLocal((::TKalTrackState*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TKalTrackState*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TKalTrack(void *p = 0);
   static void *newArray_TKalTrack(Long_t size, void *p);
   static void delete_TKalTrack(void *p);
   static void deleteArray_TKalTrack(void *p);
   static void destruct_TKalTrack(void *p);
   static Long64_t merge_TKalTrack(void *obj, TCollection *coll,TFileMergeInfo *info);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TKalTrack*)
   {
      ::TKalTrack *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TKalTrack >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TKalTrack", ::TKalTrack::Class_Version(), "TKalTrack.h", 32,
                  typeid(::TKalTrack), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TKalTrack::Dictionary, isa_proxy, 4,
                  sizeof(::TKalTrack) );
      instance.SetNew(&new_TKalTrack);
      instance.SetNewArray(&newArray_TKalTrack);
      instance.SetDelete(&delete_TKalTrack);
      instance.SetDeleteArray(&deleteArray_TKalTrack);
      instance.SetDestructor(&destruct_TKalTrack);
      instance.SetMerge(&merge_TKalTrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TKalTrack*)
   {
      return GenerateInitInstanceLocal((::TKalTrack*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TKalTrack*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TKalDetCradle(void *p = 0);
   static void *newArray_TKalDetCradle(Long_t size, void *p);
   static void delete_TKalDetCradle(void *p);
   static void deleteArray_TKalDetCradle(void *p);
   static void destruct_TKalDetCradle(void *p);
   static Long64_t merge_TKalDetCradle(void *obj, TCollection *coll,TFileMergeInfo *info);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TKalDetCradle*)
   {
      ::TKalDetCradle *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TKalDetCradle >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TKalDetCradle", ::TKalDetCradle::Class_Version(), "TKalDetCradle.h", 43,
                  typeid(::TKalDetCradle), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TKalDetCradle::Dictionary, isa_proxy, 4,
                  sizeof(::TKalDetCradle) );
      instance.SetNew(&new_TKalDetCradle);
      instance.SetNewArray(&newArray_TKalDetCradle);
      instance.SetDelete(&delete_TKalDetCradle);
      instance.SetDeleteArray(&deleteArray_TKalDetCradle);
      instance.SetDestructor(&destruct_TKalDetCradle);
      instance.SetMerge(&merge_TKalDetCradle);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TKalDetCradle*)
   {
      return GenerateInitInstanceLocal((::TKalDetCradle*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TKalDetCradle*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TVKalDetector(void *p = 0);
   static void *newArray_TVKalDetector(Long_t size, void *p);
   static void delete_TVKalDetector(void *p);
   static void deleteArray_TVKalDetector(void *p);
   static void destruct_TVKalDetector(void *p);
   static Long64_t merge_TVKalDetector(void *obj, TCollection *coll,TFileMergeInfo *info);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TVKalDetector*)
   {
      ::TVKalDetector *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TVKalDetector >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TVKalDetector", ::TVKalDetector::Class_Version(), "TVKalDetector.h", 30,
                  typeid(::TVKalDetector), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TVKalDetector::Dictionary, isa_proxy, 4,
                  sizeof(::TVKalDetector) );
      instance.SetNew(&new_TVKalDetector);
      instance.SetNewArray(&newArray_TVKalDetector);
      instance.SetDelete(&delete_TVKalDetector);
      instance.SetDeleteArray(&deleteArray_TVKalDetector);
      instance.SetDestructor(&destruct_TVKalDetector);
      instance.SetMerge(&merge_TVKalDetector);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TVKalDetector*)
   {
      return GenerateInitInstanceLocal((::TVKalDetector*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TVKalDetector*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TKalFilterCond(void *p = 0);
   static void *newArray_TKalFilterCond(Long_t size, void *p);
   static void delete_TKalFilterCond(void *p);
   static void deleteArray_TKalFilterCond(void *p);
   static void destruct_TKalFilterCond(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TKalFilterCond*)
   {
      ::TKalFilterCond *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TKalFilterCond >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TKalFilterCond", ::TKalFilterCond::Class_Version(), "TKalFilterCond.h", 25,
                  typeid(::TKalFilterCond), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TKalFilterCond::Dictionary, isa_proxy, 4,
                  sizeof(::TKalFilterCond) );
      instance.SetNew(&new_TKalFilterCond);
      instance.SetNewArray(&newArray_TKalFilterCond);
      instance.SetDelete(&delete_TKalFilterCond);
      instance.SetDeleteArray(&deleteArray_TKalFilterCond);
      instance.SetDestructor(&destruct_TKalFilterCond);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TKalFilterCond*)
   {
      return GenerateInitInstanceLocal((::TKalFilterCond*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TKalFilterCond*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TVMeasLayer::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TVMeasLayer::Class_Name()
{
   return "TVMeasLayer";
}

//______________________________________________________________________________
const char *TVMeasLayer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TVMeasLayer*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TVMeasLayer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TVMeasLayer*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TVMeasLayer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TVMeasLayer*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TVMeasLayer::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TVMeasLayer*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TVTrackHit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TVTrackHit::Class_Name()
{
   return "TVTrackHit";
}

//______________________________________________________________________________
const char *TVTrackHit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TVTrackHit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TVTrackHit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TVTrackHit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TVTrackHit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TVTrackHit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TVTrackHit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TVTrackHit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TTrackFrame::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TTrackFrame::Class_Name()
{
   return "TTrackFrame";
}

//______________________________________________________________________________
const char *TTrackFrame::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TTrackFrame*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TTrackFrame::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TTrackFrame*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TTrackFrame::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TTrackFrame*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TTrackFrame::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TTrackFrame*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TKalTrackSite::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TKalTrackSite::Class_Name()
{
   return "TKalTrackSite";
}

//______________________________________________________________________________
const char *TKalTrackSite::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TKalTrackSite*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TKalTrackSite::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TKalTrackSite*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TKalTrackSite::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TKalTrackSite*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TKalTrackSite::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TKalTrackSite*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TKalTrackState::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TKalTrackState::Class_Name()
{
   return "TKalTrackState";
}

//______________________________________________________________________________
const char *TKalTrackState::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TKalTrackState*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TKalTrackState::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TKalTrackState*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TKalTrackState::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TKalTrackState*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TKalTrackState::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TKalTrackState*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TKalTrack::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TKalTrack::Class_Name()
{
   return "TKalTrack";
}

//______________________________________________________________________________
const char *TKalTrack::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TKalTrack*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TKalTrack::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TKalTrack*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TKalTrack::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TKalTrack*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TKalTrack::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TKalTrack*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TKalDetCradle::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TKalDetCradle::Class_Name()
{
   return "TKalDetCradle";
}

//______________________________________________________________________________
const char *TKalDetCradle::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TKalDetCradle*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TKalDetCradle::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TKalDetCradle*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TKalDetCradle::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TKalDetCradle*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TKalDetCradle::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TKalDetCradle*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TVKalDetector::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TVKalDetector::Class_Name()
{
   return "TVKalDetector";
}

//______________________________________________________________________________
const char *TVKalDetector::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TVKalDetector*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TVKalDetector::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TVKalDetector*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TVKalDetector::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TVKalDetector*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TVKalDetector::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TVKalDetector*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TKalFilterCond::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TKalFilterCond::Class_Name()
{
   return "TKalFilterCond";
}

//______________________________________________________________________________
const char *TKalFilterCond::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TKalFilterCond*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TKalFilterCond::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TKalFilterCond*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TKalFilterCond::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TKalFilterCond*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TKalFilterCond::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TKalFilterCond*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TVMeasLayer::Streamer(TBuffer &R__b)
{
   // Stream an object of class TVMeasLayer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TVMeasLayer::Class(),this);
   } else {
      R__b.WriteClassBuffer(TVMeasLayer::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_TVMeasLayer(void *p) {
      delete ((::TVMeasLayer*)p);
   }
   static void deleteArray_TVMeasLayer(void *p) {
      delete [] ((::TVMeasLayer*)p);
   }
   static void destruct_TVMeasLayer(void *p) {
      typedef ::TVMeasLayer current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TVMeasLayer

//______________________________________________________________________________
void TVTrackHit::Streamer(TBuffer &R__b)
{
   // Stream an object of class TVTrackHit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TVTrackHit::Class(),this);
   } else {
      R__b.WriteClassBuffer(TVTrackHit::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_TVTrackHit(void *p) {
      delete ((::TVTrackHit*)p);
   }
   static void deleteArray_TVTrackHit(void *p) {
      delete [] ((::TVTrackHit*)p);
   }
   static void destruct_TVTrackHit(void *p) {
      typedef ::TVTrackHit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TVTrackHit

//______________________________________________________________________________
void TTrackFrame::Streamer(TBuffer &R__b)
{
   // Stream an object of class TTrackFrame.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TTrackFrame::Class(),this);
   } else {
      R__b.WriteClassBuffer(TTrackFrame::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TTrackFrame(void *p) {
      return  p ? new(p) ::TTrackFrame : new ::TTrackFrame;
   }
   static void *newArray_TTrackFrame(Long_t nElements, void *p) {
      return p ? new(p) ::TTrackFrame[nElements] : new ::TTrackFrame[nElements];
   }
   // Wrapper around operator delete
   static void delete_TTrackFrame(void *p) {
      delete ((::TTrackFrame*)p);
   }
   static void deleteArray_TTrackFrame(void *p) {
      delete [] ((::TTrackFrame*)p);
   }
   static void destruct_TTrackFrame(void *p) {
      typedef ::TTrackFrame current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TTrackFrame

//______________________________________________________________________________
void TKalTrackSite::Streamer(TBuffer &R__b)
{
   // Stream an object of class TKalTrackSite.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TKalTrackSite::Class(),this);
   } else {
      R__b.WriteClassBuffer(TKalTrackSite::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TKalTrackSite(void *p) {
      return  p ? new(p) ::TKalTrackSite : new ::TKalTrackSite;
   }
   static void *newArray_TKalTrackSite(Long_t nElements, void *p) {
      return p ? new(p) ::TKalTrackSite[nElements] : new ::TKalTrackSite[nElements];
   }
   // Wrapper around operator delete
   static void delete_TKalTrackSite(void *p) {
      delete ((::TKalTrackSite*)p);
   }
   static void deleteArray_TKalTrackSite(void *p) {
      delete [] ((::TKalTrackSite*)p);
   }
   static void destruct_TKalTrackSite(void *p) {
      typedef ::TKalTrackSite current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around the merge function.
   static Long64_t  merge_TKalTrackSite(void *obj,TCollection *coll,TFileMergeInfo *) {
      return ((::TKalTrackSite*)obj)->Merge(coll);
   }
} // end of namespace ROOT for class ::TKalTrackSite

//______________________________________________________________________________
void TKalTrackState::Streamer(TBuffer &R__b)
{
   // Stream an object of class TKalTrackState.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TKalTrackState::Class(),this);
   } else {
      R__b.WriteClassBuffer(TKalTrackState::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TKalTrackState(void *p) {
      return  p ? new(p) ::TKalTrackState : new ::TKalTrackState;
   }
   static void *newArray_TKalTrackState(Long_t nElements, void *p) {
      return p ? new(p) ::TKalTrackState[nElements] : new ::TKalTrackState[nElements];
   }
   // Wrapper around operator delete
   static void delete_TKalTrackState(void *p) {
      delete ((::TKalTrackState*)p);
   }
   static void deleteArray_TKalTrackState(void *p) {
      delete [] ((::TKalTrackState*)p);
   }
   static void destruct_TKalTrackState(void *p) {
      typedef ::TKalTrackState current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TKalTrackState

//______________________________________________________________________________
void TKalTrack::Streamer(TBuffer &R__b)
{
   // Stream an object of class TKalTrack.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TKalTrack::Class(),this);
   } else {
      R__b.WriteClassBuffer(TKalTrack::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TKalTrack(void *p) {
      return  p ? new(p) ::TKalTrack : new ::TKalTrack;
   }
   static void *newArray_TKalTrack(Long_t nElements, void *p) {
      return p ? new(p) ::TKalTrack[nElements] : new ::TKalTrack[nElements];
   }
   // Wrapper around operator delete
   static void delete_TKalTrack(void *p) {
      delete ((::TKalTrack*)p);
   }
   static void deleteArray_TKalTrack(void *p) {
      delete [] ((::TKalTrack*)p);
   }
   static void destruct_TKalTrack(void *p) {
      typedef ::TKalTrack current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around the merge function.
   static Long64_t  merge_TKalTrack(void *obj,TCollection *coll,TFileMergeInfo *) {
      return ((::TKalTrack*)obj)->Merge(coll);
   }
} // end of namespace ROOT for class ::TKalTrack

//______________________________________________________________________________
void TKalDetCradle::Streamer(TBuffer &R__b)
{
   // Stream an object of class TKalDetCradle.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TKalDetCradle::Class(),this);
   } else {
      R__b.WriteClassBuffer(TKalDetCradle::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TKalDetCradle(void *p) {
      return  p ? new(p) ::TKalDetCradle : new ::TKalDetCradle;
   }
   static void *newArray_TKalDetCradle(Long_t nElements, void *p) {
      return p ? new(p) ::TKalDetCradle[nElements] : new ::TKalDetCradle[nElements];
   }
   // Wrapper around operator delete
   static void delete_TKalDetCradle(void *p) {
      delete ((::TKalDetCradle*)p);
   }
   static void deleteArray_TKalDetCradle(void *p) {
      delete [] ((::TKalDetCradle*)p);
   }
   static void destruct_TKalDetCradle(void *p) {
      typedef ::TKalDetCradle current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around the merge function.
   static Long64_t  merge_TKalDetCradle(void *obj,TCollection *coll,TFileMergeInfo *) {
      return ((::TKalDetCradle*)obj)->Merge(coll);
   }
} // end of namespace ROOT for class ::TKalDetCradle

//______________________________________________________________________________
void TVKalDetector::Streamer(TBuffer &R__b)
{
   // Stream an object of class TVKalDetector.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TVKalDetector::Class(),this);
   } else {
      R__b.WriteClassBuffer(TVKalDetector::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TVKalDetector(void *p) {
      return  p ? new(p) ::TVKalDetector : new ::TVKalDetector;
   }
   static void *newArray_TVKalDetector(Long_t nElements, void *p) {
      return p ? new(p) ::TVKalDetector[nElements] : new ::TVKalDetector[nElements];
   }
   // Wrapper around operator delete
   static void delete_TVKalDetector(void *p) {
      delete ((::TVKalDetector*)p);
   }
   static void deleteArray_TVKalDetector(void *p) {
      delete [] ((::TVKalDetector*)p);
   }
   static void destruct_TVKalDetector(void *p) {
      typedef ::TVKalDetector current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around the merge function.
   static Long64_t  merge_TVKalDetector(void *obj,TCollection *coll,TFileMergeInfo *) {
      return ((::TVKalDetector*)obj)->Merge(coll);
   }
} // end of namespace ROOT for class ::TVKalDetector

//______________________________________________________________________________
void TKalFilterCond::Streamer(TBuffer &R__b)
{
   // Stream an object of class TKalFilterCond.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TKalFilterCond::Class(),this);
   } else {
      R__b.WriteClassBuffer(TKalFilterCond::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TKalFilterCond(void *p) {
      return  p ? new(p) ::TKalFilterCond : new ::TKalFilterCond;
   }
   static void *newArray_TKalFilterCond(Long_t nElements, void *p) {
      return p ? new(p) ::TKalFilterCond[nElements] : new ::TKalFilterCond[nElements];
   }
   // Wrapper around operator delete
   static void delete_TKalFilterCond(void *p) {
      delete ((::TKalFilterCond*)p);
   }
   static void deleteArray_TKalFilterCond(void *p) {
      delete [] ((::TKalFilterCond*)p);
   }
   static void destruct_TKalFilterCond(void *p) {
      typedef ::TKalFilterCond current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TKalFilterCond

namespace {
  void TriggerDictionaryInitialization_S4KalTrackDict_Impl() {
    static const char* headers[] = {
"TVTrackHit.h",
"TVMeasLayer.h",
"TKalTrackSite.h",
"TKalTrackState.h",
"TKalTrack.h",
"TTrackFrame.h",
"TKalDetCradle.h",
"TVKalDetector.h",
"TKalFilterCond.h",
"KalTrackDim.h",
0
    };
    static const char* includePaths[] = {
"../../include",
"/opt/root/include",
"/home/smechbal/ANALYSIS_SOFT/src/kaltracklib/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "S4KalTrackDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(Measurement layer interface class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TVMeasLayer.h")))  __attribute__((annotate("$clingAutoload$TVTrackHit.h")))  TVMeasLayer;
class __attribute__((annotate(R"ATTRDUMP(Sample hit class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TVTrackHit.h")))  TVTrackHit;
class __attribute__((annotate("$clingAutoload$TTrackFrame.h")))  __attribute__((annotate("$clingAutoload$TKalTrackSite.h")))  TTrackFrame;
class __attribute__((annotate(R"ATTRDUMP(sample measurement site class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TKalTrackSite.h")))  TKalTrackSite;
class __attribute__((annotate(R"ATTRDUMP(sample state vector class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TKalTrackState.h")))  TKalTrackState;
class __attribute__((annotate(R"ATTRDUMP(Base class for Kalman Filter)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TKalTrack.h")))  TKalTrack;
class __attribute__((annotate(R"ATTRDUMP(Base class for detector system)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TKalDetCradle.h")))  TKalDetCradle;
class __attribute__((annotate(R"ATTRDUMP(Base class for detector system)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TVKalDetector.h")))  TVKalDetector;
class __attribute__((annotate(R"ATTRDUMP(Base class for detector system)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TKalFilterCond.h")))  TKalFilterCond;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "S4KalTrackDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TVTrackHit.h"
#include "TVMeasLayer.h"
#include "TKalTrackSite.h"
#include "TKalTrackState.h"
#include "TKalTrack.h"
#include "TTrackFrame.h"
#include "TKalDetCradle.h"
#include "TVKalDetector.h"
#include "TKalFilterCond.h"
#include "KalTrackDim.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TKalDetCradle", payloadCode, "@",
"TKalFilterCond", payloadCode, "@",
"TKalTrack", payloadCode, "@",
"TKalTrackSite", payloadCode, "@",
"TKalTrackState", payloadCode, "@",
"TTrackFrame", payloadCode, "@",
"TVKalDetector", payloadCode, "@",
"TVMeasLayer", payloadCode, "@",
"TVTrackHit", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("S4KalTrackDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_S4KalTrackDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_S4KalTrackDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_S4KalTrackDict() {
  TriggerDictionaryInitialization_S4KalTrackDict_Impl();
}
