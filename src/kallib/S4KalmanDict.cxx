// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME S4KalmanDict

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
#include "TKalMatrix.h"
#include "TVKalSite.h"
#include "TVKalState.h"
#include "TVKalSystem.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TKalMatrix(void *p = 0);
   static void *newArray_TKalMatrix(Long_t size, void *p);
   static void delete_TKalMatrix(void *p);
   static void deleteArray_TKalMatrix(void *p);
   static void destruct_TKalMatrix(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TKalMatrix*)
   {
      ::TKalMatrix *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TKalMatrix >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TKalMatrix", ::TKalMatrix::Class_Version(), "TKalMatrix.h", 28,
                  typeid(::TKalMatrix), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TKalMatrix::Dictionary, isa_proxy, 4,
                  sizeof(::TKalMatrix) );
      instance.SetNew(&new_TKalMatrix);
      instance.SetNewArray(&newArray_TKalMatrix);
      instance.SetDelete(&delete_TKalMatrix);
      instance.SetDeleteArray(&deleteArray_TKalMatrix);
      instance.SetDestructor(&destruct_TKalMatrix);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TKalMatrix*)
   {
      return GenerateInitInstanceLocal((::TKalMatrix*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TKalMatrix*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_TVKalState(void *p);
   static void deleteArray_TVKalState(void *p);
   static void destruct_TVKalState(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TVKalState*)
   {
      ::TVKalState *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TVKalState >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TVKalState", ::TVKalState::Class_Version(), "TVKalState.h", 26,
                  typeid(::TVKalState), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TVKalState::Dictionary, isa_proxy, 4,
                  sizeof(::TVKalState) );
      instance.SetDelete(&delete_TVKalState);
      instance.SetDeleteArray(&deleteArray_TVKalState);
      instance.SetDestructor(&destruct_TVKalState);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TVKalState*)
   {
      return GenerateInitInstanceLocal((::TVKalState*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TVKalState*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_TVKalSite(void *p);
   static void deleteArray_TVKalSite(void *p);
   static void destruct_TVKalSite(void *p);
   static Long64_t merge_TVKalSite(void *obj, TCollection *coll,TFileMergeInfo *info);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TVKalSite*)
   {
      ::TVKalSite *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TVKalSite >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TVKalSite", ::TVKalSite::Class_Version(), "TVKalSite.h", 37,
                  typeid(::TVKalSite), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TVKalSite::Dictionary, isa_proxy, 4,
                  sizeof(::TVKalSite) );
      instance.SetDelete(&delete_TVKalSite);
      instance.SetDeleteArray(&deleteArray_TVKalSite);
      instance.SetDestructor(&destruct_TVKalSite);
      instance.SetMerge(&merge_TVKalSite);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TVKalSite*)
   {
      return GenerateInitInstanceLocal((::TVKalSite*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TVKalSite*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TVKalSystem(void *p = 0);
   static void *newArray_TVKalSystem(Long_t size, void *p);
   static void delete_TVKalSystem(void *p);
   static void deleteArray_TVKalSystem(void *p);
   static void destruct_TVKalSystem(void *p);
   static Long64_t merge_TVKalSystem(void *obj, TCollection *coll,TFileMergeInfo *info);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TVKalSystem*)
   {
      ::TVKalSystem *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TVKalSystem >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TVKalSystem", ::TVKalSystem::Class_Version(), "TVKalSystem.h", 31,
                  typeid(::TVKalSystem), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TVKalSystem::Dictionary, isa_proxy, 4,
                  sizeof(::TVKalSystem) );
      instance.SetNew(&new_TVKalSystem);
      instance.SetNewArray(&newArray_TVKalSystem);
      instance.SetDelete(&delete_TVKalSystem);
      instance.SetDeleteArray(&deleteArray_TVKalSystem);
      instance.SetDestructor(&destruct_TVKalSystem);
      instance.SetMerge(&merge_TVKalSystem);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TVKalSystem*)
   {
      return GenerateInitInstanceLocal((::TVKalSystem*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TVKalSystem*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TKalMatrix::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TKalMatrix::Class_Name()
{
   return "TKalMatrix";
}

//______________________________________________________________________________
const char *TKalMatrix::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TKalMatrix*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TKalMatrix::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TKalMatrix*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TKalMatrix::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TKalMatrix*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TKalMatrix::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TKalMatrix*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TVKalState::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TVKalState::Class_Name()
{
   return "TVKalState";
}

//______________________________________________________________________________
const char *TVKalState::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TVKalState*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TVKalState::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TVKalState*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TVKalState::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TVKalState*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TVKalState::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TVKalState*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TVKalSite::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TVKalSite::Class_Name()
{
   return "TVKalSite";
}

//______________________________________________________________________________
const char *TVKalSite::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TVKalSite*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TVKalSite::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TVKalSite*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TVKalSite::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TVKalSite*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TVKalSite::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TVKalSite*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TVKalSystem::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TVKalSystem::Class_Name()
{
   return "TVKalSystem";
}

//______________________________________________________________________________
const char *TVKalSystem::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TVKalSystem*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TVKalSystem::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TVKalSystem*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TVKalSystem::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TVKalSystem*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TVKalSystem::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TVKalSystem*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TKalMatrix::Streamer(TBuffer &R__b)
{
   // Stream an object of class TKalMatrix.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TKalMatrix::Class(),this);
   } else {
      R__b.WriteClassBuffer(TKalMatrix::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TKalMatrix(void *p) {
      return  p ? new(p) ::TKalMatrix : new ::TKalMatrix;
   }
   static void *newArray_TKalMatrix(Long_t nElements, void *p) {
      return p ? new(p) ::TKalMatrix[nElements] : new ::TKalMatrix[nElements];
   }
   // Wrapper around operator delete
   static void delete_TKalMatrix(void *p) {
      delete ((::TKalMatrix*)p);
   }
   static void deleteArray_TKalMatrix(void *p) {
      delete [] ((::TKalMatrix*)p);
   }
   static void destruct_TKalMatrix(void *p) {
      typedef ::TKalMatrix current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TKalMatrix

//______________________________________________________________________________
void TVKalState::Streamer(TBuffer &R__b)
{
   // Stream an object of class TVKalState.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TVKalState::Class(),this);
   } else {
      R__b.WriteClassBuffer(TVKalState::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_TVKalState(void *p) {
      delete ((::TVKalState*)p);
   }
   static void deleteArray_TVKalState(void *p) {
      delete [] ((::TVKalState*)p);
   }
   static void destruct_TVKalState(void *p) {
      typedef ::TVKalState current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TVKalState

//______________________________________________________________________________
void TVKalSite::Streamer(TBuffer &R__b)
{
   // Stream an object of class TVKalSite.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TVKalSite::Class(),this);
   } else {
      R__b.WriteClassBuffer(TVKalSite::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_TVKalSite(void *p) {
      delete ((::TVKalSite*)p);
   }
   static void deleteArray_TVKalSite(void *p) {
      delete [] ((::TVKalSite*)p);
   }
   static void destruct_TVKalSite(void *p) {
      typedef ::TVKalSite current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around the merge function.
   static Long64_t  merge_TVKalSite(void *obj,TCollection *coll,TFileMergeInfo *) {
      return ((::TVKalSite*)obj)->Merge(coll);
   }
} // end of namespace ROOT for class ::TVKalSite

//______________________________________________________________________________
void TVKalSystem::Streamer(TBuffer &R__b)
{
   // Stream an object of class TVKalSystem.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TVKalSystem::Class(),this);
   } else {
      R__b.WriteClassBuffer(TVKalSystem::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TVKalSystem(void *p) {
      return  p ? new(p) ::TVKalSystem : new ::TVKalSystem;
   }
   static void *newArray_TVKalSystem(Long_t nElements, void *p) {
      return p ? new(p) ::TVKalSystem[nElements] : new ::TVKalSystem[nElements];
   }
   // Wrapper around operator delete
   static void delete_TVKalSystem(void *p) {
      delete ((::TVKalSystem*)p);
   }
   static void deleteArray_TVKalSystem(void *p) {
      delete [] ((::TVKalSystem*)p);
   }
   static void destruct_TVKalSystem(void *p) {
      typedef ::TVKalSystem current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around the merge function.
   static Long64_t  merge_TVKalSystem(void *obj,TCollection *coll,TFileMergeInfo *) {
      return ((::TVKalSystem*)obj)->Merge(coll);
   }
} // end of namespace ROOT for class ::TVKalSystem

namespace {
  void TriggerDictionaryInitialization_S4KalmanDict_Impl() {
    static const char* headers[] = {
"TKalMatrix.h",
"TVKalSite.h",
"TVKalState.h",
"TVKalSystem.h",
0
    };
    static const char* includePaths[] = {
"../../include",
"/opt/root/include",
"/home/smechbal/ANALYSIS_SOFT/src/kallib/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "S4KalmanDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(Base class for Kalman matrix)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TKalMatrix.h")))  TKalMatrix;
class __attribute__((annotate(R"ATTRDUMP(Base class for state vector objects)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TVKalState.h")))  __attribute__((annotate("$clingAutoload$TVKalSite.h")))  TVKalState;
class __attribute__((annotate(R"ATTRDUMP(Base class for measurement vector objects)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TVKalSite.h")))  TVKalSite;
class __attribute__((annotate(R"ATTRDUMP(Base class for Kalman Filter)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TVKalSystem.h")))  TVKalSystem;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "S4KalmanDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TKalMatrix.h"
#include "TVKalSite.h"
#include "TVKalState.h"
#include "TVKalSystem.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TKalMatrix", payloadCode, "@",
"TVKalSite", payloadCode, "@",
"TVKalState", payloadCode, "@",
"TVKalSystem", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("S4KalmanDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_S4KalmanDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_S4KalmanDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_S4KalmanDict() {
  TriggerDictionaryInitialization_S4KalmanDict_Impl();
}
