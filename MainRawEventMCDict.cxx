// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME MainRawEventMCDict

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
#include "MainRawEventMC.h"
#include "/home/smechbal/ANALYSIS_SOFT/src/ALSim/MakeRawBPDEvent.h"
#include "/home/smechbal/ANALYSIS_SOFT/src/ALSim/MakeEventData.h"
#include "/home/smechbal/ANALYSIS_SOFT/src/ALSim/MakeRawEventMC.h"
#include "/home/smechbal/ANALYSIS_SOFT/src/ALKalman/ALKalman.h"
#include "/home/smechbal/ANALYSIS_SOFT/src/ALPatternRecognition/ALPatternRecognition.h"
#include "/home/smechbal/ANALYSIS_SOFT/src/RKFitter/RKfitter.h"
#include "/home/smechbal/ANALYSIS_SOFT/src/RKFitter/TkrData.h"
#include "/home/smechbal/ANALYSIS_SOFT/src/RKFitter/FieldMap.h"
#include "/home/smechbal/ANALYSIS_SOFT/src/RKFitter/RungeKutta4.h"

// Header files passed via #pragma extra_include

namespace {
  void TriggerDictionaryInitialization_MainRawEventMCDict_Impl() {
    static const char* headers[] = {
"MainRawEventMC.h",
"src/ALSim/MakeRawBPDEvent.h",
"src/ALSim/MakeEventData.h",
"src/ALSim/MakeRawEventMC.h",
"src/ALKalman/ALKalman.h",
"src/ALPatternRecognition/ALPatternRecognition.h",
"src/RKFitter/RKfitter.h",
"src/RKFitter/TkrData.h",
"src/RKFitter/FieldMap.h",
"src/RKFitter/RungeKutta4.h",
0
    };
    static const char* includePaths[] = {
"/home/smechbal/ANALYSIS_SOFT/include",
"/opt/root/include",
"/home/smechbal/ANALYSIS_SOFT/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "MainRawEventMCDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "MainRawEventMCDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "MainRawEventMC.h"
#include "src/ALSim/MakeRawBPDEvent.h"
#include "src/ALSim/MakeEventData.h"
#include "src/ALSim/MakeRawEventMC.h"
#include "src/ALKalman/ALKalman.h"
#include "src/ALPatternRecognition/ALPatternRecognition.h"
#include "src/RKFitter/RKfitter.h"
#include "src/RKFitter/TkrData.h"
#include "src/RKFitter/FieldMap.h"
#include "src/RKFitter/RungeKutta4.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("MainRawEventMCDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_MainRawEventMCDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_MainRawEventMCDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_MainRawEventMCDict() {
  TriggerDictionaryInitialization_MainRawEventMCDict_Impl();
}
