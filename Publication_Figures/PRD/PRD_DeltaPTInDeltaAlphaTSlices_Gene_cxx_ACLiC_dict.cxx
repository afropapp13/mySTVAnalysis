// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIhomedIafroditidIDropboxdIPhDdImyCodedI21th_assignment_CalibratedProductsdICodeRootFilesdIuboonecode_v08dImySTVAnalysisdIPublication_FiguresdIPRDdIPRD_DeltaPTInDeltaAlphaTSlices_Gene_cxx_ACLiC_dict
#define R__NO_DEPRECATION

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

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/home/afroditi/Dropbox/PhD/myCode/21th_assignment_CalibratedProducts/CodeRootFiles/uboonecode_v08/mySTVAnalysis/Publication_Figures/PRD/./PRD_DeltaPTInDeltaAlphaTSlices_Gene.cxx"

// Header files passed via #pragma extra_include

namespace {
  void TriggerDictionaryInitialization_PRD_DeltaPTInDeltaAlphaTSlices_Gene_cxx_ACLiC_dict_Impl() {
    static const char* headers[] = {
"./PRD_DeltaPTInDeltaAlphaTSlices_Gene.cxx",
0
    };
    static const char* includePaths[] = {
"/home/afroditi/root_install/include",
"/home/afroditi/root_install/etc/",
"/home/afroditi/root_install/etc//cling",
"/home/afroditi/root_install/include/",
"/home/afroditi/root_install/include",
"/home/afroditi/root_install/include/",
"/home/afroditi/Dropbox/PhD/myCode/21th_assignment_CalibratedProducts/CodeRootFiles/uboonecode_v08/mySTVAnalysis/Publication_Figures/PRD/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PRD_DeltaPTInDeltaAlphaTSlices_Gene_cxx_ACLiC_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PRD_DeltaPTInDeltaAlphaTSlices_Gene_cxx_ACLiC_dict dictionary payload"

#ifndef __ACLIC__
  #define __ACLIC__ 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "./PRD_DeltaPTInDeltaAlphaTSlices_Gene.cxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"Mapping", payloadCode, "@",
"PRD_DeltaPTInDeltaAlphaTSlices_Gene", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PRD_DeltaPTInDeltaAlphaTSlices_Gene_cxx_ACLiC_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PRD_DeltaPTInDeltaAlphaTSlices_Gene_cxx_ACLiC_dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PRD_DeltaPTInDeltaAlphaTSlices_Gene_cxx_ACLiC_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PRD_DeltaPTInDeltaAlphaTSlices_Gene_cxx_ACLiC_dict() {
  TriggerDictionaryInitialization_PRD_DeltaPTInDeltaAlphaTSlices_Gene_cxx_ACLiC_dict_Impl();
}