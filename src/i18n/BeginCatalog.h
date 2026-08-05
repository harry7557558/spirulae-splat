// Short per-language tag macros for writing catalog entries. Deliberately
// NOT #pragma once: a translation unit may include several catalogs, and each
// pairs this file with EndCatalog.h, which #undefs the tags again so EN/JA/...
// never leak into ordinary code.
//
//     #include "i18n/BeginCatalog.h"
//     namespace spirula::i18n::msg::gui {
//     SS_MSG(start_training, EN("Start Training"), JA("..."), ...);
//     }
//     #include "i18n/EndCatalog.h"
//
// The tags are written out by hand because a macro cannot define macros. The
// canary at the bottom is what keeps that list honest: add a language to
// SS_LANGUAGES without adding a tag here and the build fails HERE, naming this
// file, rather than on whichever catalog message happened to compile first.

#include "i18n/Message.h"

#define SS_I18N_TR(lang, s) ::spirula::i18n::Tr<::spirula::i18n::Lang::lang>{s}

#define EN(s)      SS_I18N_TR(en, s)
#define JA(s)      SS_I18N_TR(ja, s)
#define ZH_HANS(s) SS_I18N_TR(zh_hans, s)
#define ZH_HANT(s) SS_I18N_TR(zh_hant, s)
#define KO(s)      SS_I18N_TR(ko, s)
#define DE(s)      SS_I18N_TR(de, s)
#define FR(s)      SS_I18N_TR(fr, s)
#define ES(s)      SS_I18N_TR(es, s)
#define PT(s)      SS_I18N_TR(pt, s)
#define IT(s)      SS_I18N_TR(it, s)
#define NL(s)      SS_I18N_TR(nl, s)
#define RU(s)      SS_I18N_TR(ru, s)
#define TR(s)      SS_I18N_TR(tr, s)

#ifndef SS_I18N_TAG_CANARY
#define SS_I18N_TAG_CANARY
namespace spirula {
namespace i18n {
namespace detail {
inline constexpr Msg kTagCanary{
    "kTagCanary",
    EN("x"), JA("x"), ZH_HANS("x"), ZH_HANT("x"), KO("x"), DE("x"), FR("x"),
    ES("x"), PT("x"), IT("x"), NL("x"), RU("x"), TR("x")};
static_assert(kTagCanary.complete(),
              "i18n: src/i18n/BeginCatalog.h is missing a tag macro for a "
              "language in SS_LANGUAGES -- add it there too");
}  // namespace detail
}  // namespace i18n
}  // namespace spirula
#endif
