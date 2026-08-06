#pragma once

// The product name, as data rather than as a literal.
//
// Policy: the Latin wordmark stays the logo everywhere, and the localized name
// is in-text copy -- window titles, About, prose. This is what Blender, Krita
// and Godot do, and it avoids maintaining five logo lockups.
//
// The exception is Chinese, where a Latin-only name genuinely does not take:
// CJK markets adopt local names whether or not the vendor supplies one, and
// the failure mode this file exists to prevent is users independently
// inventing three different ones. 旋影工坊 is written with four characters that
// are IDENTICAL in Simplified and Traditional, so one name serves both
// scripts, and 旋 ("spiral, revolve") keeps the sense of *Spirula*.
//
// Japanese and Korean pin a transliteration for the same reason -- left alone,
// スピルラ / スピルーラ / スパイルラ all appear -- but keep "Spirula Studio"
// beside it, since that is what the window decoration and the download say.
//
// Russian keeps the Latin wordmark: Russian technical users do not translate
// product names, and a Cyrillic gloss belongs in prose, not in a title bar.

#include "i18n/BeginCatalog.h"

namespace spirula {
namespace i18n {
namespace msg {
namespace brand {

SS_MSG(product,
    EN("Spirula Studio"),      JA("スピルラ・スタジオ"),
    ZH_HANS("旋影工坊"),        ZH_HANT("旋影工坊"),
    KO("스피룰라 스튜디오"),      DE("Spirula Studio"),
    FR("Spirula Studio"),      ES("Spirula Studio"),
    PT("Spirula Studio"),      IT("Spirula Studio"),
    NL("Spirula Studio"),      RU("Spirula Studio"),
    TR("Spirula Studio"));

// The window decoration.
//
// Not `product` on its own: the languages that pin a local name keep the Latin
// wordmark beside it here, because the title bar is where someone matches the
// window to the thing they downloaded and to the process in a task manager.
// The languages that do not translate the name at all just repeat it.
//
// Written out per language rather than assembled from `product` and a
// separator -- a title is copy like any other, and one that reads
// "<name> - Spirula Studio" is a choice each language gets to make.
SS_MSG(window_title,
    EN("Spirula Studio"),
    JA("スピルラ・スタジオ — Spirula Studio"),
    ZH_HANS("旋影工坊 — Spirula Studio"),
    ZH_HANT("旋影工坊 — Spirula Studio"),
    KO("스피룰라 스튜디오 — Spirula Studio"),
    DE("Spirula Studio"),
    FR("Spirula Studio"),
    ES("Spirula Studio"),
    PT("Spirula Studio"),
    IT("Spirula Studio"),
    NL("Spirula Studio"),
    RU("Spirula Studio"),
    TR("Spirula Studio"));

// One line, under the wordmark on the home screen and in About.
SS_MSG(tagline,
    EN("Reconstruct 3D scenes from photos with Gaussian splatting."),
    JA("写真からガウススプラッティングで3Dシーンを再構成します。"),
    ZH_HANS("用高斯泼溅从照片重建三维场景。"),
    ZH_HANT("以高斯潑濺從相片重建三維場景。"),
    KO("사진에서 가우시안 스플래팅으로 3D 장면을 재구성합니다."),
    DE("Rekonstruiert 3D-Szenen aus Fotos mit Gaussian Splatting."),
    FR("Reconstruit des scènes 3D à partir de photos par Gaussian splatting."),
    ES("Reconstruye escenas 3D a partir de fotos con Gaussian splatting."),
    PT("Reconstrói cenas 3D a partir de fotos com Gaussian splatting."),
    IT("Ricostruisce scene 3D da fotografie con il Gaussian splatting."),
    NL("Reconstrueert 3D-scènes uit foto's met Gaussian splatting."),
    RU("Восстанавливает трёхмерные сцены из фотографий методом гауссова сплаттинга."),
    TR("Fotoğraflardan Gaussian splatting ile 3B sahneler oluşturur."));

SS_MSG(about_line,
    EN("Trains 3D Gaussian Splatting models."),
    JA("3Dガウススプラッティングのモデルを学習します。"),
    ZH_HANS("训练三维高斯泼溅模型。"),
    ZH_HANT("訓練三維高斯潑濺模型。"),
    KO("3D 가우시안 스플래팅 모델을 학습합니다."),
    DE("Trainiert 3D-Gaussian-Splatting-Modelle."),
    FR("Entraîne des modèles 3D Gaussian Splatting."),
    ES("Entrena modelos de 3D Gaussian Splatting."),
    PT("Treina modelos de 3D Gaussian Splatting."),
    IT("Addestra modelli 3D Gaussian Splatting."),
    NL("Traint 3D Gaussian Splatting-modellen."),
    RU("Обучает модели 3D Gaussian Splatting."),
    TR("3D Gaussian Splatting modelleri eğitir."));

}  // namespace brand
}  // namespace msg
}  // namespace i18n
}  // namespace spirula

#include "i18n/EndCatalog.h"
