#pragma once

// What the progress line and the log panel say while a job runs.
//
// Two kinds of text end up in that panel and only one of them is here:
//
//   OURS -- the stage names, the notes, the "here is what I did with your
//   folder" lines. They are written for the person watching, so they are
//   translated, and the stage names doubly so: the same string is the caption
//   above the progress bar.
//
//   THE CHILD PROCESS'S -- every line COLMAP or ffmpeg prints. Passed through
//   verbatim. They are English, they are what a bug report is pasted from, and
//   they are not ours to rewrite.
//
// So a log in Japanese is Japanese around English, which is honest: the
// English parts are the ones that came from somewhere else.
//
// `spirula sfm` is the exception that proves the rule -- it is a child process
// but it is OURS, so it translates itself, out of i18n/catalog/Sfm.h and
// through src/sfm/core/Log.h.

#include "i18n/BeginCatalog.h"

namespace spirula {
namespace i18n {
namespace msg {
namespace log {

// ===========================================================================
// Stages -- also the caption above the progress bar
// ===========================================================================

SS_MSG(stage_collecting_photos,
    EN("Collecting photos"),
    JA("写真を集めています"),
    ZH_HANS("正在收集照片"),
    ZH_HANT("正在收集照片"),
    KO("사진을 모으는 중"),
    DE("Fotos werden zusammengetragen"),
    FR("Collecte des photos"),
    ES("Recopilando las fotos"),
    PT("Reunindo as fotos"),
    IT("Raccolta delle foto"),
    NL("Foto's verzamelen"),
    RU("Сбор фотографий"),
    TR("Fotoğraflar toplanıyor"));

SS_MSG(stage_extract_gpu,
    EN("Extracting frames (GPU decode)"),
    JA("フレームを取り出しています（GPUデコード）"),
    ZH_HANS("正在提取帧（GPU 解码）"),
    ZH_HANT("正在擷取影格（GPU 解碼）"),
    KO("프레임을 뽑는 중(GPU 디코딩)"),
    DE("Einzelbilder werden entnommen (GPU-Dekodierung)"),
    FR("Extraction des images (décodage GPU)"),
    ES("Extrayendo fotogramas (decodificación por GPU)"),
    PT("Extraindo quadros (decodificação na GPU)"),
    IT("Estrazione dei fotogrammi (decodifica GPU)"),
    NL("Beelden uitpakken (GPU-decodering)"),
    RU("Извлечение кадров (декодирование на GPU)"),
    TR("Kareler çıkarılıyor (GPU çözümü)"));

SS_MSG(stage_extract_mask_gpu,
    EN("Extracting frames and masking (GPU)"),
    JA("フレームの取り出しとマスク作成をしています（GPU）"),
    ZH_HANS("正在提取帧并生成蒙版（GPU）"),
    ZH_HANT("正在擷取影格並產生遮罩（GPU）"),
    KO("프레임을 뽑고 마스크를 만드는 중(GPU)"),
    DE("Einzelbilder entnehmen und maskieren (GPU)"),
    FR("Extraction des images et masquage (GPU)"),
    ES("Extrayendo fotogramas y enmascarando (GPU)"),
    PT("Extraindo quadros e mascarando (GPU)"),
    IT("Estrazione dei fotogrammi e mascheratura (GPU)"),
    NL("Beelden uitpakken en maskeren (GPU)"),
    RU("Извлечение кадров и маскирование (GPU)"),
    TR("Kareler çıkarılıyor ve maskeleniyor (GPU)"));

SS_MSG(stage_extract_candidates,
    EN("Extracting candidate frames (ffmpeg)"),
    JA("候補となるフレームを取り出しています（ffmpeg）"),
    ZH_HANS("正在提取候选帧（ffmpeg）"),
    ZH_HANT("正在擷取候選影格（ffmpeg）"),
    KO("후보 프레임을 뽑는 중(ffmpeg)"),
    DE("Kandidatenbilder werden entnommen (ffmpeg)"),
    FR("Extraction des images candidates (ffmpeg)"),
    ES("Extrayendo fotogramas candidatos (ffmpeg)"),
    PT("Extraindo quadros candidatos (ffmpeg)"),
    IT("Estrazione dei fotogrammi candidati (ffmpeg)"),
    NL("Kandidaatbeelden uitpakken (ffmpeg)"),
    RU("Извлечение кадров-кандидатов (ffmpeg)"),
    TR("Aday kareler çıkarılıyor (ffmpeg)"));

SS_MSG(stage_extract_ffmpeg,
    EN("Extracting frames (ffmpeg)"),
    JA("フレームを取り出しています（ffmpeg）"),
    ZH_HANS("正在提取帧（ffmpeg）"),
    ZH_HANT("正在擷取影格（ffmpeg）"),
    KO("프레임을 뽑는 중(ffmpeg)"),
    DE("Einzelbilder werden entnommen (ffmpeg)"),
    FR("Extraction des images (ffmpeg)"),
    ES("Extrayendo fotogramas (ffmpeg)"),
    PT("Extraindo quadros (ffmpeg)"),
    IT("Estrazione dei fotogrammi (ffmpeg)"),
    NL("Beelden uitpakken (ffmpeg)"),
    RU("Извлечение кадров (ffmpeg)"),
    TR("Kareler çıkarılıyor (ffmpeg)"));

// {0} is the track number inside a multi-track file (an Insta360 .insv holds
// one per lens).
SS_MSG(stage_split_track,
    EN("Splitting video track {0}"),
    JA("動画のトラック {0} を分けています"),
    ZH_HANS("正在拆分视频轨道 {0}"),
    ZH_HANT("正在拆分影片軌道 {0}"),
    KO("영상 트랙 {0} 분리 중"),
    DE("Videospur {0} wird aufgeteilt"),
    FR("Séparation de la piste vidéo {0}"),
    ES("Separando la pista de vídeo {0}"),
    PT("Separando a trilha de vídeo {0}"),
    IT("Separazione della traccia video {0}"),
    NL("Videospoor {0} splitsen"),
    RU("Разделение видеодорожки {0}"),
    TR("{0} numaralı video izi ayrılıyor"));

SS_MSG(stage_select_sharpest,
    EN("Selecting sharpest frames (multithreaded)"),
    JA("いちばん鮮明なフレームを選んでいます（マルチスレッド）"),
    ZH_HANS("正在挑选最清晰的帧（多线程）"),
    ZH_HANT("正在挑選最清晰的影格（多執行緒）"),
    KO("가장 선명한 프레임을 고르는 중(멀티스레드)"),
    DE("Schärfste Einzelbilder werden ausgewählt (mehrere Threads)"),
    FR("Sélection des images les plus nettes (multithread)"),
    ES("Seleccionando los fotogramas más nítidos (multihilo)"),
    PT("Selecionando os quadros mais nítidos (multithread)"),
    IT("Selezione dei fotogrammi più nitidi (multithread)"),
    NL("Scherpste beelden kiezen (meerdere threads)"),
    RU("Выбор самых резких кадров (в несколько потоков)"),
    TR("En net kareler seçiliyor (çok iş parçacıklı)"));

SS_MSG(stage_masks_builtin,
    EN("Generating masks (segmentation)"),
    JA("マスクを作成しています（セグメンテーション）"),
    ZH_HANS("正在生成蒙版（分割）"),
    ZH_HANT("正在產生遮罩（分割）"),
    KO("마스크를 만드는 중(분할)"),
    DE("Masken werden erzeugt (Segmentierung)"),
    FR("Génération des masques (segmentation)"),
    ES("Generando las máscaras (segmentación)"),
    PT("Gerando as máscaras (segmentação)"),
    IT("Generazione delle maschere (segmentazione)"),
    NL("Maskers maken (segmentatie)"),
    RU("Создание масок (сегментация)"),
    TR("Maskeler oluşturuluyor (bölütleme)"));

SS_MSG(stage_masks_python,
    EN("Generating masks (external Python)"),
    JA("マスクを作成しています（外部のPython）"),
    ZH_HANS("正在生成蒙版（外部 Python）"),
    ZH_HANT("正在產生遮罩（外部 Python）"),
    KO("마스크를 만드는 중(외부 Python)"),
    DE("Masken werden erzeugt (externes Python)"),
    FR("Génération des masques (Python externe)"),
    ES("Generando las máscaras (Python externo)"),
    PT("Gerando as máscaras (Python externo)"),
    IT("Generazione delle maschere (Python esterno)"),
    NL("Maskers maken (extern Python)"),
    RU("Создание масок (внешний Python)"),
    TR("Maskeler oluşturuluyor (harici Python)"));

SS_MSG(stage_finding_features,
    EN("Finding features"),
    JA("特徴点を探しています"),
    ZH_HANS("正在寻找特征点"),
    ZH_HANT("正在尋找特徵點"),
    KO("특징점을 찾는 중"),
    DE("Merkmale werden gesucht"),
    FR("Recherche des points caractéristiques"),
    ES("Buscando puntos característicos"),
    PT("Procurando pontos característicos"),
    IT("Ricerca dei punti caratteristici"),
    NL("Kenmerken zoeken"),
    RU("Поиск особых точек"),
    TR("Öznitelikler aranıyor"));

SS_MSG(stage_matching_images,
    EN("Matching images"),
    JA("画像どうしを照合しています"),
    ZH_HANS("正在匹配图像"),
    ZH_HANT("正在比對影像"),
    KO("이미지끼리 맞춰 보는 중"),
    DE("Bilder werden einander zugeordnet"),
    FR("Mise en correspondance des images"),
    ES("Emparejando las imágenes"),
    PT("Correspondendo as imagens"),
    IT("Corrispondenza tra le immagini"),
    NL("Beelden aan elkaar koppelen"),
    RU("Сопоставление снимков"),
    TR("Görüntüler eşleştiriliyor"));

SS_MSG(stage_reconstructing,
    EN("Reconstructing cameras (the slow part)"),
    JA("カメラ位置を復元しています（時間のかかる工程）"),
    ZH_HANS("正在恢复相机位姿（最慢的一步）"),
    ZH_HANT("正在還原相機位姿（最慢的一步）"),
    KO("카메라 위치를 복원하는 중(가장 오래 걸리는 단계)"),
    DE("Kamerastandpunkte werden rekonstruiert (der langsame Teil)"),
    FR("Reconstruction des caméras (l'étape lente)"),
    ES("Reconstruyendo las cámaras (la parte lenta)"),
    PT("Reconstruindo as câmeras (a parte lenta)"),
    IT("Ricostruzione delle fotocamere (la parte lenta)"),
    NL("Camerastandpunten reconstrueren (het trage deel)"),
    RU("Восстановление положений камер (самый долгий этап)"),
    TR("Kamera konumları yeniden kuruluyor (yavaş kısım)"));

SS_MSG(stage_reconstructing_features,
    EN("Reconstructing (finding features)"),
    JA("復元しています（特徴点の抽出）"),
    ZH_HANS("正在重建（寻找特征点）"),
    ZH_HANT("正在重建（尋找特徵點）"),
    KO("복원하는 중(특징점 찾기)"),
    DE("Rekonstruktion (Merkmalssuche)"),
    FR("Reconstruction (recherche des points caractéristiques)"),
    ES("Reconstruyendo (buscando puntos característicos)"),
    PT("Reconstruindo (procurando pontos característicos)"),
    IT("Ricostruzione (ricerca dei punti caratteristici)"),
    NL("Reconstructie (kenmerken zoeken)"),
    RU("Восстановление (поиск особых точек)"),
    TR("Yeniden kuruluyor (öznitelik arama)"));

SS_MSG(stage_cleaning_up,
    EN("Cleaning up"),
    JA("後片付けをしています"),
    ZH_HANS("正在清理"),
    ZH_HANT("正在清理"),
    KO("정리하는 중"),
    DE("Aufräumen"),
    FR("Nettoyage"),
    ES("Limpiando"),
    PT("Limpando"),
    IT("Pulizia"),
    NL("Opruimen"),
    RU("Очистка"),
    TR("Temizleniyor"));

SS_MSG(stage_done,
    EN("Done"),
    JA("完了"),
    ZH_HANS("完成"),
    ZH_HANT("完成"),
    KO("완료"),
    DE("Fertig"),
    FR("Terminé"),
    ES("Listo"),
    PT("Concluído"),
    IT("Fatto"),
    NL("Klaar"),
    RU("Готово"),
    TR("Bitti"));

// ---- COLMAP's stages ------------------------------------------------------

SS_MSG(stage_vocab_download,
    EN("Downloading vocabulary tree (one-time, ~150 MB)"),
    JA("ボキャブラリツリーをダウンロードしています（初回のみ、約150 MB）"),
    ZH_HANS("正在下载词汇树（仅一次，约 150 MB）"),
    ZH_HANT("正在下載詞彙樹（僅一次，約 150 MB）"),
    KO("어휘 트리를 내려받는 중(최초 1회, 약 150 MB)"),
    DE("Vokabelbaum wird heruntergeladen (einmalig, ca. 150 MB)"),
    FR("Téléchargement de l'arbre de vocabulaire (une seule fois, ~150 Mo)"),
    ES("Descargando el árbol de vocabulario (una sola vez, ~150 MB)"),
    PT("Baixando a árvore de vocabulário (uma única vez, ~150 MB)"),
    IT("Download dell'albero di vocabolario (una sola volta, ~150 MB)"),
    NL("Vocabulaireboom downloaden (eenmalig, ~150 MB)"),
    RU("Загрузка словарного дерева (один раз, ~150 МБ)"),
    TR("Sözcük ağacı indiriliyor (tek seferlik, ~150 MB)"));

SS_MSG(stage_colmap_features,
    EN("Extracting features (colmap)"),
    JA("特徴点を抽出しています（colmap）"),
    ZH_HANS("正在提取特征点（colmap）"),
    ZH_HANT("正在擷取特徵點（colmap）"),
    KO("특징점을 뽑는 중(colmap)"),
    DE("Merkmale werden extrahiert (colmap)"),
    FR("Extraction des points caractéristiques (colmap)"),
    ES("Extrayendo puntos característicos (colmap)"),
    PT("Extraindo pontos característicos (colmap)"),
    IT("Estrazione dei punti caratteristici (colmap)"),
    NL("Kenmerken extraheren (colmap)"),
    RU("Извлечение особых точек (colmap)"),
    TR("Öznitelikler çıkarılıyor (colmap)"));

SS_MSG(stage_colmap_features_aliked,
    EN("Extracting features (colmap, ALIKED)"),
    JA("特徴点を抽出しています（colmap、ALIKED）"),
    ZH_HANS("正在提取特征点（colmap、ALIKED）"),
    ZH_HANT("正在擷取特徵點（colmap、ALIKED）"),
    KO("특징점을 뽑는 중(colmap, ALIKED)"),
    DE("Merkmale werden extrahiert (colmap, ALIKED)"),
    FR("Extraction des points caractéristiques (colmap, ALIKED)"),
    ES("Extrayendo puntos característicos (colmap, ALIKED)"),
    PT("Extraindo pontos característicos (colmap, ALIKED)"),
    IT("Estrazione dei punti caratteristici (colmap, ALIKED)"),
    NL("Kenmerken extraheren (colmap, ALIKED)"),
    RU("Извлечение особых точек (colmap, ALIKED)"),
    TR("Öznitelikler çıkarılıyor (colmap, ALIKED)"));

SS_MSG(stage_match_vocab,
    EN("Matching features (vocabulary tree)"),
    JA("特徴点を照合しています（ボキャブラリツリー）"),
    ZH_HANS("正在匹配特征点（词汇树）"),
    ZH_HANT("正在比對特徵點（詞彙樹）"),
    KO("특징점을 맞춰 보는 중(어휘 트리)"),
    DE("Merkmale werden zugeordnet (Vokabelbaum)"),
    FR("Mise en correspondance (arbre de vocabulaire)"),
    ES("Emparejando puntos característicos (árbol de vocabulario)"),
    PT("Correspondendo pontos característicos (árvore de vocabulário)"),
    IT("Corrispondenza dei punti caratteristici (albero di vocabolario)"),
    NL("Kenmerken koppelen (vocabulaireboom)"),
    RU("Сопоставление особых точек (словарное дерево)"),
    TR("Öznitelikler eşleştiriliyor (sözcük ağacı)"));

SS_MSG(stage_match_sequential,
    EN("Matching features (sequential)"),
    JA("特徴点を照合しています（連続フレーム）"),
    ZH_HANS("正在匹配特征点（顺序）"),
    ZH_HANT("正在比對特徵點（順序）"),
    KO("특징점을 맞춰 보는 중(순차)"),
    DE("Merkmale werden zugeordnet (fortlaufend)"),
    FR("Mise en correspondance (séquentielle)"),
    ES("Emparejando puntos característicos (secuencial)"),
    PT("Correspondendo pontos característicos (sequencial)"),
    IT("Corrispondenza dei punti caratteristici (sequenziale)"),
    NL("Kenmerken koppelen (opeenvolgend)"),
    RU("Сопоставление особых точек (последовательное)"),
    TR("Öznitelikler eşleştiriliyor (ardışık)"));

SS_MSG(stage_match_exhaustive,
    EN("Matching features (exhaustive)"),
    JA("特徴点を照合しています（総当たり）"),
    ZH_HANS("正在匹配特征点（穷举）"),
    ZH_HANT("正在比對特徵點（窮舉）"),
    KO("특징점을 맞춰 보는 중(전수 비교)"),
    DE("Merkmale werden zugeordnet (vollständig)"),
    FR("Mise en correspondance (exhaustive)"),
    ES("Emparejando puntos característicos (exhaustivo)"),
    PT("Correspondendo pontos característicos (exaustivo)"),
    IT("Corrispondenza dei punti caratteristici (esaustiva)"),
    NL("Kenmerken koppelen (uitputtend)"),
    RU("Сопоставление особых точек (полный перебор)"),
    TR("Öznitelikler eşleştiriliyor (tam karşılaştırma)"));

SS_MSG(stage_colmap_mapper,
    EN("Reconstructing cameras (mapper; this is the slow part)"),
    JA("カメラ位置を復元しています（mapper。時間のかかる工程です）"),
    ZH_HANS("正在恢复相机位姿（mapper，这一步最慢）"),
    ZH_HANT("正在還原相機位姿（mapper，這一步最慢）"),
    KO("카메라 위치를 복원하는 중(mapper, 가장 오래 걸리는 단계)"),
    DE("Kamerastandpunkte werden rekonstruiert (mapper; der langsame Teil)"),
    FR("Reconstruction des caméras (mapper ; c'est l'étape lente)"),
    ES("Reconstruyendo las cámaras (mapper; esta es la parte lenta)"),
    PT("Reconstruindo as câmeras (mapper; esta é a parte lenta)"),
    IT("Ricostruzione delle fotocamere (mapper; è la parte lenta)"),
    NL("Camerastandpunten reconstrueren (mapper; dit is het trage deel)"),
    RU("Восстановление положений камер (mapper; это самый долгий этап)"),
    TR("Kamera konumları yeniden kuruluyor (mapper; yavaş kısım budur)"));

SS_MSG(stage_merge_models,
    EN("Merging partial models"),
    JA("部分的なモデルをつなぎ合わせています"),
    ZH_HANS("正在合并零散的模型"),
    ZH_HANT("正在合併零散的模型"),
    KO("조각난 모델을 합치는 중"),
    DE("Teilmodelle werden zusammengeführt"),
    FR("Fusion des modèles partiels"),
    ES("Uniendo los modelos parciales"),
    PT("Juntando os modelos parciais"),
    IT("Unione dei modelli parziali"),
    NL("Deelmodellen samenvoegen"),
    RU("Объединение частичных моделей"),
    TR("Parçalı modeller birleştiriliyor"));

SS_MSG(stage_bundle_adjust,
    EN("Refining cameras (bundle adjustment)"),
    JA("カメラ位置を微調整しています（バンドル調整）"),
    ZH_HANS("正在微调相机位姿（光束法平差）"),
    ZH_HANT("正在微調相機位姿（光束法平差）"),
    KO("카메라 위치를 다듬는 중(번들 조정)"),
    DE("Kamerastandpunkte werden verfeinert (Bündelblockausgleichung)"),
    FR("Affinement des caméras (ajustement de faisceaux)"),
    ES("Afinando las cámaras (ajuste de haces)"),
    PT("Refinando as câmeras (ajuste de feixes)"),
    IT("Affinamento delle fotocamere (bundle adjustment)"),
    NL("Camerastandpunten verfijnen (bundeladjustering)"),
    RU("Уточнение положений камер (уравнивание связок)"),
    TR("Kamera konumları iyileştiriliyor (demet dengelemesi)"));

// ===========================================================================
// Notes -- what happened to the user's folder, and what to do about it
// ===========================================================================

// {0} the folder, {1} how many images.
SS_MSG(found_images,
    EN("Found {0} images in {1}"),
    JA("{1} に画像が {0} 枚見つかりました"),
    ZH_HANS("在 {1} 中找到 {0} 张图像"),
    ZH_HANT("在 {1} 中找到 {0} 張影像"),
    KO("{1}에서 이미지 {0}장을 찾았습니다"),
    DE("{0} Bilder in {1} gefunden"),
    FR("{0} images trouvées dans {1}"),
    ES("Se encontraron {0} imágenes en {1}"),
    PT("Encontradas {0} imagens em {1}"),
    IT("Trovate {0} immagini in {1}"),
    NL("{0} beelden gevonden in {1}"),
    RU("Найдено изображений: {0} (в {1})"),
    TR("{1} içinde {0} görüntü bulundu"));

SS_MSG(using_bundled_masks,
    EN("Using the masks that came with the photos: {0}"),
    JA("写真に付いていたマスクを使います: {0}"),
    ZH_HANS("使用随照片一起提供的蒙版：{0}"),
    ZH_HANT("使用隨照片一起提供的遮罩：{0}"),
    KO("사진에 딸려 온 마스크를 사용합니다: {0}"),
    DE("Die mitgelieferten Masken werden verwendet: {0}"),
    FR("Utilisation des masques fournis avec les photos : {0}"),
    ES("Se usan las máscaras que venían con las fotos: {0}"),
    PT("Usando as máscaras que vieram com as fotos: {0}"),
    IT("Si usano le maschere fornite con le foto: {0}"),
    NL("De meegeleverde maskers worden gebruikt: {0}"),
    RU("Используются маски, приложенные к фотографиям: {0}"),
    TR("Fotoğraflarla birlikte gelen maskeler kullanılıyor: {0}"));

SS_MSG(video_input,
    EN("Video: {0}"),
    JA("動画: {0}"),
    ZH_HANS("视频：{0}"),
    ZH_HANT("影片：{0}"),
    KO("영상: {0}"),
    DE("Video: {0}"),
    FR("Vidéo : {0}"),
    ES("Vídeo: {0}"),
    PT("Vídeo: {0}"),
    IT("Video: {0}"),
    NL("Video: {0}"),
    RU("Видео: {0}"),
    TR("Video: {0}"));

// {0} how many, {1} the folder.
SS_MSG(resume_keep_frames,
    EN("Resume: keeping {0} extracted frames in {1} (delete the folder to "
       "re-extract)"),
    JA("再開: {1} にある取り出し済みのフレーム {0} 枚をそのまま使います"
       "（取り直すにはフォルダーを削除してください）"),
    ZH_HANS("继续：保留 {1} 中已提取的 {0} 帧（想重新提取请删除该文件夹）"),
    ZH_HANT("繼續：保留 {1} 中已擷取的 {0} 個影格（想重新擷取請刪除該資料夾）"),
    KO("이어서 진행: {1}에 이미 뽑아 둔 프레임 {0}장을 그대로 씁니다(다시 "
       "뽑으려면 폴더를 지우세요)"),
    DE("Fortsetzen: die {0} bereits entnommenen Einzelbilder in {1} bleiben "
       "(zum erneuten Entnehmen den Ordner löschen)"),
    FR("Reprise : les {0} images déjà extraites dans {1} sont conservées "
       "(supprimez le dossier pour les réextraire)"),
    ES("Reanudar: se conservan los {0} fotogramas ya extraídos en {1} (borre "
       "la carpeta para volver a extraerlos)"),
    PT("Retomar: mantendo os {0} quadros já extraídos em {1} (apague a pasta "
       "para extrair de novo)"),
    IT("Ripresa: si tengono i {0} fotogrammi già estratti in {1} (cancelli la "
       "cartella per riestrarli)"),
    NL("Hervatten: de {0} al uitgepakte beelden in {1} blijven staan (verwijder "
       "de map om opnieuw uit te pakken)"),
    RU("Продолжение: оставляем {0} уже извлечённых кадров в {1} (чтобы извлечь "
       "заново, удалите папку)"),
    TR("Sürdürme: {1} içindeki {0} çıkarılmış kare korunuyor (yeniden çıkarmak "
       "için klasörü silin)"));

SS_MSG(resume_keep_masks,
    EN("Resume: keeping the masks in {0}"),
    JA("再開: {0} にあるマスクをそのまま使います"),
    ZH_HANS("继续：保留 {0} 中的蒙版"),
    ZH_HANT("繼續：保留 {0} 中的遮罩"),
    KO("이어서 진행: {0}에 있는 마스크를 그대로 씁니다"),
    DE("Fortsetzen: die Masken in {0} bleiben"),
    FR("Reprise : les masques dans {0} sont conservés"),
    ES("Reanudar: se conservan las máscaras en {0}"),
    PT("Retomar: mantendo as máscaras em {0}"),
    IT("Ripresa: si tengono le maschere in {0}"),
    NL("Hervatten: de maskers in {0} blijven staan"),
    RU("Продолжение: маски в {0} остаются"),
    TR("Sürdürme: {0} içindeki maskeler korunuyor"));

SS_MSG(resume_keep_frames_dir,
    EN("Resume: keeping the frames already in {0}"),
    JA("再開: すでに {0} にあるフレームをそのまま使います"),
    ZH_HANS("继续：保留 {0} 中已有的帧"),
    ZH_HANT("繼續：保留 {0} 中已有的影格"),
    KO("이어서 진행: {0}에 이미 있는 프레임을 그대로 씁니다"),
    DE("Fortsetzen: die Einzelbilder in {0} bleiben"),
    FR("Reprise : les images déjà présentes dans {0} sont conservées"),
    ES("Reanudar: se conservan los fotogramas que ya hay en {0}"),
    PT("Retomar: mantendo os quadros que já estão em {0}"),
    IT("Ripresa: si tengono i fotogrammi già presenti in {0}"),
    NL("Hervatten: de beelden die al in {0} staan blijven staan"),
    RU("Продолжение: кадры, уже лежащие в {0}, остаются"),
    TR("Sürdürme: {0} içinde zaten bulunan kareler korunuyor"));

// {0} the decoder's own message, English.
SS_MSG(decode_fallback_ffmpeg,
    EN("Built-in decoding could not handle this file ({0}); falling back to "
       "ffmpeg"),
    JA("内蔵のデコーダではこのファイルを扱えませんでした（{0}）。ffmpeg に"
       "切り替えます"),
    ZH_HANS("内置解码无法处理这个文件（{0}），改用 ffmpeg"),
    ZH_HANT("內建解碼無法處理這個檔案（{0}），改用 ffmpeg"),
    KO("내장 디코더로는 이 파일을 다룰 수 없었습니다({0}). ffmpeg으로 "
       "넘어갑니다"),
    DE("Die eingebaute Dekodierung kam mit dieser Datei nicht zurecht ({0}); es "
       "wird auf ffmpeg zurückgegriffen"),
    FR("Le décodage intégré n'a pas su traiter ce fichier ({0}) ; repli sur "
       "ffmpeg"),
    ES("La decodificación integrada no pudo con este archivo ({0}); se recurre "
       "a ffmpeg"),
    PT("A decodificação integrada não deu conta deste arquivo ({0}); recorrendo "
       "ao ffmpeg"),
    IT("La decodifica integrata non ha gestito questo file ({0}); si ripiega su "
       "ffmpeg"),
    NL("De ingebouwde decodering kon dit bestand niet aan ({0}); er wordt "
       "teruggevallen op ffmpeg"),
    RU("Встроенный декодер не справился с этим файлом ({0}); переходим на "
       "ffmpeg"),
    TR("Yerleşik çözücü bu dosyayla baş edemedi ({0}); ffmpeg'e geçiliyor"));

SS_MSG(kept_frames,
    EN("Kept {0} frames -> {1}"),
    JA("フレームを {0} 枚残しました -> {1}"),
    ZH_HANS("保留了 {0} 帧 -> {1}"),
    ZH_HANT("保留了 {0} 個影格 -> {1}"),
    KO("프레임 {0}장을 남겼습니다 -> {1}"),
    DE("{0} Einzelbilder behalten -> {1}"),
    FR("{0} images conservées -> {1}"),
    ES("Se conservaron {0} fotogramas -> {1}"),
    PT("Mantidos {0} quadros -> {1}"),
    IT("Tenuti {0} fotogrammi -> {1}"),
    NL("{0} beelden bewaard -> {1}"),
    RU("Оставлено кадров: {0} -> {1}"),
    TR("{0} kare tutuldu -> {1}"));

SS_MSG(linked_copied_kept,
    EN("  {0} linked, {1} copied, {2} already there"),
    JA("  リンク {0} 件、コピー {1} 件、既存 {2} 件"),
    ZH_HANS("  链接 {0} 个，复制 {1} 个，已有 {2} 个"),
    ZH_HANT("  連結 {0} 個，複製 {1} 個，已有 {2} 個"),
    KO("  링크 {0}개, 복사 {1}개, 이미 있던 것 {2}개"),
    DE("  {0} verknüpft, {1} kopiert, {2} schon vorhanden"),
    FR("  {0} liées, {1} copiées, {2} déjà présentes"),
    ES("  {0} enlazadas, {1} copiadas, {2} ya estaban"),
    PT("  {0} vinculadas, {1} copiadas, {2} já estavam lá"),
    IT("  {0} collegate, {1} copiate, {2} già presenti"),
    NL("  {0} gekoppeld, {1} gekopieerd, {2} stonden er al"),
    RU("  связано: {0}, скопировано: {1}, уже было: {2}"),
    TR("  {0} bağlandı, {1} kopyalandı, {2} zaten vardı"));

SS_MSG(clicks_other_input_unmasked,
    EN("Note: {0} is not the input the clicked objects were drawn on and there "
       "is no text prompt, so its frames are left unmasked. Add a prompt, or "
       "click the object on this input too."),
    JA("メモ: {0} はクリックで指定した対象を描いた入力ではなく、文字の指定も"
       "ないため、この入力のフレームはマスクされません。プロンプトを入れるか、"
       "この入力でも対象をクリックしてください。"),
    ZH_HANS("提示：{0} 不是当初点选对象所在的输入，也没有文字提示，所以它的帧"
            "不会生成蒙版。请填写提示词，或在这个输入上也点选一次对象。"),
    ZH_HANT("提示：{0} 不是當初點選對象所在的輸入，也沒有文字提示，所以它的影格"
            "不會產生遮罩。請填寫提示詞，或在這個輸入上也點選一次對象。"),
    KO("참고: {0}은(는) 클릭으로 지정한 대상을 그린 입력이 아니고 텍스트 "
       "프롬프트도 없어서, 이 입력의 프레임에는 마스크가 생기지 않습니다. "
       "프롬프트를 넣거나 이 입력에서도 대상을 클릭하세요."),
    DE("Hinweis: {0} ist nicht die Eingabe, auf der die Objekte angeklickt "
       "wurden, und es gibt keinen Texthinweis -- ihre Bilder bleiben also "
       "unmaskiert. Geben Sie einen Prompt ein, oder klicken Sie das Objekt "
       "auch auf dieser Eingabe an."),
    FR("Note : {0} n'est pas l'entrée sur laquelle les objets ont été cliqués, "
       "et il n'y a pas d'invite textuelle ; ses images restent donc non "
       "masquées. Ajoutez une invite, ou cliquez aussi l'objet sur cette "
       "entrée."),
    ES("Nota: {0} no es la entrada sobre la que se pulsaron los objetos y no "
       "hay indicación de texto, así que sus fotogramas quedan sin máscara. "
       "Añada una indicación, o pulse el objeto también en esta entrada."),
    PT("Nota: {0} não é a entrada em que os objetos foram clicados e não há "
       "texto de comando, então seus quadros ficam sem máscara. Acrescente um "
       "comando, ou clique no objeto também nesta entrada."),
    IT("Nota: {0} non è l'ingresso su cui gli oggetti sono stati cliccati e non "
       "c'è un prompt testuale, quindi i suoi fotogrammi restano senza "
       "maschera. Aggiunga un prompt, oppure clicchi l'oggetto anche su questo "
       "ingresso."),
    NL("Let op: {0} is niet de invoer waarop de objecten zijn aangeklikt en er "
       "is geen tekstprompt, dus de beelden ervan blijven ongemaskeerd. Voeg "
       "een prompt toe, of klik het object ook op deze invoer aan."),
    RU("Примечание: {0} — не тот вход, на котором отмечались объекты, и "
       "текстового запроса нет, поэтому его кадры остаются без масок. Введите "
       "запрос или отметьте объект и на этом входе."),
    TR("Not: {0}, nesnelerin tıklandığı girdi değil ve metin istemi de yok; bu "
       "yüzden kareleri maskesiz kalıyor. Bir istem ekleyin ya da nesneyi bu "
       "girdide de tıklayın."));

SS_MSG(warn_unreadable_skipped,
    EN("warning: could not read {0}; skipped"),
    JA("警告: {0} を読み込めませんでした。とばします"),
    ZH_HANS("警告：无法读取 {0}，已跳过"),
    ZH_HANT("警告：無法讀取 {0}，已略過"),
    KO("경고: {0}을(를) 읽지 못해 건너뜁니다"),
    DE("Warnung: {0} konnte nicht gelesen werden; übersprungen"),
    FR("Avertissement : {0} n'a pas pu être lu ; ignoré"),
    ES("Aviso: no se pudo leer {0}; omitido"),
    PT("Aviso: não foi possível ler {0}; ignorado"),
    IT("Avviso: non è stato possibile leggere {0}; saltato"),
    NL("Waarschuwing: {0} kon niet worden gelezen; overgeslagen"),
    RU("Предупреждение: не удалось прочитать {0}; пропущено"),
    TR("Uyarı: {0} okunamadı; atlandı"));

// ---- the built-in reconstruction's notes ----------------------------------

SS_MSG(sfm_focal_unreadable,
    EN("warning: could not read an image in {0}; leaving its focal length to "
       "be guessed"),
    JA("警告: {0} の画像を読み込めませんでした。焦点距離は推測に任せます"),
    ZH_HANS("警告：无法读取 {0} 里的图像，焦距交给自动推测"),
    ZH_HANT("警告：無法讀取 {0} 裡的影像，焦距交給自動推測"),
    KO("경고: {0}의 이미지를 읽지 못했습니다. 초점 거리는 추정에 맡깁니다"),
    DE("Warnung: In {0} konnte kein Bild gelesen werden; die Brennweite wird "
       "geschätzt"),
    FR("Avertissement : aucune image lisible dans {0} ; la focale sera devinée"),
    ES("Aviso: no se pudo leer ninguna imagen en {0}; la focal se dejará "
       "adivinar"),
    PT("Aviso: não foi possível ler nenhuma imagem em {0}; a distância focal "
       "será adivinhada"),
    IT("Avviso: non è stato possibile leggere alcuna immagine in {0}; la "
       "focale sarà indovinata"),
    NL("Waarschuwing: geen leesbaar beeld in {0}; de brandpuntsafstand wordt "
       "geraden"),
    RU("Предупреждение: не удалось прочитать ни одного изображения в {0}; "
       "фокусное расстояние будет угадано"),
    TR("Uyarı: {0} içinde okunabilir görüntü yok; odak uzaklığı tahmin "
       "edilecek"));

// {0} which input, {1} the focal in px, {2} the factor, {3} the image width.
SS_MSG(sfm_initial_focal,
    EN("Initial focal length for {0}: {1} px ({2} x {3} px wide)"),
    JA("{0} の初期焦点距離: {1} px（横幅 {3} px の {2} 倍）"),
    ZH_HANS("{0} 的初始焦距：{1} px（宽 {3} px 的 {2} 倍）"),
    ZH_HANT("{0} 的初始焦距：{1} px（寬 {3} px 的 {2} 倍）"),
    KO("{0}의 초기 초점 거리: {1} px(가로 {3} px의 {2}배)"),
    DE("Anfangsbrennweite für {0}: {1} px ({2} x {3} px Breite)"),
    FR("Focale initiale pour {0} : {1} px ({2} x {3} px de large)"),
    ES("Distancia focal inicial de {0}: {1} px ({2} x {3} px de ancho)"),
    PT("Distância focal inicial de {0}: {1} px ({2} x {3} px de largura)"),
    IT("Focale iniziale per {0}: {1} px ({2} x {3} px di larghezza)"),
    NL("Beginbrandpuntsafstand voor {0}: {1} px ({2} x {3} px breed)"),
    RU("Начальное фокусное расстояние для {0}: {1} px ({2} x {3} px ширины)"),
    TR("{0} için başlangıç odak uzaklığı: {1} px ({2} x {3} px genişlik)"));

SS_MSG(sfm_the_capture,
    EN("the capture"),
    JA("この撮影"),
    ZH_HANS("本次拍摄"),
    ZH_HANT("本次拍攝"),
    KO("이번 촬영"),
    DE("die Aufnahme"),
    FR("la prise de vue"),
    ES("la captura"),
    PT("a captura"),
    IT("la ripresa"),
    NL("de opname"),
    RU("эта съёмка"),
    TR("bu çekim"));

SS_MSG(sfm_resuming,
    EN("Resuming the previous run in {0}"),
    JA("{0} にある前回の実行を再開します"),
    ZH_HANS("从 {0} 中上一次的运行继续"),
    ZH_HANT("從 {0} 中上一次的執行繼續"),
    KO("{0}에 있는 지난번 실행을 이어서 진행합니다"),
    DE("Der vorherige Lauf in {0} wird fortgesetzt"),
    FR("Reprise de l'exécution précédente dans {0}"),
    ES("Se reanuda la ejecución anterior en {0}"),
    PT("Retomando a execução anterior em {0}"),
    IT("Si riprende l'esecuzione precedente in {0}"),
    NL("De vorige uitvoering in {0} wordt hervat"),
    RU("Продолжаем предыдущий запуск в {0}"),
    TR("{0} içindeki önceki çalıştırma sürdürülüyor"));

SS_MSG(sfm_will_overwrite,
    EN("Note: {0} already holds a reconstruction; this run writes over it."),
    JA("メモ: {0} にはすでに復元結果があります。今回の実行で上書きします。"),
    ZH_HANS("提示：{0} 中已有一份重建结果，本次运行会覆盖它。"),
    ZH_HANT("提示：{0} 中已有一份重建結果，本次執行會覆蓋它。"),
    KO("참고: {0}에 이미 복원 결과가 있습니다. 이번 실행이 덮어씁니다."),
    DE("Hinweis: {0} enthält bereits eine Rekonstruktion; dieser Lauf "
       "überschreibt sie."),
    FR("Note : {0} contient déjà une reconstruction ; cette exécution "
       "l'écrase."),
    ES("Nota: {0} ya contiene una reconstrucción; esta ejecución la sobrescribe."),
    PT("Nota: {0} já contém uma reconstrução; esta execução vai sobrescrevê-la."),
    IT("Nota: {0} contiene già una ricostruzione; questa esecuzione la "
       "sovrascrive."),
    NL("Let op: {0} bevat al een reconstructie; deze uitvoering overschrijft "
       "die."),
    RU("Примечание: в {0} уже есть реконструкция; этот запуск её перезапишет."),
    TR("Not: {0} zaten bir yeniden kurulum içeriyor; bu çalıştırma onu "
       "üzerine yazacak."));

SS_MSG(one_camera_per_folder,
    EN("images/ holds one folder per camera: switching to one camera per "
       "folder"),
    JA("images/ はカメラごとに1つのフォルダーになっています。フォルダーごとに"
       "1台のカメラとして扱います"),
    ZH_HANS("images/ 里每台相机一个文件夹：改为每个文件夹一台相机"),
    ZH_HANT("images/ 裡每台相機一個資料夾：改為每個資料夾一台相機"),
    KO("images/ 안에 카메라별로 폴더가 하나씩 있습니다: 폴더당 카메라 하나로 "
       "바꿉니다"),
    DE("images/ enthält einen Ordner je Kamera: es wird auf eine Kamera pro "
       "Ordner umgestellt"),
    FR("images/ contient un dossier par appareil : passage à un appareil par "
       "dossier"),
    ES("images/ tiene una carpeta por cámara: se cambia a una cámara por "
       "carpeta"),
    PT("images/ tem uma pasta por câmera: mudando para uma câmera por pasta"),
    IT("images/ ha una cartella per fotocamera: si passa a una fotocamera per "
       "cartella"),
    NL("images/ bevat één map per camera: er wordt overgeschakeld op één camera "
       "per map"),
    RU("В images/ по одной папке на камеру: переходим к режиму «одна камера на "
       "папку»"),
    TR("images/ her kamera için bir klasör içeriyor: klasör başına bir kameraya "
       "geçiliyor"));

SS_MSG(sfm_resume_skip_recon,
    EN("Resume: a reconstruction already exists under sparse/; skipping "
       "(delete it to reconstruct again)"),
    JA("再開: sparse/ にすでに復元結果があります。とばします"
       "（やり直すには削除してください）"),
    ZH_HANS("继续：sparse/ 下已有重建结果，跳过（想重做请先删除）"),
    ZH_HANT("繼續：sparse/ 下已有重建結果，略過（想重做請先刪除）"),
    KO("이어서 진행: sparse/ 아래에 이미 복원 결과가 있어 건너뜁니다(다시 "
       "하려면 지우세요)"),
    DE("Fortsetzen: unter sparse/ liegt bereits eine Rekonstruktion; "
       "übersprungen (zum Neuaufbau löschen)"),
    FR("Reprise : une reconstruction existe déjà sous sparse/ ; étape ignorée "
       "(supprimez-la pour recommencer)"),
    ES("Reanudar: ya hay una reconstrucción en sparse/; se omite (bórrela para "
       "reconstruir de nuevo)"),
    PT("Retomar: já existe uma reconstrução em sparse/; pulando (apague-a para "
       "reconstruir de novo)"),
    IT("Ripresa: sotto sparse/ esiste già una ricostruzione; saltata (la "
       "cancelli per rifarla)"),
    NL("Hervatten: onder sparse/ staat al een reconstructie; overgeslagen "
       "(verwijder die om opnieuw te reconstrueren)"),
    RU("Продолжение: в sparse/ уже есть реконструкция; пропускаем (удалите её, "
       "чтобы построить заново)"),
    TR("Sürdürme: sparse/ altında zaten bir yeniden kurulum var; atlanıyor "
       "(yeniden kurmak için silin)"));

SS_MSG(sfm_partial,
    EN("Note: only part of the capture reconstructed. It will still train, but "
       "expect gaps."),
    JA("メモ: 撮影の一部しか復元できませんでした。学習はできますが、"
       "欠けが出ます。"),
    ZH_HANS("提示：只重建出了拍摄内容的一部分。仍然可以训练，但会有缺口。"),
    ZH_HANT("提示：只重建出了拍攝內容的一部分。仍然可以訓練，但會有缺口。"),
    KO("참고: 촬영분의 일부만 복원되었습니다. 학습은 되지만 빈 곳이 생깁니다."),
    DE("Hinweis: Nur ein Teil der Aufnahme wurde rekonstruiert. Das Training "
       "läuft trotzdem, aber mit Lücken."),
    FR("Note : seule une partie de la prise de vue a été reconstruite. "
       "L'entraînement fonctionnera, mais avec des trous."),
    ES("Nota: solo se reconstruyó parte de la captura. Se puede entrenar igual, "
       "pero habrá huecos."),
    PT("Nota: só parte da captura foi reconstruída. Ainda dá para treinar, mas "
       "haverá lacunas."),
    IT("Nota: è stata ricostruita solo una parte della ripresa. Si può "
       "addestrare lo stesso, ma con dei buchi."),
    NL("Let op: slechts een deel van de opname is gereconstrueerd. Trainen kan "
       "nog steeds, maar met gaten."),
    RU("Примечание: восстановлена лишь часть съёмки. Обучение пойдёт, но с "
       "пробелами."),
    TR("Not: çekimin yalnızca bir bölümü yeniden kuruldu. Yine de eğitilebilir, "
       "ama boşluklar olacak."));

SS_MSG(photos_referenced_in_place,
    EN("Note: the photos are referenced where they are. If you reopen this "
       "dataset later, set image_dir to {0} under the dataset options."),
    JA("メモ: 写真は元の場所を参照しています。あとでこのデータセットを開き直す"
       "ときは、データセット設定の image_dir に {0} を指定してください。"),
    ZH_HANS("提示：照片是就地引用的。以后重新打开这个数据集时，请在数据集选项里"
            "把 image_dir 设为 {0}。"),
    ZH_HANT("提示：照片是就地引用的。以後重新開啟這個資料集時，請在資料集選項裡"
            "把 image_dir 設為 {0}。"),
    KO("참고: 사진은 원래 자리를 참조합니다. 나중에 이 데이터셋을 다시 열 때는 "
       "데이터셋 옵션의 image_dir을 {0}(으)로 지정하세요."),
    DE("Hinweis: Die Fotos werden dort referenziert, wo sie liegen. Wenn Sie "
       "diesen Datensatz später wieder öffnen, setzen Sie image_dir in den "
       "Datensatzoptionen auf {0}."),
    FR("Note : les photos sont référencées là où elles sont. Si vous rouvrez ce "
       "jeu de données plus tard, réglez image_dir sur {0} dans ses options."),
    ES("Nota: las fotos se referencian donde están. Si vuelve a abrir este "
       "conjunto de datos más adelante, ponga image_dir en {0} dentro de sus "
       "opciones."),
    PT("Nota: as fotos são referenciadas onde estão. Se reabrir este conjunto "
       "de dados mais tarde, defina image_dir como {0} nas opções dele."),
    IT("Nota: le foto sono referenziate dove si trovano. Se riapre questo set "
       "di dati più avanti, imposti image_dir su {0} tra le sue opzioni."),
    NL("Let op: de foto's worden op hun eigen plek aangehaald. Als u deze "
       "dataset later opnieuw opent, zet image_dir dan op {0} bij de "
       "datasetopties."),
    RU("Примечание: фотографии используются там, где лежат. Если позже откроете "
       "этот набор данных снова, укажите в его параметрах image_dir = {0}."),
    TR("Not: fotoğraflar bulundukları yerden kullanılıyor. Bu veri kümesini "
       "sonra yeniden açarsanız, veri kümesi seçeneklerinde image_dir değerini "
       "{0} yapın."));

// {0} the port.
SS_MSG(web_viewer_at,
    EN("Web viewer at http://localhost:{0}/"),
    JA("ウェブビューア: http://localhost:{0}/"),
    ZH_HANS("网页查看器：http://localhost:{0}/"),
    ZH_HANT("網頁檢視器：http://localhost:{0}/"),
    KO("웹 뷰어: http://localhost:{0}/"),
    DE("Web-Betrachter unter http://localhost:{0}/"),
    FR("Visionneuse web sur http://localhost:{0}/"),
    ES("Visor web en http://localhost:{0}/"),
    PT("Visualizador web em http://localhost:{0}/"),
    IT("Visualizzatore web su http://localhost:{0}/"),
    NL("Webviewer op http://localhost:{0}/"),
    RU("Веб-просмотрщик: http://localhost:{0}/"),
    TR("Web görüntüleyici: http://localhost:{0}/"));

// ---- training ------------------------------------------------------------
//
// Only the lines that report progress to whoever is watching. The warnings
// about unported flags next to these in TrainerCore.cpp stay English on
// purpose: they name command-line flags and files under docs/notes/, so they
// are addressed to someone working on this program rather than using it.

// {0} cameras parsed, {1} after splitting panoramas, {2} seed points,
// {3} the frame scale.
SS_MSG(parsed_dataset,
    EN("Cameras: {0} ({1} after splitting), seed points: {2} "
       "(train_frame_scale={3})"),
    JA("カメラ: {0}（分割後 {1}）、初期点: {2}（train_frame_scale={3}）"),
    ZH_HANS("相机：{0}（拆分后 {1}），初始点：{2}（train_frame_scale={3}）"),
    ZH_HANT("相機：{0}（拆分後 {1}），初始點：{2}（train_frame_scale={3}）"),
    KO("카메라: {0}(분할 후 {1}), 초기 점: {2}(train_frame_scale={3})"),
    DE("Kameras: {0} ({1} nach dem Aufteilen), Startpunkte: {2} "
       "(train_frame_scale={3})"),
    FR("Caméras : {0} ({1} après découpage), points initiaux : {2} "
       "(train_frame_scale={3})"),
    ES("Cámaras: {0} ({1} tras dividir), puntos iniciales: {2} "
       "(train_frame_scale={3})"),
    PT("Câmeras: {0} ({1} após dividir), pontos iniciais: {2} "
       "(train_frame_scale={3})"),
    IT("Fotocamere: {0} ({1} dopo la divisione), punti iniziali: {2} "
       "(train_frame_scale={3})"),
    NL("Camera's: {0} ({1} na het splitsen), startpunten: {2} "
       "(train_frame_scale={3})"),
    RU("Камер: {0} ({1} после разделения), начальных точек: {2} "
       "(train_frame_scale={3})"),
    TR("Kamera: {0} (bölmeden sonra {1}), başlangıç noktası: {2} "
       "(train_frame_scale={3})"));

SS_MSG(output_directory,
    EN("Output directory: {0}"),
    JA("出力先フォルダー: {0}"),
    ZH_HANS("输出文件夹：{0}"),
    ZH_HANT("輸出資料夾：{0}"),
    KO("출력 폴더: {0}"),
    DE("Ausgabeordner: {0}"),
    FR("Dossier de sortie : {0}"),
    ES("Carpeta de salida: {0}"),
    PT("Pasta de saída: {0}"),
    IT("Cartella di uscita: {0}"),
    NL("Uitvoermap: {0}"),
    RU("Папка вывода: {0}"),
    TR("Çıktı klasörü: {0}"));

SS_MSG(ckpt_adapting,
    EN("Checkpoint layout differs from this run's; adapting on the host (no "
       "extra VRAM)..."),
    JA("チェックポイントの構成が今回の実行と違います。ホスト側で合わせます"
       "（VRAMは追加で使いません）…"),
    ZH_HANS("检查点的结构与本次运行不同，正在主机端做适配（不额外占用显存）…"),
    ZH_HANT("檢查點的結構與本次執行不同，正在主機端做調整（不額外佔用顯示記憶體）…"),
    KO("체크포인트 구성이 이번 실행과 달라 호스트에서 맞추는 중입니다"
       "(VRAM을 더 쓰지 않습니다)…"),
    DE("Der Aufbau des Checkpoints weicht von diesem Lauf ab; er wird auf dem "
       "Host angepasst (kein zusätzlicher VRAM) …"),
    FR("La structure du point de reprise diffère de celle de cette exécution ; "
       "adaptation côté hôte (sans VRAM supplémentaire)…"),
    ES("La estructura del punto de control difiere de la de esta ejecución; se "
       "adapta en el anfitrión (sin VRAM adicional)…"),
    PT("A estrutura do ponto de verificação difere da desta execução; "
       "adaptando no host (sem VRAM extra)…"),
    IT("La struttura del checkpoint è diversa da quella di questa esecuzione; "
       "adattamento sull'host (senza VRAM aggiuntiva)…"),
    NL("De opbouw van het controlepunt wijkt af van die van deze uitvoering; "
       "aanpassen op de host (zonder extra VRAM)…"),
    RU("Структура контрольной точки отличается от текущего запуска; "
       "подгоняем на хосте (без дополнительной видеопамяти)…"),
    TR("Denetim noktasının yapısı bu çalıştırmadan farklı; ana makinede "
       "uyarlanıyor (ek VRAM kullanmadan)…"));

SS_MSG(resumed_from,
    EN("Resumed from {0} at step {1}"),
    JA("{0} のステップ {1} から再開しました"),
    ZH_HANS("已从 {0} 的第 {1} 步继续"),
    ZH_HANT("已從 {0} 的第 {1} 步繼續"),
    KO("{0}의 {1}단계부터 이어서 진행합니다"),
    DE("Fortgesetzt ab {0}, Schritt {1}"),
    FR("Reprise depuis {0} à l'étape {1}"),
    ES("Reanudado desde {0} en el paso {1}"),
    PT("Retomado de {0} no passo {1}"),
    IT("Ripreso da {0} al passo {1}"),
    NL("Hervat vanaf {0} bij stap {1}"),
    RU("Продолжено с {0}, шаг {1}"),
    TR("{0} konumundan {1}. adımda sürdürüldü"));

SS_MSG(checkpoint_saved,
    EN("Checkpoint saved to: {0}"),
    JA("チェックポイントを保存しました: {0}"),
    ZH_HANS("检查点已保存到：{0}"),
    ZH_HANT("檢查點已儲存到：{0}"),
    KO("체크포인트를 저장했습니다: {0}"),
    DE("Checkpoint gespeichert unter: {0}"),
    FR("Point de reprise enregistré dans : {0}"),
    ES("Punto de control guardado en: {0}"),
    PT("Ponto de verificação salvo em: {0}"),
    IT("Checkpoint salvato in: {0}"),
    NL("Controlepunt opgeslagen in: {0}"),
    RU("Контрольная точка сохранена: {0}"),
    TR("Denetim noktası şuraya kaydedildi: {0}"));

SS_MSG(eval_split_empty,
    EN("Eval: the eval split is empty; nothing to score."),
    JA("評価: 評価用の分割が空です。採点するものがありません。"),
    ZH_HANS("评估：评估集为空，没有可评分的内容。"),
    ZH_HANT("評估：評估集為空，沒有可評分的內容。"),
    KO("평가: 평가용 분할이 비어 있어 점수를 낼 것이 없습니다."),
    DE("Auswertung: Der Auswertungsanteil ist leer; nichts zu bewerten."),
    FR("Évaluation : la partie d'évaluation est vide ; rien à noter."),
    ES("Evaluación: la partición de evaluación está vacía; nada que puntuar."),
    PT("Avaliação: a divisão de avaliação está vazia; nada a pontuar."),
    IT("Valutazione: la parte di valutazione è vuota; niente da valutare."),
    NL("Evaluatie: het evaluatiedeel is leeg; niets te scoren."),
    RU("Оценка: набор для оценки пуст; оценивать нечего."),
    TR("Değerlendirme: değerlendirme bölümü boş; puanlanacak bir şey yok."));

SS_MSG(eval_views,
    EN("Eval: views to score: {0}"),
    JA("評価: 採点する視点: {0}"),
    ZH_HANS("评估：待评分的视角：{0}"),
    ZH_HANT("評估：待評分的視角：{0}"),
    KO("평가: 점수를 낼 시점: {0}"),
    DE("Auswertung: zu bewertende Ansichten: {0}"),
    FR("Évaluation : vues à noter : {0}"),
    ES("Evaluación: vistas a puntuar: {0}"),
    PT("Avaliação: vistas a pontuar: {0}"),
    IT("Valutazione: viste da valutare: {0}"),
    NL("Evaluatie: te scoren aanzichten: {0}"),
    RU("Оценка: видов для оценки: {0}"),
    TR("Değerlendirme: puanlanacak görünüm: {0}"));

SS_MSG(eval_no_views,
    EN("Eval: no views were rendered."),
    JA("評価: 描画された視点がありませんでした。"),
    ZH_HANS("评估：没有渲染出任何视角。"),
    ZH_HANT("評估：沒有算繪出任何視角。"),
    KO("평가: 렌더링된 시점이 없습니다."),
    DE("Auswertung: Es wurden keine Ansichten gerendert."),
    FR("Évaluation : aucune vue n'a été rendue."),
    ES("Evaluación: no se renderizó ninguna vista."),
    PT("Avaliação: nenhuma vista foi renderizada."),
    IT("Valutazione: nessuna vista è stata renderizzata."),
    NL("Evaluatie: er zijn geen aanzichten gerenderd."),
    RU("Оценка: ни один вид не был отрисован."),
    TR("Değerlendirme: hiçbir görünüm işlenmedi."));

SS_MSG(eval_metrics_written,
    EN("Eval metrics written to {0}"),
    JA("評価指標を {0} に書き出しました"),
    ZH_HANS("评估指标已写入 {0}"),
    ZH_HANT("評估指標已寫入 {0}"),
    KO("평가 지표를 {0}에 기록했습니다"),
    DE("Auswertungskennzahlen geschrieben nach {0}"),
    FR("Mesures d'évaluation écrites dans {0}"),
    ES("Métricas de evaluación escritas en {0}"),
    PT("Métricas de avaliação escritas em {0}"),
    IT("Metriche di valutazione scritte in {0}"),
    NL("Evaluatiecijfers weggeschreven naar {0}"),
    RU("Показатели оценки записаны в {0}"),
    TR("Değerlendirme ölçütleri {0} dosyasına yazıldı"));

// ===========================================================================
// Failures -- what the screen says when a run stops
// ===========================================================================

SS_MSG(err_cancelled,
    EN("cancelled"),
    JA("中止しました"),
    ZH_HANS("已取消"),
    ZH_HANT("已取消"),
    KO("취소했습니다"),
    DE("abgebrochen"),
    FR("annulé"),
    ES("cancelado"),
    PT("cancelado"),
    IT("annullato"),
    NL("geannuleerd"),
    RU("отменено"),
    TR("iptal edildi"));

SS_MSG(err_unfinished_run,
    EN("the output folder already holds an unfinished run (extracted frames / "
       "features / masks); tick \"Resume previous run\" to reuse it, or pick "
       "another folder"),
    JA("出力先のフォルダーには終わっていない実行の残り（取り出したフレーム、"
       "特徴点、マスク）があります。「前回の続きから」を入れて使い回すか、"
       "別のフォルダーを選んでください"),
    ZH_HANS("输出文件夹里已有一次未完成运行的残留（提取的帧／特征点／蒙版）。"
            "勾选“接着上次继续”来复用，或换一个文件夹"),
    ZH_HANT("輸出資料夾裡已有一次未完成執行的殘留（擷取的影格／特徵點／遮罩）。"
            "勾選「接著上次繼續」來重用，或換一個資料夾"),
    KO("출력 폴더에 끝나지 않은 실행의 흔적(뽑아 둔 프레임·특징점·마스크)이 "
       "있습니다. \"이전 작업 이어서\"를 켜서 재사용하거나 다른 폴더를 "
       "고르세요"),
    DE("Der Ausgabeordner enthält bereits einen unfertigen Lauf (entnommene "
       "Einzelbilder / Merkmale / Masken); haken Sie \"Vorherigen Lauf "
       "fortsetzen\" an, um ihn weiterzuverwenden, oder wählen Sie einen "
       "anderen Ordner"),
    FR("Le dossier de sortie contient déjà une exécution inachevée (images "
       "extraites, points caractéristiques, masques) ; cochez « Reprendre "
       "l'exécution précédente » pour la réutiliser, ou choisissez un autre "
       "dossier"),
    ES("La carpeta de salida ya contiene una ejecución sin terminar "
       "(fotogramas extraídos, puntos característicos, máscaras); marque "
       "\"Reanudar la ejecución anterior\" para aprovecharla, o elija otra "
       "carpeta"),
    PT("A pasta de saída já contém uma execução inacabada (quadros extraídos, "
       "pontos característicos, máscaras); marque \"Retomar a execução "
       "anterior\" para reaproveitá-la, ou escolha outra pasta"),
    IT("La cartella di uscita contiene già un'esecuzione incompiuta "
       "(fotogrammi estratti, punti caratteristici, maschere); spunti "
       "\"Riprendi l'esecuzione precedente\" per riusarla, oppure scelga "
       "un'altra cartella"),
    NL("De uitvoermap bevat al een onafgemaakte uitvoering (uitgepakte beelden "
       "/ kenmerken / maskers); vink \"Vorige uitvoering hervatten\" aan om die "
       "te hergebruiken, of kies een andere map"),
    RU("В папке вывода уже есть незавершённый запуск (извлечённые кадры, особые "
       "точки, маски); отметьте «Продолжить предыдущий запуск», чтобы "
       "использовать его, или выберите другую папку"),
    TR("Çıktı klasöründe yarım kalmış bir çalıştırma var (çıkarılmış kareler / "
       "öznitelikler / maskeler); yeniden kullanmak için \"Önceki çalıştırmayı "
       "sürdür\" seçeneğini işaretleyin ya da başka bir klasör seçin"));

// {0} the program that would not start.
SS_MSG(err_spawn_recon,
    EN("could not start the reconstruction ({0})"),
    JA("復元処理を起動できませんでした（{0}）"),
    ZH_HANS("无法启动重建程序（{0}）"),
    ZH_HANT("無法啟動重建程式（{0}）"),
    KO("복원 프로그램을 시작하지 못했습니다({0})"),
    DE("Die Rekonstruktion konnte nicht gestartet werden ({0})"),
    FR("Impossible de lancer la reconstruction ({0})"),
    ES("No se pudo iniciar la reconstrucción ({0})"),
    PT("Não foi possível iniciar a reconstrução ({0})"),
    IT("Non è stato possibile avviare la ricostruzione ({0})"),
    NL("De reconstructie kon niet worden gestart ({0})"),
    RU("Не удалось запустить реконструкцию ({0})"),
    TR("Yeniden kurulum başlatılamadı ({0})"));

SS_MSG(err_recon_failed,
    EN("reconstruction failed (see the log). Common causes: too few "
       "overlapping images, not enough overlap between them, or the wrong "
       "camera model for the lens."),
    JA("復元に失敗しました（ログを見てください）。よくある原因は、重なりのある"
       "画像が少なすぎる、画像どうしの重なりが足りない、レンズに合わない"
       "カメラモデルを選んでいる、などです。"),
    ZH_HANS("重建失败（请看日志）。常见原因：有重叠的图像太少、图像之间重叠不够、"
            "或者选的相机模型与镜头不匹配。"),
    ZH_HANT("重建失敗（請看日誌）。常見原因：有重疊的影像太少、影像之間重疊不夠、"
            "或者選的相機模型與鏡頭不符。"),
    KO("복원에 실패했습니다(로그를 보세요). 흔한 원인은 겹치는 이미지가 너무 "
       "적음, 이미지끼리 겹침이 부족함, 렌즈에 맞지 않는 카메라 모델 선택 "
       "등입니다."),
    DE("Die Rekonstruktion ist fehlgeschlagen (siehe Protokoll). Häufige "
       "Ursachen: zu wenige überlappende Bilder, zu geringe Überlappung, oder "
       "das falsche Kameramodell für das Objektiv."),
    FR("La reconstruction a échoué (voir le journal). Causes fréquentes : trop "
       "peu d'images qui se recouvrent, recouvrement insuffisant entre elles, "
       "ou un modèle de caméra qui ne correspond pas à l'objectif."),
    ES("La reconstrucción falló (vea el registro). Causas frecuentes: "
       "demasiadas pocas imágenes solapadas, solape insuficiente entre ellas, "
       "o un modelo de cámara que no corresponde al objetivo."),
    PT("A reconstrução falhou (veja o registro). Causas comuns: imagens "
       "sobrepostas de menos, sobreposição insuficiente entre elas, ou o modelo "
       "de câmera errado para a lente."),
    IT("La ricostruzione è fallita (veda il registro). Cause frequenti: troppe "
       "poche immagini sovrapposte, sovrapposizione insufficiente fra loro, "
       "oppure il modello di fotocamera sbagliato per l'obiettivo."),
    NL("De reconstructie is mislukt (zie het logboek). Veelvoorkomende "
       "oorzaken: te weinig overlappende beelden, te weinig overlap ertussen, "
       "of het verkeerde cameramodel voor de lens."),
    RU("Реконструкция не удалась (смотрите журнал). Обычные причины: слишком "
       "мало перекрывающихся снимков, недостаточное перекрытие между ними или "
       "неподходящая модель камеры для объектива."),
    TR("Yeniden kurulum başarısız oldu (günlüğe bakın). Sık görülen nedenler: "
       "örtüşen görüntünün çok az olması, aralarındaki örtüşmenin yetersizliği "
       "veya objektife uymayan kamera modeli."));

SS_MSG(err_no_reconstruction,
    EN("no reconstruction was produced -- the images may not overlap enough. "
       "Try more photos around the subject, or a higher quality setting."),
    JA("復元結果ができませんでした。画像どうしの重なりが足りない可能性が"
       "あります。被写体のまわりで写真を増やすか、品質を上げてみてください。"),
    ZH_HANS("没有产生任何重建结果——图像之间的重叠可能不够。试着围着拍摄对象多拍"
            "一些照片，或者把质量调高。"),
    ZH_HANT("沒有產生任何重建結果——影像之間的重疊可能不夠。試著繞著拍攝對象多拍"
            "一些照片，或者把品質調高。"),
    KO("복원 결과가 나오지 않았습니다. 이미지끼리 겹침이 부족할 수 있습니다. "
       "피사체 주위에서 사진을 더 찍거나 품질을 높여 보세요."),
    DE("Es wurde keine Rekonstruktion erzeugt -- die Bilder überlappen "
       "vermutlich zu wenig. Nehmen Sie mehr Fotos rund um das Motiv auf, oder "
       "wählen Sie eine höhere Qualitätsstufe."),
    FR("Aucune reconstruction n'a été produite -- les images ne se recouvrent "
       "peut-être pas assez. Prenez plus de photos autour du sujet, ou montez "
       "le réglage de qualité."),
    ES("No se produjo ninguna reconstrucción: puede que las imágenes no se "
       "solapen lo bastante. Haga más fotos alrededor del motivo, o suba el "
       "ajuste de calidad."),
    PT("Nenhuma reconstrução foi produzida -- talvez as imagens não se "
       "sobreponham o bastante. Tire mais fotos ao redor do objeto, ou aumente "
       "a qualidade."),
    IT("Non è stata prodotta alcuna ricostruzione: forse le immagini non si "
       "sovrappongono abbastanza. Scatti più foto attorno al soggetto, oppure "
       "alzi il livello di qualità."),
    NL("Er is geen reconstructie ontstaan -- de beelden overlappen wellicht te "
       "weinig. Maak meer foto's rond het onderwerp, of zet de kwaliteit "
       "hoger."),
    RU("Реконструкция не получена — вероятно, снимки перекрываются слишком "
       "мало. Снимите больше кадров вокруг объекта или поднимите уровень "
       "качества."),
    TR("Hiçbir yeniden kurulum üretilmedi -- görüntüler yeterince örtüşmüyor "
       "olabilir. Nesnenin çevresinde daha çok fotoğraf çekin ya da kalite "
       "ayarını yükseltin."));

}  // namespace log
}  // namespace msg
}  // namespace i18n
}  // namespace spirula

#include "i18n/EndCatalog.h"
