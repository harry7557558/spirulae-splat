#pragma once

// The New Dataset screen, the mask preview, and the model licence prompt --
// everything between "here are my photos" and "here is a dataset to train on".
//
// This is where most of the beginner-facing copy lives, so it is also where a
// bad translation costs the most. The messages attached to something
// irreversible or expensive (overwriting a reconstruction, a 700 MB download,
// accepting a licence) are marked with a comment; those get human review in
// every language before shipping, machine translation is not acceptable for
// them.
//
// Same two rules as Gui.h: no sentence built from fragments, no
// plural-sensitive counting.

#include "i18n/BeginCatalog.h"

namespace spirula {
namespace i18n {
namespace msg {
namespace dataset {

// ===========================================================================
// Screen chrome
// ===========================================================================

SS_MSG(title_from_video,
    EN("Create Dataset from Video"),
    JA("動画からデータセットを作成"),
    ZH_HANS("从视频创建数据集"),
    ZH_HANT("從影片建立資料集"),
    KO("동영상으로 데이터셋 만들기"),
    DE("Datensatz aus Video erstellen"),
    FR("Créer un jeu de données à partir d'une vidéo"),
    ES("Crear un conjunto de datos a partir de un vídeo"),
    PT("Criar um conjunto de dados a partir de um vídeo"),
    IT("Crea un set di dati da un video"),
    NL("Dataset maken uit video"),
    RU("Создание набора данных из видео"),
    TR("Videodan veri kümesi oluştur"));

SS_MSG(title_from_photos,
    EN("Create Dataset from Photos"),
    JA("写真からデータセットを作成"),
    ZH_HANS("从照片创建数据集"),
    ZH_HANT("從相片建立資料集"),
    KO("사진으로 데이터셋 만들기"),
    DE("Datensatz aus Fotos erstellen"),
    FR("Créer un jeu de données à partir de photos"),
    ES("Crear un conjunto de datos a partir de fotos"),
    PT("Criar um conjunto de dados a partir de fotos"),
    IT("Crea un set di dati da fotografie"),
    NL("Dataset maken uit foto's"),
    RU("Создание набора данных из фотографий"),
    TR("Fotoğraflardan veri kümesi oluştur"));

SS_MSG(section_settings,
    EN("Settings"),      JA("設定"),          ZH_HANS("设置"),     ZH_HANT("設定"),
    KO("설정"),           DE("Einstellungen"), FR("Réglages"),    ES("Ajustes"),
    PT("Configurações"), IT("Impostazioni"), NL("Instellingen"), RU("Настройки"),
    TR("Ayarlar"));

SS_MSG(section_advanced,
    EN("Advanced"),      JA("詳細設定"),      ZH_HANS("高级"),     ZH_HANT("進階"),
    KO("고급"),           DE("Erweitert"),    FR("Avancé"),       ES("Avanzado"),
    PT("Avançado"),      IT("Avanzate"),     NL("Geavanceerd"),  RU("Дополнительно"),
    TR("Gelişmiş"));

SS_MSG(create_dataset,
    EN("Create Dataset"), JA("データセットを作成"), ZH_HANS("创建数据集"),
    ZH_HANT("建立資料集"), KO("데이터셋 만들기"), DE("Datensatz erstellen"),
    FR("Créer le jeu de données"), ES("Crear el conjunto de datos"),
    PT("Criar o conjunto de dados"), IT("Crea il set di dati"),
    NL("Dataset maken"), RU("Создать набор данных"), TR("Veri kümesini oluştur"));

SS_MSG(pick_input_first,
    EN("pick the input and the output folder first"),
    JA("先に入力と出力フォルダを選んでください"),
    ZH_HANS("请先选好输入和输出文件夹"),
    ZH_HANT("請先選好輸入和輸出資料夾"),
    KO("먼저 입력과 출력 폴더를 고르세요"),
    DE("zuerst die Eingabe und den Ausgabeordner wählen"),
    FR("choisissez d'abord l'entrée et le dossier de sortie"),
    ES("elija primero la entrada y la carpeta de salida"),
    PT("escolha primeiro a entrada e a pasta de saída"),
    IT("scelga prima l'ingresso e la cartella di destinazione"),
    NL("kies eerst de invoer en de uitvoermap"),
    RU("сначала выберите вход и папку результатов"),
    TR("önce girdiyi ve çıktı klasörünü seçin"));

SS_MSG(cancel,
    EN("Cancel"),        JA("キャンセル"),    ZH_HANS("取消"),     ZH_HANT("取消"),
    KO("취소"),           DE("Abbrechen"),    FR("Annuler"),      ES("Cancelar"),
    PT("Cancelar"),      IT("Annulla"),      NL("Annuleren"),    RU("Отмена"),
    TR("İptal"));

// {0} is the stage the runner is on, in English -- it is a diagnostic.
SS_MSG(stage_running,
    EN("{0} ..."),       JA("{0} …"),        ZH_HANS("{0} …"),   ZH_HANT("{0} …"),
    KO("{0} …"),          DE("{0} …"),        FR("{0}…"),         ES("{0}…"),
    PT("{0}…"),          IT("{0}…"),         NL("{0}…"),         RU("{0}…"),
    TR("{0}…"));

SS_MSG(open_in_trainer,
    EN("Open in Trainer"),
    JA("トレーナーで開く"),
    ZH_HANS("在训练器中打开"),
    ZH_HANT("在訓練器中開啟"),
    KO("트레이너에서 열기"),
    DE("Im Trainer öffnen"),
    FR("Ouvrir dans l'atelier"),
    ES("Abrir en el entrenador"),
    PT("Abrir no treinador"),
    IT("Apri nell'addestratore"),
    NL("Openen in de trainer"),
    RU("Открыть в тренажёре"),
    TR("Eğiticide aç"));

SS_MSG(done_at,
    EN("Done: {0}"),     JA("完了: {0}"),     ZH_HANS("完成：{0}"), ZH_HANT("完成：{0}"),
    KO("완료: {0}"),      DE("Fertig: {0}"),  FR("Terminé : {0}"), ES("Listo: {0}"),
    PT("Pronto: {0}"),   IT("Fatto: {0}"),   NL("Klaar: {0}"),   RU("Готово: {0}"),
    TR("Bitti: {0}"));

SS_MSG(failed,
    EN("Failed: {0}"),   JA("失敗: {0}"),     ZH_HANS("失败：{0}"), ZH_HANT("失敗：{0}"),
    KO("실패: {0}"),      DE("Fehlgeschlagen: {0}"), FR("Échec : {0}"),
    ES("Error: {0}"),    PT("Falhou: {0}"),  IT("Non riuscito: {0}"),
    NL("Mislukt: {0}"),  RU("Ошибка: {0}"),  TR("Başarısız: {0}"));

SS_MSG(cancelled,
    EN("Cancelled."),    JA("中止しました。"), ZH_HANS("已取消。"),  ZH_HANT("已取消。"),
    KO("취소했습니다."),  DE("Abgebrochen."), FR("Annulé."),      ES("Cancelado."),
    PT("Cancelado."),    IT("Annullato."),   NL("Geannuleerd."), RU("Отменено."),
    TR("İptal edildi."));

SS_MSG(partial_reconstruction,
    EN("Only part of the capture reconstructed -- it will train, but expect "
       "gaps. More overlap between shots, or a higher quality setting, "
       "usually fixes it."),
    JA("撮影の一部しか再構成できませんでした。学習は行えますが、欠けが出ます。"
       "撮影どうしの重なりを増やすか、品質設定を上げると直ることが多いです。"),
    ZH_HANS("只重建出了拍摄内容的一部分——仍然可以训练，但会有缺口。"
            "增加拍摄之间的重叠，或者提高质量设置，通常就能解决。"),
    ZH_HANT("只重建出了拍攝內容的一部分——仍然可以訓練，但會有缺口。"
            "增加拍攝之間的重疊，或者提高品質設定，通常就能解決。"),
    KO("촬영분의 일부만 재구성되었습니다. 학습은 되지만 빈 곳이 생깁니다. "
       "촬영끼리 더 많이 겹치게 하거나 품질 설정을 올리면 대개 해결됩니다."),
    DE("Nur ein Teil der Aufnahme wurde rekonstruiert -- trainieren lässt es "
       "sich, aber mit Lücken. Mehr Überlappung zwischen den Aufnahmen oder "
       "eine höhere Qualitätsstufe hilft meist."),
    FR("Seule une partie de la prise a été reconstruite : elle s'entraînera, "
       "mais avec des trous. Plus de recouvrement entre les prises, ou un "
       "réglage de qualité plus élevé, corrige généralement le problème."),
    ES("Solo se reconstruyó parte de la captura: se puede entrenar, pero "
       "habrá huecos. Más solapamiento entre tomas, o una calidad más alta, "
       "suele arreglarlo."),
    PT("Só parte da captura foi reconstruída: dá para treinar, mas haverá "
       "falhas. Mais sobreposição entre as tomadas, ou uma qualidade mais "
       "alta, costuma resolver."),
    IT("È stata ricostruita solo una parte della ripresa: si può addestrare, "
       "ma con dei vuoti. Più sovrapposizione tra gli scatti, o una qualità "
       "più alta, di solito risolve."),
    NL("Slechts een deel van de opname is gereconstrueerd -- trainen kan, "
       "maar met gaten. Meer overlap tussen de opnamen, of een hogere "
       "kwaliteitsinstelling, helpt meestal."),
    RU("Восстановилась лишь часть съёмки — обучать можно, но с пробелами. "
       "Обычно помогает большее перекрытие между кадрами или более высокая "
       "настройка качества."),
    TR("Çekimin yalnızca bir bölümü yeniden oluşturuldu -- eğitilebilir ama "
       "boşluklar olacak. Çekimler arasında daha çok örtüşme ya da daha "
       "yüksek bir kalite ayarı genelde sorunu çözer."));

// ===========================================================================
// Inputs and output folder
// ===========================================================================

SS_MSG(browse,
    EN("Browse..."),     JA("参照…"),         ZH_HANS("浏览…"),    ZH_HANT("瀏覽…"),
    KO("찾아보기…"),      DE("Durchsuchen …"), FR("Parcourir…"),  ES("Examinar…"),
    PT("Procurar…"),     IT("Sfoglia…"),     NL("Bladeren…"),    RU("Обзор…"),
    TR("Gözat…"));

SS_MSG(remove,
    EN("Remove"),        JA("削除"),          ZH_HANS("移除"),     ZH_HANT("移除"),
    KO("제거"),           DE("Entfernen"),    FR("Retirer"),      ES("Quitar"),
    PT("Remover"),       IT("Rimuovi"),      NL("Verwijderen"),  RU("Убрать"),
    TR("Kaldır"));

SS_MSG(add_video,
    EN("Add video..."),  JA("動画を追加…"),   ZH_HANS("添加视频…"), ZH_HANT("新增影片…"),
    KO("동영상 추가…"),   DE("Video hinzufügen …"), FR("Ajouter une vidéo…"),
    ES("Añadir un vídeo…"), PT("Adicionar um vídeo…"), IT("Aggiungi un video…"),
    NL("Video toevoegen…"), RU("Добавить видео…"), TR("Video ekle…"));

SS_MSG(add_video_help,
    EN("Add another clip to this dataset. Several videos reconstruct together "
       "as one scene: each gets its own folder of frames, and its own camera, "
       "so they may come from different lenses. Click several files in the "
       "dialog to take them all, or drop them onto this window."),
    JA("このデータセットにクリップをもう1本追加します。複数の動画は1つの"
       "シーンとしてまとめて再構成されます。それぞれにフレーム用のフォルダと"
       "カメラが割り当てられるので、レンズが違っていてもかまいません。"
       "ダイアログで複数のファイルを選ぶか、このウィンドウにドロップして"
       "ください。"),
    ZH_HANS("再往这个数据集里加一段视频。多段视频会作为同一个场景一起重建："
            "每段有自己的帧文件夹和自己的相机，所以它们可以来自不同镜头。"
            "在对话框里点选多个文件，或者把它们拖到这个窗口。"),
    ZH_HANT("再往這個資料集裡加一段影片。多段影片會作為同一個場景一起重建："
            "每段有自己的影格資料夾和自己的相機，所以它們可以來自不同鏡頭。"
            "在對話框裡點選多個檔案，或者把它們拖到這個視窗。"),
    KO("이 데이터셋에 클립을 하나 더 추가합니다. 여러 동영상은 하나의 장면으로 "
       "함께 재구성됩니다. 각각 프레임 폴더와 카메라를 따로 가지므로 서로 다른 "
       "렌즈여도 괜찮습니다. 대화 상자에서 여러 파일을 클릭하거나 이 창에 "
       "끌어다 놓으세요."),
    DE("Diesem Datensatz einen weiteren Clip hinzufügen. Mehrere Videos werden "
       "gemeinsam als eine Szene rekonstruiert: jedes bekommt einen eigenen "
       "Bilderordner und eine eigene Kamera, sie dürfen also von "
       "verschiedenen Objektiven stammen. Im Dialog mehrere Dateien anklicken "
       "oder sie in dieses Fenster ziehen."),
    FR("Ajouter une autre séquence à ce jeu de données. Plusieurs vidéos sont "
       "reconstruites ensemble comme une seule scène : chacune a son dossier "
       "d'images et sa propre caméra, elles peuvent donc venir d'objectifs "
       "différents. Cliquez plusieurs fichiers dans la boîte de dialogue, ou "
       "déposez-les sur cette fenêtre."),
    ES("Añadir otro clip a este conjunto de datos. Varios vídeos se "
       "reconstruyen juntos como una sola escena: cada uno tiene su carpeta "
       "de fotogramas y su propia cámara, así que pueden venir de objetivos "
       "distintos. Haga clic en varios archivos en el diálogo, o arrástrelos "
       "a esta ventana."),
    PT("Adicionar outro clipe a este conjunto de dados. Vários vídeos são "
       "reconstruídos juntos como uma única cena: cada um ganha sua pasta de "
       "quadros e sua própria câmera, então podem vir de lentes diferentes. "
       "Clique em vários arquivos na caixa de diálogo, ou arraste-os para "
       "esta janela."),
    IT("Aggiunge un'altra clip a questo set di dati. Più video vengono "
       "ricostruiti insieme come un'unica scena: ciascuno ha la sua cartella "
       "di fotogrammi e la sua fotocamera, quindi possono venire da obiettivi "
       "diversi. Clicchi più file nella finestra di dialogo, oppure li "
       "trascini su questa finestra."),
    NL("Nog een clip aan deze dataset toevoegen. Meerdere video's worden samen "
       "als één scène gereconstrueerd: elke krijgt een eigen map met beelden "
       "en een eigen camera, dus ze mogen van verschillende lenzen komen. "
       "Klik meerdere bestanden aan in het dialoogvenster, of sleep ze naar "
       "dit venster."),
    RU("Добавить в этот набор ещё один ролик. Несколько видео восстанавливаются "
       "вместе как одна сцена: у каждого своя папка кадров и своя камера, так "
       "что объективы могут быть разными. Выберите в диалоге несколько файлов "
       "или перетащите их в это окно."),
    TR("Bu veri kümesine bir klip daha ekleyin. Birden çok video tek bir sahne "
       "olarak birlikte yeniden oluşturulur: her birinin kendi kare klasörü "
       "ve kendi kamerası olur, dolayısıyla farklı objektiflerden gelebilirler. "
       "İletişim kutusunda birkaç dosyaya tıklayın ya da onları bu pencereye "
       "bırakın."));

SS_MSG(add_photos,
    EN("Add photos..."), JA("写真を追加…"),   ZH_HANS("添加照片…"), ZH_HANT("新增相片…"),
    KO("사진 추가…"),     DE("Fotos hinzufügen …"), FR("Ajouter des photos…"),
    ES("Añadir fotos…"), PT("Adicionar fotos…"), IT("Aggiungi fotografie…"),
    NL("Foto's toevoegen…"), RU("Добавить фотографии…"), TR("Fotoğraf ekle…"));

SS_MSG(add_photos_help,
    EN("Add a folder of photos. On its own it is read where it is; alongside "
       "another input its images are linked into the dataset, because the "
       "reconstruction reads one folder tree."),
    JA("写真のフォルダを追加します。単独ならその場所のまま読み込みます。"
       "ほかの入力と一緒の場合は、再構成が1つのフォルダツリーを読むため、"
       "画像はデータセットにリンクされます。"),
    ZH_HANS("添加一个照片文件夹。只有它时就地读取；和其他输入放在一起时，"
            "它的图像会被链接进数据集，因为重建只读取一棵文件夹树。"),
    ZH_HANT("新增一個相片資料夾。只有它時就地讀取；和其他輸入放在一起時，"
            "它的影像會被連結進資料集，因為重建只讀取一棵資料夾樹。"),
    KO("사진 폴더를 추가합니다. 하나뿐이면 있는 자리에서 그대로 읽고, 다른 "
       "입력과 함께라면 재구성이 하나의 폴더 트리만 읽기 때문에 이미지가 "
       "데이터셋으로 링크됩니다."),
    DE("Einen Fotoordner hinzufügen. Allein wird er dort gelesen, wo er liegt; "
       "neben einer anderen Eingabe werden seine Bilder in den Datensatz "
       "verlinkt, weil die Rekonstruktion einen einzigen Ordnerbaum liest."),
    FR("Ajouter un dossier de photos. Seul, il est lu là où il se trouve ; à "
       "côté d'une autre entrée, ses images sont liées dans le jeu de "
       "données, car la reconstruction ne lit qu'une seule arborescence."),
    ES("Añadir una carpeta de fotos. Sola, se lee donde está; junto a otra "
       "entrada, sus imágenes se enlazan al conjunto de datos, porque la "
       "reconstrucción lee un único árbol de carpetas."),
    PT("Adicionar uma pasta de fotos. Sozinha, ela é lida onde está; ao lado "
       "de outra entrada, suas imagens são vinculadas ao conjunto de dados, "
       "porque a reconstrução lê uma única árvore de pastas."),
    IT("Aggiunge una cartella di fotografie. Da sola viene letta dov'è; "
       "insieme a un altro ingresso le sue immagini vengono collegate nel set "
       "di dati, perché la ricostruzione legge un solo albero di cartelle."),
    NL("Een fotomap toevoegen. Alleen wordt die gelezen waar hij staat; naast "
       "een andere invoer worden de beelden in de dataset gelinkt, omdat de "
       "reconstructie één mappenboom leest."),
    RU("Добавить папку с фотографиями. В одиночку она читается там, где лежит; "
       "рядом с другим входом её снимки связываются с набором данных, потому "
       "что реконструкция читает одно дерево папок."),
    TR("Bir fotoğraf klasörü ekleyin. Tek başınaysa bulunduğu yerde okunur; "
       "başka bir girdiyle birlikteyse görüntüleri veri kümesine bağlanır, "
       "çünkü yeniden oluşturma tek bir klasör ağacı okur."));

SS_MSG(no_input_yet,
    EN("no input picked yet"),
    JA("入力がまだ選ばれていません"),
    ZH_HANS("还没有选择输入"),
    ZH_HANT("還沒有選擇輸入"),
    KO("아직 입력을 고르지 않았습니다"),
    DE("noch keine Eingabe gewählt"),
    FR("aucune entrée choisie pour l'instant"),
    ES("aún no se ha elegido ninguna entrada"),
    PT("nenhuma entrada escolhida ainda"),
    IT("nessun ingresso scelto finora"),
    NL("nog geen invoer gekozen"),
    RU("вход ещё не выбран"),
    TR("henüz girdi seçilmedi"));

SS_MSG(kind_video_file,
    EN("video file"),    JA("動画ファイル"),   ZH_HANS("视频文件"),  ZH_HANT("影片檔"),
    KO("동영상 파일"),    DE("Videodatei"),   FR("fichier vidéo"), ES("archivo de vídeo"),
    PT("arquivo de vídeo"), IT("file video"), NL("videobestand"),
    RU("видеофайл"),     TR("video dosyası"));

SS_MSG(kind_photo_folder,
    EN("photo folder"),  JA("写真フォルダ"),   ZH_HANS("照片文件夹"), ZH_HANT("相片資料夾"),
    KO("사진 폴더"),      DE("Fotoordner"),   FR("dossier de photos"),
    ES("carpeta de fotos"), PT("pasta de fotos"), IT("cartella di fotografie"),
    NL("fotomap"),       RU("папка с фотографиями"), TR("fotoğraf klasörü"));

SS_MSG(kind_video_file_masks,
    EN("video file + masks"),
    JA("動画ファイル＋マスク"),
    ZH_HANS("视频文件 + 蒙版"),
    ZH_HANT("影片檔 + 遮罩"),
    KO("동영상 파일 + 마스크"),
    DE("Videodatei + Masken"),
    FR("fichier vidéo + masques"),
    ES("archivo de vídeo + máscaras"),
    PT("arquivo de vídeo + máscaras"),
    IT("file video + maschere"),
    NL("videobestand + maskers"),
    RU("видеофайл + маски"),
    TR("video dosyası + maskeler"));

SS_MSG(kind_photo_folder_masks,
    EN("photo folder + masks"),
    JA("写真フォルダ＋マスク"),
    ZH_HANS("照片文件夹 + 蒙版"),
    ZH_HANT("相片資料夾 + 遮罩"),
    KO("사진 폴더 + 마스크"),
    DE("Fotoordner + Masken"),
    FR("dossier de photos + masques"),
    ES("carpeta de fotos + máscaras"),
    PT("pasta de fotos + máscaras"),
    IT("cartella di fotografie + maschere"),
    NL("fotomap + maskers"),
    RU("папка с фотографиями + маски"),
    TR("fotoğraf klasörü + maskeler"));

// {0} is a folder name under images/.
SS_MSG(row_video_to,
    EN("video -> images/{0}"),       JA("動画 → images/{0}"),
    ZH_HANS("视频 → images/{0}"),     ZH_HANT("影片 → images/{0}"),
    KO("동영상 → images/{0}"),        DE("Video -> images/{0}"),
    FR("vidéo -> images/{0}"),       ES("vídeo -> images/{0}"),
    PT("vídeo -> images/{0}"),       IT("video -> images/{0}"),
    NL("video -> images/{0}"),       RU("видео -> images/{0}"),
    TR("video -> images/{0}"));

SS_MSG(row_photos_to,
    EN("photos -> images/{0}"),      JA("写真 → images/{0}"),
    ZH_HANS("照片 → images/{0}"),     ZH_HANT("相片 → images/{0}"),
    KO("사진 → images/{0}"),          DE("Fotos -> images/{0}"),
    FR("photos -> images/{0}"),      ES("fotos -> images/{0}"),
    PT("fotos -> images/{0}"),       IT("fotografie -> images/{0}"),
    NL("foto's -> images/{0}"),      RU("фотографии -> images/{0}"),
    TR("fotoğraflar -> images/{0}"));

SS_MSG(row_video_masks_to,
    EN("video + masks -> images/{0}"),   JA("動画＋マスク → images/{0}"),
    ZH_HANS("视频 + 蒙版 → images/{0}"),  ZH_HANT("影片 + 遮罩 → images/{0}"),
    KO("동영상 + 마스크 → images/{0}"),   DE("Video + Masken -> images/{0}"),
    FR("vidéo + masques -> images/{0}"), ES("vídeo + máscaras -> images/{0}"),
    PT("vídeo + máscaras -> images/{0}"), IT("video + maschere -> images/{0}"),
    NL("video + maskers -> images/{0}"), RU("видео + маски -> images/{0}"),
    TR("video + maskeler -> images/{0}"));

SS_MSG(row_photos_masks_to,
    EN("photos + masks -> images/{0}"),  JA("写真＋マスク → images/{0}"),
    ZH_HANS("照片 + 蒙版 → images/{0}"),  ZH_HANT("相片 + 遮罩 → images/{0}"),
    KO("사진 + 마스크 → images/{0}"),     DE("Fotos + Masken -> images/{0}"),
    FR("photos + masques -> images/{0}"), ES("fotos + máscaras -> images/{0}"),
    PT("fotos + máscaras -> images/{0}"), IT("fotografie + maschere -> images/{0}"),
    NL("foto's + maskers -> images/{0}"), RU("фотографии + маски -> images/{0}"),
    TR("fotoğraflar + maskeler -> images/{0}"));

SS_MSG(existing_masks_tooltip,
    EN("Masks already made for these images:\n{0}\n\nThey are used as they "
       "are -- nothing is segmented for this input."),
    JA("これらの画像用にすでに用意されているマスクです:\n{0}\n\n"
       "そのまま使われ、この入力に対してセグメンテーションは行いません。"),
    ZH_HANS("已经为这些图像准备好的蒙版：\n{0}\n\n它们会被原样使用——"
            "不会对这个输入再做分割。"),
    ZH_HANT("已經為這些影像準備好的遮罩：\n{0}\n\n它們會被原樣使用——"
            "不會對這個輸入再做分割。"),
    KO("이 이미지들에 이미 만들어진 마스크입니다:\n{0}\n\n그대로 쓰이며, 이 "
       "입력에 대해서는 분할을 하지 않습니다."),
    DE("Für diese Bilder liegen bereits Masken vor:\n{0}\n\nSie werden "
       "unverändert genutzt -- für diese Eingabe wird nichts segmentiert."),
    FR("Des masques existent déjà pour ces images :\n{0}\n\nIls sont utilisés "
       "tels quels : rien n'est segmenté pour cette entrée."),
    ES("Ya hay máscaras hechas para estas imágenes:\n{0}\n\nSe usan tal cual: "
       "no se segmenta nada para esta entrada."),
    PT("Já existem máscaras feitas para estas imagens:\n{0}\n\nElas são usadas "
       "como estão -- nada é segmentado para esta entrada."),
    IT("Per queste immagini esistono già delle maschere:\n{0}\n\nVengono usate "
       "così come sono: per questo ingresso non si segmenta nulla."),
    NL("Er zijn al maskers voor deze beelden:\n{0}\n\nZe worden gebruikt zoals "
       "ze zijn -- voor deze invoer wordt niets gesegmenteerd."),
    RU("Для этих изображений уже есть маски:\n{0}\n\nОни используются как "
       "есть — для этого входа ничего не сегментируется."),
    TR("Bu görüntüler için hazır maskeler var:\n{0}\n\nOlduğu gibi "
       "kullanılırlar -- bu girdi için bir şey bölütlenmez."));

SS_MSG(output_folder,
    EN("output dataset folder"),
    JA("出力先のデータセットフォルダ"),
    ZH_HANS("输出数据集文件夹"),
    ZH_HANT("輸出資料集資料夾"),
    KO("출력 데이터셋 폴더"),
    DE("Ausgabeordner des Datensatzes"),
    FR("dossier du jeu de données de sortie"),
    ES("carpeta del conjunto de datos de salida"),
    PT("pasta do conjunto de dados de saída"),
    IT("cartella del set di dati in uscita"),
    NL("uitvoermap van de dataset"),
    RU("папка выходного набора данных"),
    TR("çıktı veri kümesi klasörü"));

SS_MSG(resume_previous,
    EN("Resume previous run"),
    JA("前回の実行を再開する"),
    ZH_HANS("继续上次的运行"),
    ZH_HANT("繼續上次的執行"),
    KO("이전 실행 이어서 하기"),
    DE("Vorherigen Lauf fortsetzen"),
    FR("Reprendre l'exécution précédente"),
    ES("Reanudar la ejecución anterior"),
    PT("Retomar a execução anterior"),
    IT("Riprendi l'esecuzione precedente"),
    NL("Vorige run hervatten"),
    RU("Продолжить прошлый запуск"),
    TR("Önceki çalıştırmayı sürdür"));

SS_MSG(resume_previous_help,
    EN("This folder holds a previous (possibly interrupted) run. Checked, the "
       "finished parts are reused -- extracted frames, masks, features and "
       "matches -- and only what is missing runs. Unchecked, a folder with "
       "none of that is required; nothing is ever deleted automatically."),
    JA("このフォルダには前回の（途中で止まったかもしれない）実行が残っています。"
       "チェックすると、切り出したフレーム、マスク、特徴、マッチなど完了済みの"
       "部分を再利用し、足りないところだけを実行します。外した場合は、それらが"
       "何もないフォルダが必要です。自動で削除されるものはありません。"),
    ZH_HANS("这个文件夹里有上一次（可能被中断的）运行。勾选后会复用已完成的部分"
            "——抽出的帧、蒙版、特征和匹配——只跑缺的那些。不勾选则需要一个"
            "不含这些内容的文件夹；任何东西都不会被自动删除。"),
    ZH_HANT("這個資料夾裡有上一次（可能被中斷的）執行。勾選後會重用已完成的部分"
            "——抽出的影格、遮罩、特徵和配對——只跑缺的那些。不勾選則需要一個"
            "不含這些內容的資料夾；任何東西都不會被自動刪除。"),
    KO("이 폴더에는 지난번(중간에 멈췄을 수도 있는) 실행이 남아 있습니다. "
       "체크하면 추출한 프레임, 마스크, 특징점, 매칭 등 끝난 부분을 재사용하고 "
       "빠진 것만 실행합니다. 체크를 해제하면 그런 것이 하나도 없는 폴더가 "
       "필요합니다. 무엇도 자동으로 지워지지 않습니다."),
    DE("In diesem Ordner liegt ein früherer (womöglich abgebrochener) Lauf. "
       "Angehakt werden die fertigen Teile weiterverwendet -- extrahierte "
       "Bilder, Masken, Merkmale und Zuordnungen -- und nur das Fehlende läuft. "
       "Nicht angehakt wird ein Ordner ohne all das verlangt; gelöscht wird "
       "nie etwas von selbst."),
    FR("Ce dossier contient une exécution précédente (peut-être interrompue). "
       "Coché, les parties terminées sont réutilisées -- images extraites, "
       "masques, points caractéristiques et appariements -- et seul ce qui "
       "manque est calculé. Décoché, un dossier vierge est exigé ; rien n'est "
       "jamais supprimé automatiquement."),
    ES("Esta carpeta contiene una ejecución anterior (quizá interrumpida). "
       "Marcado, se reutilizan las partes terminadas -- fotogramas extraídos, "
       "máscaras, características y correspondencias -- y solo se calcula lo "
       "que falta. Sin marcar, se exige una carpeta sin nada de eso; nunca se "
       "borra nada automáticamente."),
    PT("Esta pasta contém uma execução anterior (talvez interrompida). "
       "Marcado, as partes concluídas são reaproveitadas -- quadros "
       "extraídos, máscaras, características e correspondências -- e só o que "
       "falta é calculado. Desmarcado, exige-se uma pasta sem nada disso; "
       "nada é apagado automaticamente."),
    IT("Questa cartella contiene un'esecuzione precedente (forse "
       "interrotta). Selezionato, le parti già finite vengono riusate -- "
       "fotogrammi estratti, maschere, caratteristiche e corrispondenze -- e "
       "si calcola solo ciò che manca. Deselezionato, serve una cartella "
       "priva di tutto questo; nulla viene mai cancellato da solo."),
    NL("Deze map bevat een eerdere (mogelijk afgebroken) run. Aangevinkt "
       "worden de afgeronde delen hergebruikt -- uitgehaalde beelden, "
       "maskers, kenmerken en overeenkomsten -- en draait alleen wat "
       "ontbreekt. Uitgevinkt is een map zonder dat alles vereist; er wordt "
       "nooit iets automatisch verwijderd."),
    RU("В этой папке лежит прошлый (возможно, прерванный) запуск. С флажком "
       "готовые части используются повторно — извлечённые кадры, маски, "
       "особые точки и соответствия — и считается только недостающее. Без "
       "флажка требуется папка, где ничего этого нет; автоматически ничего "
       "никогда не удаляется."),
    TR("Bu klasörde önceki bir (belki yarıda kalmış) çalıştırma var. "
       "İşaretliyken bitmiş parçalar yeniden kullanılır -- çıkarılan kareler, "
       "maskeler, öznitelikler ve eşleşmeler -- ve yalnızca eksik olan "
       "çalışır. İşaretsizken bunların hiçbirini içermeyen bir klasör gerekir; "
       "hiçbir şey kendiliğinden silinmez."));

SS_MSG(unfinished_run_detected,
    EN("(unfinished run detected in this folder)"),
    JA("（このフォルダに未完了の実行があります）"),
    ZH_HANS("（在这个文件夹里发现了未完成的运行）"),
    ZH_HANT("（在這個資料夾裡發現了未完成的執行）"),
    KO("(이 폴더에서 끝나지 않은 실행을 찾았습니다)"),
    DE("(unvollständiger Lauf in diesem Ordner gefunden)"),
    FR("(exécution inachevée trouvée dans ce dossier)"),
    ES("(se encontró una ejecución sin terminar en esta carpeta)"),
    PT("(execução inacabada encontrada nesta pasta)"),
    IT("(trovata un'esecuzione incompleta in questa cartella)"),
    NL("(onafgemaakte run in deze map gevonden)"),
    RU("(в этой папке найден незавершённый запуск)"),
    TR("(bu klasörde yarım kalmış bir çalıştırma bulundu)"));

// IRREVERSIBLE -- human review in every language.
SS_MSG(would_overwrite_model,
    EN("This folder already holds a reconstruction (sparse/). Creating the "
       "dataset writes a new one over it -- point the output somewhere else "
       "to keep the old one."),
    JA("このフォルダにはすでに再構成結果（sparse/）があります。データセットを"
       "作成すると新しい結果で上書きされます。古いものを残したい場合は、"
       "出力先を別の場所にしてください。"),
    ZH_HANS("这个文件夹里已经有一份重建结果（sparse/）。创建数据集会用新的覆盖"
            "掉它——想保留旧的，请把输出指到别处。"),
    ZH_HANT("這個資料夾裡已經有一份重建結果（sparse/）。建立資料集會用新的覆蓋"
            "掉它——想保留舊的，請把輸出指到別處。"),
    KO("이 폴더에는 이미 재구성 결과(sparse/)가 있습니다. 데이터셋을 만들면 그 "
       "위에 새로 덮어씁니다. 예전 것을 남기려면 출력을 다른 곳으로 지정하세요."),
    DE("In diesem Ordner liegt bereits eine Rekonstruktion (sparse/). Beim "
       "Erstellen des Datensatzes wird sie überschrieben -- die Ausgabe "
       "woandershin legen, um die alte zu behalten."),
    FR("Ce dossier contient déjà une reconstruction (sparse/). Créer le jeu de "
       "données l'écrasera : pointez la sortie ailleurs pour conserver "
       "l'ancienne."),
    ES("Esta carpeta ya contiene una reconstrucción (sparse/). Crear el "
       "conjunto de datos la sobrescribirá: apunte la salida a otro sitio "
       "para conservar la anterior."),
    PT("Esta pasta já contém uma reconstrução (sparse/). Criar o conjunto de "
       "dados vai sobrescrevê-la -- aponte a saída para outro lugar para "
       "manter a antiga."),
    IT("Questa cartella contiene già una ricostruzione (sparse/). Creare il "
       "set di dati la sovrascriverà: indichi un'altra destinazione per "
       "conservare quella vecchia."),
    NL("Deze map bevat al een reconstructie (sparse/). Het maken van de "
       "dataset schrijft er een nieuwe overheen -- wijs de uitvoer ergens "
       "anders aan om de oude te behouden."),
    RU("В этой папке уже есть реконструкция (sparse/). Создание набора данных "
       "перезапишет её — укажите другую папку результатов, чтобы сохранить "
       "прежнюю."),
    TR("Bu klasörde zaten bir yeniden oluşturma var (sparse/). Veri kümesini "
       "oluşturmak onun üzerine yenisini yazar -- eskisini korumak için "
       "çıktıyı başka bir yere yönlendirin."));

// ===========================================================================
// Reconstruction engine
// ===========================================================================

SS_MSG(reconstruction,
    EN("Reconstruction:"), JA("再構成:"),   ZH_HANS("重建："),   ZH_HANT("重建："),
    KO("재구성:"),         DE("Rekonstruktion:"), FR("Reconstruction :"),
    ES("Reconstrucción:"), PT("Reconstrução:"), IT("Ricostruzione:"),
    NL("Reconstructie:"), RU("Реконструкция:"), TR("Yeniden oluşturma:"));

SS_MSG(engine_builtin,
    EN("Built-in (GPU)"), JA("内蔵（GPU）"),  ZH_HANS("内置（GPU）"), ZH_HANT("內建（GPU）"),
    KO("내장(GPU)"),      DE("Eingebaut (GPU)"), FR("Intégrée (GPU)"),
    ES("Integrada (GPU)"), PT("Integrada (GPU)"), IT("Integrata (GPU)"),
    NL("Ingebouwd (GPU)"), RU("Встроенная (GPU)"), TR("Yerleşik (GPU)"));

SS_MSG(engine_builtin_help,
    EN("This program's own structure-from-motion. Nothing to install, runs on "
       "the GPU."),
    JA("このプログラム自身の Structure from Motion です。インストール不要で、"
       "GPU 上で動きます。"),
    ZH_HANS("本程序自带的运动恢复结构。无需安装，在 GPU 上运行。"),
    ZH_HANT("本程式自帶的運動恢復結構。無需安裝，在 GPU 上執行。"),
    KO("이 프로그램에 들어 있는 Structure from Motion입니다. 따로 설치할 것 "
       "없이 GPU에서 돌아갑니다."),
    DE("Die eigene Structure-from-Motion dieses Programms. Nichts zu "
       "installieren, läuft auf der GPU."),
    FR("Le structure-from-motion intégré à ce programme. Rien à installer, "
       "tourne sur le GPU."),
    ES("El structure-from-motion propio de este programa. Nada que instalar, "
       "funciona en la GPU."),
    PT("O structure-from-motion do próprio programa. Nada a instalar, roda na "
       "GPU."),
    IT("Lo structure-from-motion integrato in questo programma. Niente da "
       "installare, gira sulla GPU."),
    NL("De eigen structure-from-motion van dit programma. Niets te "
       "installeren, draait op de GPU."),
    RU("Собственная реализация structure-from-motion в этой программе. Ставить "
       "ничего не нужно, работает на GPU."),
    TR("Bu programın kendi structure-from-motion'ı. Kurulacak bir şey yok, "
       "GPU üzerinde çalışır."));

SS_MSG(engine_colmap,
    EN("COLMAP (installed separately)"),
    JA("COLMAP（別途インストール）"),
    ZH_HANS("COLMAP（需另行安装）"),
    ZH_HANT("COLMAP（需另行安裝）"),
    KO("COLMAP(따로 설치)"),
    DE("COLMAP (separat installiert)"),
    FR("COLMAP (installé séparément)"),
    ES("COLMAP (instalado aparte)"),
    PT("COLMAP (instalado à parte)"),
    IT("COLMAP (installato a parte)"),
    NL("COLMAP (apart geïnstalleerd)"),
    RU("COLMAP (устанавливается отдельно)"),
    TR("COLMAP (ayrıca kurulur)"));

SS_MSG(engine_colmap_help,
    EN("Drive an external COLMAP instead. Worth having for comparison, and "
       "for the features it has that the built-in engine does not yet -- "
       "neural features and matching, in particular."),
    JA("代わりに外部の COLMAP を動かします。比較用として、また内蔵エンジンに"
       "まだない機能、とくにニューラル特徴とマッチングのために役立ちます。"),
    ZH_HANS("改为驱动外部的 COLMAP。留着它有用：可以做对比，也能用上内置引擎"
            "还没有的功能，尤其是神经网络特征和匹配。"),
    ZH_HANT("改為驅動外部的 COLMAP。留著它有用：可以做對比，也能用上內建引擎"
            "還沒有的功能，尤其是神經網路特徵和配對。"),
    KO("대신 외부 COLMAP을 실행합니다. 비교용으로도, 내장 엔진에 아직 없는 "
       "기능—특히 신경망 특징점과 매칭—을 쓰기 위해서도 쓸모가 있습니다."),
    DE("Stattdessen ein externes COLMAP steuern. Nützlich zum Vergleich und "
       "wegen der Funktionen, die die eingebaute Engine noch nicht hat -- "
       "neuronale Merkmale und Zuordnung vor allem."),
    FR("Piloter un COLMAP externe à la place. Utile pour comparer, et pour ce "
       "qu'il sait faire et que le moteur intégré ne fait pas encore -- les "
       "descripteurs et l'appariement neuronaux, notamment."),
    ES("Controlar un COLMAP externo en su lugar. Útil para comparar y por lo "
       "que sabe hacer y el motor integrado todavía no: en particular, "
       "características y emparejamiento neuronales."),
    PT("Controlar um COLMAP externo em vez disso. Útil para comparar e pelo "
       "que ele faz e o motor integrado ainda não: em especial, "
       "características e correspondência neurais."),
    IT("Pilotare invece un COLMAP esterno. Utile per confrontare e per ciò "
       "che sa fare e il motore integrato ancora no: in particolare "
       "caratteristiche e corrispondenze neurali."),
    NL("In plaats daarvan een externe COLMAP aansturen. Nuttig om te "
       "vergelijken en voor wat het kan en de ingebouwde engine nog niet -- "
       "neurale kenmerken en matching in het bijzonder."),
    RU("Вместо этого запускать внешний COLMAP. Пригодится для сравнения и "
       "ради того, чего встроенный движок пока не умеет, — прежде всего "
       "нейросетевых признаков и сопоставления."),
    TR("Onun yerine harici bir COLMAP çalıştırın. Karşılaştırma için ve "
       "yerleşik motorda henüz olmayan yetenekler için -- özellikle sinir "
       "ağı tabanlı öznitelikler ve eşleştirme için -- işe yarar."));

// ===========================================================================
// The basics
// ===========================================================================

SS_MSG(quality,
    EN("Quality"),       JA("品質"),          ZH_HANS("质量"),     ZH_HANT("品質"),
    KO("품질"),           DE("Qualität"),     FR("Qualité"),      ES("Calidad"),
    PT("Qualidade"),     IT("Qualità"),      NL("Kwaliteit"),    RU("Качество"),
    TR("Kalite"));

SS_MSG(quality_fast,
    EN("Fast"),          JA("高速"),          ZH_HANS("快速"),     ZH_HANT("快速"),
    KO("빠름"),           DE("Schnell"),      FR("Rapide"),       ES("Rápida"),
    PT("Rápida"),        IT("Veloce"),       NL("Snel"),         RU("Быстро"),
    TR("Hızlı"));

SS_MSG(quality_balanced,
    EN("Balanced"),      JA("バランス"),      ZH_HANS("均衡"),     ZH_HANT("均衡"),
    KO("균형"),           DE("Ausgewogen"),   FR("Équilibrée"),   ES("Equilibrada"),
    PT("Equilibrada"),   IT("Bilanciata"),   NL("Gebalanceerd"), RU("Сбалансированно"),
    TR("Dengeli"));

SS_MSG(quality_high_recommended,
    EN("High (recommended)"),
    JA("高（推奨）"),     ZH_HANS("高（推荐）"), ZH_HANT("高（建議）"),
    KO("높음(권장)"),     DE("Hoch (empfohlen)"), FR("Élevée (recommandée)"),
    ES("Alta (recomendada)"), PT("Alta (recomendada)"), IT("Alta (consigliata)"),
    NL("Hoog (aanbevolen)"), RU("Высокое (рекомендуется)"), TR("Yüksek (önerilen)"));

SS_MSG(quality_maximum,
    EN("Maximum"),       JA("最大"),          ZH_HANS("最高"),     ZH_HANT("最高"),
    KO("최대"),           DE("Maximal"),      FR("Maximale"),     ES("Máxima"),
    PT("Máxima"),        IT("Massima"),      NL("Maximaal"),     RU("Максимальное"),
    TR("En yüksek"));

SS_MSG(quality_high,
    EN("High quality"),  JA("高品質"),        ZH_HANS("高质量"),   ZH_HANT("高品質"),
    KO("고품질"),         DE("Hohe Qualität"), FR("Haute qualité"),
    ES("Alta calidad"),  PT("Alta qualidade"), IT("Alta qualità"),
    NL("Hoge kwaliteit"), RU("Высокое качество"), TR("Yüksek kalite"));

SS_MSG(quality_help_builtin,
    EN("Working resolution, how many features are found per image, and how "
       "many image pairs are compared. Higher finds more cameras in difficult "
       "scenes and takes longer -- roughly quadratically."),
    JA("作業解像度、1枚あたりに検出する特徴の数、比較する画像ペアの数をまとめて"
       "決めます。上げるほど難しいシーンでもカメラが見つかりますが、時間は"
       "おおよそ二乗で伸びます。"),
    ZH_HANS("工作分辨率、每张图像找多少特征，以及比较多少对图像。调高在困难场景"
            "里能找到更多相机，但耗时大致按平方增长。"),
    ZH_HANT("工作解析度、每張影像找多少特徵，以及比較多少對影像。調高在困難場景"
            "裡能找到更多相機，但耗時大致按平方增長。"),
    KO("작업 해상도, 이미지당 찾는 특징점 수, 비교하는 이미지 쌍의 수를 함께 "
       "정합니다. 높일수록 어려운 장면에서도 카메라를 더 많이 찾지만 시간은 "
       "대략 제곱으로 늘어납니다."),
    DE("Arbeitsauflösung, wie viele Merkmale je Bild gefunden werden und wie "
       "viele Bildpaare verglichen werden. Höher findet in schwierigen Szenen "
       "mehr Kameras und dauert länger -- ungefähr quadratisch."),
    FR("Résolution de travail, nombre de points caractéristiques par image et "
       "nombre de paires d'images comparées. Plus haut trouve plus de caméras "
       "dans les scènes difficiles et prend plus de temps -- à peu près au "
       "carré."),
    ES("Resolución de trabajo, cuántas características se buscan por imagen y "
       "cuántos pares de imágenes se comparan. Más alto encuentra más cámaras "
       "en escenas difíciles y tarda más, aproximadamente al cuadrado."),
    PT("Resolução de trabalho, quantas características são encontradas por "
       "imagem e quantos pares de imagens são comparados. Mais alto encontra "
       "mais câmeras em cenas difíceis e demora mais, aproximadamente ao "
       "quadrado."),
    IT("Risoluzione di lavoro, quante caratteristiche si cercano per immagine "
       "e quante coppie di immagini si confrontano. Più alto trova più "
       "fotocamere nelle scene difficili e richiede più tempo, all'incirca in "
       "modo quadratico."),
    NL("Werkresolutie, hoeveel kenmerken per beeld worden gevonden en hoeveel "
       "beeldparen worden vergeleken. Hoger vindt meer camera's in lastige "
       "scènes en duurt langer -- ongeveer kwadratisch."),
    RU("Рабочее разрешение, сколько особых точек ищется на снимке и сколько "
       "пар снимков сравнивается. Выше — больше найденных камер в сложных "
       "сценах и дольше, примерно квадратично."),
    TR("Çalışma çözünürlüğü, görüntü başına kaç öznitelik bulunacağı ve kaç "
       "görüntü çiftinin karşılaştırılacağı. Yüksek olan zor sahnelerde daha "
       "çok kamera bulur ve kabaca karesel olarak uzar."));

SS_MSG(quality_help_colmap,
    EN("Feature count used for matching (4k / 8k / 16k). Higher finds more "
       "cameras in difficult scenes but matching is O(n^2) in feature count."),
    JA("マッチングに使う特徴の数です（4k / 8k / 16k）。増やすほど難しいシーンでも"
       "カメラが見つかりますが、マッチングの計算量は特徴数の2乗です。"),
    ZH_HANS("用于匹配的特征数量（4k / 8k / 16k）。调高在困难场景里能找到更多相机，"
            "但匹配的复杂度是特征数的平方。"),
    ZH_HANT("用於配對的特徵數量（4k / 8k / 16k）。調高在困難場景裡能找到更多相機，"
            "但配對的複雜度是特徵數的平方。"),
    KO("매칭에 쓰는 특징점 수입니다(4k / 8k / 16k). 높이면 어려운 장면에서 "
       "카메라를 더 찾지만 매칭 비용은 특징점 수의 제곱입니다."),
    DE("Zahl der Merkmale für die Zuordnung (4k / 8k / 16k). Höher findet in "
       "schwierigen Szenen mehr Kameras, doch die Zuordnung wächst quadratisch "
       "mit der Merkmalszahl."),
    FR("Nombre de points caractéristiques utilisés pour l'appariement "
       "(4k / 8k / 16k). Plus haut trouve plus de caméras dans les scènes "
       "difficiles, mais l'appariement est en O(n²) du nombre de points."),
    ES("Número de características usadas para el emparejamiento "
       "(4k / 8k / 16k). Más alto encuentra más cámaras en escenas difíciles, "
       "pero el emparejamiento es O(n²) en el número de características."),
    PT("Número de características usadas na correspondência (4k / 8k / 16k). "
       "Mais alto encontra mais câmeras em cenas difíceis, mas a "
       "correspondência é O(n²) no número de características."),
    IT("Numero di caratteristiche usate per la corrispondenza "
       "(4k / 8k / 16k). Più alto trova più fotocamere nelle scene difficili, "
       "ma la corrispondenza è O(n²) nel numero di caratteristiche."),
    NL("Aantal kenmerken voor het matchen (4k / 8k / 16k). Hoger vindt meer "
       "camera's in lastige scènes, maar matchen is O(n²) in het aantal "
       "kenmerken."),
    RU("Число особых точек для сопоставления (4k / 8k / 16k). Больше — больше "
       "найденных камер в сложных сценах, но сопоставление растёт как квадрат "
       "числа точек."),
    TR("Eşleştirmede kullanılan öznitelik sayısı (4k / 8k / 16k). Yüksek olan "
       "zor sahnelerde daha çok kamera bulur ama eşleştirme öznitelik "
       "sayısına göre O(n²)'dir."));

SS_MSG(camera_lens,
    EN("Camera / lens"), JA("カメラ／レンズ"), ZH_HANS("相机／镜头"), ZH_HANT("相機／鏡頭"),
    KO("카메라 / 렌즈"), DE("Kamera / Objektiv"), FR("Caméra / objectif"),
    ES("Cámara / objetivo"), PT("Câmera / lente"), IT("Fotocamera / obiettivo"),
    NL("Camera / lens"), RU("Камера / объектив"), TR("Kamera / objektif"));

SS_MSG(camera_lens_help,
    EN("The lens distortion the reconstruction fits. OpenCV suits nearly "
       "every phone and camera. Pick a fisheye model for an action camera or "
       "a 360 rig -- a fisheye reconstructed as a normal lens comes out "
       "badly, and nothing detects that for you. Pinhole is only for images "
       "that are already undistorted."),
    JA("再構成が当てはめるレンズ歪みモデルです。OpenCV はほぼすべてのスマホ・"
       "カメラに合います。アクションカメラや360度リグでは魚眼モデルを選んで"
       "ください。魚眼を通常レンズとして再構成すると結果は悪くなり、それを"
       "自動で検出する仕組みはありません。Pinhole は歪み補正済みの画像専用です。"),
    ZH_HANS("重建要拟合的镜头畸变模型。OpenCV 几乎适用于所有手机和相机。"
            "运动相机或 360 相机组请选鱼眼模型——把鱼眼当普通镜头重建结果会很差，"
            "而且没有任何机制会替你发现。Pinhole 只用于已经去畸变的图像。"),
    ZH_HANT("重建要擬合的鏡頭變形模型。OpenCV 幾乎適用於所有手機和相機。"
            "運動相機或 360 相機組請選魚眼模型——把魚眼當普通鏡頭重建結果會很差，"
            "而且沒有任何機制會替你發現。Pinhole 只用於已經去變形的影像。"),
    KO("재구성이 맞출 렌즈 왜곡 모델입니다. OpenCV는 거의 모든 휴대폰과 카메라에 "
       "맞습니다. 액션캠이나 360 리그라면 어안 모델을 고르세요. 어안을 일반 "
       "렌즈로 재구성하면 결과가 나빠지는데, 그걸 대신 알아채 주는 장치는 "
       "없습니다. Pinhole은 이미 왜곡을 보정한 이미지에만 씁니다."),
    DE("Das Verzeichnungsmodell, das die Rekonstruktion anpasst. OpenCV passt "
       "zu fast jedem Telefon und jeder Kamera. Für eine Actionkamera oder ein "
       "360-Rig ein Fischaugenmodell wählen -- ein als Normalobjektiv "
       "rekonstruiertes Fischauge wird schlecht, und niemand merkt das für "
       "Sie. Pinhole ist nur für bereits entzerrte Bilder."),
    FR("Le modèle de distorsion que la reconstruction ajuste. OpenCV convient "
       "à presque tous les téléphones et appareils. Choisissez un modèle "
       "fisheye pour une caméra d'action ou un rig 360 : un fisheye "
       "reconstruit comme un objectif normal donne un mauvais résultat, et "
       "rien ne le détecte pour vous. Pinhole ne sert qu'aux images déjà "
       "corrigées."),
    ES("El modelo de distorsión que ajusta la reconstrucción. OpenCV vale "
       "para casi todos los teléfonos y cámaras. Elija un modelo de ojo de "
       "pez para una cámara de acción o un equipo 360: un ojo de pez "
       "reconstruido como objetivo normal sale mal, y nada lo detecta por "
       "usted. Pinhole solo sirve para imágenes ya corregidas."),
    PT("O modelo de distorção que a reconstrução ajusta. OpenCV serve para "
       "quase todo telefone e câmera. Escolha um modelo olho de peixe para "
       "uma câmera de ação ou um conjunto 360: um olho de peixe reconstruído "
       "como lente normal sai ruim, e nada detecta isso por você. Pinhole só "
       "serve para imagens já corrigidas."),
    IT("Il modello di distorsione che la ricostruzione adatta. OpenCV va bene "
       "per quasi ogni telefono e fotocamera. Per una action cam o un rig 360 "
       "scelga un modello fisheye: un fisheye ricostruito come obiettivo "
       "normale viene male, e nulla se ne accorge al posto suo. Pinhole serve "
       "solo per immagini già corrette."),
    NL("Het vervormingsmodel dat de reconstructie past. OpenCV past bij bijna "
       "elke telefoon en camera. Kies een fisheye-model voor een actiecamera "
       "of een 360-rig -- een fisheye die als gewone lens wordt "
       "gereconstrueerd komt er slecht uit, en niets merkt dat voor u op. "
       "Pinhole is alleen voor al ontvormde beelden."),
    RU("Модель искажений объектива, которую подгоняет реконструкция. OpenCV "
       "подходит почти любому телефону и фотоаппарату. Для экшн-камеры или "
       "360-риг выберите модель фишай: фишай, восстановленный как обычный "
       "объектив, выходит плохо, и заметить это за вас некому. Pinhole — "
       "только для уже исправленных изображений."),
    TR("Yeniden oluşturmanın uyduracağı objektif bozulma modeli. OpenCV neredeyse "
       "her telefona ve kameraya uyar. Aksiyon kamerası veya 360 düzeneği için "
       "balıkgözü modeli seçin -- normal objektif gibi yeniden oluşturulan bir "
       "balıkgözü kötü çıkar ve bunu sizin yerinize fark eden bir şey yoktur. "
       "Pinhole yalnızca bozulması giderilmiş görüntüler içindir."));

SS_MSG(colmap_one_lens_warning,
    EN("COLMAP fits this one lens model to every input. Switch to the "
       "built-in reconstruction to give each input its own."),
    JA("COLMAP はこの1つのレンズモデルをすべての入力に当てはめます。入力ごとに"
       "別のモデルを使うには、内蔵の再構成に切り替えてください。"),
    ZH_HANS("COLMAP 会把这一个镜头模型套用到所有输入上。想让每个输入各用各的，"
            "请切换到内置重建。"),
    ZH_HANT("COLMAP 會把這一個鏡頭模型套用到所有輸入上。想讓每個輸入各用各的，"
            "請切換到內建重建。"),
    KO("COLMAP은 이 렌즈 모델 하나를 모든 입력에 적용합니다. 입력마다 따로 "
       "쓰려면 내장 재구성으로 바꾸세요."),
    DE("COLMAP legt dieses eine Objektivmodell über jede Eingabe. Für ein "
       "eigenes Modell je Eingabe zur eingebauten Rekonstruktion wechseln."),
    FR("COLMAP applique ce seul modèle d'objectif à toutes les entrées. "
       "Passez à la reconstruction intégrée pour en donner un à chacune."),
    ES("COLMAP aplica este único modelo de objetivo a todas las entradas. "
       "Cambie a la reconstrucción integrada para dar uno propio a cada una."),
    PT("O COLMAP aplica este único modelo de lente a todas as entradas. Mude "
       "para a reconstrução integrada para dar um a cada uma."),
    IT("COLMAP applica questo unico modello di obiettivo a tutti gli "
       "ingressi. Passi alla ricostruzione integrata per darne uno a "
       "ciascuno."),
    NL("COLMAP past dit ene lensmodel op elke invoer toe. Schakel over op de "
       "ingebouwde reconstructie om elke invoer een eigen model te geven."),
    RU("COLMAP применяет одну эту модель объектива ко всем входам. Чтобы у "
       "каждого входа была своя, переключитесь на встроенную реконструкцию."),
    TR("COLMAP bu tek objektif modelini bütün girdilere uygular. Her girdiye "
       "kendi modelini vermek için yerleşik yeniden oluşturmaya geçin."));

SS_MSG(camera_lens_per_input,
    EN("Camera / lens per input"),
    JA("入力ごとのカメラ／レンズ"),
    ZH_HANS("每个输入的相机／镜头"),
    ZH_HANT("每個輸入的相機／鏡頭"),
    KO("입력별 카메라 / 렌즈"),
    DE("Kamera / Objektiv je Eingabe"),
    FR("Caméra / objectif par entrée"),
    ES("Cámara / objetivo por entrada"),
    PT("Câmera / lente por entrada"),
    IT("Fotocamera / obiettivo per ingresso"),
    NL("Camera / lens per invoer"),
    RU("Камера / объектив для каждого входа"),
    TR("Girdi başına kamera / objektif"));

SS_MSG(camera_lens_per_input_help,
    EN("Each input's images go into their own folder and get their own "
       "camera, so each can have its own lens model and starting focal "
       "length. A 360 file's two lens tracks share the row -- they are the "
       "same lens twice -- but are still solved as two cameras."),
    JA("入力ごとに画像は専用のフォルダに入り、カメラも別々になります。そのため"
       "レンズモデルと初期焦点距離を入力ごとに設定できます。360度ファイルの"
       "2つのレンズトラックは同じレンズが2つなので1行にまとまりますが、"
       "解かれるときは2台のカメラとして扱われます。"),
    ZH_HANS("每个输入的图像各进各的文件夹，也各有各的相机，所以镜头模型和初始"
            "焦距都可以分别设置。360 文件的两条镜头轨道共用一行——它们是同一"
            "只镜头的两份——但求解时仍算两台相机。"),
    ZH_HANT("每個輸入的影像各進各的資料夾，也各有各的相機，所以鏡頭模型和初始"
            "焦距都可以分別設定。360 檔案的兩條鏡頭軌道共用一列——它們是同一"
            "顆鏡頭的兩份——但求解時仍算兩台相機。"),
    KO("입력마다 이미지가 각자의 폴더에 들어가고 카메라도 따로 생기므로, 렌즈 "
       "모델과 시작 초점거리를 각각 정할 수 있습니다. 360 파일의 두 렌즈 트랙은 "
       "같은 렌즈가 둘이라 한 줄을 함께 쓰지만, 풀 때는 두 대의 카메라로 "
       "다룹니다."),
    DE("Die Bilder jeder Eingabe kommen in einen eigenen Ordner und bekommen "
       "eine eigene Kamera, also darf jede ihr eigenes Objektivmodell und ihre "
       "eigene Startbrennweite haben. Die zwei Objektivspuren einer "
       "360-Datei teilen sich die Zeile -- es ist zweimal dasselbe Objektiv "
       "-- werden aber trotzdem als zwei Kameras gelöst."),
    FR("Les images de chaque entrée vont dans leur propre dossier et "
       "reçoivent leur propre caméra ; chacune peut donc avoir son modèle "
       "d'objectif et sa focale de départ. Les deux pistes d'un fichier 360 "
       "partagent la ligne -- c'est deux fois le même objectif -- mais sont "
       "quand même résolues comme deux caméras."),
    ES("Las imágenes de cada entrada van a su propia carpeta y reciben su "
       "propia cámara, así que cada una puede tener su modelo de objetivo y "
       "su focal inicial. Las dos pistas de un archivo 360 comparten la fila "
       "-- es el mismo objetivo dos veces -- pero se resuelven igualmente "
       "como dos cámaras."),
    PT("As imagens de cada entrada vão para a própria pasta e ganham a "
       "própria câmera, então cada uma pode ter seu modelo de lente e sua "
       "distância focal inicial. As duas trilhas de um arquivo 360 dividem a "
       "linha -- é a mesma lente duas vezes -- mas mesmo assim são resolvidas "
       "como duas câmeras."),
    IT("Le immagini di ogni ingresso vanno in una cartella propria e ricevono "
       "una fotocamera propria, così ciascuna può avere il suo modello di "
       "obiettivo e la sua focale iniziale. Le due tracce di un file 360 "
       "condividono la riga -- è lo stesso obiettivo due volte -- ma vengono "
       "comunque risolte come due fotocamere."),
    NL("De beelden van elke invoer gaan in een eigen map en krijgen een eigen "
       "camera, dus elk mag zijn eigen lensmodel en beginbrandpuntsafstand "
       "hebben. De twee lenssporen van een 360-bestand delen de regel -- het "
       "is tweemaal dezelfde lens -- maar worden toch als twee camera's "
       "opgelost."),
    RU("Снимки каждого входа попадают в свою папку и получают свою камеру, так "
       "что у каждого может быть своя модель объектива и своё начальное "
       "фокусное расстояние. Две дорожки 360-файла делят одну строку — это "
       "один и тот же объектив дважды, — но решаются всё равно как две камеры."),
    TR("Her girdinin görüntüleri kendi klasörüne gider ve kendi kamerasını "
       "alır, dolayısıyla her biri kendi objektif modeline ve başlangıç odak "
       "uzaklığına sahip olabilir. Bir 360 dosyasının iki objektif izi aynı "
       "satırı paylaşır -- aynı objektifin iki kopyasıdır -- ama yine de iki "
       "kamera olarak çözülür."));

SS_MSG(focal_x_width,
    EN("x width"),       JA("×幅"),          ZH_HANS("× 宽度"),   ZH_HANT("× 寬度"),
    KO("× 너비"),         DE("× Breite"),     FR("× largeur"),    ES("× ancho"),
    PT("× largura"),     IT("× larghezza"),  NL("× breedte"),    RU("× ширина"),
    TR("× genişlik"));

SS_MSG(focal_x_width_help,
    EN("Starting focal length for this input, as a fraction of its image "
       "width (fx = fy = factor x width) -- the width is only known once the "
       "frames exist, which is why it is not in pixels here. 0 reads EXIF and "
       "falls back to a guess from the image size. Worth setting for a "
       "fisheye, where a bad guess can stop the reconstruction from starting "
       "at all; an Insta360 X5 is ~0.269, which .insv files are filled in "
       "with."),
    JA("この入力の初期焦点距離を、画像幅に対する割合で指定します"
       "（fx = fy = 係数 × 幅）。幅はフレームができて初めて分かるため、ここでは"
       "ピクセルで指定できません。0 なら EXIF を読み、なければ画像サイズから"
       "推定します。魚眼では推定が外れると再構成がそもそも始まらないことがある"
       "ので、指定する価値があります。Insta360 X5 はおよそ 0.269 で、.insv では"
       "自動で入ります。"),
    ZH_HANS("这个输入的初始焦距，按图像宽度的比例给出（fx = fy = 系数 × 宽度）。"
            "宽度要等帧生成后才知道，所以这里不用像素。填 0 会读 EXIF，读不到"
            "就按图像尺寸估计。鱼眼值得设一下：估错可能让重建根本起不来。"
            "Insta360 X5 约为 0.269，.insv 文件会自动填上。"),
    ZH_HANT("這個輸入的初始焦距，按影像寬度的比例給出（fx = fy = 係數 × 寬度）。"
            "寬度要等影格產生後才知道，所以這裡不用像素。填 0 會讀 EXIF，讀不到"
            "就按影像尺寸估計。魚眼值得設一下：估錯可能讓重建根本起不來。"
            "Insta360 X5 約為 0.269，.insv 檔會自動填上。"),
    KO("이 입력의 시작 초점거리를 이미지 너비에 대한 비율로 지정합니다"
       "(fx = fy = 계수 × 너비). 너비는 프레임이 생겨야 알 수 있어서 여기서는 "
       "픽셀로 지정하지 않습니다. 0이면 EXIF를 읽고, 없으면 이미지 크기로 "
       "추정합니다. 어안에서는 추정이 어긋나면 재구성이 아예 시작되지 않을 수 "
       "있어 설정할 값어치가 있습니다. Insta360 X5는 약 0.269이며 .insv 파일에는 "
       "자동으로 채워집니다."),
    DE("Startbrennweite für diese Eingabe, als Bruchteil ihrer Bildbreite "
       "(fx = fy = Faktor × Breite) -- die Breite steht erst fest, wenn die "
       "Bilder da sind, deshalb hier nicht in Pixeln. 0 liest EXIF und fällt "
       "auf eine Schätzung aus der Bildgröße zurück. Bei einem Fischauge "
       "lohnt es sich, denn eine schlechte Schätzung kann die Rekonstruktion "
       "ganz verhindern; eine Insta360 X5 liegt bei etwa 0,269, womit "
       ".insv-Dateien gefüllt werden."),
    FR("Focale de départ pour cette entrée, en fraction de la largeur d'image "
       "(fx = fy = facteur × largeur) -- la largeur n'est connue qu'une fois "
       "les images extraites, d'où l'absence de pixels ici. 0 lit l'EXIF et "
       "retombe sur une estimation d'après la taille d'image. Utile pour un "
       "fisheye, où une mauvaise estimation peut empêcher la reconstruction "
       "de démarrer ; une Insta360 X5 vaut environ 0,269, valeur inscrite "
       "d'office pour les fichiers .insv."),
    ES("Focal inicial de esta entrada, como fracción del ancho de imagen "
       "(fx = fy = factor × ancho): el ancho no se conoce hasta que existen "
       "los fotogramas, por eso aquí no va en píxeles. 0 lee el EXIF y "
       "recurre a una estimación por el tamaño de imagen. Vale la pena "
       "fijarla en un ojo de pez, donde una mala estimación puede impedir que "
       "la reconstrucción arranque; una Insta360 X5 ronda 0,269, valor que se "
       "rellena solo para archivos .insv."),
    PT("Distância focal inicial desta entrada, como fração da largura da "
       "imagem (fx = fy = fator × largura) -- a largura só é conhecida depois "
       "que os quadros existem, por isso aqui não é em pixels. 0 lê o EXIF e "
       "recorre a uma estimativa pelo tamanho da imagem. Vale definir num "
       "olho de peixe, onde um palpite ruim pode impedir a reconstrução de "
       "começar; uma Insta360 X5 fica em ~0,269, valor preenchido "
       "automaticamente para arquivos .insv."),
    IT("Focale iniziale per questo ingresso, come frazione della larghezza "
       "dell'immagine (fx = fy = fattore × larghezza): la larghezza si conosce "
       "solo quando i fotogrammi esistono, ecco perché qui non è in pixel. 0 "
       "legge l'EXIF e ripiega su una stima dalla dimensione dell'immagine. "
       "Conviene impostarla per un fisheye, dove una stima sbagliata può "
       "impedire del tutto l'avvio della ricostruzione; una Insta360 X5 sta "
       "attorno a 0,269, valore che i file .insv ricevono da soli."),
    NL("Beginbrandpuntsafstand voor deze invoer, als fractie van de "
       "beeldbreedte (fx = fy = factor × breedte) -- de breedte is pas bekend "
       "als de beelden er zijn, vandaar geen pixels hier. 0 leest EXIF en valt "
       "terug op een schatting uit de beeldgrootte. De moeite waard bij een "
       "fisheye, waar een slechte schatting de reconstructie helemaal kan "
       "blokkeren; een Insta360 X5 zit rond 0,269, waarmee .insv-bestanden "
       "worden ingevuld."),
    RU("Начальное фокусное расстояние для этого входа — как доля ширины "
       "изображения (fx = fy = коэффициент × ширина). Ширина известна только "
       "после появления кадров, поэтому здесь не пиксели. 0 читает EXIF, а при "
       "его отсутствии оценивает по размеру изображения. Для фишая задать "
       "стоит: плохая догадка может вовсе не дать реконструкции начаться; у "
       "Insta360 X5 это примерно 0,269, и для .insv значение подставляется "
       "само."),
    TR("Bu girdi için başlangıç odak uzaklığı, görüntü genişliğinin bir kesri "
       "olarak (fx = fy = katsayı × genişlik) -- genişlik ancak kareler "
       "oluştuğunda bilindiğinden burada piksel cinsinden verilmez. 0, EXIF'i "
       "okur ve bulamazsa görüntü boyutundan tahmin eder. Balıkgözünde "
       "ayarlamaya değer: kötü bir tahmin yeniden oluşturmanın hiç "
       "başlamamasına yol açabilir. Insta360 X5 için yaklaşık 0,269'dur ve "
       ".insv dosyalarına kendiliğinden yazılır."));

SS_MSG(camera_sharing,
    EN("Camera sharing"), JA("カメラの共有"),  ZH_HANS("相机共享"),  ZH_HANT("相機共用"),
    KO("카메라 공유"),    DE("Kamera teilen"), FR("Partage de caméra"),
    ES("Cámara compartida"), PT("Compartilhamento de câmera"),
    IT("Condivisione fotocamera"), NL("Camera delen"), RU("Общая камера"),
    TR("Kamera paylaşımı"));

SS_MSG(camera_sharing_one,
    EN("one shared camera"),
    JA("共有カメラ1台"), ZH_HANS("共用一台相机"), ZH_HANT("共用一台相機"),
    KO("공유 카메라 하나"), DE("eine gemeinsame Kamera"), FR("une caméra partagée"),
    ES("una cámara compartida"), PT("uma câmera compartilhada"),
    IT("una fotocamera condivisa"), NL("één gedeelde camera"),
    RU("одна общая камера"), TR("tek ortak kamera"));

SS_MSG(camera_sharing_folder,
    EN("one camera per folder"),
    JA("フォルダごとに1台"), ZH_HANS("每个文件夹一台相机"), ZH_HANT("每個資料夾一台相機"),
    KO("폴더마다 카메라 하나"), DE("eine Kamera je Ordner"),
    FR("une caméra par dossier"), ES("una cámara por carpeta"),
    PT("uma câmera por pasta"), IT("una fotocamera per cartella"),
    NL("één camera per map"), RU("по камере на папку"),
    TR("klasör başına bir kamera"));

SS_MSG(camera_sharing_image,
    EN("one camera per image"),
    JA("画像ごとに1台"), ZH_HANS("每张图像一台相机"), ZH_HANT("每張影像一台相機"),
    KO("이미지마다 카메라 하나"), DE("eine Kamera je Bild"),
    FR("une caméra par image"), ES("una cámara por imagen"),
    PT("uma câmera por imagem"), IT("una fotocamera per immagine"),
    NL("één camera per beeld"), RU("по камере на снимок"),
    TR("görüntü başına bir kamera"));

SS_MSG(camera_sharing_help,
    EN("How lens parameters are shared. \"Shared\" when everything was shot "
       "with one camera at one zoom. \"Per folder\" for a multi-camera rig "
       "organized one subfolder per camera -- a multi-track 360 video "
       "switches to this on its own. \"Per image\" when zoom or focus changed "
       "between shots."),
    JA("レンズパラメータをどう共有するかです。すべて同じカメラ・同じズームで"
       "撮ったなら「共有カメラ1台」。カメラごとにサブフォルダを分けたマルチ"
       "カメラ装置なら「フォルダごと」で、複数トラックの360度動画では自動的に"
       "これになります。ショットごとにズームやピントが変わったなら「画像ごと」。"),
    ZH_HANS("镜头参数如何共享。全部用同一台相机、同一焦段拍的选“共用一台相机”。"
            "按相机分子文件夹的多相机装置选“每个文件夹一台”——多轨 360 视频会"
            "自动切到这一项。若各次拍摄之间变过焦距或对焦，选“每张图像一台”。"),
    ZH_HANT("鏡頭參數如何共用。全部用同一台相機、同一焦段拍的選「共用一台相機」。"
            "按相機分子資料夾的多相機裝置選「每個資料夾一台」——多軌 360 影片會"
            "自動切到這一項。若各次拍攝之間變過焦距或對焦，選「每張影像一台」。"),
    KO("렌즈 파라미터를 어떻게 공유할지입니다. 전부 같은 카메라·같은 줌으로 "
       "찍었다면 '공유 카메라 하나'. 카메라마다 하위 폴더를 나눈 다중 카메라 "
       "장치라면 '폴더마다 하나'이고, 다중 트랙 360 동영상은 알아서 이쪽으로 "
       "바뀝니다. 촬영 사이에 줌이나 초점이 바뀌었다면 '이미지마다 하나'."),
    DE("Wie Objektivparameter geteilt werden. „Eine gemeinsame Kamera“, wenn "
       "alles mit einer Kamera bei einer Brennweite aufgenommen wurde. „Je "
       "Ordner“ für ein Mehrkamera-Rig mit einem Unterordner pro Kamera -- "
       "ein mehrspuriges 360-Video schaltet von selbst darauf um. „Je Bild“, "
       "wenn Zoom oder Fokus zwischen den Aufnahmen wechselten."),
    FR("Comment les paramètres d'objectif sont partagés. « Une caméra "
       "partagée » si tout a été pris avec un seul appareil à une seule "
       "focale. « Par dossier » pour un rig multi-caméras avec un sous-dossier "
       "par caméra -- une vidéo 360 multipiste bascule là-dessus toute seule. "
       "« Par image » si le zoom ou la mise au point a changé entre les "
       "prises."),
    ES("Cómo se comparten los parámetros de objetivo. «Una cámara compartida» "
       "si todo se tomó con una cámara a un mismo zoom. «Por carpeta» para un "
       "equipo multicámara con una subcarpeta por cámara: un vídeo 360 "
       "multipista cambia solo a esta opción. «Por imagen» si el zoom o el "
       "enfoque cambiaron entre tomas."),
    PT("Como os parâmetros de lente são compartilhados. “Uma câmera "
       "compartilhada” se tudo foi feito com uma câmera num mesmo zoom. “Por "
       "pasta” para um conjunto multicâmera com uma subpasta por câmera -- um "
       "vídeo 360 multipista muda sozinho para isso. “Por imagem” se o zoom "
       "ou o foco mudou entre as tomadas."),
    IT("Come si condividono i parametri dell'obiettivo. «Una fotocamera "
       "condivisa» se è stato ripreso tutto con una fotocamera a uno stesso "
       "zoom. «Per cartella» per un rig multi-fotocamera con una sottocartella "
       "per fotocamera: un video 360 multitraccia passa da solo a questa "
       "voce. «Per immagine» se zoom o messa a fuoco sono cambiati tra gli "
       "scatti."),
    NL("Hoe lensparameters worden gedeeld. ‘Eén gedeelde camera’ als alles met "
       "één camera bij één zoomstand is opgenomen. ‘Per map’ voor een "
       "meercamera-opstelling met een submap per camera -- een 360-video met "
       "meerdere sporen schakelt hier vanzelf op over. ‘Per beeld’ als zoom of "
       "scherpstelling tussen opnamen veranderde."),
    RU("Как разделяются параметры объектива. «Одна общая камера» — если всё "
       "снято одним аппаратом на одном зуме. «По папке» — для многокамерной "
       "установки, где каждой камере отведена подпапка; многодорожечное "
       "360-видео переключается на это само. «По снимку» — если между кадрами "
       "менялись зум или фокус."),
    TR("Objektif parametrelerinin nasıl paylaşılacağı. Her şey tek kamerayla "
       "ve tek yakınlaştırmayla çekildiyse “tek ortak kamera”. Kamera başına "
       "bir alt klasöre ayrılmış çoklu kamera düzeneği için “klasör başına” -- "
       "çok izli bir 360 video kendiliğinden buna geçer. Çekimler arasında "
       "yakınlaştırma veya odak değiştiyse “görüntü başına”."));

SS_MSG(image_matching,
    EN("Image matching"), JA("画像のマッチング"), ZH_HANS("图像匹配"),  ZH_HANT("影像配對"),
    KO("이미지 매칭"),    DE("Bildzuordnung"), FR("Appariement d'images"),
    ES("Emparejamiento de imágenes"), PT("Correspondência de imagens"),
    IT("Corrispondenza tra immagini"), NL("Beeldmatching"),
    RU("Сопоставление снимков"), TR("Görüntü eşleştirme"));

SS_MSG(matching_automatic,
    EN("Automatic"),     JA("自動"),          ZH_HANS("自动"),     ZH_HANT("自動"),
    KO("자동"),           DE("Automatisch"),  FR("Automatique"),  ES("Automático"),
    PT("Automático"),    IT("Automatico"),   NL("Automatisch"),  RU("Автоматически"),
    TR("Otomatik"));

SS_MSG(matching_every_pair,
    EN("Every pair (best, slowest)"),
    JA("すべての組み合わせ（最良・最遅）"),
    ZH_HANS("所有配对（最好、最慢）"),
    ZH_HANT("所有配對（最好、最慢）"),
    KO("모든 쌍(가장 좋고 가장 느림)"),
    DE("Jedes Paar (bestes, langsamstes)"),
    FR("Toutes les paires (meilleur, plus lent)"),
    ES("Todos los pares (mejor, más lento)"),
    PT("Todos os pares (melhor, mais lento)"),
    IT("Tutte le coppie (migliore, più lento)"),
    NL("Elk paar (beste, traagste)"),
    RU("Все пары (лучше всего, медленнее всего)"),
    TR("Her çift (en iyi, en yavaş)"));

SS_MSG(matching_neighbours,
    EN("Neighbouring frames (video)"),
    JA("隣り合うフレーム（動画）"),
    ZH_HANS("相邻帧（视频）"),
    ZH_HANT("相鄰影格（影片）"),
    KO("이웃 프레임(동영상)"),
    DE("Benachbarte Bilder (Video)"),
    FR("Images voisines (vidéo)"),
    ES("Fotogramas vecinos (vídeo)"),
    PT("Quadros vizinhos (vídeo)"),
    IT("Fotogrammi vicini (video)"),
    NL("Naburige beelden (video)"),
    RU("Соседние кадры (видео)"),
    TR("Komşu kareler (video)"));

SS_MSG(matching_gpu_preselect,
    EN("GPU pre-selection (large captures)"),
    JA("GPU による事前選択（大規模な撮影）"),
    ZH_HANS("GPU 预筛选（大规模拍摄）"),
    ZH_HANT("GPU 預篩選（大規模拍攝）"),
    KO("GPU 사전 선별(대규모 촬영)"),
    DE("GPU-Vorauswahl (große Aufnahmen)"),
    FR("Présélection GPU (grandes prises)"),
    ES("Preselección en GPU (capturas grandes)"),
    PT("Pré-seleção na GPU (capturas grandes)"),
    IT("Preselezione su GPU (riprese grandi)"),
    NL("GPU-voorselectie (grote opnamen)"),
    RU("Предварительный отбор на GPU (крупные съёмки)"),
    TR("GPU ön seçimi (büyük çekimler)"));

SS_MSG(matching_help_builtin,
    EN("Which pairs of images are compared. Automatic is right almost always: "
       "neighbouring frames for video, every pair below a hundred images, GPU "
       "pre-selection above that."),
    JA("どの画像の組み合わせを比較するかです。ほとんどの場合「自動」で正しく、"
       "動画なら隣接フレーム、100 枚未満ならすべての組み合わせ、それ以上なら "
       "GPU による事前選択が使われます。"),
    ZH_HANS("比较哪些图像配对。“自动”几乎总是对的：视频用相邻帧，不足一百张时"
            "比较所有配对，超过则用 GPU 预筛选。"),
    ZH_HANT("比較哪些影像配對。「自動」幾乎總是對的：影片用相鄰影格，不足一百張時"
            "比較所有配對，超過則用 GPU 預篩選。"),
    KO("어떤 이미지 쌍을 비교할지입니다. 거의 언제나 '자동'이 맞습니다. 동영상은 "
       "이웃 프레임, 100장 미만이면 모든 쌍, 그보다 많으면 GPU 사전 선별을 "
       "씁니다."),
    DE("Welche Bildpaare verglichen werden. Automatisch ist fast immer "
       "richtig: benachbarte Bilder bei Video, jedes Paar unter hundert "
       "Bildern, darüber GPU-Vorauswahl."),
    FR("Quelles paires d'images sont comparées. Automatique a presque "
       "toujours raison : images voisines pour la vidéo, toutes les paires en "
       "dessous de cent images, présélection GPU au-delà."),
    ES("Qué pares de imágenes se comparan. Automático acierta casi siempre: "
       "fotogramas vecinos en vídeo, todos los pares por debajo de cien "
       "imágenes y preselección en GPU por encima."),
    PT("Quais pares de imagens são comparados. Automático acerta quase "
       "sempre: quadros vizinhos em vídeo, todos os pares abaixo de cem "
       "imagens e pré-seleção na GPU acima disso."),
    IT("Quali coppie di immagini vengono confrontate. Automatico è quasi "
       "sempre giusto: fotogrammi vicini per il video, tutte le coppie sotto "
       "le cento immagini, preselezione su GPU oltre."),
    NL("Welke beeldparen worden vergeleken. Automatisch klopt bijna altijd: "
       "naburige beelden bij video, elk paar onder de honderd beelden, "
       "GPU-voorselectie daarboven."),
    RU("Какие пары снимков сравниваются. «Автоматически» почти всегда верно: "
       "соседние кадры для видео, все пары при менее чем ста снимках и "
       "предварительный отбор на GPU сверх того."),
    TR("Hangi görüntü çiftlerinin karşılaştırılacağı. Otomatik neredeyse her "
       "zaman doğrudur: videoda komşu kareler, yüz görüntünün altında her "
       "çift, üstünde GPU ön seçimi."));

SS_MSG(matching_exhaustive,
    EN("Exhaustive"),    JA("総当たり"),      ZH_HANS("穷举"),     ZH_HANT("窮舉"),
    KO("전수"),           DE("Vollständig"),  FR("Exhaustif"),    ES("Exhaustivo"),
    PT("Exaustivo"),     IT("Esaustivo"),    NL("Uitputtend"),   RU("Полный перебор"),
    TR("Kapsamlı"));

SS_MSG(matching_sequential,
    EN("Sequential"),    JA("逐次"),          ZH_HANS("顺序"),     ZH_HANT("循序"),
    KO("순차"),           DE("Sequenziell"),  FR("Séquentiel"),   ES("Secuencial"),
    PT("Sequencial"),    IT("Sequenziale"),  NL("Sequentieel"),  RU("Последовательное"),
    TR("Sıralı"));

SS_MSG(matching_vocab_tree,
    EN("Vocabulary tree"), JA("ボキャブラリツリー"), ZH_HANS("词汇树"), ZH_HANT("詞彙樹"),
    KO("어휘 트리"),      DE("Vokabularbaum"), FR("Arbre de vocabulaire"),
    ES("Árbol de vocabulario"), PT("Árvore de vocabulário"),
    IT("Albero di vocabolario"), NL("Vocabulaireboom"),
    RU("Словарное дерево"), TR("Sözcük ağacı"));

SS_MSG(matching_help_colmap,
    EN("How image pairs are matched. Exhaustive tries every pair (best "
       "quality, fine up to a few hundred images). Sequential matches temporal "
       "neighbors (video). Vocabulary tree scales to thousands of unordered "
       "photos."),
    JA("画像ペアのマッチング方法です。総当たりはすべての組み合わせを試します"
       "（品質は最良で、数百枚までなら問題ありません）。逐次は時間的に隣り合う"
       "フレームを照合します（動画向け）。ボキャブラリツリーは順不同の数千枚"
       "規模まで対応できます。"),
    ZH_HANS("图像配对的匹配方式。穷举会试遍所有配对（质量最好，几百张以内没问题）。"
            "顺序匹配时间上相邻的帧（视频用）。词汇树可以扩展到数千张无序照片。"),
    ZH_HANT("影像配對的比對方式。窮舉會試遍所有配對（品質最好，幾百張以內沒問題）。"
            "循序比對時間上相鄰的影格（影片用）。詞彙樹可以擴展到數千張無序相片。"),
    KO("이미지 쌍을 어떻게 매칭할지입니다. 전수는 모든 쌍을 시도합니다(품질이 "
       "가장 좋고 수백 장까지는 괜찮습니다). 순차는 시간적으로 이웃한 프레임을 "
       "맞춥니다(동영상용). 어휘 트리는 순서 없는 수천 장까지 감당합니다."),
    DE("Wie Bildpaare zugeordnet werden. Vollständig probiert jedes Paar "
       "(beste Qualität, bis zu einigen hundert Bildern unproblematisch). "
       "Sequenziell ordnet zeitliche Nachbarn zu (Video). Der Vokabularbaum "
       "skaliert auf Tausende ungeordneter Fotos."),
    FR("Comment les paires d'images sont appariées. Exhaustif essaie chaque "
       "paire (meilleure qualité, sans problème jusqu'à quelques centaines "
       "d'images). Séquentiel apparie les voisins temporels (vidéo). L'arbre "
       "de vocabulaire monte à des milliers de photos sans ordre."),
    ES("Cómo se emparejan los pares de imágenes. Exhaustivo prueba todos los "
       "pares (mejor calidad, sin problema hasta unos cientos de imágenes). "
       "Secuencial empareja vecinos temporales (vídeo). El árbol de "
       "vocabulario escala a miles de fotos sin orden."),
    PT("Como os pares de imagens são correspondidos. Exaustivo tenta todos os "
       "pares (melhor qualidade, tranquilo até algumas centenas de imagens). "
       "Sequencial casa vizinhos temporais (vídeo). A árvore de vocabulário "
       "escala para milhares de fotos sem ordem."),
    IT("Come vengono messe in corrispondenza le coppie di immagini. Esaustivo "
       "prova ogni coppia (qualità migliore, tranquillo fino a qualche "
       "centinaio di immagini). Sequenziale abbina i vicini temporali "
       "(video). L'albero di vocabolario arriva a migliaia di foto senza "
       "ordine."),
    NL("Hoe beeldparen worden gematcht. Uitputtend probeert elk paar (beste "
       "kwaliteit, tot enkele honderden beelden prima). Sequentieel matcht "
       "tijdelijke buren (video). De vocabulaireboom schaalt naar duizenden "
       "ongeordende foto's."),
    RU("Как сопоставляются пары снимков. Полный перебор пробует все пары "
       "(лучшее качество, до нескольких сотен снимков вполне нормально). "
       "Последовательное сопоставляет соседей по времени (видео). Словарное "
       "дерево тянет тысячи неупорядоченных фотографий."),
    TR("Görüntü çiftlerinin nasıl eşleştirileceği. Kapsamlı her çifti dener "
       "(en iyi kalite, birkaç yüz görüntüye kadar sorunsuz). Sıralı zamansal "
       "komşuları eşleştirir (video). Sözcük ağacı sırasız binlerce fotoğrafa "
       "ölçeklenir."));

SS_MSG(loop_closure,
    EN("Loop closure detection"),
    JA("ループ閉じ込みの検出"),
    ZH_HANS("回环检测"),
    ZH_HANT("迴環偵測"),
    KO("루프 클로저 검출"),
    DE("Schleifenschluss erkennen"),
    FR("Détection de fermeture de boucle"),
    ES("Detección de cierre de bucle"),
    PT("Detecção de fechamento de laço"),
    IT("Rilevamento della chiusura del giro"),
    NL("Lusdetectie"),
    RU("Обнаружение замыкания петли"),
    TR("Çevrim kapanışı algılama"));

SS_MSG(loop_closure_help_colmap,
    EN("Retrieve visually similar non-neighbour frames through the vocabulary "
       "tree and match them, closing loops when the camera revisits a spot. "
       "SIFT features only."),
    JA("ボキャブラリツリーで見た目の似た非隣接フレームを探して照合し、カメラが"
       "同じ場所に戻ったときにループを閉じます。SIFT 特徴でのみ使えます。"),
    ZH_HANS("通过词汇树找出视觉上相似但不相邻的帧并加以匹配，在相机回到同一位置时"
            "闭合回环。仅适用于 SIFT 特征。"),
    ZH_HANT("透過詞彙樹找出視覺上相似但不相鄰的影格並加以比對，在相機回到同一位置時"
            "閉合迴環。僅適用於 SIFT 特徵。"),
    KO("어휘 트리로 시각적으로 비슷한 비이웃 프레임을 찾아 매칭해, 카메라가 같은 "
       "자리로 돌아왔을 때 루프를 닫습니다. SIFT 특징점에서만 됩니다."),
    DE("Über den Vokabularbaum optisch ähnliche, nicht benachbarte Bilder "
       "finden und zuordnen, damit sich Schleifen schließen, wenn die Kamera "
       "an einen Ort zurückkehrt. Nur mit SIFT-Merkmalen."),
    FR("Retrouver via l'arbre de vocabulaire des images non voisines mais "
       "visuellement proches et les apparier, ce qui ferme la boucle quand la "
       "caméra repasse au même endroit. Uniquement avec les points SIFT."),
    ES("Recuperar mediante el árbol de vocabulario fotogramas no vecinos pero "
       "visualmente parecidos y emparejarlos, cerrando el bucle cuando la "
       "cámara vuelve a un mismo punto. Solo con características SIFT."),
    PT("Recuperar pela árvore de vocabulário quadros não vizinhos mas "
       "visualmente parecidos e casá-los, fechando o laço quando a câmera "
       "volta ao mesmo ponto. Só com características SIFT."),
    IT("Recuperare tramite l'albero di vocabolario fotogrammi non vicini ma "
       "visivamente simili e abbinarli, chiudendo il giro quando la "
       "fotocamera ripassa da uno stesso punto. Solo con caratteristiche "
       "SIFT."),
    NL("Via de vocabulaireboom visueel gelijkende niet-naburige beelden "
       "opzoeken en matchen, zodat de lus sluit als de camera een plek weer "
       "aandoet. Alleen met SIFT-kenmerken."),
    RU("Находить через словарное дерево внешне похожие несоседние кадры и "
       "сопоставлять их, замыкая петлю, когда камера возвращается на то же "
       "место. Только с признаками SIFT."),
    TR("Sözcük ağacı üzerinden görsel olarak benzeyen komşu olmayan kareleri "
       "bulup eşleştirir; böylece kamera aynı noktaya döndüğünde çevrim "
       "kapanır. Yalnızca SIFT öznitelikleriyle."));

SS_MSG(loop_closure_help_builtin,
    EN("Sequential matching only pairs each frame with its temporal "
       "neighbours, so a capture that walks around a subject and returns has "
       "nothing joining the two ends -- one weak step then splits the "
       "reconstruction into pieces. This also matches frames that look alike "
       "wherever they fall in the video, which closes the loop. Costs a "
       "pair-selection pass and roughly twice the matching time; enabled by "
       "default. Under \"Automatic\" it only applies below 100 frames -- "
       "above that, matching is content-based already."),
    JA("逐次マッチングは各フレームを時間的な隣とだけ組にするため、被写体の"
       "まわりを一周して戻ってくる撮影では両端をつなぐものがありません。"
       "弱いつなぎ目が1か所あるだけで再構成はばらばらに割れてしまいます。"
       "これを有効にすると、動画のどこにあっても見た目の似たフレームどうしを"
       "照合するので、ループが閉じます。ペア選択のパスが1回増え、マッチング"
       "時間はおよそ2倍になります。既定で有効です。「自動」では 100 フレーム"
       "未満のときだけ適用されます。それ以上ではマッチングがすでに内容ベース"
       "だからです。"),
    ZH_HANS("顺序匹配只把每一帧和它时间上的邻居配对，所以绕着被摄物走一圈再回来的"
            "拍摄，两端之间没有任何连接——只要有一处连接较弱，重建就会碎成几块。"
            "开启后还会匹配视频中任何位置上看起来相似的帧，从而闭合回环。代价是"
            "多一遍配对筛选和大约两倍的匹配时间；默认开启。在“自动”下只在不足"
            "100 帧时生效——超过之后匹配本来就是基于内容的。"),
    ZH_HANT("循序比對只把每一影格和它時間上的鄰居配對，所以繞著被攝物走一圈再回來的"
            "拍攝，兩端之間沒有任何連接——只要有一處連接較弱，重建就會碎成幾塊。"
            "開啟後還會比對影片中任何位置上看起來相似的影格，從而閉合迴環。代價是"
            "多一遍配對篩選和大約兩倍的比對時間；預設開啟。在「自動」下只在不足"
            "100 影格時生效——超過之後比對本來就是基於內容的。"),
    KO("순차 매칭은 각 프레임을 시간적으로 이웃한 것과만 짝지으므로, 피사체 "
       "주위를 한 바퀴 돌아 돌아오는 촬영에서는 양 끝을 잇는 것이 없습니다. "
       "약한 연결이 한 군데만 있어도 재구성이 조각으로 갈라집니다. 이 옵션은 "
       "동영상 어디에 있든 비슷해 보이는 프레임끼리도 매칭해 루프를 닫습니다. "
       "쌍 선별 패스 한 번과 대략 두 배의 매칭 시간이 들며 기본으로 켜져 "
       "있습니다. '자동'에서는 100프레임 미만일 때만 적용됩니다. 그보다 많으면 "
       "매칭이 이미 내용 기반입니다."),
    DE("Sequenzielle Zuordnung paart jedes Bild nur mit seinen zeitlichen "
       "Nachbarn; bei einer Aufnahme, die um ein Motiv herumgeht und "
       "zurückkommt, verbindet also nichts die beiden Enden -- eine einzige "
       "schwache Stelle zerlegt die Rekonstruktion dann in Teile. Dies ordnet "
       "zusätzlich ähnlich aussehende Bilder zu, egal wo sie im Video liegen, "
       "und schließt so die Schleife. Kostet einen Durchgang zur Paarauswahl "
       "und etwa die doppelte Zuordnungszeit; standardmäßig an. Unter "
       "„Automatisch“ greift es nur unter 100 Bildern -- darüber ist die "
       "Zuordnung ohnehin inhaltsbasiert."),
    FR("L'appariement séquentiel ne relie chaque image qu'à ses voisines dans "
       "le temps ; une prise qui fait le tour d'un sujet et revient n'a donc "
       "rien qui joigne les deux extrémités -- un seul maillon faible et la "
       "reconstruction se scinde. Ceci apparie en plus les images qui se "
       "ressemblent où qu'elles soient dans la vidéo, ce qui ferme la boucle. "
       "Coûte une passe de sélection de paires et environ le double du temps "
       "d'appariement ; activé par défaut. En « Automatique », ne s'applique "
       "qu'en dessous de 100 images : au-delà, l'appariement est déjà fondé "
       "sur le contenu."),
    ES("El emparejamiento secuencial solo une cada fotograma con sus vecinos "
       "temporales, así que una captura que rodea un sujeto y vuelve no tiene "
       "nada que una los dos extremos: un solo eslabón débil parte la "
       "reconstrucción en trozos. Esto empareja además fotogramas parecidos "
       "estén donde estén en el vídeo, cerrando el bucle. Cuesta una pasada "
       "de selección de pares y aproximadamente el doble de tiempo de "
       "emparejamiento; activado por defecto. En «Automático» solo se aplica "
       "por debajo de 100 fotogramas: por encima, el emparejamiento ya es por "
       "contenido."),
    PT("A correspondência sequencial só liga cada quadro aos seus vizinhos no "
       "tempo, então uma captura que dá a volta num objeto e retorna não tem "
       "nada unindo as duas pontas -- um único elo fraco parte a reconstrução "
       "em pedaços. Isto casa também quadros parecidos onde quer que estejam "
       "no vídeo, fechando o laço. Custa uma passagem de seleção de pares e "
       "cerca do dobro do tempo de correspondência; ligado por padrão. Em "
       "“Automático” só vale abaixo de 100 quadros: acima disso, a "
       "correspondência já é por conteúdo."),
    IT("La corrispondenza sequenziale accoppia ogni fotogramma solo con i "
       "vicini nel tempo, quindi una ripresa che gira attorno a un soggetto e "
       "torna indietro non ha nulla che unisca i due capi: basta un anello "
       "debole e la ricostruzione si spezza. Questo abbina anche fotogrammi "
       "simili ovunque cadano nel video, chiudendo il giro. Costa una passata "
       "di selezione delle coppie e circa il doppio del tempo di "
       "corrispondenza; attivo di default. In «Automatico» vale solo sotto i "
       "100 fotogrammi: oltre, la corrispondenza è già basata sul contenuto."),
    NL("Sequentieel matchen koppelt elk beeld alleen aan zijn buren in de "
       "tijd, dus bij een opname die om een onderwerp heen loopt en terugkomt "
       "verbindt niets de twee uiteinden -- één zwakke schakel en de "
       "reconstructie valt uiteen. Dit matcht daarnaast beelden die op elkaar "
       "lijken, waar ze ook in de video vallen, en sluit zo de lus. Kost een "
       "extra selectieronde en ongeveer het dubbele van de matchtijd; staat "
       "standaard aan. Onder ‘Automatisch’ geldt het alleen onder 100 "
       "beelden -- daarboven is het matchen al inhoudsgebaseerd."),
    RU("Последовательное сопоставление связывает каждый кадр только с "
       "соседями по времени, поэтому у съёмки, которая обходит объект и "
       "возвращается, два конца ничем не соединены — достаточно одного "
       "слабого звена, и реконструкция распадается на куски. Здесь вдобавок "
       "сопоставляются похожие кадры, где бы они ни оказались в видео, и петля "
       "замыкается. Стоит одного прохода отбора пар и примерно вдвое большего "
       "времени сопоставления; включено по умолчанию. При «Автоматически» "
       "применяется только до 100 кадров — выше сопоставление и так идёт по "
       "содержимому."),
    TR("Sıralı eşleştirme her kareyi yalnızca zamandaki komşularıyla "
       "eşleştirir; bir öznenin çevresini dolaşıp geri dönen bir çekimde iki "
       "ucu birleştiren hiçbir şey olmaz -- tek bir zayıf halka bile yeniden "
       "oluşturmayı parçalara böler. Bu seçenek ayrıca videoda nerede olursa "
       "olsun birbirine benzeyen kareleri de eşleştirerek çevrimi kapatır. Bir "
       "çift seçme geçişi ve kabaca iki katı eşleştirme süresi gerektirir; "
       "varsayılan olarak açıktır. “Otomatik”te yalnızca 100 karenin altında "
       "geçerlidir -- üstünde eşleştirme zaten içerik temellidir."));

SS_MSG(frames_per_second,
    EN("Frames per second"),
    JA("1秒あたりのフレーム数"),
    ZH_HANS("每秒帧数"),
    ZH_HANT("每秒影格數"),
    KO("초당 프레임 수"),
    DE("Bilder pro Sekunde"),
    FR("Images par seconde"),
    ES("Fotogramas por segundo"),
    PT("Quadros por segundo"),
    IT("Fotogrammi al secondo"),
    NL("Beelden per seconde"),
    RU("Кадров в секунду"),
    TR("Saniyedeki kare"));

SS_MSG(frames_per_second_help,
    EN("How many frames to keep per second of video. 1-3 is right for a slow "
       "walkthrough; more only helps if the camera moved fast. Applies to "
       "every video in the list."),
    JA("動画1秒あたり何フレーム残すかです。ゆっくり歩いて撮ったなら 1〜3 が"
       "適切で、それ以上が効くのはカメラが速く動いたときだけです。リスト内の"
       "すべての動画に適用されます。"),
    ZH_HANS("每秒视频保留多少帧。慢慢走着拍的话 1-3 就合适；更高只有在相机移动"
            "很快时才有用。对列表中的所有视频都生效。"),
    ZH_HANT("每秒影片保留多少影格。慢慢走著拍的話 1-3 就合適；更高只有在相機移動"
            "很快時才有用。對清單中的所有影片都生效。"),
    KO("동영상 1초당 몇 프레임을 남길지입니다. 천천히 걸으며 찍었다면 1~3이 "
       "알맞고, 그보다 높이는 건 카메라가 빠르게 움직였을 때만 도움이 됩니다. "
       "목록의 모든 동영상에 적용됩니다."),
    DE("Wie viele Bilder je Sekunde Video behalten werden. 1-3 passt für "
       "einen langsamen Rundgang; mehr hilft nur, wenn die Kamera schnell "
       "bewegt wurde. Gilt für jedes Video in der Liste."),
    FR("Combien d'images conserver par seconde de vidéo. 1 à 3 convient à une "
       "déambulation lente ; davantage n'aide que si la caméra bougeait vite. "
       "S'applique à toutes les vidéos de la liste."),
    ES("Cuántos fotogramas conservar por segundo de vídeo. De 1 a 3 va bien "
       "para un recorrido lento; más solo ayuda si la cámara se movía rápido. "
       "Se aplica a todos los vídeos de la lista."),
    PT("Quantos quadros manter por segundo de vídeo. De 1 a 3 serve para um "
       "percurso lento; mais só ajuda se a câmera se moveu rápido. Vale para "
       "todos os vídeos da lista."),
    IT("Quanti fotogrammi tenere per ogni secondo di video. Da 1 a 3 va bene "
       "per una camminata lenta; di più serve solo se la fotocamera si "
       "muoveva in fretta. Vale per tutti i video dell'elenco."),
    NL("Hoeveel beelden per seconde video bewaard blijven. 1-3 past bij een "
       "rustige rondgang; meer helpt alleen als de camera snel bewoog. Geldt "
       "voor elke video in de lijst."),
    RU("Сколько кадров оставлять на секунду видео. 1-3 подходит для "
       "неторопливого обхода; больше помогает, только если камера двигалась "
       "быстро. Применяется ко всем видео в списке."),
    TR("Videonun her saniyesinden kaç karenin tutulacağı. Yavaş bir gezinti "
       "için 1-3 uygundur; daha fazlası yalnızca kamera hızlı hareket ettiyse "
       "işe yarar. Listedeki bütün videolara uygulanır."));

SS_MSG(sharpness_window,
    EN("Sharpness window"),
    JA("シャープさの判定窓"),
    ZH_HANS("清晰度窗口"),
    ZH_HANT("清晰度視窗"),
    KO("선명도 창"),
    DE("Schärfefenster"),
    FR("Fenêtre de netteté"),
    ES("Ventana de nitidez"),
    PT("Janela de nitidez"),
    IT("Finestra di nitidezza"),
    NL("Scherptevenster"),
    RU("Окно резкости"),
    TR("Keskinlik penceresi"));

SS_MSG(sharpness_window_help,
    EN("Look at this many candidate frames for each one kept and keep the "
       "least motion-blurred. 1 turns the selection off."),
    JA("1枚残すごとにこの数だけ候補フレームを見て、いちばんブレの少ないものを"
       "残します。1 にすると選別しません。"),
    ZH_HANS("每保留一帧就看这么多候选帧，从中挑运动模糊最少的。设为 1 就不做筛选。"),
    ZH_HANT("每保留一影格就看這麼多候選影格，從中挑運動模糊最少的。設為 1 就不做篩選。"),
    KO("한 장을 남길 때마다 이만큼의 후보 프레임을 보고 흔들림이 가장 적은 것을 "
       "남깁니다. 1이면 선별하지 않습니다."),
    DE("Für jedes behaltene Bild so viele Kandidaten ansehen und das am "
       "wenigsten verwackelte nehmen. 1 schaltet die Auswahl ab."),
    FR("Examiner ce nombre d'images candidates pour chaque image conservée et "
       "garder la moins floue. 1 désactive la sélection."),
    ES("Mirar esta cantidad de fotogramas candidatos por cada uno conservado "
       "y quedarse con el menos movido. 1 desactiva la selección."),
    PT("Olhar esta quantidade de quadros candidatos para cada um mantido e "
       "ficar com o menos tremido. 1 desliga a seleção."),
    IT("Guardare questo numero di fotogrammi candidati per ognuno tenuto e "
       "conservare il meno mosso. 1 disattiva la selezione."),
    NL("Zoveel kandidaatbeelden bekijken voor elk bewaard beeld en het minst "
       "bewogen exemplaar houden. 1 zet de selectie uit."),
    RU("Просматривать столько кадров-кандидатов на каждый оставленный и брать "
       "наименее смазанный. 1 отключает отбор."),
    TR("Tutulan her kare için bu kadar aday kareye bakıp en az bulanık olanı "
       "tutar. 1, seçimi kapatır."));

// {0} is a build-configuration note from the backend, in English.
SS_MSG(build_note,
    EN("Note: {0}."),    JA("注意: {0}。"),   ZH_HANS("注意：{0}。"), ZH_HANT("注意：{0}。"),
    KO("참고: {0}."),     DE("Hinweis: {0}."), FR("Remarque : {0}."),
    ES("Nota: {0}."),    PT("Observação: {0}."), IT("Nota: {0}."),
    NL("Let op: {0}."),  RU("Примечание: {0}."), TR("Not: {0}."));


// ===========================================================================
// Masking
// ===========================================================================

SS_MSG(masks_found_all,
    EN("Masks were found beside the images; they are used as they are."),
    JA("画像の隣にマスクが見つかりました。そのまま使われます。"),
    ZH_HANS("在图像旁边找到了蒙版，将按原样使用。"),
    ZH_HANT("在影像旁邊找到了遮罩，將按原樣使用。"),
    KO("이미지 옆에서 마스크를 찾았습니다. 그대로 사용합니다."),
    DE("Neben den Bildern wurden Masken gefunden; sie werden unverändert "
       "genutzt."),
    FR("Des masques ont été trouvés à côté des images ; ils sont utilisés tels "
       "quels."),
    ES("Se encontraron máscaras junto a las imágenes; se usan tal cual."),
    PT("Foram encontradas máscaras ao lado das imagens; elas são usadas como "
       "estão."),
    IT("Accanto alle immagini sono state trovate delle maschere; vengono usate "
       "così come sono."),
    NL("Naast de beelden zijn maskers gevonden; ze worden gebruikt zoals ze "
       "zijn."),
    RU("Рядом со снимками найдены маски; они используются как есть."),
    TR("Görüntülerin yanında maskeler bulundu; oldukları gibi kullanılıyor."));

// {0} inputs with masks, {1} inputs in total.
SS_MSG(masks_found_some,
    EN("{0} of the {1} inputs came with masks; those are used as they are."),
    JA("{1} 件の入力のうち {0} 件にマスクが付いていました。それらはそのまま"
       "使われます。"),
    ZH_HANS("{1} 个输入中有 {0} 个自带蒙版；这些会按原样使用。"),
    ZH_HANT("{1} 個輸入中有 {0} 個自帶遮罩；這些會按原樣使用。"),
    KO("입력 {1}개 가운데 {0}개에 마스크가 딸려 있습니다. 그것들은 그대로 "
       "사용합니다."),
    DE("{0} der {1} Eingaben brachten Masken mit; diese werden unverändert "
       "genutzt."),
    FR("{0} des {1} entrées sont fournies avec des masques ; ceux-là sont "
       "utilisés tels quels."),
    ES("{0} de las {1} entradas traían máscaras; esas se usan tal cual."),
    PT("{0} das {1} entradas vieram com máscaras; essas são usadas como estão."),
    IT("{0} dei {1} ingressi hanno già delle maschere; quelle vengono usate "
       "così come sono."),
    NL("{0} van de {1} invoeren kwamen met maskers; die worden gebruikt zoals "
       "ze zijn."),
    RU("Маски были у {0} из {1} входов; они используются как есть."),
    TR("{1} girdiden {0} tanesi maskeyle geldi; onlar oldukları gibi "
       "kullanılıyor."));

SS_MSG(mask_enable,
    EN("Remove moving or unwanted objects"),
    JA("動くものや不要なものを取り除く"),
    ZH_HANS("移除移动或不需要的物体"),
    ZH_HANT("移除移動或不需要的物體"),
    KO("움직이거나 원하지 않는 물체 제거"),
    DE("Bewegte oder unerwünschte Objekte entfernen"),
    FR("Retirer les objets mobiles ou indésirables"),
    ES("Quitar objetos en movimiento o no deseados"),
    PT("Remover objetos em movimento ou indesejados"),
    IT("Rimuovere oggetti in movimento o indesiderati"),
    NL("Bewegende of ongewenste objecten verwijderen"),
    RU("Убрать движущиеся или лишние объекты"),
    TR("Hareketli veya istenmeyen nesneleri kaldır"));

SS_MSG(mask_enable_help,
    EN("Describe what should not be part of the scene -- people walking "
       "through, parked cars, your own shadow -- and it is cut out of both "
       "the camera solve and the training. This is the single biggest quality "
       "win on a capture with anything moving in it. Inputs that arrived with "
       "masks of their own keep them; this is for the rest."),
    JA("シーンに含めたくないもの、たとえば通りがかりの人、駐車中の車、自分の影"
       "などを書いてください。カメラの推定と学習の両方から取り除かれます。"
       "動くものが写っている撮影では、これがいちばん効く品質改善です。"
       "すでにマスクが付いている入力はそのまま使われ、これはそれ以外に"
       "適用されます。"),
    ZH_HANS("写出不该属于场景的东西——路过的行人、停着的车、你自己的影子——"
            "它们会同时从相机求解和训练中被剔除。对于画面里有任何移动物体的"
            "拍摄，这是提升质量最有效的一招。自带蒙版的输入会保留自己的蒙版，"
            "这里针对的是其余的。"),
    ZH_HANT("寫出不該屬於場景的東西——路過的行人、停著的車、你自己的影子——"
            "它們會同時從相機求解和訓練中被剔除。對於畫面裡有任何移動物體的"
            "拍攝，這是提升品質最有效的一招。自帶遮罩的輸入會保留自己的遮罩，"
            "這裡針對的是其餘的。"),
    KO("장면에 들어가면 안 되는 것을 적으세요. 지나가는 사람, 주차된 차, 자기 "
       "그림자 같은 것들이며, 카메라 계산과 학습 양쪽에서 잘려 나갑니다. 움직이는 "
       "것이 있는 촬영에서는 이것이 품질을 가장 크게 끌어올립니다. 이미 마스크가 "
       "딸려 온 입력은 그것을 그대로 쓰고, 이 설정은 나머지에 적용됩니다."),
    DE("Beschreiben Sie, was nicht zur Szene gehören soll -- durchlaufende "
       "Personen, parkende Autos, der eigene Schatten -- und es fällt sowohl "
       "aus der Kameraberechnung als auch aus dem Training heraus. Bei einer "
       "Aufnahme mit irgendetwas Bewegtem ist das der größte einzelne "
       "Qualitätsgewinn. Eingaben, die eigene Masken mitbrachten, behalten "
       "sie; dies gilt dem Rest."),
    FR("Décrivez ce qui ne doit pas faire partie de la scène -- passants, "
       "voitures garées, votre propre ombre -- et cela est retiré à la fois du "
       "calcul des caméras et de l'entraînement. Sur une prise où quelque "
       "chose bouge, c'est le plus gros gain de qualité possible. Les entrées "
       "arrivées avec leurs propres masques les gardent ; ceci vaut pour les "
       "autres."),
    ES("Describa lo que no debe formar parte de la escena -- gente que pasa, "
       "coches aparcados, su propia sombra -- y se elimina tanto del cálculo "
       "de cámaras como del entrenamiento. En una captura con algo en "
       "movimiento, es la mayor mejora de calidad posible. Las entradas que "
       "llegaron con sus propias máscaras las conservan; esto vale para el "
       "resto."),
    PT("Descreva o que não deve fazer parte da cena -- pessoas passando, "
       "carros estacionados, sua própria sombra -- e isso é retirado tanto do "
       "cálculo das câmeras quanto do treinamento. Numa captura com qualquer "
       "coisa em movimento, é o maior ganho de qualidade que existe. Entradas "
       "que chegaram com máscaras próprias as mantêm; isto vale para as "
       "demais."),
    IT("Descriva ciò che non deve far parte della scena -- passanti, "
       "automobili in sosta, la propria ombra -- e verrà tolto sia dal calcolo "
       "delle fotocamere sia dall'addestramento. In una ripresa con qualcosa "
       "in movimento è il singolo guadagno di qualità più grande. Gli ingressi "
       "arrivati con maschere proprie le conservano; questo vale per gli "
       "altri."),
    NL("Beschrijf wat geen deel van de scène hoort te zijn -- voorbijgangers, "
       "geparkeerde auto's, je eigen schaduw -- en het valt zowel uit de "
       "cameraberekening als uit de training weg. Bij een opname waarin iets "
       "beweegt is dit de grootste kwaliteitswinst die er is. Invoeren die met "
       "eigen maskers kwamen, houden die; dit geldt voor de rest."),
    RU("Опишите, чего в сцене быть не должно, — прохожих, припаркованных машин, "
       "собственной тени, — и это будет исключено и из расчёта камер, и из "
       "обучения. Для съёмки, где хоть что-то движется, это самый большой "
       "выигрыш в качестве. Входы, пришедшие со своими масками, их сохраняют; "
       "это для остальных."),
    TR("Sahnenin parçası olmaması gerekeni yazın -- geçen insanlar, park etmiş "
       "arabalar, kendi gölgeniz -- ve bu hem kamera çözümünden hem de "
       "eğitimden çıkarılır. İçinde hareketli bir şey olan çekimlerde kaliteyi "
       "en çok artıran tek şey budur. Kendi maskeleriyle gelen girdiler "
       "onları korur; bu ayar geri kalanlar içindir."));

SS_MSG(mask_model,
    EN("Model"),         JA("モデル"),        ZH_HANS("模型"),     ZH_HANT("模型"),
    KO("모델"),           DE("Modell"),       FR("Modèle"),       ES("Modelo"),
    PT("Modelo"),        IT("Modello"),      NL("Model"),        RU("Модель"),
    TR("Model"));

SS_MSG(mask_model_needs_download,
    EN("{0}  (download)"),
    JA("{0}（ダウンロード）"),
    ZH_HANS("{0}（需下载）"),
    ZH_HANT("{0}（需下載）"),
    KO("{0}(내려받기)"),
    DE("{0}  (Download)"),
    FR("{0}  (à télécharger)"),
    ES("{0}  (descarga)"),
    PT("{0}  (baixar)"),
    IT("{0}  (da scaricare)"),
    NL("{0}  (downloaden)"),
    RU("{0}  (загрузка)"),
    TR("{0}  (indirilecek)"));

SS_MSG(mask_get_model,
    EN("Get the model"), JA("モデルを取得"),   ZH_HANS("获取模型"),  ZH_HANT("取得模型"),
    KO("모델 받기"),      DE("Modell holen"), FR("Obtenir le modèle"),
    ES("Obtener el modelo"), PT("Obter o modelo"), IT("Ottieni il modello"),
    NL("Model ophalen"), RU("Получить модель"), TR("Modeli getir"));

SS_MSG(mask_one_time_download,
    EN("one-time download, kept for next time"),
    JA("初回だけのダウンロードで、次回以降は再利用します"),
    ZH_HANS("只需下载一次，之后会一直保留"),
    ZH_HANT("只需下載一次，之後會一直保留"),
    KO("한 번만 내려받고 다음부터는 그대로 씁니다"),
    DE("einmaliger Download, bleibt für das nächste Mal erhalten"),
    FR("téléchargement unique, conservé pour la prochaine fois"),
    ES("descarga única, se conserva para la próxima vez"),
    PT("download único, guardado para a próxima vez"),
    IT("scaricamento una tantum, resta per la prossima volta"),
    NL("eenmalige download, blijft bewaard voor de volgende keer"),
    RU("загрузка один раз, дальше используется сохранённая"),
    TR("bir kez indirilir, bir dahaki sefere saklanır"));

SS_MSG(stop,
    EN("Stop"),          JA("停止"),          ZH_HANS("停止"),     ZH_HANT("停止"),
    KO("중지"),           DE("Stopp"),        FR("Arrêter"),      ES("Detener"),
    PT("Parar"),         IT("Ferma"),        NL("Stoppen"),      RU("Стоп"),
    TR("Durdur"));

SS_MSG(mask_model_ready,
    EN("Model ready."),  JA("モデルの準備ができました。"), ZH_HANS("模型已就绪。"),
    ZH_HANT("模型已就緒。"), KO("모델이 준비되었습니다."), DE("Modell bereit."),
    FR("Modèle prêt."),  ES("Modelo listo."), PT("Modelo pronto."),
    IT("Modello pronto."), NL("Model gereed."), RU("Модель готова."),
    TR("Model hazır."));

SS_MSG(mask_try,
    EN("Try the mask..."),
    JA("マスクを試す…"),  ZH_HANS("试一下蒙版…"), ZH_HANT("試一下遮罩…"),
    KO("마스크 시험해 보기…"), DE("Maske ausprobieren …"), FR("Essayer le masque…"),
    ES("Probar la máscara…"), PT("Testar a máscara…"), IT("Prova la maschera…"),
    NL("Masker uitproberen…"), RU("Проверить маску…"), TR("Maskeyi dene…"));

SS_MSG(mask_try_help,
    EN("Run the prompt on one real frame and see exactly what would be cut "
       "out, before committing to the whole capture."),
    JA("撮影全体に適用する前に、実際のフレーム1枚でプロンプトを試し、何が"
       "取り除かれるのかを確かめられます。"),
    ZH_HANS("在对整段拍摄下手之前，先在一张真实帧上跑一遍提示词，看清到底会剔除什么。"),
    ZH_HANT("在對整段拍攝下手之前，先在一張真實影格上跑一遍提示詞，看清到底會剔除什麼。"),
    KO("전체 촬영에 적용하기 전에, 실제 프레임 한 장에 프롬프트를 돌려 무엇이 "
       "잘려 나가는지 정확히 확인합니다."),
    DE("Die Eingabe an einem echten Bild ausprobieren und genau sehen, was "
       "wegfiele, bevor man sich für die ganze Aufnahme entscheidet."),
    FR("Essayer l'invite sur une vraie image et voir exactement ce qui serait "
       "retiré, avant de l'appliquer à toute la prise."),
    ES("Probar la indicación en un fotograma real y ver exactamente qué se "
       "eliminaría, antes de aplicarlo a toda la captura."),
    PT("Testar o comando num quadro real e ver exatamente o que seria "
       "removido, antes de aplicar à captura inteira."),
    IT("Provare il testo su un fotogramma reale e vedere esattamente che cosa "
       "verrebbe tolto, prima di applicarlo a tutta la ripresa."),
    NL("De prompt op één echt beeld draaien en precies zien wat eruit zou "
       "gaan, voordat je de hele opname vastlegt."),
    RU("Проверить запрос на одном настоящем кадре и увидеть, что именно будет "
       "вырезано, прежде чем применять ко всей съёмке."),
    TR("Bütün çekime uygulamadan önce istemi gerçek bir kare üzerinde "
       "çalıştırıp neyin çıkarılacağını tam olarak görün."));

SS_MSG(mask_on_input,
    EN("on input"),      JA("対象の入力"),    ZH_HANS("在输入"),   ZH_HANT("在輸入"),
    KO("대상 입력"),      DE("an Eingabe"),   FR("sur l'entrée"), ES("en la entrada"),
    PT("na entrada"),    IT("sull'ingresso"), NL("op invoer"),   RU("на входе"),
    TR("girdi üzerinde"));

SS_MSG(mask_on_input_help,
    EN("Which input \"Try the mask\" runs on. The text prompt applies to every "
       "input; clicked objects only to this one."),
    JA("「マスクを試す」をどの入力で実行するかです。テキストのプロンプトは"
       "すべての入力に適用されますが、クリックした物体はこの入力だけに"
       "適用されます。"),
    ZH_HANS("“试一下蒙版”在哪个输入上运行。文本提示对所有输入生效；点选的物体"
            "只对这一个输入生效。"),
    ZH_HANT("「試一下遮罩」在哪個輸入上執行。文字提示對所有輸入生效；點選的物體"
            "只對這一個輸入生效。"),
    KO("[마스크 시험해 보기]를 어느 입력에서 실행할지입니다. 텍스트 프롬프트는 "
       "모든 입력에 적용되지만, 클릭한 물체는 이 입력에만 적용됩니다."),
    DE("An welcher Eingabe „Maske ausprobieren“ läuft. Der Textbefehl gilt für "
       "jede Eingabe; angeklickte Objekte nur für diese."),
    FR("Sur quelle entrée « Essayer le masque » s'exécute. L'invite textuelle "
       "vaut pour toutes les entrées ; les objets cliqués seulement pour "
       "celle-ci."),
    ES("En qué entrada se ejecuta «Probar la máscara». La indicación de texto "
       "vale para todas las entradas; los objetos marcados, solo para esta."),
    PT("Em qual entrada “Testar a máscara” roda. O comando de texto vale para "
       "todas as entradas; os objetos clicados só para esta."),
    IT("Su quale ingresso viene eseguito «Prova la maschera». Il testo vale "
       "per tutti gli ingressi; gli oggetti cliccati solo per questo."),
    NL("Op welke invoer ‘Masker uitproberen’ draait. De tekstprompt geldt voor "
       "elke invoer; aangeklikte objecten alleen voor deze."),
    RU("На каком входе выполняется «Проверить маску». Текстовый запрос "
       "действует на все входы, а отмеченные объекты — только на этот."),
    TR("“Maskeyi dene”nin hangi girdide çalışacağı. Metin istemi bütün "
       "girdilere uygulanır; tıklanan nesneler yalnızca buna."));

SS_MSG(mask_no_text_prompts,
    EN("This model has no text understanding -- use \"Try the mask\" and click "
       "the object instead."),
    JA("このモデルはテキストを理解しません。「マスクを試す」から対象を"
       "クリックして指定してください。"),
    ZH_HANS("这个模型不理解文字——请改用“试一下蒙版”，直接点选目标物体。"),
    ZH_HANT("這個模型不理解文字——請改用「試一下遮罩」，直接點選目標物體。"),
    KO("이 모델은 텍스트를 이해하지 못합니다. [마스크 시험해 보기]에서 물체를 "
       "직접 클릭하세요."),
    DE("Dieses Modell versteht keinen Text -- stattdessen „Maske ausprobieren“ "
       "öffnen und das Objekt anklicken."),
    FR("Ce modèle ne comprend pas le texte : utilisez plutôt « Essayer le "
       "masque » et cliquez sur l'objet."),
    ES("Este modelo no entiende texto: use «Probar la máscara» y haga clic en "
       "el objeto."),
    PT("Este modelo não entende texto -- use “Testar a máscara” e clique no "
       "objeto."),
    IT("Questo modello non capisce il testo: usi «Prova la maschera» e clicchi "
       "sull'oggetto."),
    NL("Dit model begrijpt geen tekst -- gebruik ‘Masker uitproberen’ en klik "
       "het object aan."),
    RU("Эта модель не понимает текст — откройте «Проверить маску» и укажите "
       "объект щелчком."),
    TR("Bu model metni anlamaz -- bunun yerine “Maskeyi dene”yi açıp nesneye "
       "tıklayın."));

// {0} is a count, labelled rather than inflected.
SS_MSG(mask_clicked_objects,
    EN("Clicked objects: {0}. They are tracked through the capture."),
    JA("クリックした物体: {0} 件。撮影全体を通して追跡されます。"),
    ZH_HANS("已点选的物体：{0} 个。它们会在整段拍摄中被跟踪。"),
    ZH_HANT("已點選的物體：{0} 個。它們會在整段拍攝中被追蹤。"),
    KO("클릭한 물체: {0}개. 촬영 전체에 걸쳐 추적됩니다."),
    DE("Angeklickte Objekte: {0}. Sie werden durch die Aufnahme hindurch "
       "verfolgt."),
    FR("Objets cliqués : {0}. Ils sont suivis tout au long de la prise."),
    ES("Objetos marcados: {0}. Se siguen a lo largo de toda la captura."),
    PT("Objetos clicados: {0}. Eles são rastreados por toda a captura."),
    IT("Oggetti cliccati: {0}. Vengono seguiti per tutta la ripresa."),
    NL("Aangeklikte objecten: {0}. Ze worden door de hele opname gevolgd."),
    RU("Отмечено объектов: {0}. Они отслеживаются по всей съёмке."),
    TR("Tıklanan nesneler: {0}. Çekim boyunca izlenirler."));

SS_MSG(mask_forget_clicks,
    EN("Forget them"),   JA("忘れる"),        ZH_HANS("忘掉它们"),  ZH_HANT("忘掉它們"),
    KO("잊기"),           DE("Verwerfen"),    FR("Les oublier"),  ES("Olvidarlos"),
    PT("Esquecê-los"),   IT("Dimenticali"),  NL("Vergeten"),     RU("Забыть их"),
    TR("Unut"));

SS_MSG(mask_forget_clicks_help,
    EN("Objects you pointed at in \"Try the mask\". Each is followed from the "
       "frame you clicked it on, using the model's video memory, so you do not "
       "have to click every frame."),
    JA("「マスクを試す」で指し示した物体です。クリックしたフレームを起点に、"
       "モデルの動画メモリを使って追跡されるので、毎フレームでクリックする"
       "必要はありません。"),
    ZH_HANS("你在“试一下蒙版”里点选的物体。每个都从你点它的那一帧开始，"
            "借助模型的视频记忆一路跟踪，不必逐帧点击。"),
    ZH_HANT("你在「試一下遮罩」裡點選的物體。每個都從你點它的那一影格開始，"
            "藉助模型的影片記憶一路追蹤，不必逐格點擊。"),
    KO("[마스크 시험해 보기]에서 지목한 물체들입니다. 각각 클릭한 프레임에서부터 "
       "모델의 비디오 메모리로 따라가므로 프레임마다 클릭할 필요가 없습니다."),
    DE("Objekte, auf die Sie in „Maske ausprobieren“ gezeigt haben. Jedes wird "
       "ab dem angeklickten Bild über das Videogedächtnis des Modells "
       "weiterverfolgt, Sie müssen also nicht jedes Bild anklicken."),
    FR("Les objets désignés dans « Essayer le masque ». Chacun est suivi à "
       "partir de l'image où vous l'avez cliqué, grâce à la mémoire vidéo du "
       "modèle : inutile de cliquer sur chaque image."),
    ES("Los objetos que señaló en «Probar la máscara». Cada uno se sigue desde "
       "el fotograma donde hizo clic, usando la memoria de vídeo del modelo, "
       "así que no hace falta pulsar en cada fotograma."),
    PT("Os objetos que você apontou em “Testar a máscara”. Cada um é seguido a "
       "partir do quadro em que você clicou, usando a memória de vídeo do "
       "modelo, então não é preciso clicar em cada quadro."),
    IT("Gli oggetti indicati in «Prova la maschera». Ciascuno viene seguito a "
       "partire dal fotogramma su cui ha cliccato, grazie alla memoria video "
       "del modello: non serve cliccare su ogni fotogramma."),
    NL("De objecten die je in ‘Masker uitproberen’ hebt aangewezen. Elk wordt "
       "vanaf het aangeklikte beeld gevolgd met het videogeheugen van het "
       "model, dus je hoeft niet elk beeld aan te klikken."),
    RU("Объекты, которые вы указали в «Проверить маску». Каждый отслеживается "
       "от кадра, где вы по нему щёлкнули, через видеопамять модели, так что "
       "щёлкать на каждом кадре не нужно."),
    TR("“Maskeyi dene”de işaret ettiğiniz nesneler. Her biri tıkladığınız "
       "kareden başlayarak modelin video belleğiyle izlenir, yani her kareye "
       "tıklamanız gerekmez."));

// {0} is the input the clicks were drawn on.
SS_MSG(mask_clicks_one_input_only,
    EN("They describe one frame of {0}, so only that input is masked. Add a "
       "text prompt to cover the others."),
    JA("これらは {0} の1フレームを指しているため、マスクされるのはその入力だけ"
       "です。ほかもカバーするにはテキストのプロンプトを追加してください。"),
    ZH_HANS("它们描述的是 {0} 的某一帧，所以只有这个输入会被蒙版。"
            "要覆盖其他输入，请加上文本提示。"),
    ZH_HANT("它們描述的是 {0} 的某一影格，所以只有這個輸入會被遮罩。"
            "要涵蓋其他輸入，請加上文字提示。"),
    KO("이것들은 {0}의 한 프레임을 가리키므로 그 입력만 마스킹됩니다. 나머지도 "
       "덮으려면 텍스트 프롬프트를 추가하세요."),
    DE("Sie beschreiben ein Bild von {0}, also wird nur diese Eingabe "
       "maskiert. Für die übrigen einen Textbefehl ergänzen."),
    FR("Ils décrivent une image de {0} : seule cette entrée est masquée. "
       "Ajoutez une invite textuelle pour couvrir les autres."),
    ES("Describen un fotograma de {0}, así que solo se enmascara esa entrada. "
       "Añada una indicación de texto para cubrir las demás."),
    PT("Eles descrevem um quadro de {0}, então só essa entrada é mascarada. "
       "Adicione um comando de texto para cobrir as outras."),
    IT("Descrivono un fotogramma di {0}, quindi viene mascherato solo quell'"
       "ingresso. Aggiunga un testo per coprire gli altri."),
    NL("Ze beschrijven één beeld van {0}, dus alleen die invoer wordt "
       "gemaskeerd. Voeg een tekstprompt toe om de rest te dekken."),
    RU("Они описывают один кадр входа {0}, поэтому маскируется только он. "
       "Добавьте текстовый запрос, чтобы охватить остальные."),
    TR("Bunlar {0} girdisinin bir karesini betimliyor, dolayısıyla yalnızca o "
       "girdi maskelenir. Diğerlerini de kapsamak için bir metin istemi "
       "ekleyin."));

SS_MSG(mask_pick_input_first,
    EN("Pick the photos or video first."),
    JA("先に写真か動画を選んでください。"),
    ZH_HANS("请先选择照片或视频。"),
    ZH_HANT("請先選擇相片或影片。"),
    KO("먼저 사진이나 동영상을 고르세요."),
    DE("Zuerst die Fotos oder das Video wählen."),
    FR("Choisissez d'abord les photos ou la vidéo."),
    ES("Elija primero las fotos o el vídeo."),
    PT("Escolha primeiro as fotos ou o vídeo."),
    IT("Scelga prima le fotografie o il video."),
    NL("Kies eerst de foto's of de video."),
    RU("Сначала выберите фотографии или видео."),
    TR("Önce fotoğrafları veya videoyu seçin."));

SS_MSG(mask_remove_named,
    EN("Remove what I name"),
    JA("指定したものを取り除く"),
    ZH_HANS("移除我指名的东西"),
    ZH_HANT("移除我指名的東西"),
    KO("내가 말한 것을 제거"),
    DE("Entfernen, was ich nenne"),
    FR("Retirer ce que je nomme"),
    ES("Quitar lo que yo nombre"),
    PT("Remover o que eu nomear"),
    IT("Rimuovere ciò che indico"),
    NL("Verwijderen wat ik noem"),
    RU("Убрать то, что я назову"),
    TR("Adını verdiğimi kaldır"));

SS_MSG(mask_keep_named,
    EN("Keep only what I name"),
    JA("指定したものだけを残す"),
    ZH_HANS("只保留我指名的东西"),
    ZH_HANT("只保留我指名的東西"),
    KO("내가 말한 것만 남기기"),
    DE("Nur behalten, was ich nenne"),
    FR("Ne garder que ce que je nomme"),
    ES("Conservar solo lo que yo nombre"),
    PT("Manter só o que eu nomear"),
    IT("Tenere solo ciò che indico"),
    NL("Alleen houden wat ik noem"),
    RU("Оставить только то, что я назову"),
    TR("Yalnızca adını verdiğimi tut"));

SS_MSG(mask_polarity_help,
    EN("\"Remove\" is for distractors. \"Keep only\" is for capturing a single "
       "object, where everything around it should be ignored."),
    JA("「取り除く」は邪魔物向けです。「だけを残す」は単体の被写体を撮るとき、"
       "まわりのすべてを無視したい場合に使います。"),
    ZH_HANS("“移除”用于干扰物。“只保留”用于拍摄单个物体，此时周围的一切都该被忽略。"),
    ZH_HANT("「移除」用於干擾物。「只保留」用於拍攝單個物體，此時周圍的一切都該被忽略。"),
    KO("'제거'는 방해물용입니다. '만 남기기'는 물체 하나를 촬영할 때, 그 주위의 "
       "모든 것을 무시하고 싶을 때 씁니다."),
    DE("„Entfernen“ ist für Störendes. „Nur behalten“ ist für die Aufnahme "
       "eines einzelnen Objekts, bei der alles ringsum ignoriert werden soll."),
    FR("« Retirer » sert pour les gêneurs. « Ne garder que » sert à "
       "photographier un objet unique, où tout ce qui l'entoure doit être "
       "ignoré."),
    ES("«Quitar» es para elementos molestos. «Conservar solo» es para capturar "
       "un único objeto, cuando todo lo que lo rodea debe ignorarse."),
    PT("“Remover” é para elementos indesejados. “Manter só” é para capturar um "
       "único objeto, quando tudo à volta deve ser ignorado."),
    IT("«Rimuovere» serve per gli elementi di disturbo. «Tenere solo» serve a "
       "riprendere un singolo oggetto, quando tutto ciò che lo circonda va "
       "ignorato."),
    NL("‘Verwijderen’ is voor stoorelementen. ‘Alleen houden’ is voor het "
       "opnemen van één object, waarbij alles eromheen genegeerd moet worden."),
    RU("«Убрать» — для помех. «Оставить только» — для съёмки одного предмета, "
       "когда всё вокруг него нужно игнорировать."),
    TR("“Kaldır” istenmeyen ögeler içindir. “Yalnızca tut” tek bir nesneyi "
       "çekerken, çevresindeki her şeyin yok sayılması gerektiğinde "
       "kullanılır."));

SS_MSG(mask_what_to_keep,
    EN("What to keep"),  JA("残すもの"),      ZH_HANS("要保留什么"), ZH_HANT("要保留什麼"),
    KO("남길 것"),        DE("Was bleiben soll"), FR("Ce qu'il faut garder"),
    ES("Qué conservar"), PT("O que manter"), IT("Che cosa tenere"),
    NL("Wat te houden"), RU("Что оставить"), TR("Ne tutulacak"));

SS_MSG(mask_what_to_remove,
    EN("What to remove"), JA("取り除くもの"),  ZH_HANS("要移除什么"), ZH_HANT("要移除什麼"),
    KO("제거할 것"),      DE("Was weg soll"), FR("Ce qu'il faut retirer"),
    ES("Qué quitar"),    PT("O que remover"), IT("Che cosa rimuovere"),
    NL("Wat te verwijderen"), RU("Что убрать"), TR("Ne kaldırılacak"));

SS_MSG(mask_prompt_help_keep,
    EN("Plain words, separated by semicolons. Everything NOT matching them is "
       "cut out of the reconstruction."),
    JA("ふつうの言葉をセミコロンで区切って書きます。それに当てはまらないものは"
       "すべて再構成から取り除かれます。"),
    ZH_HANS("用平常的词语，以分号分隔。凡是不匹配的都会从重建中剔除。"),
    ZH_HANT("用平常的詞語，以分號分隔。凡是不符合的都會從重建中剔除。"),
    KO("평범한 낱말을 세미콜론으로 구분해 적으세요. 거기 해당하지 않는 것은 모두 "
       "재구성에서 제외됩니다."),
    DE("Einfache Wörter, durch Semikolon getrennt. Alles, was NICHT darauf "
       "passt, fällt aus der Rekonstruktion heraus."),
    FR("Des mots ordinaires, séparés par des points-virgules. Tout ce qui n'y "
       "correspond PAS est retiré de la reconstruction."),
    ES("Palabras corrientes, separadas por punto y coma. Todo lo que NO "
       "coincida se elimina de la reconstrucción."),
    PT("Palavras simples, separadas por ponto e vírgula. Tudo o que NÃO "
       "corresponder é retirado da reconstrução."),
    IT("Parole comuni, separate da punto e virgola. Tutto ciò che NON "
       "corrisponde viene tolto dalla ricostruzione."),
    NL("Gewone woorden, gescheiden door puntkomma's. Alles wat er NIET bij "
       "past valt uit de reconstructie."),
    RU("Обычные слова через точку с запятой. Всё, что им НЕ соответствует, "
       "исключается из реконструкции."),
    TR("Noktalı virgülle ayrılmış sıradan sözcükler. Bunlara uymayan her şey "
       "yeniden oluşturmadan çıkarılır."));

SS_MSG(mask_prompt_help_remove,
    EN("Plain words, separated by semicolons. Everything matching them is cut "
       "out of the reconstruction."),
    JA("ふつうの言葉をセミコロンで区切って書きます。当てはまるものは再構成から"
       "取り除かれます。"),
    ZH_HANS("用平常的词语，以分号分隔。凡是匹配的都会从重建中剔除。"),
    ZH_HANT("用平常的詞語，以分號分隔。凡是符合的都會從重建中剔除。"),
    KO("평범한 낱말을 세미콜론으로 구분해 적으세요. 거기 해당하는 것은 재구성에서 "
       "제외됩니다."),
    DE("Einfache Wörter, durch Semikolon getrennt. Alles, was darauf passt, "
       "fällt aus der Rekonstruktion heraus."),
    FR("Des mots ordinaires, séparés par des points-virgules. Tout ce qui y "
       "correspond est retiré de la reconstruction."),
    ES("Palabras corrientes, separadas por punto y coma. Todo lo que coincida "
       "se elimina de la reconstrucción."),
    PT("Palavras simples, separadas por ponto e vírgula. Tudo o que "
       "corresponder é retirado da reconstrução."),
    IT("Parole comuni, separate da punto e virgola. Tutto ciò che corrisponde "
       "viene tolto dalla ricostruzione."),
    NL("Gewone woorden, gescheiden door puntkomma's. Alles wat erbij past valt "
       "uit de reconstructie."),
    RU("Обычные слова через точку с запятой. Всё, что им соответствует, "
       "исключается из реконструкции."),
    TR("Noktalı virgülle ayrılmış sıradan sözcükler. Bunlara uyan her şey "
       "yeniden oluşturmadan çıkarılır."));

SS_MSG(mask_but_remove,
    EN("...but remove"), JA("…ただし取り除く"), ZH_HANS("…但要移除"), ZH_HANT("…但要移除"),
    KO("…단, 제거할 것"), DE("… aber entfernen"), FR("… mais retirer"),
    ES("… pero quitar"), PT("… mas remover"),  IT("… ma rimuovere"),
    NL("… maar verwijderen"), RU("…но убрать"), TR("…ama kaldır"));

SS_MSG(mask_but_keep,
    EN("...but keep"),   JA("…ただし残す"),   ZH_HANS("…但要保留"), ZH_HANT("…但要保留"),
    KO("…단, 남길 것"),   DE("… aber behalten"), FR("… mais garder"),
    ES("… pero conservar"), PT("… mas manter"), IT("… ma tenere"),
    NL("… maar houden"), RU("…но оставить"), TR("…ama tut"));

SS_MSG(mask_negative_help_keep,
    EN("Exceptions that go even though they match the line above. Optional."),
    JA("上の行に当てはまっても取り除きたい例外です。省略できます。"),
    ZH_HANS("即使符合上一行也要剔除的例外。可以留空。"),
    ZH_HANT("即使符合上一行也要剔除的例外。可以留空。"),
    KO("위 줄에 해당하더라도 제거할 예외입니다. 비워 둬도 됩니다."),
    DE("Ausnahmen, die trotz Treffer in der Zeile darüber wegfallen. "
       "Optional."),
    FR("Exceptions qui partent quand même bien qu'elles correspondent à la "
       "ligne du dessus. Facultatif."),
    ES("Excepciones que se van aunque coincidan con la línea de arriba. "
       "Opcional."),
    PT("Exceções que saem mesmo correspondendo à linha acima. Opcional."),
    IT("Eccezioni che vanno via pur corrispondendo alla riga sopra. "
       "Facoltativo."),
    NL("Uitzonderingen die toch weggaan hoewel ze bij de regel hierboven "
       "passen. Optioneel."),
    RU("Исключения, которые всё же убираются, хотя и подходят под строку выше. "
       "Необязательно."),
    TR("Yukarıdaki satıra uysa bile yine de çıkarılacak istisnalar. İsteğe "
       "bağlı."));

SS_MSG(mask_negative_help_remove,
    EN("Exceptions that stay even though they match the line above. "
       "Optional."),
    JA("上の行に当てはまっても残したい例外です。省略できます。"),
    ZH_HANS("即使符合上一行也要保留的例外。可以留空。"),
    ZH_HANT("即使符合上一行也要保留的例外。可以留空。"),
    KO("위 줄에 해당하더라도 남겨 둘 예외입니다. 비워 둬도 됩니다."),
    DE("Ausnahmen, die trotz Treffer in der Zeile darüber bleiben. Optional."),
    FR("Exceptions qui restent bien qu'elles correspondent à la ligne du "
       "dessus. Facultatif."),
    ES("Excepciones que se quedan aunque coincidan con la línea de arriba. "
       "Opcional."),
    PT("Exceções que ficam mesmo correspondendo à linha acima. Opcional."),
    IT("Eccezioni che restano pur corrispondendo alla riga sopra. "
       "Facoltativo."),
    NL("Uitzonderingen die blijven hoewel ze bij de regel hierboven passen. "
       "Optioneel."),
    RU("Исключения, которые остаются, хотя и подходят под строку выше. "
       "Необязательно."),
    TR("Yukarıdaki satıra uysa bile kalacak istisnalar. İsteğe bağlı."));

SS_MSG(use_found_masks,
    EN("Use the masks found next to the photos"),
    JA("写真のとなりで見つかったマスクを使う"),
    ZH_HANS("使用照片旁边找到的蒙版"),
    ZH_HANT("使用照片旁邊找到的遮罩"),
    KO("사진 옆에서 찾은 마스크 사용"),
    DE("Neben den Fotos gefundene Masken verwenden"),
    FR("Utiliser les masques trouvés à côté des photos"),
    ES("Usar las máscaras encontradas junto a las fotos"),
    PT("Usar as máscaras encontradas ao lado das fotos"),
    IT("Usare le maschere trovate accanto alle foto"),
    NL("De maskers gebruiken die naast de foto's staan"),
    RU("Использовать маски, найденные рядом с фотографиями"),
    TR("Fotoğrafların yanında bulunan maskeleri kullan"));

SS_MSG(use_found_masks_help,
    EN("A `masks` folder beside or under the photos is picked up on its own, "
       "because that is where a prepared capture keeps them. Turn this off "
       "for a dataset whose `masks` folder belongs to something else -- a "
       "different set of views, or an experiment you are not reproducing."),
    JA("写真のとなりや下にある `masks` フォルダーは自動で拾います。用意済みの"
       "データセットはそこに置くからです。その `masks` が別のもの――別の視点の"
       "集まりや、いま再現しようとしていない実験――に属している場合は、"
       "これを外してください。"),
    ZH_HANS("照片旁边或下面的 `masks` 文件夹会被自动采用，因为准备好的数据集"
            "就放在那里。如果那个 `masks` 属于别的东西——另一组视角，或者你"
            "并不打算复现的实验——请关掉这一项。"),
    ZH_HANT("照片旁邊或下面的 `masks` 資料夾會被自動採用，因為準備好的資料集"
            "就放在那裡。如果那個 `masks` 屬於別的東西——另一組視角，或者你"
            "並不打算重現的實驗——請關掉這一項。"),
    KO("사진 옆이나 아래의 `masks` 폴더는 자동으로 사용됩니다. 준비된 촬영본이 "
       "마스크를 두는 자리이기 때문입니다. 그 `masks`가 다른 것에 속한다면"
       "(다른 시점 묶음이거나, 지금 재현하려는 것이 아닌 실험이라면) 이 "
       "항목을 끄세요."),
    DE("Ein `masks`-Ordner neben oder unter den Fotos wird von selbst "
       "übernommen -- dort legt eine vorbereitete Aufnahme sie ab. Schalten "
       "Sie das aus, wenn dieser Ordner zu etwas anderem gehört: zu einem "
       "anderen Satz Ansichten oder zu einem Versuch, den Sie nicht "
       "nachstellen."),
    FR("Un dossier `masks` à côté ou sous les photos est repris tout seul, "
       "car c'est là qu'une prise de vue préparée les garde. Décochez pour un "
       "jeu de données dont le dossier `masks` appartient à autre chose : un "
       "autre ensemble de vues, ou une expérience que vous ne reproduisez "
       "pas."),
    ES("Una carpeta `masks` junto a las fotos o debajo de ellas se toma sola, "
       "porque ahí es donde las guarda una captura preparada. Desactívelo si "
       "esa carpeta pertenece a otra cosa: a otro conjunto de vistas, o a un "
       "experimento que no está reproduciendo."),
    PT("Uma pasta `masks` ao lado das fotos ou abaixo delas é adotada "
       "sozinha, porque é ali que uma captura preparada as guarda. Desligue "
       "isto se essa pasta pertencer a outra coisa: a outro conjunto de "
       "vistas, ou a um experimento que você não está reproduzindo."),
    IT("Una cartella `masks` accanto alle foto o sotto di esse viene presa da "
       "sola, perché è lì che una ripresa preparata le tiene. Lo disattivi se "
       "quella cartella appartiene ad altro: a un altro insieme di viste, o a "
       "un esperimento che non sta riproducendo."),
    NL("Een map `masks` naast of onder de foto's wordt vanzelf overgenomen -- "
       "daar bewaart een voorbereide opname ze. Zet dit uit voor een dataset "
       "waarvan die map bij iets anders hoort: een andere reeks aanzichten, "
       "of een proef die u niet nabootst."),
    RU("Папка `masks` рядом с фотографиями или под ними подхватывается сама — "
       "именно там подготовленная съёмка их держит. Снимите галочку, если эта "
       "папка относится к чему-то другому: к другому набору видов или к "
       "опыту, который вы не повторяете."),
    TR("Fotoğrafların yanındaki ya da altındaki bir `masks` klasörü kendinden "
       "alınır; hazırlanmış bir çekim maskeleri orada tutar. O klasör başka "
       "bir şeye aitse -- başka bir görünüm kümesine ya da yeniden "
       "üretmediğiniz bir denemeye -- bunu kapatın."));

// ---------------------------------------------------------------------------
// Lens models offered by the built-in reconstruction (SfmRunner.h).
//
// The model NAMES -- OpenCV, Kannala-Brandt, thin prism -- are the names of
// the distortion models themselves and stay put; what is translated is the
// parenthetical that says which one to pick.
// ---------------------------------------------------------------------------

SS_MSG(lens_opencv,
    EN("OpenCV (most cameras)"),
    JA("OpenCV（ほとんどのカメラ）"),
    ZH_HANS("OpenCV（多数相机）"),
    ZH_HANT("OpenCV（多數相機）"),
    KO("OpenCV(대부분의 카메라)"),
    DE("OpenCV (die meisten Kameras)"),
    FR("OpenCV (la plupart des appareils)"),
    ES("OpenCV (la mayoría de cámaras)"),
    PT("OpenCV (a maioria das câmeras)"),
    IT("OpenCV (quasi tutte le fotocamere)"),
    NL("OpenCV (de meeste camera's)"),
    RU("OpenCV (большинство камер)"),
    TR("OpenCV (çoğu kamera)"));

SS_MSG(lens_pinhole,
    EN("Pinhole (no distortion)"),
    JA("ピンホール（歪みなし）"),
    ZH_HANS("针孔（无畸变）"),
    ZH_HANT("針孔（無變形）"),
    KO("핀홀(왜곡 없음)"),
    DE("Lochkamera (ohne Verzeichnung)"),
    FR("Sténopé (sans distorsion)"),
    ES("Estenopeica (sin distorsión)"),
    PT("Estenopeica (sem distorção)"),
    IT("Stenopeica (senza distorsione)"),
    NL("Gaatjescamera (zonder vervorming)"),
    RU("Точечная камера (без дисторсии)"),
    TR("İğne deliği (bozulmasız)"));

SS_MSG(lens_simple_pinhole,
    EN("Simple pinhole"),
    JA("簡易ピンホール"),
    ZH_HANS("简化针孔"),
    ZH_HANT("簡化針孔"),
    KO("간단 핀홀"),
    DE("Einfache Lochkamera"),
    FR("Sténopé simple"),
    ES("Estenopeica simple"),
    PT("Estenopeica simples"),
    IT("Stenopeica semplice"),
    NL("Eenvoudige gaatjescamera"),
    RU("Простая точечная камера"),
    TR("Basit iğne deliği"));

SS_MSG(lens_radial,
    EN("Radial"),
    JA("放射方向のみ"),
    ZH_HANS("仅径向畸变"),
    ZH_HANT("僅徑向變形"),
    KO("방사 왜곡만"),
    DE("Radial"),
    FR("Radiale"),
    ES("Radial"),
    PT("Radial"),
    IT("Radiale"),
    NL("Radiaal"),
    RU("Радиальная"),
    TR("Işınsal"));

SS_MSG(lens_full_opencv,
    EN("Full OpenCV"),
    JA("OpenCV（全パラメータ）"),
    ZH_HANS("OpenCV（完整参数）"),
    ZH_HANT("OpenCV（完整參數）"),
    KO("OpenCV(전체 계수)"),
    DE("OpenCV (vollständig)"),
    FR("OpenCV complet"),
    ES("OpenCV completo"),
    PT("OpenCV completo"),
    IT("OpenCV completo"),
    NL("Volledige OpenCV"),
    RU("OpenCV (полная модель)"),
    TR("Tam OpenCV"));

SS_MSG(lens_fisheye_kb,
    EN("Fisheye (Kannala-Brandt)"),
    JA("魚眼（Kannala-Brandt）"),
    ZH_HANS("鱼眼（Kannala-Brandt）"),
    ZH_HANT("魚眼（Kannala-Brandt）"),
    KO("어안(Kannala-Brandt)"),
    DE("Fisheye (Kannala-Brandt)"),
    FR("Fisheye (Kannala-Brandt)"),
    ES("Ojo de pez (Kannala-Brandt)"),
    PT("Olho de peixe (Kannala-Brandt)"),
    IT("Fisheye (Kannala-Brandt)"),
    NL("Fisheye (Kannala-Brandt)"),
    RU("Фишай (Каннала — Брандт)"),
    TR("Balıkgözü (Kannala-Brandt)"));

SS_MSG(lens_fisheye_thin_prism,
    EN("Fisheye (thin prism)"),
    JA("魚眼（薄プリズム）"),
    ZH_HANS("鱼眼（薄棱镜）"),
    ZH_HANT("魚眼（薄稜鏡）"),
    KO("어안(얇은 프리즘)"),
    DE("Fisheye (dünnes Prisma)"),
    FR("Fisheye (prisme mince)"),
    ES("Ojo de pez (prisma delgado)"),
    PT("Olho de peixe (prisma fino)"),
    IT("Fisheye (prisma sottile)"),
    NL("Fisheye (dun prisma)"),
    RU("Фишай (тонкая призма)"),
    TR("Balıkgözü (ince prizma)"));

SS_MSG(lens_equirectangular,
    EN("Equirectangular (360\xc2\xb0 panorama)"),
    JA("正距円筒（360度パノラマ）"),
    ZH_HANS("等距柱状（360 度全景）"),
    ZH_HANT("等距柱狀（360 度全景）"),
    KO("정거원통(360도 파노라마)"),
    DE("Äquirektangulär (360\xc2\xb0-Panorama)"),
    FR("Équirectangulaire (panorama 360\xc2\xb0)"),
    ES("Equirrectangular (panorama 360\xc2\xb0)"),
    PT("Equirretangular (panorama 360\xc2\xb0)"),
    IT("Equirettangolare (panorama 360\xc2\xb0)"),
    NL("Equirectangulair (360\xc2\xb0-panorama)"),
    RU("Эквиректангулярная (панорама 360\xc2\xb0)"),
    TR("Eş dikdörtgen (360\xc2\xb0 panorama)"));

// ---------------------------------------------------------------------------
// The common-subject palette (src/app/gui/MaskPrompt.h)
//
// The chip LABELS below are translated; the words they put in the box are
// not, and must not be. SAM 3's text encoder is trained on English, and a
// prompt in anything else finds noticeably less of what it names -- so the
// box stays English however the interface is set, and this palette is how a
// user who does not write English still gets a good prompt.
//
// These are label translations, not dictionary entries: each one should be
// what a speaker would call the thing IN A PHOTOGRAPH, not the closest
// dictionary word to the English.
// ---------------------------------------------------------------------------

SS_MSG(mask_english_only,
    EN("The model reads English only, so this box stays English whatever the "
       "interface language is. Pick from the list below, or type English "
       "words."),
    JA("モデルは英語しか読み取れないため、この欄は表示言語に関わらず英語のまま"
       "です。下の一覧から選ぶか、英語で入力してください。"),
    ZH_HANS("模型只能读英文，所以无论界面用哪种语言，这个输入框都保持英文。"
            "可以从下面的列表里选，或者直接输入英文词。"),
    ZH_HANT("模型只能讀英文，所以無論介面用哪種語言，這個輸入框都維持英文。"
            "可以從下面的清單挑選，或直接輸入英文詞。"),
    KO("모델은 영어만 읽으므로 인터페이스 언어와 관계없이 이 칸은 영어로 "
       "유지됩니다. 아래 목록에서 고르거나 영어 낱말을 입력하세요."),
    DE("Das Modell liest nur Englisch, deshalb bleibt dieses Feld englisch, "
       "unabhängig von der Sprache der Oberfläche. Wählen Sie aus der Liste "
       "unten, oder tippen Sie englische Wörter."),
    FR("Le modèle ne lit que l'anglais : ce champ reste donc en anglais quelle "
       "que soit la langue de l'interface. Choisissez dans la liste ci-dessous, "
       "ou tapez des mots anglais."),
    ES("El modelo solo lee inglés, así que este campo sigue en inglés sea cual "
       "sea el idioma de la interfaz. Elija de la lista de abajo, o escriba "
       "palabras en inglés."),
    PT("O modelo só lê inglês, por isso este campo continua em inglês seja qual "
       "for o idioma da interface. Escolha na lista abaixo, ou digite palavras "
       "em inglês."),
    IT("Il modello legge solo l'inglese, quindi questo campo resta in inglese "
       "qualunque sia la lingua dell'interfaccia. Scelga dall'elenco qui sotto, "
       "oppure scriva parole inglesi."),
    NL("Het model leest alleen Engels, dus dit veld blijft Engels ongeacht de "
       "taal van de interface. Kies uit de lijst hieronder, of typ Engelse "
       "woorden."),
    RU("Модель читает только по-английски, поэтому это поле остаётся "
       "английским при любом языке интерфейса. Выберите из списка ниже или "
       "введите английские слова."),
    TR("Model yalnızca İngilizce okur, bu yüzden arayüz dili ne olursa olsun "
       "bu kutu İngilizce kalır. Aşağıdaki listeden seçin veya İngilizce "
       "sözcük yazın."));

SS_MSG(mask_subjects,
    EN("Common subjects"),
    JA("よく使う対象"),
    ZH_HANS("常见对象"),
    ZH_HANT("常見對象"),
    KO("자주 쓰는 대상"),
    DE("Häufige Motive"),
    FR("Sujets courants"),
    ES("Sujetos habituales"),
    PT("Assuntos comuns"),
    IT("Soggetti comuni"),
    NL("Veelgebruikte onderwerpen"),
    RU("Частые объекты"),
    TR("Sık kullanılan konular"));

SS_MSG(mask_subjects_help,
    EN("Click to put the English word in the box above; click again to take it "
       "out. A highlighted chip is already in the box."),
    JA("押すと上の欄に英語の語が入り、もう一度押すと外れます。色が付いている"
       "ものは既に入っています。"),
    ZH_HANS("点一下把这个英文词放进上面的框，再点一下取出。高亮的表示已经在框里。"),
    ZH_HANT("按一下把這個英文詞放進上面的框，再按一下取出。標亮的表示已經在框裡。"),
    KO("누르면 위 칸에 영어 낱말이 들어가고, 다시 누르면 빠집니다. 강조된 것은 "
       "이미 들어 있는 것입니다."),
    DE("Klicken setzt das englische Wort in das Feld oben, nochmals klicken "
       "nimmt es wieder heraus. Hervorgehobene stehen bereits darin."),
    FR("Cliquez pour mettre le mot anglais dans le champ ci-dessus ; cliquez à "
       "nouveau pour l'enlever. Les éléments en surbrillance y sont déjà."),
    ES("Pulse para poner la palabra inglesa en el campo de arriba; púlsela otra "
       "vez para quitarla. Las resaltadas ya están dentro."),
    PT("Clique para pôr a palavra em inglês no campo acima; clique de novo para "
       "tirá-la. As destacadas já estão lá."),
    IT("Clicchi per mettere la parola inglese nel campo qui sopra; clicchi di "
       "nuovo per toglierla. Quelle evidenziate ci sono già."),
    NL("Klik om het Engelse woord in het veld hierboven te zetten; klik "
       "nogmaals om het eruit te halen. Gemarkeerde staan er al in."),
    RU("Нажмите, чтобы поставить английское слово в поле выше; нажмите ещё "
       "раз, чтобы убрать. Выделенные уже стоят в поле."),
    TR("Tıklayınca İngilizce sözcük yukarıdaki kutuya girer, yeniden "
       "tıklayınca çıkar. Vurgulu olanlar zaten kutudadır."));

SS_MSG(subj_person,
    EN("Person"),       JA("人物"),        ZH_HANS("人"),      ZH_HANT("人"),
    KO("사람"),          DE("Person"),     FR("Personne"),    ES("Persona"),
    PT("Pessoa"),       IT("Persona"),    NL("Persoon"),     RU("Человек"),
    TR("İnsan"));

SS_MSG(subj_hand,
    EN("Hand"),         JA("手"),          ZH_HANS("手"),      ZH_HANT("手"),
    KO("손"),            DE("Hand"),       FR("Main"),        ES("Mano"),
    PT("Mão"),          IT("Mano"),       NL("Hand"),        RU("Рука"),
    TR("El"));

SS_MSG(subj_dog,
    EN("Dog"),          JA("犬"),          ZH_HANS("狗"),      ZH_HANT("狗"),
    KO("개"),            DE("Hund"),       FR("Chien"),       ES("Perro"),
    PT("Cachorro"),     IT("Cane"),       NL("Hond"),        RU("Собака"),
    TR("Köpek"));

SS_MSG(subj_animal,
    EN("Animal"),       JA("動物"),        ZH_HANS("动物"),     ZH_HANT("動物"),
    KO("동물"),          DE("Tier"),       FR("Animal"),      ES("Animal"),
    PT("Animal"),       IT("Animale"),    NL("Dier"),        RU("Животное"),
    TR("Hayvan"));

SS_MSG(subj_car,
    EN("Car"),          JA("車"),          ZH_HANS("汽车"),     ZH_HANT("汽車"),
    KO("자동차"),        DE("Auto"),       FR("Voiture"),     ES("Coche"),
    PT("Carro"),        IT("Automobile"), NL("Auto"),        RU("Машина"),
    TR("Araba"));

SS_MSG(subj_bicycle,
    EN("Bicycle"),      JA("自転車"),      ZH_HANS("自行车"),   ZH_HANT("腳踏車"),
    KO("자전거"),        DE("Fahrrad"),    FR("Vélo"),        ES("Bicicleta"),
    PT("Bicicleta"),    IT("Bicicletta"), NL("Fiets"),       RU("Велосипед"),
    TR("Bisiklet"));

SS_MSG(subj_vehicle,
    EN("Vehicle"),      JA("乗り物"),      ZH_HANS("车辆"),     ZH_HANT("車輛"),
    KO("차량"),          DE("Fahrzeug"),   FR("Véhicule"),    ES("Vehículo"),
    PT("Veículo"),      IT("Veicolo"),    NL("Voertuig"),    RU("Транспорт"),
    TR("Araç"));

SS_MSG(subj_sky,
    EN("Sky"),          JA("空"),          ZH_HANS("天空"),     ZH_HANT("天空"),
    KO("하늘"),          DE("Himmel"),     FR("Ciel"),        ES("Cielo"),
    PT("Céu"),          IT("Cielo"),      NL("Lucht"),       RU("Небо"),
    TR("Gökyüzü"));

SS_MSG(subj_shadow,
    EN("Shadow"),       JA("影"),          ZH_HANS("影子"),     ZH_HANT("影子"),
    KO("그림자"),        DE("Schatten"),   FR("Ombre"),       ES("Sombra"),
    PT("Sombra"),       IT("Ombra"),      NL("Schaduw"),     RU("Тень"),
    TR("Gölge"));

SS_MSG(subj_water,
    EN("Water"),        JA("水面"),        ZH_HANS("水面"),     ZH_HANT("水面"),
    KO("물"),            DE("Wasser"),     FR("Eau"),         ES("Agua"),
    PT("Água"),         IT("Acqua"),      NL("Water"),       RU("Вода"),
    TR("Su"));

SS_MSG(subj_reflection,
    EN("Reflection"),   JA("映り込み"),    ZH_HANS("倒影"),     ZH_HANT("倒影"),
    KO("반사"),          DE("Spiegelung"), FR("Reflet"),      ES("Reflejo"),
    PT("Reflexo"),      IT("Riflesso"),   NL("Weerspiegeling"), RU("Отражение"),
    TR("Yansıma"));

SS_MSG(subj_tripod,
    EN("Tripod"),       JA("三脚"),        ZH_HANS("三脚架"),   ZH_HANT("三腳架"),
    KO("삼각대"),        DE("Stativ"),     FR("Trépied"),     ES("Trípode"),
    PT("Tripé"),        IT("Treppiede"),  NL("Statief"),     RU("Штатив"),
    TR("Tripod"));

// The black area outside the image circle of a 360 / fisheye camera. Worth a
// chip of its own: it is the one thing on this list that is not part of the
// scene at all, and every Insta360 or fisheye capture has it.
SS_MSG(subj_fisheye_border,
    EN("Fisheye border"),
    JA("魚眼の黒枠"),
    ZH_HANS("鱼眼黑边"),
    ZH_HANT("魚眼黑邊"),
    KO("어안 검은 테두리"),
    DE("Fisheye-Rand"),
    FR("Bord du fisheye"),
    ES("Borde de ojo de pez"),
    PT("Borda olho-de-peixe"),
    IT("Bordo fisheye"),
    NL("Fisheye-rand"),
    RU("Чёрный край фишая"),
    TR("Balıkgözü kenarı"));

SS_MSG(subj_watermark,
    EN("Watermark or timestamp"),
    JA("透かし・日時表示"),
    ZH_HANS("水印或时间戳"),
    ZH_HANT("浮水印或時間戳"),
    KO("워터마크나 날짜 표시"),
    DE("Wasserzeichen oder Zeitstempel"),
    FR("Filigrane ou horodatage"),
    ES("Marca de agua o fecha"),
    PT("Marca d'água ou data"),
    IT("Filigrana o data"),
    NL("Watermerk of tijdstempel"),
    RU("Водяной знак или дата"),
    TR("Filigran veya zaman damgası"));

SS_MSG(subj_person_painting,
    EN("Person in a painting"),
    JA("絵の中の人物"),
    ZH_HANS("画里的人"),
    ZH_HANT("畫裡的人"),
    KO("그림 속 인물"),
    DE("Person auf einem Gemälde"),
    FR("Personnage dans un tableau"),
    ES("Persona en un cuadro"),
    PT("Pessoa num quadro"),
    IT("Persona in un dipinto"),
    NL("Persoon op een schilderij"),
    RU("Человек на картине"),
    TR("Tablodaki insan"));

SS_MSG(subj_statue,
    EN("Statue of a person"),
    JA("人物の彫像"),
    ZH_HANS("人像雕塑"),
    ZH_HANT("人像雕塑"),
    KO("인물 조각상"),
    DE("Statue eines Menschen"),
    FR("Statue de personne"),
    ES("Estatua de una persona"),
    PT("Estátua de pessoa"),
    IT("Statua di persona"),
    NL("Standbeeld van een persoon"),
    RU("Статуя человека"),
    TR("İnsan heykeli"));

SS_MSG(subj_mannequin,
    EN("Mannequin"),
    JA("マネキン"),
    ZH_HANS("人体模特"),
    ZH_HANT("人體模特"),
    KO("마네킹"),
    DE("Schaufensterpuppe"),
    FR("Mannequin"),
    ES("Maniquí"),
    PT("Manequim"),
    IT("Manichino"),
    NL("Etalagepop"),
    RU("Манекен"),
    TR("Manken"));

// ===========================================================================
// Mask preview window
// ===========================================================================

SS_MSG(preview_title,
    EN("Try the mask"),  JA("マスクを試す"),   ZH_HANS("试一下蒙版"), ZH_HANT("試一下遮罩"),
    KO("마스크 시험해 보기"), DE("Maske ausprobieren"), FR("Essayer le masque"),
    ES("Probar la máscara"), PT("Testar a máscara"), IT("Prova la maschera"),
    NL("Masker uitproberen"), RU("Проверить маску"), TR("Maskeyi dene"));

SS_MSG(preview_legend,
    EN("Red = removed from the reconstruction. Adjust the prompt until only "
       "what you want gone is red."),
    JA("赤は再構成から取り除かれる部分です。消したいものだけが赤くなるまで"
       "プロンプトを調整してください。"),
    ZH_HANS("红色 = 会从重建中剔除。请调整提示词，直到只有你想去掉的东西是红色。"),
    ZH_HANT("紅色 = 會從重建中剔除。請調整提示詞，直到只有你想去掉的東西是紅色。"),
    KO("빨강 = 재구성에서 제외됩니다. 없애고 싶은 것만 빨갛게 될 때까지 프롬프트를 "
       "다듬으세요."),
    DE("Rot = fällt aus der Rekonstruktion heraus. Den Text so lange anpassen, "
       "bis nur noch rot ist, was verschwinden soll."),
    FR("Rouge = retiré de la reconstruction. Ajustez l'invite jusqu'à ce que "
       "seul ce que vous voulez supprimer soit rouge."),
    ES("Rojo = se elimina de la reconstrucción. Ajuste la indicación hasta que "
       "solo esté en rojo lo que quiere quitar."),
    PT("Vermelho = removido da reconstrução. Ajuste o comando até que só o que "
       "você quer tirar fique vermelho."),
    IT("Rosso = tolto dalla ricostruzione. Regoli il testo finché è rosso solo "
       "ciò che vuole eliminare."),
    NL("Rood = valt uit de reconstructie. Pas de prompt aan tot alleen rood is "
       "wat je weg wilt hebben."),
    RU("Красное — то, что уйдёт из реконструкции. Правьте запрос, пока красным "
       "не останется только лишнее."),
    TR("Kırmızı = yeniden oluşturmadan çıkarılır. Yalnızca gitmesini "
       "istedikleriniz kırmızı olana dek istemi ayarlayın."));

SS_MSG(preview_polarity_help,
    EN("\"Remove\" is for distractors -- people, cars, the photographer's "
       "shadow. \"Keep only\" is for object captures, where everything but the "
       "subject should be ignored."),
    JA("「取り除く」は邪魔物、たとえば人、車、撮影者の影に使います。"
       "「だけを残す」は物体の撮影用で、被写体以外をすべて無視したいときに"
       "使います。"),
    ZH_HANS("“移除”用于干扰物——行人、汽车、摄影者的影子。“只保留”用于物体拍摄，"
            "此时除主体外的一切都该被忽略。"),
    ZH_HANT("「移除」用於干擾物——行人、汽車、攝影者的影子。「只保留」用於物體拍攝，"
            "此時除主體外的一切都該被忽略。"),
    KO("'제거'는 지나가는 사람, 자동차, 촬영자의 그림자 같은 방해물용입니다. "
       "'만 남기기'는 물체 촬영용으로, 피사체 말고는 전부 무시하고 싶을 때 "
       "씁니다."),
    DE("„Entfernen“ ist für Störendes -- Passanten, Autos, der eigene "
       "Schatten. „Nur behalten“ ist für Objektaufnahmen, bei denen alles "
       "außer dem Motiv ignoriert werden soll."),
    FR("« Retirer » sert pour les gêneurs -- passants, voitures, l'ombre du "
       "photographe. « Ne garder que » sert aux prises d'objet, où tout sauf "
       "le sujet doit être ignoré."),
    ES("«Quitar» es para elementos molestos: transeúntes, coches, la sombra "
       "del fotógrafo. «Conservar solo» es para capturas de objetos, donde "
       "todo salvo el sujeto debe ignorarse."),
    PT("“Remover” é para elementos indesejados: pessoas passando, carros, a "
       "sombra do fotógrafo. “Manter só” é para capturas de objetos, em que "
       "tudo além do sujeito deve ser ignorado."),
    IT("«Rimuovere» serve per gli elementi di disturbo: passanti, automobili, "
       "l'ombra del fotografo. «Tenere solo» serve alle riprese di oggetti, "
       "dove tutto tranne il soggetto va ignorato."),
    NL("‘Verwijderen’ is voor stoorelementen -- voorbijgangers, auto's, de "
       "schaduw van de fotograaf. ‘Alleen houden’ is voor objectopnamen, "
       "waarbij alles behalve het onderwerp genegeerd moet worden."),
    RU("«Убрать» — для помех: прохожих, машин, тени фотографа. «Оставить "
       "только» — для съёмки предметов, когда всё, кроме объекта, нужно "
       "игнорировать."),
    TR("“Kaldır” istenmeyenler içindir -- geçen insanlar, arabalar, "
       "fotoğrafçının gölgesi. “Yalnızca tut” nesne çekimleri içindir; orada "
       "özne dışındaki her şey yok sayılmalıdır."));

SS_MSG(preview_what_kept,
    EN("What should be kept?"),
    JA("何を残しますか？"),
    ZH_HANS("要保留什么？"),
    ZH_HANT("要保留什麼？"),
    KO("무엇을 남길까요?"),
    DE("Was soll bleiben?"),
    FR("Que faut-il garder ?"),
    ES("¿Qué se debe conservar?"),
    PT("O que deve ser mantido?"),
    IT("Che cosa va tenuto?"),
    NL("Wat moet blijven?"),
    RU("Что оставить?"),
    TR("Ne tutulsun?"));

SS_MSG(preview_what_removed,
    EN("What should be removed?"),
    JA("何を取り除きますか？"),
    ZH_HANS("要移除什么？"),
    ZH_HANT("要移除什麼？"),
    KO("무엇을 제거할까요?"),
    DE("Was soll weg?"),
    FR("Que faut-il retirer ?"),
    ES("¿Qué se debe quitar?"),
    PT("O que deve ser removido?"),
    IT("Che cosa va rimosso?"),
    NL("Wat moet weg?"),
    RU("Что убрать?"),
    TR("Ne kaldırılsın?"));

SS_MSG(preview_prompt_help_keep,
    EN("Plain words for the subject of the capture, separated by semicolons. "
       "Everything else is cut out of the reconstruction."),
    JA("撮影の被写体をふつうの言葉でセミコロン区切りに書きます。それ以外は"
       "すべて再構成から取り除かれます。"),
    ZH_HANS("用平常的词语写出拍摄的主体，以分号分隔。其余的一切都会从重建中剔除。"),
    ZH_HANT("用平常的詞語寫出拍攝的主體，以分號分隔。其餘的一切都會從重建中剔除。"),
    KO("촬영의 피사체를 평범한 낱말로 세미콜론으로 구분해 적으세요. 그 밖의 "
       "모든 것은 재구성에서 제외됩니다."),
    DE("Einfache Wörter für das Motiv der Aufnahme, durch Semikolon getrennt. "
       "Alles andere fällt aus der Rekonstruktion heraus."),
    FR("Des mots ordinaires pour le sujet de la prise, séparés par des "
       "points-virgules. Tout le reste est retiré de la reconstruction."),
    ES("Palabras corrientes para el sujeto de la captura, separadas por punto "
       "y coma. Todo lo demás se elimina de la reconstrucción."),
    PT("Palavras simples para o sujeito da captura, separadas por ponto e "
       "vírgula. Todo o resto é retirado da reconstrução."),
    IT("Parole comuni per il soggetto della ripresa, separate da punto e "
       "virgola. Tutto il resto viene tolto dalla ricostruzione."),
    NL("Gewone woorden voor het onderwerp van de opname, gescheiden door "
       "puntkomma's. Al het andere valt uit de reconstructie."),
    RU("Обычные слова для объекта съёмки через точку с запятой. Всё остальное "
       "исключается из реконструкции."),
    TR("Çekimin öznesi için noktalı virgülle ayrılmış sıradan sözcükler. Geri "
       "kalan her şey yeniden oluşturmadan çıkarılır."));

SS_MSG(preview_prompt_help_remove,
    EN("Plain words for the things to take out of the reconstruction, "
       "separated by semicolons. Anything that moved, reflected, or was not "
       "part of the scene is a good candidate."),
    JA("再構成から外したいものをふつうの言葉でセミコロン区切りに書きます。"
       "動いたもの、映り込んだもの、シーンの一部ではなかったものが候補です。"),
    ZH_HANS("用平常的词语写出要从重建中拿掉的东西，以分号分隔。凡是移动过的、"
            "有反射的、本来就不属于场景的，都是合适的候选。"),
    ZH_HANT("用平常的詞語寫出要從重建中拿掉的東西，以分號分隔。凡是移動過的、"
            "有反射的、本來就不屬於場景的，都是合適的候選。"),
    KO("재구성에서 빼고 싶은 것들을 평범한 낱말로 세미콜론으로 구분해 적으세요. "
       "움직였던 것, 비쳤던 것, 원래 장면의 일부가 아니었던 것이 좋은 "
       "후보입니다."),
    DE("Einfache Wörter für das, was aus der Rekonstruktion heraus soll, durch "
       "Semikolon getrennt. Alles, was sich bewegte, spiegelte oder nicht zur "
       "Szene gehörte, ist ein guter Kandidat."),
    FR("Des mots ordinaires pour ce qu'il faut sortir de la reconstruction, "
       "séparés par des points-virgules. Tout ce qui bougeait, se reflétait ou "
       "ne faisait pas partie de la scène est un bon candidat."),
    ES("Palabras corrientes para lo que hay que sacar de la reconstrucción, "
       "separadas por punto y coma. Todo lo que se movía, se reflejaba o no "
       "formaba parte de la escena es buen candidato."),
    PT("Palavras simples para o que deve sair da reconstrução, separadas por "
       "ponto e vírgula. Tudo o que se movia, refletia ou não fazia parte da "
       "cena é bom candidato."),
    IT("Parole comuni per ciò che va tolto dalla ricostruzione, separate da "
       "punto e virgola. Tutto ciò che si muoveva, si rifletteva o non faceva "
       "parte della scena è un buon candidato."),
    NL("Gewone woorden voor wat uit de reconstructie moet, gescheiden door "
       "puntkomma's. Alles wat bewoog, weerspiegelde of geen deel van de scène "
       "was, is een goede kandidaat."),
    RU("Обычные слова для того, что нужно вынести из реконструкции, через "
       "точку с запятой. Хорошие кандидаты — всё, что двигалось, отражалось "
       "или не относилось к сцене."),
    TR("Yeniden oluşturmadan çıkarılacak şeyler için noktalı virgülle ayrılmış "
       "sıradan sözcükler. Hareket eden, yansıyan ya da sahnenin parçası "
       "olmayan her şey iyi bir adaydır."));

SS_MSG(preview_but_remove_these,
    EN("...but remove these"),
    JA("…ただしこれらは取り除く"),
    ZH_HANS("…但要移除这些"),
    ZH_HANT("…但要移除這些"),
    KO("…단, 이것들은 제거"),
    DE("… aber diese entfernen"),
    FR("… mais retirer ceci"),
    ES("… pero quitar estos"),
    PT("… mas remover estes"),
    IT("… ma rimuovere questi"),
    NL("… maar deze verwijderen"),
    RU("…но эти убрать"),
    TR("…ama bunları kaldır"));

SS_MSG(preview_but_keep_these,
    EN("...but keep these"),
    JA("…ただしこれらは残す"),
    ZH_HANS("…但要保留这些"),
    ZH_HANT("…但要保留這些"),
    KO("…단, 이것들은 남기기"),
    DE("… aber diese behalten"),
    FR("… mais garder ceci"),
    ES("… pero conservar estos"),
    PT("… mas manter estes"),
    IT("… ma tenere questi"),
    NL("… maar deze houden"),
    RU("…но эти оставить"),
    TR("…ama bunları tut"));

SS_MSG(preview_negative_help_keep,
    EN("Exceptions: things that match the line above but should still go. "
       "Optional."),
    JA("例外です。上の行に当てはまるけれど、それでも取り除きたいものを書きます。"
       "省略できます。"),
    ZH_HANS("例外：符合上一行、但仍然应当剔除的东西。可以留空。"),
    ZH_HANT("例外：符合上一行、但仍然應當剔除的東西。可以留空。"),
    KO("예외입니다. 위 줄에 해당하지만 그래도 없애야 할 것들. 비워 둬도 됩니다."),
    DE("Ausnahmen: Dinge, die auf die Zeile darüber passen, aber trotzdem weg "
       "sollen. Optional."),
    FR("Exceptions : ce qui correspond à la ligne du dessus mais doit quand "
       "même partir. Facultatif."),
    ES("Excepciones: cosas que coinciden con la línea de arriba pero deben "
       "irse igual. Opcional."),
    PT("Exceções: coisas que correspondem à linha acima mas devem sair mesmo "
       "assim. Opcional."),
    IT("Eccezioni: cose che corrispondono alla riga sopra ma devono comunque "
       "sparire. Facoltativo."),
    NL("Uitzonderingen: dingen die bij de regel hierboven passen maar toch weg "
       "moeten. Optioneel."),
    RU("Исключения: то, что подходит под строку выше, но всё же должно уйти. "
       "Необязательно."),
    TR("İstisnalar: yukarıdaki satıra uyan ama yine de gitmesi gerekenler. "
       "İsteğe bağlı."));

SS_MSG(preview_negative_help_remove,
    EN("Exceptions: things that match the line above but should stay. "
       "Optional."),
    JA("例外です。上の行に当てはまるけれど、残したいものを書きます。"
       "省略できます。"),
    ZH_HANS("例外：符合上一行、但应当保留的东西。可以留空。"),
    ZH_HANT("例外：符合上一行、但應當保留的東西。可以留空。"),
    KO("예외입니다. 위 줄에 해당하지만 남겨야 할 것들. 비워 둬도 됩니다."),
    DE("Ausnahmen: Dinge, die auf die Zeile darüber passen, aber bleiben "
       "sollen. Optional."),
    FR("Exceptions : ce qui correspond à la ligne du dessus mais doit rester. "
       "Facultatif."),
    ES("Excepciones: cosas que coinciden con la línea de arriba pero deben "
       "quedarse. Opcional."),
    PT("Exceções: coisas que correspondem à linha acima mas devem ficar. "
       "Opcional."),
    IT("Eccezioni: cose che corrispondono alla riga sopra ma devono restare. "
       "Facoltativo."),
    NL("Uitzonderingen: dingen die bij de regel hierboven passen maar moeten "
       "blijven. Optioneel."),
    RU("Исключения: то, что подходит под строку выше, но должно остаться. "
       "Необязательно."),
    TR("İstisnalar: yukarıdaki satıra uyan ama kalması gerekenler. İsteğe "
       "bağlı."));

SS_MSG(objects_to_click,
    EN("Objects to click on"),
    JA("クリックして指定する物体"),
    ZH_HANS("要点选的物体"),
    ZH_HANT("要點選的物體"),
    KO("클릭할 물체"),
    DE("Objekte zum Anklicken"),
    FR("Objets à cliquer"),
    ES("Objetos que marcar"),
    PT("Objetos para clicar"),
    IT("Oggetti su cui cliccare"),
    NL("Objecten om aan te klikken"),
    RU("Объекты для указания"),
    TR("Tıklanacak nesneler"));

SS_MSG(objects_to_click_help,
    EN("One object per thing you want. SAM finds a single object per prompt, "
       "so clicking a person and then a car with the same object selected "
       "gives one mask that fits neither -- open a second object instead. "
       "Clicks belong to the frame you made them on: scrub to a later frame "
       "and click again to correct an object that has drifted."),
    JA("欲しいものごとに物体を1つ用意します。SAM は1つのプロンプトにつき1つの"
       "物体しか見つけないので、同じ物体を選んだまま人と車をクリックすると、"
       "どちらにも合わない1つのマスクになります。代わりに物体を追加して"
       "ください。クリックはそれを行ったフレームに属します。追跡がずれた"
       "物体は、後のフレームまで送ってからもう一度クリックすると直せます。"),
    ZH_HANS("你想要的每样东西各占一个物体。SAM 每条提示只找一个物体，所以在同一"
            "个物体上先点人再点车，会得到一个两边都不合的蒙版——请另开一个物体。"
            "点击属于你点它的那一帧：如果某个物体跟丢了，拖到后面的帧再点一次"
            "就能纠正。"),
    ZH_HANT("你想要的每樣東西各佔一個物體。SAM 每條提示只找一個物體，所以在同一"
            "個物體上先點人再點車，會得到一個兩邊都不合的遮罩——請另開一個物體。"
            "點擊屬於你點它的那一影格：如果某個物體跟丟了，拖到後面的影格再點一次"
            "就能糾正。"),
    KO("원하는 것마다 물체를 하나씩 두세요. SAM은 프롬프트 하나당 물체 하나만 "
       "찾으므로, 같은 물체를 고른 채 사람과 자동차를 차례로 클릭하면 둘 다 맞지 "
       "않는 마스크 하나가 나옵니다. 대신 물체를 하나 더 여세요. 클릭은 그것을 "
       "찍은 프레임에 속합니다. 추적이 어긋난 물체는 뒤쪽 프레임으로 옮겨 다시 "
       "클릭하면 바로잡을 수 있습니다."),
    DE("Ein Objekt je Sache, die Sie wollen. SAM findet pro Eingabe genau ein "
       "Objekt; klickt man also bei demselben ausgewählten Objekt erst eine "
       "Person und dann ein Auto an, entsteht eine Maske, die zu keinem von "
       "beiden passt -- stattdessen ein zweites Objekt anlegen. Klicks gehören "
       "zu dem Bild, auf dem sie gemacht wurden: zu einem späteren Bild "
       "spulen und erneut klicken, um ein abgedriftetes Objekt zu korrigieren."),
    FR("Un objet par chose voulue. SAM ne trouve qu'un seul objet par invite : "
       "cliquer une personne puis une voiture avec le même objet sélectionné "
       "donne un masque qui ne convient ni à l'une ni à l'autre -- ouvrez "
       "plutôt un deuxième objet. Les clics appartiennent à l'image où ils "
       "ont été faits : avancez à une image ultérieure et recliquez pour "
       "corriger un objet qui a dérivé."),
    ES("Un objeto por cada cosa que quiera. SAM encuentra un solo objeto por "
       "indicación, así que marcar una persona y luego un coche con el mismo "
       "objeto seleccionado da una máscara que no encaja con ninguno: abra un "
       "segundo objeto. Los clics pertenecen al fotograma en que se hicieron: "
       "avance a un fotograma posterior y vuelva a marcar para corregir un "
       "objeto que se ha desviado."),
    PT("Um objeto para cada coisa que você quer. O SAM encontra um único "
       "objeto por comando, então clicar numa pessoa e depois num carro com o "
       "mesmo objeto selecionado dá uma máscara que não serve para nenhum dos "
       "dois -- abra um segundo objeto. Os cliques pertencem ao quadro em que "
       "foram feitos: avance para um quadro posterior e clique de novo para "
       "corrigir um objeto que se desviou."),
    IT("Un oggetto per ogni cosa che le serve. SAM trova un solo oggetto per "
       "richiesta, quindi cliccare una persona e poi un'automobile con lo "
       "stesso oggetto selezionato dà una maschera che non va bene per "
       "nessuno dei due: apra piuttosto un secondo oggetto. I clic "
       "appartengono al fotogramma su cui sono stati fatti: vada a un "
       "fotogramma successivo e clicchi di nuovo per correggere un oggetto che "
       "è andato alla deriva."),
    NL("Eén object per ding dat je wilt. SAM vindt per prompt maar één object, "
       "dus eerst een persoon en dan een auto aanklikken met hetzelfde object "
       "geselecteerd geeft één masker dat bij geen van beide past -- open in "
       "plaats daarvan een tweede object. Klikken horen bij het beeld waarop "
       "je ze zette: spoel naar een later beeld en klik opnieuw om een object "
       "te corrigeren dat is afgedwaald."),
    RU("По одному объекту на каждую нужную вещь. SAM находит по одному объекту "
       "на запрос, так что если при выбранном объекте щёлкнуть сначала "
       "человека, а потом машину, выйдет одна маска, не подходящая ни тому, ни "
       "другому — заведите второй объект. Щелчки принадлежат тому кадру, где "
       "сделаны: перемотайте на более поздний кадр и щёлкните снова, чтобы "
       "поправить сбившийся объект."),
    TR("İstediğiniz her şey için bir nesne. SAM istem başına tek bir nesne "
       "bulur; aynı nesne seçiliyken önce bir insana sonra bir arabaya "
       "tıklarsanız ikisine de uymayan tek bir maske çıkar -- bunun yerine "
       "ikinci bir nesne açın. Tıklamalar yapıldıkları kareye aittir: kayan "
       "bir nesneyi düzeltmek için ileri bir kareye gidip yeniden tıklayın."));

// {0} object number, {1} clicks on this frame, {2} clicks on other frames.
SS_MSG(object_with_clicks,
    EN("Object {0} ({1} here, {2} elsewhere)"),
    JA("物体 {0}（このフレーム {1}、他 {2}）"),
    ZH_HANS("物体 {0}（本帧 {1}，其他 {2}）"),
    ZH_HANT("物體 {0}（本影格 {1}，其他 {2}）"),
    KO("물체 {0}(여기 {1}, 다른 곳 {2})"),
    DE("Objekt {0} ({1} hier, {2} anderswo)"),
    FR("Objet {0} ({1} ici, {2} ailleurs)"),
    ES("Objeto {0} ({1} aquí, {2} en otros)"),
    PT("Objeto {0} ({1} aqui, {2} em outros)"),
    IT("Oggetto {0} ({1} qui, {2} altrove)"),
    NL("Object {0} ({1} hier, {2} elders)"),
    RU("Объект {0} ({1} здесь, {2} в других кадрах)"),
    TR("Nesne {0} ({1} burada, {2} başka yerde)"));

SS_MSG(object_no_clicks,
    EN("Object {0} (no clicks yet)"),
    JA("物体 {0}（まだクリックなし）"),
    ZH_HANS("物体 {0}（还没有点击）"),
    ZH_HANT("物體 {0}（還沒有點擊）"),
    KO("물체 {0}(아직 클릭 없음)"),
    DE("Objekt {0} (noch keine Klicks)"),
    FR("Objet {0} (aucun clic pour l'instant)"),
    ES("Objeto {0} (aún sin clics)"),
    PT("Objeto {0} (ainda sem cliques)"),
    IT("Oggetto {0} (ancora nessun clic)"),
    NL("Object {0} (nog geen klikken)"),
    RU("Объект {0} (пока без щелчков)"),
    TR("Nesne {0} (henüz tıklama yok)"));

SS_MSG(object_clear,
    EN("Clear"),         JA("消去"),          ZH_HANS("清除"),     ZH_HANT("清除"),
    KO("지우기"),         DE("Löschen"),      FR("Effacer"),      ES("Borrar"),
    PT("Limpar"),        IT("Cancella"),     NL("Wissen"),       RU("Очистить"),
    TR("Temizle"));

SS_MSG(object_another,
    EN("Another object"), JA("物体を追加"),   ZH_HANS("再加一个物体"), ZH_HANT("再加一個物體"),
    KO("물체 하나 더"),   DE("Weiteres Objekt"), FR("Autre objet"),
    ES("Otro objeto"),   PT("Outro objeto"), IT("Un altro oggetto"),
    NL("Nog een object"), RU("Ещё объект"),  TR("Başka bir nesne"));

SS_MSG(object_another_help,
    EN("Adds an object for the next thing you click on."),
    JA("次にクリックするもののために物体を1つ追加します。"),
    ZH_HANS("为你接下来要点选的东西新增一个物体。"),
    ZH_HANT("為你接下來要點選的東西新增一個物體。"),
    KO("다음에 클릭할 것을 위해 물체를 하나 추가합니다."),
    DE("Legt ein Objekt für das an, was Sie als Nächstes anklicken."),
    FR("Ajoute un objet pour ce que vous cliquerez ensuite."),
    ES("Añade un objeto para lo próximo que marque."),
    PT("Adiciona um objeto para o próximo item em que você clicar."),
    IT("Aggiunge un oggetto per la prossima cosa su cui cliccherà."),
    NL("Voegt een object toe voor wat je hierna aanklikt."),
    RU("Добавляет объект для того, на что вы щёлкнете следующим."),
    TR("Bir sonraki tıklayacağınız şey için bir nesne ekler."));

SS_MSG(object_clear_all,
    EN("Clear all"),     JA("すべて消去"),    ZH_HANS("全部清除"),  ZH_HANT("全部清除"),
    KO("모두 지우기"),    DE("Alle löschen"), FR("Tout effacer"), ES("Borrar todo"),
    PT("Limpar tudo"),   IT("Cancella tutto"), NL("Alles wissen"),
    RU("Очистить всё"),  TR("Hepsini temizle"));

// {0} is the object number under the cursor.
SS_MSG(click_tooltip,
    EN("Left-click: this is object {0}.  Right-click: not this."),
    JA("左クリック: これが物体 {0} です。  右クリック: これは違います。"),
    ZH_HANS("左键：这是物体 {0}。  右键：不是这个。"),
    ZH_HANT("左鍵：這是物體 {0}。  右鍵：不是這個。"),
    KO("왼쪽 클릭: 이것이 물체 {0}입니다.  오른쪽 클릭: 이건 아닙니다."),
    DE("Linksklick: das ist Objekt {0}.  Rechtsklick: das nicht."),
    FR("Clic gauche : c'est l'objet {0}.  Clic droit : pas ça."),
    ES("Clic izquierdo: esto es el objeto {0}.  Clic derecho: esto no."),
    PT("Clique esquerdo: isto é o objeto {0}.  Clique direito: isto não."),
    IT("Clic sinistro: questo è l'oggetto {0}.  Clic destro: questo no."),
    NL("Linkerklik: dit is object {0}.  Rechterklik: dit niet."),
    RU("Левый щелчок: это объект {0}.  Правый щелчок: это не он."),
    TR("Sol tık: bu, nesne {0}.  Sağ tık: bu değil."));

SS_MSG(preview_frame,
    EN("frame {0}"),     JA("フレーム {0}"),  ZH_HANS("第 {0} 帧"), ZH_HANT("第 {0} 影格"),
    KO("{0}번 프레임"),   DE("Bild {0}"),     FR("image {0}"),    ES("fotograma {0}"),
    PT("quadro {0}"),    IT("fotogramma {0}"), NL("beeld {0}"),   RU("кадр {0}"),
    TR("kare {0}"));

SS_MSG(preview_frame_help,
    EN("A few frames from across the capture. Check a prompt on more than one "
       "before running -- and, for a video, click here to correct an object "
       "part way through: what you draw is used from this frame on."),
    JA("撮影全体から取った数枚のフレームです。実行する前に複数のフレームで"
       "プロンプトを確かめてください。動画では、途中で物体を修正したいときに"
       "ここでクリックします。描いた内容はこのフレーム以降に使われます。"),
    ZH_HANS("从整段拍摄中取出的几帧。运行前请在不止一帧上检查提示词——"
            "对视频而言，还可以在这里点选来修正中途的物体：你画的内容会从这一帧"
            "起生效。"),
    ZH_HANT("從整段拍攝中取出的幾格。執行前請在不止一格上檢查提示詞——"
            "對影片而言，還可以在這裡點選來修正中途的物體：你畫的內容會從這一格"
            "起生效。"),
    KO("촬영 전체에서 뽑은 몇 장의 프레임입니다. 실행하기 전에 두 장 이상에서 "
       "프롬프트를 확인하세요. 동영상이라면 중간에 물체를 바로잡을 때 여기서 "
       "클릭하면 됩니다. 표시한 내용은 이 프레임부터 적용됩니다."),
    DE("Ein paar Bilder aus der ganzen Aufnahme. Prüfen Sie eine Eingabe an "
       "mehr als einem, bevor Sie starten -- und bei einem Video hier klicken, "
       "um ein Objekt mittendrin zu korrigieren: was Sie einzeichnen, gilt ab "
       "diesem Bild."),
    FR("Quelques images prises tout au long de la prise. Vérifiez une invite "
       "sur plus d'une avant de lancer -- et, pour une vidéo, cliquez ici pour "
       "corriger un objet en cours de route : ce que vous tracez s'applique à "
       "partir de cette image."),
    ES("Unos cuantos fotogramas de toda la captura. Compruebe una indicación "
       "en más de uno antes de lanzar; y, en un vídeo, marque aquí para "
       "corregir un objeto a mitad de camino: lo que dibuje se aplica desde "
       "este fotograma."),
    PT("Alguns quadros de toda a captura. Confira um comando em mais de um "
       "antes de rodar -- e, num vídeo, clique aqui para corrigir um objeto no "
       "meio do caminho: o que você marcar vale a partir deste quadro."),
    IT("Alcuni fotogrammi presi da tutta la ripresa. Verifichi un testo su più "
       "d'uno prima di avviare; e, in un video, clicchi qui per correggere un "
       "oggetto a metà strada: quello che traccia vale da questo fotogramma in "
       "poi."),
    NL("Een paar beelden uit de hele opname. Controleer een prompt op meer dan "
       "één voordat je start -- en klik bij een video hier om halverwege een "
       "object te corrigeren: wat je aanwijst geldt vanaf dit beeld."),
    RU("Несколько кадров со всей съёмки. Перед запуском проверьте запрос "
       "больше чем на одном; а в видео щёлкните здесь, чтобы поправить объект "
       "по ходу: то, что вы отметите, действует начиная с этого кадра."),
    TR("Çekimin tamamından alınmış birkaç kare. Çalıştırmadan önce istemi "
       "birden çok karede sınayın -- videoda ise yolun ortasında bir nesneyi "
       "düzeltmek için buraya tıklayın: çizdikleriniz bu kareden itibaren "
       "geçerlidir."));

SS_MSG(preview_try_it,
    EN("Try it"),        JA("試す"),          ZH_HANS("试一下"),   ZH_HANT("試一下"),
    KO("해 보기"),        DE("Ausprobieren"), FR("Essayer"),      ES("Probar"),
    PT("Testar"),        IT("Prova"),        NL("Uitproberen"),  RU("Проверить"),
    TR("Dene"));

// Shown in the preview when this machine has no in-process video decoding, so
// the frame has to come from ffmpeg, and ffmpeg is not there either.
// {0} the command that was looked for.
SS_MSG(preview_needs_ffmpeg,
    EN("This video can only be read here with ffmpeg, which was not found "
       "('{0}'). Install it, or set its path under Tool locations."),
    JA("この動画をここで読むには ffmpeg が必要ですが、見つかりませんでした"
       "（{0}）。ffmpeg を入れるか、「外部ツールの場所」でパスを指定して"
       "ください。"),
    ZH_HANS("这里只能用 ffmpeg 读取该视频，但没有找到它（{0}）。"
            "请安装 ffmpeg，或在“工具位置”中指定它的路径。"),
    ZH_HANT("這裡只能用 ffmpeg 讀取該影片，但沒有找到它（{0}）。"
            "請安裝 ffmpeg，或在「工具位置」中指定它的路徑。"),
    KO("이 동영상은 여기서 ffmpeg으로만 읽을 수 있는데 ffmpeg을 찾지 "
       "못했습니다({0}). ffmpeg을 설치하거나 '도구 위치'에서 경로를 "
       "지정하세요."),
    DE("Dieses Video lässt sich hier nur mit ffmpeg lesen, und ffmpeg wurde "
       "nicht gefunden ('{0}'). Installieren Sie es, oder tragen Sie seinen "
       "Pfad unter \"Speicherorte der Werkzeuge\" ein."),
    FR("Cette vidéo ne peut être lue ici qu'avec ffmpeg, qui est introuvable "
       "(« {0} »). Installez-le, ou indiquez son chemin sous « Emplacement des "
       "outils »."),
    ES("Aquí este vídeo solo se puede leer con ffmpeg, que no se ha encontrado "
       "(«{0}»). Instálelo, o indique su ruta en «Ubicación de las "
       "herramientas»."),
    PT("Aqui este vídeo só pode ser lido com o ffmpeg, que não foi encontrado "
       "(\"{0}\"). Instale-o, ou informe o caminho dele em \"Local das "
       "ferramentas\"."),
    IT("Qui questo video si può leggere solo con ffmpeg, che non è stato "
       "trovato (\"{0}\"). Lo installi, oppure indichi il suo percorso in "
       "\"Percorsi degli strumenti\"."),
    NL("Deze video kan hier alleen met ffmpeg worden gelezen, en ffmpeg is "
       "niet gevonden ('{0}'). Installeer het, of geef het pad op onder "
       "\"Locatie van hulpprogramma's\"."),
    RU("Здесь это видео читается только через ffmpeg, а он не найден "
       "(«{0}»). Установите его или укажите путь в разделе «Расположение "
       "инструментов»."),
    TR("Bu video burada yalnızca ffmpeg ile okunabilir, ffmpeg ise bulunamadı "
       "('{0}'). Kurun ya da yolunu \"Araç konumları\" altında belirtin."));

// Neither the built-in decoder nor ffmpeg produced a picture.
SS_MSG(preview_frame_unreadable,
    EN("No frame could be read from this video."),
    JA("この動画からフレームを読み取れませんでした。"),
    ZH_HANS("无法从这个视频中读取任何一帧。"),
    ZH_HANT("無法從這個影片中讀取任何一格。"),
    KO("이 동영상에서 프레임을 읽지 못했습니다."),
    DE("Aus diesem Video ließ sich kein Bild lesen."),
    FR("Aucune image n'a pu être lue dans cette vidéo."),
    ES("No se ha podido leer ningún fotograma de este vídeo."),
    PT("Não foi possível ler nenhum quadro deste vídeo."),
    IT("Non è stato possibile leggere alcun fotogramma da questo video."),
    NL("Er kon geen beeld uit deze video worden gelezen."),
    RU("Из этого видео не удалось прочитать ни одного кадра."),
    TR("Bu videodan hiçbir kare okunamadı."));

SS_MSG(preview_kept_fraction,
    EN("{0}% of the frame is kept"),
    JA("フレームの {0}% が残ります"),
    ZH_HANS("这一帧保留了 {0}%"),
    ZH_HANT("這一影格保留了 {0}%"),
    KO("프레임의 {0}%가 남습니다"),
    DE("{0} % des Bildes bleiben erhalten"),
    FR("{0} % de l'image est conservé"),
    ES("se conserva el {0} % del fotograma"),
    PT("{0}% do quadro é mantido"),
    IT("resta il {0}% del fotogramma"),
    NL("{0}% van het beeld blijft over"),
    RU("остаётся {0}% кадра"),
    TR("karenin %{0}'i tutuluyor"));

SS_MSG(preview_almost_nothing_kept,
    EN("Almost nothing is left -- the prompt matched very little of the "
       "frame."),
    JA("ほとんど何も残っていません。プロンプトがフレームのごく一部にしか"
       "当てはまりませんでした。"),
    ZH_HANS("几乎什么都没剩下——提示词只匹配到画面里很小的一部分。"),
    ZH_HANT("幾乎什麼都沒剩下——提示詞只匹配到畫面裡很小的一部分。"),
    KO("거의 아무것도 남지 않았습니다. 프롬프트가 프레임의 아주 일부에만 "
       "맞았습니다."),
    DE("Es bleibt fast nichts übrig -- die Eingabe traf nur sehr wenig des "
       "Bildes."),
    FR("Il ne reste presque rien : l'invite n'a couvert qu'une toute petite "
       "partie de l'image."),
    ES("Casi no queda nada: la indicación coincidió con muy poco del "
       "fotograma."),
    PT("Quase nada sobrou -- o comando correspondeu a muito pouco do quadro."),
    IT("Non resta quasi nulla: il testo ha corrisposto a pochissimo del "
       "fotogramma."),
    NL("Er blijft bijna niets over -- de prompt paste bij maar heel weinig van "
       "het beeld."),
    RU("Почти ничего не осталось — запрос совпал лишь с малой частью кадра."),
    TR("Neredeyse hiçbir şey kalmadı -- istem karenin çok azıyla eşleşti."));

SS_MSG(preview_almost_all_masked,
    EN("Almost everything is masked out -- did you mean \"Keep only what I "
       "name\"?"),
    JA("ほとんどすべてがマスクされています。「指定したものだけを残す」の"
       "つもりではありませんか？"),
    ZH_HANS("几乎所有内容都被蒙掉了——你是不是想选“只保留我指名的东西”？"),
    ZH_HANT("幾乎所有內容都被遮掉了——你是不是想選「只保留我指名的東西」？"),
    KO("거의 모든 것이 마스킹되었습니다. '내가 말한 것만 남기기'를 뜻하신 "
       "건가요?"),
    DE("Fast alles ist maskiert -- war „Nur behalten, was ich nenne“ gemeint?"),
    FR("Presque tout est masqué -- vouliez-vous dire « Ne garder que ce que je "
       "nomme » ?"),
    ES("Casi todo está enmascarado: ¿quería decir «Conservar solo lo que yo "
       "nombre»?"),
    PT("Quase tudo está mascarado -- você queria dizer “Manter só o que eu "
       "nomear”?"),
    IT("È mascherato quasi tutto: intendeva «Tenere solo ciò che indico»?"),
    NL("Bijna alles is gemaskeerd -- bedoelde je ‘Alleen houden wat ik noem’?"),
    RU("Замаскировано почти всё — вы имели в виду «Оставить только то, что я "
       "назову»?"),
    TR("Neredeyse her şey maskelendi -- “Yalnızca adını verdiğimi tut”u mu "
       "demek istediniz?"));

// ===========================================================================
// Advanced: built-in reconstruction
// ===========================================================================

SS_MSG(capture_type,
    EN("Capture type"),  JA("撮影の種類"),    ZH_HANS("拍摄类型"),  ZH_HANT("拍攝類型"),
    KO("촬영 종류"),      DE("Aufnahmeart"),  FR("Type de prise"), ES("Tipo de captura"),
    PT("Tipo de captura"), IT("Tipo di ripresa"), NL("Soort opname"),
    RU("Тип съёмки"),    TR("Çekim türü"));

SS_MSG(capture_photos,
    EN("Individual photos"),
    JA("個別の写真"),    ZH_HANS("单张照片"),  ZH_HANT("單張相片"),
    KO("낱장 사진"),     DE("Einzelfotos"),   FR("Photos individuelles"),
    ES("Fotos sueltas"), PT("Fotos avulsas"), IT("Fotografie singole"),
    NL("Losse foto's"),  RU("Отдельные фотографии"), TR("Tek tek fotoğraflar"));

SS_MSG(capture_video,
    EN("Video frames"),  JA("動画のフレーム"), ZH_HANS("视频帧"),   ZH_HANT("影片影格"),
    KO("동영상 프레임"), DE("Videobilder"),   FR("Images de vidéo"),
    ES("Fotogramas de vídeo"), PT("Quadros de vídeo"), IT("Fotogrammi video"),
    NL("Videobeelden"),  RU("Кадры видео"),  TR("Video kareleri"));

SS_MSG(capture_internet,
    EN("Unordered internet collection"),
    JA("順不同のインターネット写真集"),
    ZH_HANS("无序的网络图片集"),
    ZH_HANT("無序的網路圖片集"),
    KO("순서 없는 인터넷 사진 모음"),
    DE("Ungeordnete Internetsammlung"),
    FR("Collection internet sans ordre"),
    ES("Colección de internet sin orden"),
    PT("Coleção da internet sem ordem"),
    IT("Raccolta internet senza ordine"),
    NL("Ongeordende internetverzameling"),
    RU("Неупорядоченная подборка из интернета"),
    TR("Sırasız internet derlemesi"));

SS_MSG(capture_type_help,
    EN("What the input is, which sets the pairing strategy and how forgiving "
       "the mapper is. Set from the input type when you picked it."),
    JA("入力が何であるかです。ペアの作り方と、マッパーがどれだけ寛容かが"
       "これで決まります。入力を選んだ時点で種類から設定されます。"),
    ZH_HANS("输入是什么，它决定了配对策略以及建图器有多宽容。选择输入时会按"
            "输入类型自动设定。"),
    ZH_HANT("輸入是什麼，它決定了配對策略以及建圖器有多寬容。選擇輸入時會按"
            "輸入類型自動設定。"),
    KO("입력이 무엇인지입니다. 짝짓기 전략과 매퍼가 얼마나 관대한지가 여기서 "
       "정해집니다. 입력을 고를 때 그 종류에서 자동으로 설정됩니다."),
    DE("Was die Eingabe ist; davon hängen die Paarungsstrategie und die "
       "Nachsicht des Mappers ab. Wird beim Auswählen aus der Eingabeart "
       "gesetzt."),
    FR("Ce qu'est l'entrée, ce qui fixe la stratégie d'appariement et la "
       "tolérance du mapper. Défini d'après le type d'entrée au moment du "
       "choix."),
    ES("Qué es la entrada, lo que fija la estrategia de emparejamiento y la "
       "tolerancia del mapeador. Se establece según el tipo de entrada al "
       "elegirla."),
    PT("O que é a entrada, o que define a estratégia de pareamento e a "
       "tolerância do mapeador. Definido pelo tipo de entrada ao escolhê-la."),
    IT("Che cosa è l'ingresso: da questo dipendono la strategia di "
       "accoppiamento e quanto è tollerante il mapper. Impostato dal tipo di "
       "ingresso al momento della scelta."),
    NL("Wat de invoer is; dat bepaalt de koppelstrategie en hoe mild de mapper "
       "is. Wordt bij het kiezen uit het invoertype ingesteld."),
    RU("Чем является вход; от этого зависят стратегия составления пар и "
       "снисходительность маппера. Задаётся по типу входа при его выборе."),
    TR("Girdinin ne olduğu; eşleştirme stratejisini ve haritalayıcının ne "
       "kadar hoşgörülü olduğunu bu belirler. Girdiyi seçerken türünden "
       "ayarlanır."));

SS_MSG(features,
    EN("Features"),      JA("特徴"),          ZH_HANS("特征"),     ZH_HANT("特徵"),
    KO("특징점"),         DE("Merkmale"),     FR("Points caractéristiques"),
    ES("Características"), PT("Características"), IT("Caratteristiche"),
    NL("Kenmerken"),     RU("Особые точки"), TR("Öznitelikler"));

SS_MSG(features_sift,
    EN("SIFT (classic)"), JA("SIFT（古典的）"), ZH_HANS("SIFT（经典）"),
    ZH_HANT("SIFT（經典）"), KO("SIFT(고전적)"), DE("SIFT (klassisch)"),
    FR("SIFT (classique)"), ES("SIFT (clásico)"), PT("SIFT (clássico)"),
    IT("SIFT (classico)"), NL("SIFT (klassiek)"), RU("SIFT (классический)"),
    TR("SIFT (klasik)"));

SS_MSG(features_aliked_n16,
    EN("ALIKED N16-rot (learned)"),
    JA("ALIKED N16-rot（学習型）"),
    ZH_HANS("ALIKED N16-rot（学习型）"),
    ZH_HANT("ALIKED N16-rot（學習型）"),
    KO("ALIKED N16-rot(학습형)"),
    DE("ALIKED N16-rot (gelernt)"),
    FR("ALIKED N16-rot (appris)"),
    ES("ALIKED N16-rot (aprendido)"),
    PT("ALIKED N16-rot (aprendido)"),
    IT("ALIKED N16-rot (appreso)"),
    NL("ALIKED N16-rot (geleerd)"),
    RU("ALIKED N16-rot (обученный)"),
    TR("ALIKED N16-rot (öğrenilmiş)"));

SS_MSG(features_aliked_n32,
    EN("ALIKED N32 (learned, wider)"),
    JA("ALIKED N32（学習型・広め）"),
    ZH_HANS("ALIKED N32（学习型，更宽）"),
    ZH_HANT("ALIKED N32（學習型，更寬）"),
    KO("ALIKED N32(학습형, 더 넓음)"),
    DE("ALIKED N32 (gelernt, breiter)"),
    FR("ALIKED N32 (appris, plus large)"),
    ES("ALIKED N32 (aprendido, más amplio)"),
    PT("ALIKED N32 (aprendido, mais amplo)"),
    IT("ALIKED N32 (appreso, più ampio)"),
    NL("ALIKED N32 (geleerd, breder)"),
    RU("ALIKED N32 (обученный, шире)"),
    TR("ALIKED N32 (öğrenilmiş, daha geniş)"));

SS_MSG(features_help,
    EN("Which detector and descriptor. SIFT is the classic one and needs "
       "nothing downloaded. The ALIKED options are a learned frontend: they "
       "fetch a small checkpoint (3-4 MB) on first use, find fewer but "
       "better-localized keypoints, and match markedly more image pairs on "
       "hard captures. N32 samples more positions per descriptor -- slower, "
       "slightly stronger."),
    JA("どの検出器と記述子を使うかです。SIFT は古典的な選択で、何もダウンロード"
       "しません。ALIKED は学習型のフロントエンドで、初回に小さな"
       "チェックポイント（3〜4 MB）を取得します。キーポイントの数は少ない"
       "ものの位置が正確で、難しい撮影ではマッチする画像ペアが目に見えて"
       "増えます。N32 は記述子あたりのサンプル位置が多く、遅いぶん少し強力です。"),
    ZH_HANS("使用哪种检测器和描述子。SIFT 是经典选择，什么都不用下载。"
            "ALIKED 是学习型前端：首次使用会取一个小的检查点（3-4 MB），"
            "找到的关键点更少但定位更准，在困难拍摄上匹配上的图像对明显更多。"
            "N32 每个描述子采样更多位置——更慢，略强。"),
    ZH_HANT("使用哪種偵測器和描述子。SIFT 是經典選擇，什麼都不用下載。"
            "ALIKED 是學習型前端：首次使用會取一個小的檢查點（3-4 MB），"
            "找到的關鍵點更少但定位更準，在困難拍攝上配對上的影像對明顯更多。"
            "N32 每個描述子取樣更多位置——更慢，略強。"),
    KO("어떤 검출기와 기술자를 쓸지입니다. SIFT는 고전적인 선택이며 내려받을 것이 "
       "없습니다. ALIKED는 학습형 프런트엔드로, 처음 쓸 때 작은 체크포인트"
       "(3~4 MB)를 받아옵니다. 키포인트 수는 적지만 위치가 더 정확하고, 어려운 "
       "촬영에서 매칭되는 이미지 쌍이 눈에 띄게 늘어납니다. N32는 기술자당 더 "
       "많은 위치를 표본화합니다 — 더 느리고 조금 더 강합니다."),
    DE("Welcher Detektor und Deskriptor. SIFT ist der klassische und braucht "
       "keinen Download. Die ALIKED-Varianten sind ein gelerntes Frontend: sie "
       "holen beim ersten Gebrauch einen kleinen Prüfpunkt (3-4 MB), finden "
       "weniger, aber genauer verortete Schlüsselpunkte und ordnen bei "
       "schwierigen Aufnahmen deutlich mehr Bildpaare zu. N32 tastet je "
       "Deskriptor mehr Stellen ab -- langsamer, etwas stärker."),
    FR("Quel détecteur et quel descripteur. SIFT est le classique et ne "
       "demande aucun téléchargement. Les options ALIKED forment un frontal "
       "appris : elles récupèrent un petit point de sauvegarde (3-4 Mo) au "
       "premier usage, trouvent moins de points mais mieux localisés, et "
       "apparient nettement plus de paires sur les prises difficiles. N32 "
       "échantillonne plus de positions par descripteur -- plus lent, un peu "
       "plus solide."),
    ES("Qué detector y descriptor. SIFT es el clásico y no necesita descargar "
       "nada. Las opciones ALIKED son un frontal aprendido: bajan un pequeño "
       "punto de control (3-4 MB) la primera vez, encuentran menos puntos "
       "pero mejor localizados, y emparejan bastantes más pares de imágenes "
       "en capturas difíciles. N32 muestrea más posiciones por descriptor: "
       "más lento, algo más robusto."),
    PT("Qual detector e descritor. O SIFT é o clássico e não precisa baixar "
       "nada. As opções ALIKED são um front-end aprendido: buscam um pequeno "
       "ponto de verificação (3-4 MB) no primeiro uso, encontram menos pontos "
       "mas melhor localizados, e casam bem mais pares de imagens em capturas "
       "difíceis. O N32 amostra mais posições por descritor -- mais lento, um "
       "pouco mais forte."),
    IT("Quale rilevatore e descrittore. SIFT è quello classico e non richiede "
       "scaricamenti. Le opzioni ALIKED sono un frontend appreso: al primo "
       "uso prelevano un piccolo punto di controllo (3-4 MB), trovano meno "
       "punti chiave ma meglio localizzati e mettono in corrispondenza molte "
       "più coppie nelle riprese difficili. N32 campiona più posizioni per "
       "descrittore: più lento, un po' più robusto."),
    NL("Welke detector en descriptor. SIFT is de klassieke en hoeft niets te "
       "downloaden. De ALIKED-opties zijn een geleerde frontend: ze halen bij "
       "eerste gebruik een klein controlepunt (3-4 MB) op, vinden minder maar "
       "beter geplaatste sleutelpunten en matchen bij lastige opnamen "
       "merkbaar meer beeldparen. N32 bemonstert meer posities per descriptor "
       "-- trager, iets sterker."),
    RU("Какой детектор и дескриптор. SIFT — классический, ничего скачивать не "
       "нужно. Варианты ALIKED — обученный фронтенд: при первом использовании "
       "они получают небольшую контрольную точку (3-4 МБ), находят меньше "
       "ключевых точек, но точнее локализованных, и на сложных съёмках "
       "сопоставляют заметно больше пар. N32 берёт больше позиций на "
       "дескриптор — медленнее, чуть надёжнее."),
    TR("Hangi bulucu ve betimleyici. SIFT klasik olanıdır ve indirme "
       "gerektirmez. ALIKED seçenekleri öğrenilmiş bir ön uçtur: ilk "
       "kullanımda küçük bir denetim noktası (3-4 MB) indirir, daha az ama "
       "daha iyi konumlanmış anahtar nokta bulur ve zor çekimlerde belirgin "
       "biçimde daha çok görüntü çiftini eşleştirir. N32 betimleyici başına "
       "daha çok konum örnekler -- daha yavaş, biraz daha güçlü."));

SS_MSG(matcher,
    EN("Matcher"),       JA("マッチャー"),    ZH_HANS("匹配器"),   ZH_HANT("配對器"),
    KO("매처"),           DE("Zuordner"),     FR("Appariement"),  ES("Emparejador"),
    PT("Correspondedor"), IT("Abbinatore"),  NL("Matcher"),      RU("Сопоставитель"),
    TR("Eşleştirici"));

SS_MSG(matcher_brute_force,
    EN("Brute force"),   JA("総当たり"),      ZH_HANS("暴力匹配"),  ZH_HANT("暴力比對"),
    KO("완전 탐색"),      DE("Brute Force"),  FR("Force brute"),  ES("Fuerza bruta"),
    PT("Força bruta"),   IT("Forza bruta"),  NL("Brute kracht"), RU("Полный перебор"),
    TR("Kaba kuvvet"));

SS_MSG(matcher_lightglue,
    EN("LightGlue (learned)"),
    JA("LightGlue（学習型）"), ZH_HANS("LightGlue（学习型）"),
    ZH_HANT("LightGlue（學習型）"), KO("LightGlue(학습형)"),
    DE("LightGlue (gelernt)"), FR("LightGlue (appris)"),
    ES("LightGlue (aprendido)"), PT("LightGlue (aprendido)"),
    IT("LightGlue (appreso)"), NL("LightGlue (geleerd)"),
    RU("LightGlue (обученный)"), TR("LightGlue (öğrenilmiş)"));

SS_MSG(matcher_help,
    EN("How descriptors are matched. LightGlue is a learned matcher: it finds "
       "far more correct correspondences on hard pairs, and costs tens of "
       "milliseconds per pair instead of a few, so it runs behind pair "
       "selection."),
    JA("記述子をどう照合するかです。LightGlue は学習型のマッチャーで、難しい"
       "ペアでも正しい対応を格段に多く見つけます。1 ペアあたり数ミリ秒ではなく"
       "数十ミリ秒かかるため、ペア選択の後段で動きます。"),
    ZH_HANS("描述子如何匹配。LightGlue 是学习型匹配器：在困难配对上能找到多得多"
            "的正确对应，代价是每对要几十毫秒而不是几毫秒，所以它跑在配对筛选之后。"),
    ZH_HANT("描述子如何比對。LightGlue 是學習型配對器：在困難配對上能找到多得多"
            "的正確對應，代價是每對要幾十毫秒而不是幾毫秒，所以它跑在配對篩選之後。"),
    KO("기술자를 어떻게 맞출지입니다. LightGlue는 학습형 매처로, 어려운 쌍에서 "
       "올바른 대응을 훨씬 많이 찾아냅니다. 쌍당 몇 밀리초가 아니라 수십 밀리초가 "
       "들기 때문에 쌍 선별 뒤에 실행됩니다."),
    DE("Wie Deskriptoren zugeordnet werden. LightGlue ist ein gelernter "
       "Zuordner: er findet bei schwierigen Paaren weit mehr richtige "
       "Entsprechungen und kostet je Paar zig statt weniger Millisekunden, "
       "läuft also hinter der Paarauswahl."),
    FR("Comment les descripteurs sont appariés. LightGlue est un apparieur "
       "appris : il trouve bien plus de correspondances justes sur les paires "
       "difficiles, et coûte des dizaines de millisecondes par paire au lieu "
       "de quelques-unes, donc il s'exécute après la sélection des paires."),
    ES("Cómo se emparejan los descriptores. LightGlue es un emparejador "
       "aprendido: encuentra muchas más correspondencias correctas en pares "
       "difíciles y cuesta decenas de milisegundos por par en vez de unos "
       "pocos, así que se ejecuta tras la selección de pares."),
    PT("Como os descritores são casados. O LightGlue é um correspondedor "
       "aprendido: encontra muito mais correspondências corretas em pares "
       "difíceis e custa dezenas de milissegundos por par em vez de poucos, "
       "então roda depois da seleção de pares."),
    IT("Come vengono abbinati i descrittori. LightGlue è un abbinatore "
       "appreso: trova molte più corrispondenze corrette sulle coppie "
       "difficili e costa decine di millisecondi per coppia invece di pochi, "
       "perciò gira dopo la selezione delle coppie."),
    NL("Hoe descriptors worden gematcht. LightGlue is een geleerde matcher: "
       "die vindt bij lastige paren veel meer juiste overeenkomsten en kost "
       "tientallen milliseconden per paar in plaats van enkele, dus draait hij "
       "na de paarselectie."),
    RU("Как сопоставляются дескрипторы. LightGlue — обученный сопоставитель: "
       "на трудных парах он находит гораздо больше верных соответствий и "
       "стоит десятков миллисекунд на пару вместо единиц, поэтому работает "
       "после отбора пар."),
    TR("Betimleyicilerin nasıl eşleştirileceği. LightGlue öğrenilmiş bir "
       "eşleştiricidir: zor çiftlerde çok daha fazla doğru karşılık bulur ve "
       "çift başına birkaç yerine onlarca milisaniye harcar, bu yüzden çift "
       "seçiminin ardından çalışır."));

SS_MSG(matcher_needs_learned,
    EN("LightGlue needs the learned descriptors -- pick an ALIKED frontend "
       "above to enable it."),
    JA("LightGlue には学習型の記述子が必要です。上で ALIKED のフロントエンドを"
       "選ぶと使えるようになります。"),
    ZH_HANS("LightGlue 需要学习型描述子——请在上面选一个 ALIKED 前端来启用它。"),
    ZH_HANT("LightGlue 需要學習型描述子——請在上面選一個 ALIKED 前端來啟用它。"),
    KO("LightGlue에는 학습형 기술자가 필요합니다. 위에서 ALIKED 프런트엔드를 "
       "고르면 켜집니다."),
    DE("LightGlue braucht die gelernten Deskriptoren -- oben ein "
       "ALIKED-Frontend wählen, um es freizuschalten."),
    FR("LightGlue exige les descripteurs appris -- choisissez un frontal "
       "ALIKED ci-dessus pour l'activer."),
    ES("LightGlue necesita los descriptores aprendidos: elija arriba un "
       "frontal ALIKED para activarlo."),
    PT("O LightGlue precisa dos descritores aprendidos -- escolha acima um "
       "front-end ALIKED para habilitá-lo."),
    IT("LightGlue richiede i descrittori appresi: scelga sopra un frontend "
       "ALIKED per abilitarlo."),
    NL("LightGlue heeft de geleerde descriptors nodig -- kies hierboven een "
       "ALIKED-frontend om het aan te zetten."),
    RU("LightGlue нужны обученные дескрипторы — выберите выше фронтенд ALIKED, "
       "чтобы он стал доступен."),
    TR("LightGlue öğrenilmiş betimleyicilere ihtiyaç duyar -- etkinleştirmek "
       "için yukarıdan bir ALIKED ön ucu seçin."));

SS_MSG(mapper_schedule,
    EN("Mapper schedule"), JA("マッパーの進め方"), ZH_HANS("建图策略"), ZH_HANT("建圖策略"),
    KO("매퍼 진행 방식"), DE("Mapper-Ablauf"), FR("Ordonnancement du mapper"),
    ES("Estrategia del mapeador"), PT("Estratégia do mapeador"),
    IT("Strategia del mapper"), NL("Mapper-schema"),
    RU("Схема работы маппера"), TR("Haritalayıcı planı"));

SS_MSG(mapper_flat,
    EN("Flat (one reconstruction)"),
    JA("フラット（1つの再構成）"),
    ZH_HANS("扁平（单个重建）"),
    ZH_HANT("扁平（單個重建）"),
    KO("평면(재구성 하나)"),
    DE("Flach (eine Rekonstruktion)"),
    FR("Plat (une seule reconstruction)"),
    ES("Plano (una sola reconstrucción)"),
    PT("Plano (uma única reconstrução)"),
    IT("Piatto (una sola ricostruzione)"),
    NL("Vlak (één reconstructie)"),
    RU("Плоская (одна реконструкция)"),
    TR("Düz (tek yeniden oluşturma)"));

SS_MSG(mapper_bottom_up,
    EN("Bottom-up (atoms, merged upwards)"),
    JA("ボトムアップ（小さな塊を作って統合）"),
    ZH_HANS("自下而上（先建小块再向上合并）"),
    ZH_HANT("由下而上（先建小塊再向上合併）"),
    KO("상향식(작은 덩어리를 만들어 위로 병합)"),
    DE("Von unten (Atome, nach oben zusammengeführt)"),
    FR("Ascendant (atomes, fusionnés vers le haut)"),
    ES("Ascendente (átomos, fusionados hacia arriba)"),
    PT("Ascendente (átomos, mesclados para cima)"),
    IT("Dal basso (atomi, uniti verso l'alto)"),
    NL("Bottom-up (atomen, naar boven samengevoegd)"),
    RU("Снизу вверх (атомы, объединяемые вверх)"),
    TR("Aşağıdan yukarı (atomlar, yukarı doğru birleştirilir)"));

SS_MSG(mapper_schedule_help,
    EN("How the scene is built. Flat grows one reconstruction image by image, "
       "and is the default for any capture. Bottom-up cuts the view graph "
       "into small groups, reconstructs them independently and merges upwards "
       "-- worth trying on a large capture, where the flat schedule's "
       "whole-model passes start to dominate."),
    JA("シーンをどう組み立てるかです。フラットは1つの再構成を画像ごとに"
       "育てていく方式で、どの撮影でも既定です。ボトムアップはビューグラフを"
       "小さなグループに切り、それぞれを独立に再構成してから上へ統合します。"
       "大規模な撮影では、フラット方式のモデル全体を走査するパスが支配的に"
       "なってくるので、試す価値があります。"),
    ZH_HANS("场景如何构建。扁平方式逐张图像地扩展同一个重建，是任何拍摄的默认选择。"
            "自下而上把视图图切成小组，各自独立重建后再向上合并——在大规模拍摄上"
            "值得一试，因为那时扁平方式的全模型遍历会开始主导耗时。"),
    ZH_HANT("場景如何建構。扁平方式逐張影像地擴展同一個重建，是任何拍攝的預設選擇。"
            "由下而上把視圖圖切成小組，各自獨立重建後再向上合併——在大規模拍攝上"
            "值得一試，因為那時扁平方式的全模型走訪會開始主導耗時。"),
    KO("장면을 어떻게 만들지입니다. 평면 방식은 재구성 하나를 이미지마다 키워 "
       "가며, 어떤 촬영에서든 기본값입니다. 상향식은 뷰 그래프를 작은 묶음으로 "
       "자르고 각각 독립적으로 재구성한 뒤 위로 병합합니다. 규모가 큰 촬영에서는 "
       "평면 방식의 전체 모델 패스가 시간을 지배하기 시작하므로 시도해 볼 만합니다."),
    DE("Wie die Szene aufgebaut wird. Flach lässt eine Rekonstruktion Bild für "
       "Bild wachsen und ist die Vorgabe für jede Aufnahme. Von unten "
       "zerlegt den Sichtgraphen in kleine Gruppen, rekonstruiert sie "
       "unabhängig und führt sie nach oben zusammen -- bei einer großen "
       "Aufnahme einen Versuch wert, wo die Ganzmodell-Durchgänge des flachen "
       "Ablaufs zu dominieren beginnen."),
    FR("Comment la scène est construite. Le mode plat fait croître une seule "
       "reconstruction image par image, et c'est la valeur par défaut pour "
       "toute prise. Le mode ascendant découpe le graphe de vues en petits "
       "groupes, les reconstruit indépendamment et fusionne vers le haut -- à "
       "essayer sur une grande prise, où les passes sur le modèle entier du "
       "mode plat finissent par dominer."),
    ES("Cómo se construye la escena. El modo plano hace crecer una sola "
       "reconstrucción imagen a imagen, y es lo predeterminado en cualquier "
       "captura. El ascendente corta el grafo de vistas en grupos pequeños, "
       "los reconstruye por separado y los fusiona hacia arriba: vale la pena "
       "probarlo en capturas grandes, donde las pasadas sobre el modelo "
       "entero del modo plano empiezan a dominar."),
    PT("Como a cena é construída. O modo plano faz uma única reconstrução "
       "crescer imagem a imagem, e é o padrão para qualquer captura. O "
       "ascendente corta o grafo de vistas em grupos pequenos, reconstrói "
       "cada um em separado e mescla para cima -- vale tentar numa captura "
       "grande, em que as passagens sobre o modelo inteiro do modo plano "
       "começam a dominar."),
    IT("Come viene costruita la scena. Il modo piatto fa crescere una sola "
       "ricostruzione immagine dopo immagine ed è l'impostazione predefinita "
       "per qualsiasi ripresa. Quello dal basso taglia il grafo delle viste in "
       "piccoli gruppi, li ricostruisce in modo indipendente e li unisce verso "
       "l'alto: vale la pena provarlo su una ripresa grande, dove le passate "
       "sull'intero modello del modo piatto cominciano a dominare."),
    NL("Hoe de scène wordt opgebouwd. Vlak laat één reconstructie beeld voor "
       "beeld groeien en is de standaard voor elke opname. Bottom-up knipt de "
       "beeldgraaf in kleine groepen, reconstrueert die apart en voegt ze naar "
       "boven samen -- het proberen waard bij een grote opname, waar de "
       "passes over het hele model van het vlakke schema gaan overheersen."),
    RU("Как строится сцена. Плоская схема наращивает одну реконструкцию снимок "
       "за снимком и подходит любой съёмке. Схема снизу вверх режет граф видов "
       "на небольшие группы, восстанавливает их независимо и объединяет "
       "вверх — стоит попробовать на крупной съёмке, где проходы плоской схемы "
       "по всей модели начинают преобладать."),
    TR("Sahnenin nasıl kurulacağı. Düz plan tek bir yeniden oluşturmayı "
       "görüntü görüntü büyütür ve her çekim için varsayılandır. Aşağıdan "
       "yukarı plan görünüm çizgesini küçük öbeklere böler, her birini ayrı "
       "yeniden oluşturur ve yukarı doğru birleştirir -- düz planın tüm model "
       "geçişlerinin ağır basmaya başladığı büyük çekimlerde denemeye "
       "değer."));

SS_MSG(sequential_overlap,
    EN("Sequential overlap"),
    JA("逐次マッチングの重なり"),
    ZH_HANS("顺序匹配的重叠数"),
    ZH_HANT("循序比對的重疊數"),
    KO("순차 겹침"),
    DE("Sequenzielle Überlappung"),
    FR("Recouvrement séquentiel"),
    ES("Solapamiento secuencial"),
    PT("Sobreposição sequencial"),
    IT("Sovrapposizione sequenziale"),
    NL("Sequentiële overlap"),
    RU("Перекрытие при последовательном сопоставлении"),
    TR("Sıralı örtüşme"));

SS_MSG(sequential_overlap_help,
    EN("How many neighbouring frames each frame is matched against."),
    JA("各フレームを、隣り合う何フレームと照合するかです。"),
    ZH_HANS("每一帧要和多少相邻帧做匹配。"),
    ZH_HANT("每一影格要和多少相鄰影格做比對。"),
    KO("각 프레임을 이웃한 몇 개의 프레임과 매칭할지입니다."),
    DE("Mit wie vielen benachbarten Bildern jedes Bild verglichen wird."),
    FR("Nombre d'images voisines auxquelles chaque image est appariée."),
    ES("Con cuántos fotogramas vecinos se empareja cada fotograma."),
    PT("Com quantos quadros vizinhos cada quadro é comparado."),
    IT("Con quanti fotogrammi vicini viene confrontato ogni fotogramma."),
    NL("Met hoeveel naburige beelden elk beeld wordt gematcht."),
    RU("Со сколькими соседними кадрами сопоставляется каждый кадр."),
    TR("Her karenin kaç komşu kareyle eşleştirileceği."));

SS_MSG(initial_focal_px,
    EN("Initial focal length (px, 0 = auto)"),
    JA("初期焦点距離（px、0 で自動）"),
    ZH_HANS("初始焦距（像素，0 = 自动）"),
    ZH_HANT("初始焦距（像素，0 = 自動）"),
    KO("초기 초점거리(px, 0 = 자동)"),
    DE("Anfangsbrennweite (px, 0 = automatisch)"),
    FR("Focale initiale (px, 0 = auto)"),
    ES("Focal inicial (px, 0 = automático)"),
    PT("Distância focal inicial (px, 0 = automático)"),
    IT("Focale iniziale (px, 0 = automatico)"),
    NL("Beginbrandpuntsafstand (px, 0 = automatisch)"),
    RU("Начальное фокусное расстояние (пикс., 0 — авто)"),
    TR("Başlangıç odak uzaklığı (px, 0 = otomatik)"));

SS_MSG(initial_focal_px_help,
    EN("Starting guess for the focal length, in pixels of the source image. 0 "
       "reads EXIF and falls back to a guess from the image size. Worth "
       "setting for a fisheye, where a bad initial guess can stop the "
       "reconstruction from starting at all."),
    JA("焦点距離の初期推定値を、元画像のピクセル単位で指定します。0 なら EXIF を"
       "読み、なければ画像サイズから推定します。魚眼では初期値が悪いと再構成が"
       "そもそも始まらないことがあるので、指定する価値があります。"),
    ZH_HANS("焦距的初始猜测值，以源图像的像素为单位。填 0 会读 EXIF，读不到就按"
            "图像尺寸估计。鱼眼值得设一下：初值不好可能让重建根本起不来。"),
    ZH_HANT("焦距的初始猜測值，以來源影像的像素為單位。填 0 會讀 EXIF，讀不到就按"
            "影像尺寸估計。魚眼值得設一下：初值不好可能讓重建根本起不來。"),
    KO("초점거리의 초기 추정값을 원본 이미지의 픽셀 단위로 지정합니다. 0이면 "
       "EXIF를 읽고, 없으면 이미지 크기로 추정합니다. 어안에서는 초기 추정이 "
       "나쁘면 재구성이 아예 시작되지 않을 수 있어 설정할 값어치가 있습니다."),
    DE("Anfangsschätzung der Brennweite, in Pixeln des Ausgangsbildes. 0 liest "
       "EXIF und fällt auf eine Schätzung aus der Bildgröße zurück. Bei einem "
       "Fischauge lohnt es sich, denn eine schlechte Anfangsschätzung kann die "
       "Rekonstruktion ganz verhindern."),
    FR("Estimation de départ de la focale, en pixels de l'image source. 0 lit "
       "l'EXIF et retombe sur une estimation d'après la taille d'image. Utile "
       "pour un fisheye, où une mauvaise valeur de départ peut empêcher la "
       "reconstruction de démarrer."),
    ES("Estimación inicial de la focal, en píxeles de la imagen de origen. 0 "
       "lee el EXIF y recurre a una estimación por el tamaño de imagen. Vale "
       "la pena fijarla en un ojo de pez, donde un mal valor inicial puede "
       "impedir que la reconstrucción arranque."),
    PT("Palpite inicial da distância focal, em pixels da imagem de origem. 0 "
       "lê o EXIF e recorre a uma estimativa pelo tamanho da imagem. Vale "
       "definir num olho de peixe, onde um palpite inicial ruim pode impedir a "
       "reconstrução de começar."),
    IT("Stima iniziale della focale, in pixel dell'immagine di partenza. 0 "
       "legge l'EXIF e ripiega su una stima dalla dimensione dell'immagine. "
       "Conviene impostarla per un fisheye, dove una cattiva stima iniziale "
       "può impedire del tutto l'avvio della ricostruzione."),
    NL("Beginschatting van de brandpuntsafstand, in pixels van het "
       "bronbeeld. 0 leest EXIF en valt terug op een schatting uit de "
       "beeldgrootte. De moeite waard bij een fisheye, waar een slechte "
       "beginschatting de reconstructie helemaal kan blokkeren."),
    RU("Начальная оценка фокусного расстояния в пикселях исходного "
       "изображения. 0 читает EXIF, а при его отсутствии оценивает по размеру "
       "изображения. Для фишая задать стоит: плохая начальная оценка может "
       "вовсе не дать реконструкции начаться."),
    TR("Odak uzaklığı için başlangıç tahmini, kaynak görüntünün pikselleri "
       "cinsinden. 0, EXIF'i okur ve bulamazsa görüntü boyutundan tahmin "
       "eder. Balıkgözünde ayarlamaya değer: kötü bir başlangıç tahmini "
       "yeniden oluşturmanın hiç başlamamasına yol açabilir."));

SS_MSG(max_features_auto,
    EN("Max features per image (0 = auto)"),
    JA("画像あたりの特徴数の上限（0 で自動）"),
    ZH_HANS("每张图像的最大特征数（0 = 自动）"),
    ZH_HANT("每張影像的最大特徵數（0 = 自動）"),
    KO("이미지당 최대 특징점 수(0 = 자동)"),
    DE("Höchstzahl Merkmale je Bild (0 = automatisch)"),
    FR("Points max par image (0 = auto)"),
    ES("Características máximas por imagen (0 = automático)"),
    PT("Máximo de características por imagem (0 = automático)"),
    IT("Caratteristiche massime per immagine (0 = automatico)"),
    NL("Max. kenmerken per beeld (0 = automatisch)"),
    RU("Предел особых точек на снимок (0 — авто)"),
    TR("Görüntü başına en çok öznitelik (0 = otomatik)"));

SS_MSG(max_features_auto_help,
    EN("Keypoints kept per image -- largest scales first for SIFT, highest "
       "detection scores for a learned frontend. Overrides the quality preset "
       "when non-zero. The two are not comparable: SIFT wants tens of "
       "thousands, ALIKED a few thousand."),
    JA("1枚あたりに残すキーポイントの数です。SIFT ではスケールの大きい順、"
       "学習型フロントエンドでは検出スコアの高い順に残します。0 以外なら"
       "品質プリセットより優先されます。両者の数は比べられません。SIFT は"
       "数万、ALIKED は数千を求めます。"),
    ZH_HANS("每张图像保留的关键点数量——SIFT 按尺度从大到小，学习型前端按检测"
            "分数从高到低。非零时会覆盖质量预设。两者的数值不可比：SIFT 要几万个，"
            "ALIKED 只要几千个。"),
    ZH_HANT("每張影像保留的關鍵點數量——SIFT 按尺度從大到小，學習型前端按偵測"
            "分數從高到低。非零時會覆蓋品質預設。兩者的數值不可比：SIFT 要幾萬個，"
            "ALIKED 只要幾千個。"),
    KO("이미지당 남기는 키포인트 수입니다. SIFT는 스케일이 큰 것부터, 학습형 "
       "프런트엔드는 검출 점수가 높은 것부터 남깁니다. 0이 아니면 품질 프리셋보다 "
       "우선합니다. 두 값은 서로 비교할 수 없습니다. SIFT는 수만 개, ALIKED는 "
       "수천 개를 원합니다."),
    DE("Schlüsselpunkte je Bild -- bei SIFT die größten Skalen zuerst, bei "
       "einem gelernten Frontend die höchsten Erkennungswerte. Übergeht die "
       "Qualitätsvorgabe, wenn ungleich null. Die Zahlen sind nicht "
       "vergleichbar: SIFT will Zehntausende, ALIKED ein paar Tausend."),
    FR("Points clés conservés par image -- les plus grandes échelles d'abord "
       "pour SIFT, les meilleurs scores de détection pour un frontal appris. "
       "Prend le pas sur le préréglage de qualité s'il est non nul. Les deux "
       "ne sont pas comparables : SIFT en veut des dizaines de milliers, "
       "ALIKED quelques milliers."),
    ES("Puntos clave conservados por imagen: las escalas mayores primero en "
       "SIFT, las puntuaciones de detección más altas en un frontal "
       "aprendido. Si no es cero, prevalece sobre el ajuste de calidad. Los "
       "dos no son comparables: SIFT quiere decenas de miles, ALIKED unos "
       "pocos miles."),
    PT("Pontos-chave mantidos por imagem -- as maiores escalas primeiro no "
       "SIFT, as maiores pontuações de detecção num front-end aprendido. "
       "Prevalece sobre a predefinição de qualidade quando não é zero. Os dois "
       "não são comparáveis: o SIFT quer dezenas de milhares, o ALIKED alguns "
       "milhares."),
    IT("Punti chiave tenuti per immagine: le scale maggiori prima per SIFT, i "
       "punteggi di rilevamento più alti per un frontend appreso. Se diverso "
       "da zero, prevale sulla preimpostazione di qualità. I due numeri non "
       "sono confrontabili: SIFT ne vuole decine di migliaia, ALIKED qualche "
       "migliaio."),
    NL("Sleutelpunten per beeld -- bij SIFT de grootste schalen eerst, bij een "
       "geleerde frontend de hoogste detectiescores. Gaat boven de "
       "kwaliteitsvoorinstelling als het niet nul is. De twee zijn niet "
       "vergelijkbaar: SIFT wil er tienduizenden, ALIKED een paar duizend."),
    RU("Сколько ключевых точек оставлять на снимок — у SIFT сначала самые "
       "крупные масштабы, у обученного фронтенда самые высокие оценки "
       "детекции. Ненулевое значение важнее пресета качества. Числа несравнимы: "
       "SIFT хочет десятки тысяч, ALIKED — несколько тысяч."),
    TR("Görüntü başına tutulan anahtar nokta sayısı -- SIFT'te önce en büyük "
       "ölçekler, öğrenilmiş bir ön uçta en yüksek bulma puanları. Sıfır "
       "değilse kalite hazır ayarını geçersiz kılar. İkisi kıyaslanabilir "
       "değildir: SIFT on binlerce, ALIKED birkaç bin ister."));

SS_MSG(max_image_size_auto,
    EN("Max image size (0 = auto)"),
    JA("画像サイズの上限（0 で自動）"),
    ZH_HANS("最大图像尺寸（0 = 自动）"),
    ZH_HANT("最大影像尺寸（0 = 自動）"),
    KO("최대 이미지 크기(0 = 자동)"),
    DE("Maximale Bildgröße (0 = automatisch)"),
    FR("Taille d'image max (0 = auto)"),
    ES("Tamaño máximo de imagen (0 = automático)"),
    PT("Tamanho máximo da imagem (0 = automático)"),
    IT("Dimensione massima immagine (0 = automatico)"),
    NL("Max. beeldgrootte (0 = automatisch)"),
    RU("Предел размера изображения (0 — авто)"),
    TR("En büyük görüntü boyutu (0 = otomatik)"));

SS_MSG(max_image_size_auto_help,
    EN("Longest edge the feature extractor runs on; bigger images are "
       "downscaled first. Keypoints are still reported in the source image's "
       "pixels."),
    JA("特徴抽出を行う長辺の長さです。これより大きい画像は先に縮小されます。"
       "キーポイントの座標は元画像のピクセルで報告されます。"),
    ZH_HANS("特征提取所用的最长边长度；更大的图像会先缩小。关键点坐标仍以源图像"
            "的像素给出。"),
    ZH_HANT("特徵擷取所用的最長邊長度；更大的影像會先縮小。關鍵點座標仍以來源影像"
            "的像素給出。"),
    KO("특징 추출을 수행하는 긴 변의 길이입니다. 그보다 큰 이미지는 먼저 "
       "축소됩니다. 키포인트 좌표는 여전히 원본 이미지의 픽셀로 보고됩니다."),
    DE("Längste Kante, auf der die Merkmalsextraktion läuft; größere Bilder "
       "werden zuvor verkleinert. Schlüsselpunkte werden weiterhin in Pixeln "
       "des Ausgangsbildes angegeben."),
    FR("Plus grand côté sur lequel l'extraction de points s'exécute ; les "
       "images plus grandes sont d'abord réduites. Les points clés restent "
       "exprimés en pixels de l'image source."),
    ES("Lado más largo sobre el que se ejecuta la extracción de "
       "características; las imágenes mayores se reducen antes. Los puntos "
       "clave se siguen dando en píxeles de la imagen de origen."),
    PT("Maior lado sobre o qual a extração de características roda; imagens "
       "maiores são reduzidas antes. Os pontos-chave continuam em pixels da "
       "imagem de origem."),
    IT("Lato più lungo su cui gira l'estrazione delle caratteristiche; le "
       "immagini più grandi vengono prima ridotte. I punti chiave restano "
       "espressi in pixel dell'immagine di partenza."),
    NL("Langste zijde waarop de kenmerkextractie draait; grotere beelden "
       "worden eerst verkleind. Sleutelpunten worden nog steeds in pixels van "
       "het bronbeeld gegeven."),
    RU("Наибольшая сторона, на которой работает выделение особых точек; "
       "изображения крупнее сначала уменьшаются. Координаты точек всё равно "
       "даются в пикселях исходного изображения."),
    TR("Öznitelik çıkarımının çalıştığı en uzun kenar; daha büyük görüntüler "
       "önce küçültülür. Anahtar noktalar yine kaynak görüntünün pikselleri "
       "cinsinden bildirilir."));

SS_MSG(keep_intermediate,
    EN("Keep intermediate files"),
    JA("中間ファイルを残す"),
    ZH_HANS("保留中间文件"),
    ZH_HANT("保留中間檔案"),
    KO("중간 파일 남기기"),
    DE("Zwischendateien behalten"),
    FR("Conserver les fichiers intermédiaires"),
    ES("Conservar los archivos intermedios"),
    PT("Manter os arquivos intermediários"),
    IT("Conservare i file intermedi"),
    NL("Tussenbestanden bewaren"),
    RU("Сохранять промежуточные файлы"),
    TR("Ara dosyaları sakla"));

SS_MSG(keep_intermediate_help,
    EN("Keep features/ and matches.bin in the output folder after a "
       "successful run. They are large, and only useful for re-running the "
       "mapper by hand with spirula-sfm."),
    JA("実行が成功したあとも、出力フォルダに features/ と matches.bin を"
       "残します。サイズが大きく、spirula-sfm で手動でマッパーを再実行する"
       "とき以外は使いません。"),
    ZH_HANS("运行成功后仍在输出文件夹里保留 features/ 和 matches.bin。它们体积很大，"
            "只有在用 spirula-sfm 手动重跑建图时才有用。"),
    ZH_HANT("執行成功後仍在輸出資料夾裡保留 features/ 和 matches.bin。它們體積很大，"
            "只有在用 spirula-sfm 手動重跑建圖時才有用。"),
    KO("실행이 성공한 뒤에도 출력 폴더에 features/와 matches.bin을 남깁니다. "
       "크기가 크고, spirula-sfm으로 매퍼를 손수 다시 돌릴 때만 쓸모가 있습니다."),
    DE("features/ und matches.bin nach einem erfolgreichen Lauf im "
       "Ausgabeordner behalten. Sie sind groß und nur nützlich, um den Mapper "
       "von Hand mit spirula-sfm erneut laufen zu lassen."),
    FR("Conserver features/ et matches.bin dans le dossier de sortie après "
       "une exécution réussie. Ils sont volumineux et ne servent qu'à relancer "
       "le mapper à la main avec spirula-sfm."),
    ES("Conservar features/ y matches.bin en la carpeta de salida tras una "
       "ejecución correcta. Son grandes y solo sirven para volver a lanzar el "
       "mapeador a mano con spirula-sfm."),
    PT("Manter features/ e matches.bin na pasta de saída após uma execução "
       "bem-sucedida. São grandes e só servem para rodar o mapeador à mão com "
       "o spirula-sfm."),
    IT("Conservare features/ e matches.bin nella cartella di destinazione dopo "
       "un'esecuzione riuscita. Sono grandi e servono solo per rilanciare a "
       "mano il mapper con spirula-sfm."),
    NL("features/ en matches.bin na een geslaagde run in de uitvoermap "
       "bewaren. Ze zijn groot en alleen nuttig om de mapper met de hand "
       "opnieuw te draaien met spirula-sfm."),
    RU("Оставлять features/ и matches.bin в папке результатов после успешного "
       "запуска. Они большие и нужны, только чтобы вручную перезапустить "
       "маппер через spirula-sfm."),
    TR("Başarılı bir çalıştırmadan sonra features/ ve matches.bin dosyalarını "
       "çıktı klasöründe tutar. Büyüktürler ve yalnızca haritalayıcıyı "
       "spirula-sfm ile elle yeniden çalıştırmak için işe yararlar."));

SS_MSG(extra_sfm_flags_hint,
    EN("extra spirula-sfm flags, e.g. --max-error 2"),
    JA("spirula-sfm への追加オプション（例: --max-error 2）"),
    ZH_HANS("额外的 spirula-sfm 参数，例如 --max-error 2"),
    ZH_HANT("額外的 spirula-sfm 參數，例如 --max-error 2"),
    KO("추가 spirula-sfm 옵션, 예: --max-error 2"),
    DE("zusätzliche spirula-sfm-Optionen, z. B. --max-error 2"),
    FR("options spirula-sfm supplémentaires, p. ex. --max-error 2"),
    ES("opciones adicionales de spirula-sfm, p. ej. --max-error 2"),
    PT("opções adicionais do spirula-sfm, por exemplo --max-error 2"),
    IT("opzioni aggiuntive per spirula-sfm, ad es. --max-error 2"),
    NL("extra spirula-sfm-opties, bijv. --max-error 2"),
    RU("дополнительные ключи spirula-sfm, например --max-error 2"),
    TR("ek spirula-sfm seçenekleri, örn. --max-error 2"));

SS_MSG(extra_sfm_flags_help,
    EN("Passed to `spirula-sfm auto` verbatim. Everything this panel does not "
       "show is reachable here; run `spirula-sfm auto --help` for the list."),
    JA("`spirula-sfm auto` にそのまま渡されます。このパネルに出ていない設定は"
       "すべてここから指定できます。一覧は `spirula-sfm auto --help` で"
       "確認できます。"),
    ZH_HANS("原样传给 `spirula-sfm auto`。这个面板没有列出的一切都可以在这里指定；"
            "运行 `spirula-sfm auto --help` 查看完整列表。"),
    ZH_HANT("原樣傳給 `spirula-sfm auto`。這個面板沒有列出的一切都可以在這裡指定；"
            "執行 `spirula-sfm auto --help` 查看完整清單。"),
    KO("`spirula-sfm auto`에 그대로 전달됩니다. 이 패널에 없는 것은 모두 여기서 "
       "지정할 수 있습니다. 목록은 `spirula-sfm auto --help`로 확인하세요."),
    DE("Wird unverändert an `spirula-sfm auto` weitergereicht. Alles, was "
       "dieses Fenster nicht zeigt, ist hier erreichbar; die Liste liefert "
       "`spirula-sfm auto --help`."),
    FR("Transmis tel quel à `spirula-sfm auto`. Tout ce que ce panneau "
       "n'affiche pas est accessible ici ; lancez `spirula-sfm auto --help` "
       "pour la liste."),
    ES("Se pasa tal cual a `spirula-sfm auto`. Todo lo que este panel no "
       "muestra se alcanza desde aquí; ejecute `spirula-sfm auto --help` para "
       "ver la lista."),
    PT("Repassado tal e qual para `spirula-sfm auto`. Tudo o que este painel "
       "não mostra é alcançável aqui; rode `spirula-sfm auto --help` para ver "
       "a lista."),
    IT("Passato così com'è a `spirula-sfm auto`. Tutto ciò che questo pannello "
       "non mostra è raggiungibile qui; per l'elenco esegua `spirula-sfm auto "
       "--help`."),
    NL("Wordt letterlijk doorgegeven aan `spirula-sfm auto`. Alles wat dit "
       "paneel niet toont, is hier bereikbaar; draai `spirula-sfm auto --help` "
       "voor de lijst."),
    RU("Передаётся в `spirula-sfm auto` как есть. Всё, чего нет на этой панели, "
       "доступно отсюда; список выдаёт `spirula-sfm auto --help`."),
    TR("`spirula-sfm auto` komutuna olduğu gibi aktarılır. Bu panelin "
       "göstermediği her şeye buradan ulaşılır; liste için `spirula-sfm auto "
       "--help` çalıştırın."));

SS_MSG(section_fallbacks,
    EN("Fallbacks"),     JA("フォールバック"), ZH_HANS("回退方式"),  ZH_HANT("回退方式"),
    KO("대체 수단"),      DE("Ausweichwege"), FR("Solutions de repli"),
    ES("Alternativas"),  PT("Alternativas"), IT("Ripieghi"),
    NL("Terugvalopties"), RU("Запасные пути"), TR("Yedek yollar"));

SS_MSG(use_ffmpeg,
    EN("Extract frames with ffmpeg"),
    JA("フレームの切り出しに ffmpeg を使う"),
    ZH_HANS("用 ffmpeg 抽取帧"),
    ZH_HANT("用 ffmpeg 抽取影格"),
    KO("ffmpeg으로 프레임 추출"),
    DE("Bilder mit ffmpeg extrahieren"),
    FR("Extraire les images avec ffmpeg"),
    ES("Extraer los fotogramas con ffmpeg"),
    PT("Extrair os quadros com ffmpeg"),
    IT("Estrarre i fotogrammi con ffmpeg"),
    NL("Beelden met ffmpeg uithalen"),
    RU("Извлекать кадры через ffmpeg"),
    TR("Kareleri ffmpeg ile çıkar"));

SS_MSG(use_ffmpeg_help,
    EN("Use an external ffmpeg instead of decoding on the GPU. Worth trying "
       "for a codec or colour transfer the driver mishandles."),
    JA("GPU でデコードする代わりに外部の ffmpeg を使います。ドライバの扱いが"
       "おかしいコーデックや色変換特性のときに試す価値があります。"),
    ZH_HANS("改用外部 ffmpeg，而不在 GPU 上解码。遇到驱动处理不当的编解码器或"
            "色彩传递特性时值得一试。"),
    ZH_HANT("改用外部 ffmpeg，而不在 GPU 上解碼。遇到驅動處理不當的編解碼器或"
            "色彩傳遞特性時值得一試。"),
    KO("GPU에서 디코딩하는 대신 외부 ffmpeg을 씁니다. 드라이버가 잘못 다루는 "
       "코덱이나 색 전달 특성일 때 시도해 볼 만합니다."),
    DE("Ein externes ffmpeg statt der GPU-Dekodierung benutzen. Einen Versuch "
       "wert bei einem Codec oder einer Farbübertragung, mit denen der Treiber "
       "nicht zurechtkommt."),
    FR("Utiliser un ffmpeg externe plutôt que le décodage sur le GPU. À "
       "essayer pour un codec ou une fonction de transfert que le pilote gère "
       "mal."),
    ES("Usar un ffmpeg externo en lugar de decodificar en la GPU. Vale la pena "
       "probarlo con un códec o una transferencia de color que el controlador "
       "maneje mal."),
    PT("Usar um ffmpeg externo em vez de decodificar na GPU. Vale tentar com "
       "um codec ou uma transferência de cor que o driver trate mal."),
    IT("Usare un ffmpeg esterno invece della decodifica su GPU. Vale la pena "
       "provarlo per un codec o una funzione di trasferimento che il driver "
       "gestisce male."),
    NL("Een externe ffmpeg gebruiken in plaats van decoderen op de GPU. Het "
       "proberen waard bij een codec of kleuroverdracht die het "
       "stuurprogramma verkeerd aanpakt."),
    RU("Использовать внешний ffmpeg вместо декодирования на GPU. Стоит "
       "попробовать при кодеке или передаточной функции цвета, с которыми "
       "драйвер справляется плохо."),
    TR("GPU'da çözmek yerine harici bir ffmpeg kullanın. Sürücünün yanlış "
       "işlediği bir kodek veya renk aktarımı için denemeye değer."));

SS_MSG(use_ffmpeg_always,
    EN("This build always uses ffmpeg for video."),
    JA("このビルドは動画に常に ffmpeg を使います。"),
    ZH_HANS("这个版本处理视频时始终使用 ffmpeg。"),
    ZH_HANT("這個版本處理影片時始終使用 ffmpeg。"),
    KO("이 빌드는 동영상에 항상 ffmpeg을 씁니다."),
    DE("Dieser Build benutzt für Video immer ffmpeg."),
    FR("Cette version utilise toujours ffmpeg pour la vidéo."),
    ES("Esta compilación siempre usa ffmpeg para el vídeo."),
    PT("Esta compilação sempre usa ffmpeg para vídeo."),
    IT("Questa build usa sempre ffmpeg per il video."),
    NL("Deze build gebruikt voor video altijd ffmpeg."),
    RU("Эта сборка всегда использует ffmpeg для видео."),
    TR("Bu sürüm video için her zaman ffmpeg kullanır."));

SS_MSG(use_python_masking,
    EN("Mask with the external Python script"),
    JA("外部の Python スクリプトでマスクする"),
    ZH_HANS("用外部 Python 脚本做蒙版"),
    ZH_HANT("用外部 Python 指令稿做遮罩"),
    KO("외부 Python 스크립트로 마스킹"),
    DE("Mit dem externen Python-Skript maskieren"),
    FR("Masquer avec le script Python externe"),
    ES("Enmascarar con el script externo de Python"),
    PT("Mascarar com o script externo em Python"),
    IT("Mascherare con lo script Python esterno"),
    NL("Maskeren met het externe Python-script"),
    RU("Маскировать внешним скриптом Python"),
    TR("Harici Python betiğiyle maskele"));

SS_MSG(use_python_masking_help,
    EN("Use scripts/mask.py through an external Python with "
       "lang-segment-anything, instead of the built-in segmentation."),
    JA("内蔵のセグメンテーションの代わりに、lang-segment-anything を入れた"
       "外部の Python で scripts/mask.py を使います。"),
    ZH_HANS("不用内置分割，而是通过装有 lang-segment-anything 的外部 Python "
            "运行 scripts/mask.py。"),
    ZH_HANT("不用內建分割，而是透過裝有 lang-segment-anything 的外部 Python "
            "執行 scripts/mask.py。"),
    KO("내장 분할 대신, lang-segment-anything이 설치된 외부 Python으로 "
       "scripts/mask.py를 실행합니다."),
    DE("scripts/mask.py über ein externes Python mit lang-segment-anything "
       "benutzen statt der eingebauten Segmentierung."),
    FR("Utiliser scripts/mask.py via un Python externe doté de "
       "lang-segment-anything, au lieu de la segmentation intégrée."),
    ES("Usar scripts/mask.py mediante un Python externo con "
       "lang-segment-anything, en lugar de la segmentación integrada."),
    PT("Usar scripts/mask.py por meio de um Python externo com "
       "lang-segment-anything, em vez da segmentação integrada."),
    IT("Usare scripts/mask.py tramite un Python esterno con "
       "lang-segment-anything, invece della segmentazione integrata."),
    NL("scripts/mask.py gebruiken via een externe Python met "
       "lang-segment-anything, in plaats van de ingebouwde segmentatie."),
    RU("Использовать scripts/mask.py через внешний Python с "
       "lang-segment-anything вместо встроенной сегментации."),
    TR("Yerleşik bölütleme yerine, lang-segment-anything kurulu harici bir "
       "Python üzerinden scripts/mask.py kullanır."));

// ===========================================================================
// Advanced: external COLMAP
//
// The one block still English-only, and deliberately last: every entry names
// a COLMAP parameter (Mapper.abs_pose_min_num_inliers, SiftMatching.max_ratio)
// and is read alongside COLMAP's own documentation, which exists in English
// only. A translated sentence wrapped around an untranslatable identifier
// helps less here than anywhere else in the app, and this panel is reached by
// someone who has already chosen to drive an external COLMAP by hand.
//
// The parameter names must stay verbatim in every language when these are
// translated; only the prose around them changes.
// ===========================================================================

SS_MSG_EN(colmap_initial_focal, "Initial focal length (x width, 0 = unknown)");
SS_MSG_EN(colmap_initial_focal_help,
    "Seed COLMAP with fx = fy = factor * image width (principal point "
    "centered, zero distortion) instead of its generic guess. A known focal "
    "length stabilizes mapper initialization a lot, especially for fisheye "
    "lenses. Insta360 X5: ~0.269 (set automatically for .insv input).");
SS_MSG_EN(colmap_camera_params, "Initial camera params");
SS_MSG_EN(colmap_camera_params_hint, "fx,fy,cx,cy,... (overrides focal length)");
SS_MSG_EN(colmap_camera_params_help,
    "Raw ImageReader.camera_params for the selected camera model (full "
    "calibration prior). Leave empty to use the focal-length factor above, or "
    "both empty for COLMAP's default initialization.");
SS_MSG_EN(colmap_max_features, "Max features (0 = auto)");
SS_MSG_EN(colmap_max_features_help,
    "SiftExtraction / AlikedExtraction .max_num_features; overrides the "
    "Quality preset when non-zero.");
SS_MSG_EN(colmap_max_image_size, "Max image size (0 = off)");
SS_MSG_EN(colmap_max_image_size_help,
    "FeatureExtraction.max_image_size: downscale images beyond this for "
    "feature extraction.");
SS_MSG_EN(colmap_seq_overlap_help,
    "How many neighboring frames each frame is matched against (sequential "
    "matcher).");
SS_MSG_EN(colmap_quadratic_overlap, "Quadratic overlap");
SS_MSG_EN(colmap_quadratic_overlap_help,
    "Additionally match frame i against frames i +- 2^k (sequential matcher). "
    "Helps close loops in longer captures; enabled by default.");
SS_MSG_EN(colmap_lightglue, "LightGlue matching");
SS_MSG_EN(colmap_lightglue_help,
    "Neural feature matcher (FeatureMatching.type *_LIGHTGLUE): more matches "
    "on hard pairs than brute-force descriptor distance. Default for ALIKED "
    "features; also works with SIFT.");
SS_MSG_EN(colmap_affine_sift, "Affine SIFT + guided matching");
SS_MSG_EN(colmap_affine_sift_help,
    "SiftExtraction.estimate_affine_shape + FeatureMatching.guided_matching: "
    "slower but more robust matching.");
SS_MSG_EN(colmap_distortion_refinement, "Distortion refinement");
SS_MSG_EN(colmap_extra_auto, "Auto");
SS_MSG_EN(colmap_extra_during, "During mapping");
SS_MSG_EN(colmap_extra_final, "Final pass only");
SS_MSG_EN(colmap_distortion_refinement_help,
    "When distortion coefficients are optimized. \"Final pass only\" holds "
    "them fixed during mapping (Mapper.ba_refine_extra_params 0) -- more "
    "stable for low-distortion perspective lenses -- and recovers them in the "
    "final refinement pass. Auto: final-pass-only for perspective models, "
    "during mapping for fisheye.");
SS_MSG_EN(colmap_min_matches, "Min matches per pair (0 = default)");
SS_MSG_EN(colmap_min_matches_help,
    "Mapper.min_num_matches (default 15): image pairs with fewer inlier "
    "matches are ignored by the mapper. Raise to suppress spurious "
    "registrations, lower for sparse overlap.");

SS_MSG_EN(colmap_repetitive, "Repetitive scenes");
SS_MSG_EN(colmap_repetitive_help,
    "Large scenes with repeating structure (several similar rooms, tiled "
    "facades) often weld physically different but similar-looking parts "
    "together. These make matching and registration stricter to suppress "
    "that; 0 = COLMAP default.");
SS_MSG_EN(colmap_repetitive_level, "Repetitive level");
SS_MSG_EN(colmap_rep_off, "Off (COLMAP defaults)");
SS_MSG_EN(colmap_rep_low, "Low");
SS_MSG_EN(colmap_rep_medium, "Medium");
SS_MSG_EN(colmap_rep_high, "High");
SS_MSG_EN(colmap_rep_custom, "Custom");
SS_MSG_EN(colmap_repetitive_level_help,
    "How aggressively wrong matches are suppressed; fills the fields below. "
    "Low: mild tightening, keeps registration rate. Medium: good first "
    "attempt for multi-room indoor captures. High: for heavy repetition "
    "(identical rooms/facades) -- expect fewer registered images if overlap "
    "is thin.");
SS_MSG_EN(colmap_match_ratio, "Match ratio test (0 = default 0.8)");
SS_MSG_EN(colmap_match_ratio_help,
    "SiftMatching.max_ratio, the Lowe ratio test: a feature match is kept "
    "only when its best match is this much better than the second best. LOWER "
    "is stricter -- try 0.6-0.7 when repetitive texture creates false "
    "matches. SIFT only.");
SS_MSG_EN(colmap_min_inliers_pair, "Min inliers per pair (0 = default 15)");
SS_MSG_EN(colmap_min_inliers_pair_help,
    "TwoViewGeometry.min_num_inliers: image pairs whose geometric "
    "verification finds fewer inliers are discarded outright. Raise to 50-100 "
    "so weakly-supported (usually false) links between similar-looking areas "
    "never enter the database.");
SS_MSG_EN(colmap_min_inliers_reg, "Min inliers to register (0 = default 30)");
SS_MSG_EN(colmap_min_inliers_reg_help,
    "Mapper.abs_pose_min_num_inliers: minimum absolute-pose inliers to "
    "register an image into the model. Raise to 50-100 to stop images from "
    "registering onto the wrong (similar-looking) part of the scene.");
SS_MSG_EN(colmap_min_inlier_ratio, "Min inlier ratio to register (0 = default 0.25)");
SS_MSG_EN(colmap_min_inlier_ratio_help,
    "Mapper.abs_pose_min_inlier_ratio: minimum fraction of 2D-3D "
    "correspondences that must be pose inliers. Try 0.35-0.5 for stricter "
    "registration.");
SS_MSG_EN(colmap_max_reg_error, "Max registration error px (0 = default 12)");
SS_MSG_EN(colmap_max_reg_error_help,
    "Mapper.abs_pose_max_error: reprojection error threshold (px) for "
    "absolute-pose RANSAC when registering images. Lower (6-8) = stricter; "
    "combine with the inlier thresholds above.");
SS_MSG_EN(colmap_gpu_ba, "GPU bundle adjustment");
SS_MSG_EN(colmap_gpu_ba_help, "Mapper.ba_use_gpu.");
SS_MSG_EN(colmap_gpu_ba_fisheye,
    "Mapper.ba_use_gpu -- unavailable: COLMAP's GPU bundle adjustment does "
    "not support fisheye camera models yet.");
SS_MSG_EN(colmap_merge_models, "Merge partial models");
SS_MSG_EN(colmap_merge_models_help,
    "When the mapper splits the scene into several partial models, try colmap "
    "model_merger to fuse them (kept only when the merged model registers "
    "more images). The trainer otherwise auto-picks the largest partial "
    "model.");
SS_MSG_EN(colmap_final_ba, "Final refinement pass");
SS_MSG_EN(colmap_final_ba_help,
    "Run bundle_adjuster after mapping on the largest (or merged) model, "
    "refining focal length, principal point, and distortion.");
SS_MSG_EN(colmap_vocab_tree_hint, "vocabulary tree (auto find/download)");
SS_MSG_EN(colmap_vocab_tree, "vocab tree");

// ===========================================================================
// Tool locations
// ===========================================================================

SS_MSG(section_tool_locations,
    EN("Tool locations"), JA("外部ツールの場所"), ZH_HANS("工具位置"), ZH_HANT("工具位置"),
    KO("도구 위치"),      DE("Speicherorte der Werkzeuge"), FR("Emplacement des outils"),
    ES("Ubicación de las herramientas"), PT("Local das ferramentas"),
    IT("Percorsi degli strumenti"), NL("Locatie van hulpprogramma's"),
    RU("Расположение инструментов"), TR("Araç konumları"));

SS_MSG(colmap_executable,
    EN("colmap executable"),
    JA("colmap の実行ファイル"),
    ZH_HANS("colmap 可执行文件"),
    ZH_HANT("colmap 執行檔"),
    KO("colmap 실행 파일"),
    DE("colmap-Programmdatei"),
    FR("exécutable colmap"),
    ES("ejecutable de colmap"),
    PT("executável do colmap"),
    IT("eseguibile colmap"),
    NL("colmap-programma"),
    RU("исполняемый файл colmap"),
    TR("colmap çalıştırılabiliri"));

SS_MSG(ffmpeg_executable,
    EN("ffmpeg executable"),
    JA("ffmpeg の実行ファイル"),
    ZH_HANS("ffmpeg 可执行文件"),
    ZH_HANT("ffmpeg 執行檔"),
    KO("ffmpeg 실행 파일"),
    DE("ffmpeg-Programmdatei"),
    FR("exécutable ffmpeg"),
    ES("ejecutable de ffmpeg"),
    PT("executável do ffmpeg"),
    IT("eseguibile ffmpeg"),
    NL("ffmpeg-programma"),
    RU("исполняемый файл ffmpeg"),
    TR("ffmpeg çalıştırılabiliri"));

SS_MSG(ffmpeg_executable_help_fallback,
    EN("Only used when frame extraction falls back to ffmpeg."),
    JA("フレームの切り出しが ffmpeg にフォールバックしたときだけ使われます。"),
    ZH_HANS("只有当抽帧回退到 ffmpeg 时才会用到。"),
    ZH_HANT("只有當抽格回退到 ffmpeg 時才會用到。"),
    KO("프레임 추출이 ffmpeg으로 넘어갈 때만 쓰입니다."),
    DE("Wird nur benutzt, wenn die Bildextraktion auf ffmpeg zurückfällt."),
    FR("Utilisé seulement quand l'extraction des images retombe sur ffmpeg."),
    ES("Solo se usa cuando la extracción de fotogramas recurre a ffmpeg."),
    PT("Só é usado quando a extração de quadros recorre ao ffmpeg."),
    IT("Usato solo quando l'estrazione dei fotogrammi ripiega su ffmpeg."),
    NL("Wordt alleen gebruikt als het uithalen van beelden terugvalt op "
       "ffmpeg."),
    RU("Используется, только если извлечение кадров переходит на ffmpeg."),
    TR("Yalnızca kare çıkarma ffmpeg'e düştüğünde kullanılır."));

SS_MSG(ffmpeg_executable_help_always,
    EN("Used to extract frames from video."),
    JA("動画からフレームを切り出すのに使われます。"),
    ZH_HANS("用于从视频中抽取帧。"),
    ZH_HANT("用於從影片中抽取影格。"),
    KO("동영상에서 프레임을 뽑는 데 쓰입니다."),
    DE("Wird zum Extrahieren der Bilder aus dem Video benutzt."),
    FR("Sert à extraire les images de la vidéo."),
    ES("Se usa para extraer fotogramas del vídeo."),
    PT("Usado para extrair quadros do vídeo."),
    IT("Serve a estrarre i fotogrammi dal video."),
    NL("Wordt gebruikt om beelden uit de video te halen."),
    RU("Используется для извлечения кадров из видео."),
    TR("Videodan kare çıkarmak için kullanılır."));

SS_MSG(python_executable,
    EN("python executable"),
    JA("python の実行ファイル"),
    ZH_HANS("python 可执行文件"),
    ZH_HANT("python 執行檔"),
    KO("python 실행 파일"),
    DE("python-Programmdatei"),
    FR("exécutable python"),
    ES("ejecutable de python"),
    PT("executável do python"),
    IT("eseguibile python"),
    NL("python-programma"),
    RU("исполняемый файл python"),
    TR("python çalıştırılabiliri"));

SS_MSG(python_executable_help,
    EN("Only used by the external masking script."),
    JA("外部のマスキングスクリプトだけが使います。"),
    ZH_HANS("只有外部的蒙版脚本会用到。"),
    ZH_HANT("只有外部的遮罩指令稿會用到。"),
    KO("외부 마스킹 스크립트만 사용합니다."),
    DE("Wird nur vom externen Maskierungsskript benutzt."),
    FR("Utilisé uniquement par le script de masquage externe."),
    ES("Solo lo usa el script externo de enmascarado."),
    PT("Usado apenas pelo script externo de mascaramento."),
    IT("Usato solo dallo script esterno di mascheratura."),
    NL("Wordt alleen gebruikt door het externe maskeerscript."),
    RU("Используется только внешним скриптом маскирования."),
    TR("Yalnızca harici maskeleme betiği kullanır."));

// ===========================================================================
// The segmentation checkpoints (src/app/gui/ModelCache.cpp)
//
// The blurbs all quote the same measurement -- one instance, 1080p frames,
// laptop GPU -- so that they can actually be compared. Keep that parallel
// structure when translating; it is the only thing that makes the choice
// makeable.
// ===========================================================================

SS_MSG(model_sam3_label,
    EN("SAM 3 (recommended)"),
    JA("SAM 3（推奨）"),  ZH_HANS("SAM 3（推荐）"), ZH_HANT("SAM 3（建議）"),
    KO("SAM 3(권장)"),   DE("SAM 3 (empfohlen)"), FR("SAM 3 (recommandé)"),
    ES("SAM 3 (recomendado)"), PT("SAM 3 (recomendado)"), IT("SAM 3 (consigliato)"),
    NL("SAM 3 (aanbevolen)"), RU("SAM 3 (рекомендуется)"), TR("SAM 3 (önerilen)"));

SS_MSG(model_sam3_blurb,
    EN("Understands text prompts -- type what to mask out. 707 MB, ~2 GB "
       "VRAM, about 1 s per frame on a laptop GPU -- 3x slower than any "
       "SAM 2.1 below."),
    JA("テキストのプロンプトを理解します。消したいものを入力してください。"
       "707 MB、VRAM 約 2 GB、ノート PC の GPU で 1 フレームおよそ 1 秒。"
       "下の SAM 2.1 のどれよりも 3 倍遅くなります。"),
    ZH_HANS("能理解文字提示——直接输入要遮掉什么。707 MB，约 2 GB 显存，"
            "笔记本 GPU 上每帧约 1 秒——比下面任何一个 SAM 2.1 都慢 3 倍。"),
    ZH_HANT("能理解文字提示——直接輸入要遮掉什麼。707 MB，約 2 GB 顯示記憶體，"
            "筆電 GPU 上每格約 1 秒——比下面任何一個 SAM 2.1 都慢 3 倍。"),
    KO("텍스트 프롬프트를 이해합니다. 가릴 것을 입력하세요. 707 MB, VRAM 약 "
       "2 GB, 노트북 GPU에서 프레임당 약 1초 — 아래 어떤 SAM 2.1보다도 3배 "
       "느립니다."),
    DE("Versteht Texteingaben -- schreiben Sie hinein, was maskiert werden "
       "soll. 707 MB, etwa 2 GB VRAM, rund 1 s je Bild auf einer Laptop-GPU "
       "-- dreimal langsamer als jedes SAM 2.1 darunter."),
    FR("Comprend les invites textuelles -- écrivez ce qu'il faut masquer. "
       "707 Mo, environ 2 Go de VRAM, à peu près 1 s par image sur un GPU "
       "d'ordinateur portable -- 3 fois plus lent que n'importe quel SAM 2.1 "
       "ci-dessous."),
    ES("Entiende indicaciones de texto: escriba qué enmascarar. 707 MB, unos "
       "2 GB de VRAM, alrededor de 1 s por fotograma en una GPU de portátil: "
       "3 veces más lento que cualquier SAM 2.1 de abajo."),
    PT("Entende comandos de texto -- escreva o que mascarar. 707 MB, cerca de "
       "2 GB de VRAM, aproximadamente 1 s por quadro numa GPU de notebook -- "
       "3 vezes mais lento que qualquer SAM 2.1 abaixo."),
    IT("Capisce il testo: scriva che cosa mascherare. 707 MB, circa 2 GB di "
       "VRAM, all'incirca 1 s per fotogramma su una GPU da portatile: 3 volte "
       "più lento di qualsiasi SAM 2.1 qui sotto."),
    NL("Begrijpt tekstprompts -- typ wat gemaskeerd moet worden. 707 MB, "
       "ongeveer 2 GB VRAM, zo'n 1 s per beeld op een laptop-GPU -- 3 keer "
       "trager dan elke SAM 2.1 hieronder."),
    RU("Понимает текстовые запросы — напишите, что замаскировать. 707 МБ, "
       "около 2 ГБ видеопамяти, примерно 1 с на кадр на ноутбучной "
       "видеокарте — втрое медленнее любой SAM 2.1 ниже."),
    TR("Metin istemlerini anlar -- neyin maskeleneceğini yazın. 707 MB, ~2 GB "
       "VRAM, dizüstü GPU'da kare başına yaklaşık 1 sn -- aşağıdaki her "
       "SAM 2.1'den 3 kat yavaş."));

SS_MSG(model_sam3_f16_label,
    EN("SAM 3, full precision"),
    JA("SAM 3、フル精度"), ZH_HANS("SAM 3，全精度"), ZH_HANT("SAM 3，全精度"),
    KO("SAM 3, 전체 정밀도"), DE("SAM 3, volle Genauigkeit"),
    FR("SAM 3, pleine précision"), ES("SAM 3, precisión completa"),
    PT("SAM 3, precisão total"), IT("SAM 3, piena precisione"),
    NL("SAM 3, volledige precisie"), RU("SAM 3, полная точность"),
    TR("SAM 3, tam duyarlık"));

SS_MSG(model_sam3_f16_blurb,
    EN("The same model without file quantization. Slightly better masks, much "
       "bigger download, same speed."),
    JA("同じモデルのファイル量子化なし版です。マスクはわずかに良くなり、"
       "ダウンロードはずっと大きく、速度は同じです。"),
    ZH_HANS("同一个模型，不做文件量化。蒙版略好，下载大得多，速度相同。"),
    ZH_HANT("同一個模型，不做檔案量化。遮罩略好，下載大得多，速度相同。"),
    KO("같은 모델의 파일 양자화를 하지 않은 판입니다. 마스크가 조금 낫고, "
       "내려받기는 훨씬 크며, 속도는 같습니다."),
    DE("Dasselbe Modell ohne Dateiquantisierung. Etwas bessere Masken, viel "
       "größerer Download, gleiche Geschwindigkeit."),
    FR("Le même modèle sans quantification du fichier. Masques légèrement "
       "meilleurs, téléchargement bien plus lourd, même vitesse."),
    ES("El mismo modelo sin cuantización del archivo. Máscaras algo mejores, "
       "descarga mucho mayor, misma velocidad."),
    PT("O mesmo modelo sem quantização do arquivo. Máscaras um pouco "
       "melhores, download bem maior, mesma velocidade."),
    IT("Lo stesso modello senza quantizzazione del file. Maschere un po' "
       "migliori, scaricamento molto più grande, stessa velocità."),
    NL("Hetzelfde model zonder bestandskwantisatie. Iets betere maskers, veel "
       "grotere download, dezelfde snelheid."),
    RU("Та же модель без квантования файла. Маски чуть лучше, загрузка "
       "гораздо больше, скорость та же."),
    TR("Aynı modelin dosya nicemlemesi olmayan hâli. Maskeler biraz daha iyi, "
       "indirme çok daha büyük, hız aynı."));

SS_MSG(model_sam21_large_label,
    EN("SAM 2.1 Large"),  JA("SAM 2.1 Large"), ZH_HANS("SAM 2.1 Large"),
    ZH_HANT("SAM 2.1 Large"), KO("SAM 2.1 Large"), DE("SAM 2.1 Large"),
    FR("SAM 2.1 Large"),  ES("SAM 2.1 Large"), PT("SAM 2.1 Large"),
    IT("SAM 2.1 Large"),  NL("SAM 2.1 Large"), RU("SAM 2.1 Large"),
    TR("SAM 2.1 Large"));

SS_MSG(model_sam21_large_blurb,
    EN("Click or draw a box to select an object; no text prompts. The most "
       "accurate of the four and the one to pick for thin structure -- hair, "
       "railings, foliage. ~470 ms per frame. Apache-2.0."),
    JA("クリックまたは矩形で対象を選びます。テキストのプロンプトはありません。"
       "4つの中でいちばん正確で、髪、手すり、葉のような細い構造にはこれを"
       "選んでください。1 フレームおよそ 470 ms。Apache-2.0。"),
    ZH_HANS("用点击或拉框来选对象；不支持文字提示。四者中最准确，头发、栏杆、"
            "枝叶这类细结构就选它。每帧约 470 毫秒。Apache-2.0。"),
    ZH_HANT("用點擊或拉框來選物件；不支援文字提示。四者中最準確，頭髮、欄杆、"
            "枝葉這類細結構就選它。每格約 470 毫秒。Apache-2.0。"),
    KO("클릭하거나 상자를 그려 물체를 고릅니다. 텍스트 프롬프트는 없습니다. "
       "넷 중 가장 정확하며 머리카락, 난간, 잎사귀 같은 가느다란 구조에는 이것을 "
       "고르세요. 프레임당 약 470 ms. Apache-2.0."),
    DE("Zum Auswählen anklicken oder einen Rahmen ziehen; keine Texteingaben. "
       "Das genaueste der vier und die Wahl für feine Strukturen -- Haare, "
       "Geländer, Laub. Etwa 470 ms je Bild. Apache-2.0."),
    FR("Cliquez ou tracez un cadre pour sélectionner un objet ; pas d'invite "
       "textuelle. Le plus précis des quatre et celui à prendre pour les "
       "structures fines -- cheveux, garde-corps, feuillage. Environ 470 ms "
       "par image. Apache-2.0."),
    ES("Haga clic o dibuje un recuadro para elegir un objeto; sin "
       "indicaciones de texto. El más preciso de los cuatro y el indicado "
       "para estructuras finas: pelo, barandillas, follaje. Unos 470 ms por "
       "fotograma. Apache-2.0."),
    PT("Clique ou desenhe uma caixa para escolher um objeto; sem comandos de "
       "texto. O mais preciso dos quatro e o indicado para estruturas finas: "
       "cabelo, corrimãos, folhagem. Cerca de 470 ms por quadro. Apache-2.0."),
    IT("Clicchi o tracci un rettangolo per scegliere un oggetto; niente "
       "testo. Il più preciso dei quattro e quello da prendere per le "
       "strutture sottili: capelli, ringhiere, fogliame. Circa 470 ms per "
       "fotogramma. Apache-2.0."),
    NL("Klik of trek een kader om een object te kiezen; geen tekstprompts. Het "
       "nauwkeurigste van de vier en de keuze voor fijne structuur -- haar, "
       "leuningen, gebladerte. Ongeveer 470 ms per beeld. Apache-2.0."),
    RU("Щелчок или рамка выбирают объект; текстовых запросов нет. Самая точная "
       "из четырёх и та, что нужна для тонких структур — волос, перил, "
       "листвы. Около 470 мс на кадр. Apache-2.0."),
    TR("Nesne seçmek için tıklayın veya kutu çizin; metin istemi yok. Dördü "
       "arasında en doğru olanı ve ince yapılar için seçilecek olanı -- saç, "
       "korkuluk, yaprak. Kare başına ~470 ms. Apache-2.0."));

SS_MSG(model_sam21_baseplus_label,
    EN("SAM 2.1 Base+"),  JA("SAM 2.1 Base+"), ZH_HANS("SAM 2.1 Base+"),
    ZH_HANT("SAM 2.1 Base+"), KO("SAM 2.1 Base+"), DE("SAM 2.1 Base+"),
    FR("SAM 2.1 Base+"),  ES("SAM 2.1 Base+"), PT("SAM 2.1 Base+"),
    IT("SAM 2.1 Base+"),  NL("SAM 2.1 Base+"), RU("SAM 2.1 Base+"),
    TR("SAM 2.1 Base+"));

SS_MSG(model_sam21_baseplus_blurb,
    EN("Clicks and boxes only. Close to Large on most subjects at two thirds "
       "the time, ~320 ms per frame. Apache-2.0."),
    JA("クリックと矩形のみです。ほとんどの被写体で Large に近い結果を3分の2の"
       "時間で出します。1 フレームおよそ 320 ms。Apache-2.0。"),
    ZH_HANS("只支持点击和拉框。在多数对象上接近 Large，耗时只有三分之二，"
            "每帧约 320 毫秒。Apache-2.0。"),
    ZH_HANT("只支援點擊和拉框。在多數物件上接近 Large，耗時只有三分之二，"
            "每格約 320 毫秒。Apache-2.0。"),
    KO("클릭과 상자만 됩니다. 대부분의 피사체에서 Large에 가까운 결과를 3분의 2 "
       "시간에 냅니다. 프레임당 약 320 ms. Apache-2.0."),
    DE("Nur Klicks und Rahmen. Bei den meisten Motiven nah an Large, in zwei "
       "Dritteln der Zeit, etwa 320 ms je Bild. Apache-2.0."),
    FR("Clics et cadres seulement. Proche de Large sur la plupart des sujets "
       "en deux tiers du temps, environ 320 ms par image. Apache-2.0."),
    ES("Solo clics y recuadros. Cerca de Large en la mayoría de sujetos en dos "
       "tercios del tiempo, unos 320 ms por fotograma. Apache-2.0."),
    PT("Só cliques e caixas. Perto do Large na maioria dos sujeitos em dois "
       "terços do tempo, cerca de 320 ms por quadro. Apache-2.0."),
    IT("Solo clic e rettangoli. Vicino a Large su quasi tutti i soggetti in "
       "due terzi del tempo, circa 320 ms per fotogramma. Apache-2.0."),
    NL("Alleen klikken en kaders. Dicht bij Large op de meeste onderwerpen in "
       "twee derde van de tijd, ongeveer 320 ms per beeld. Apache-2.0."),
    RU("Только щелчки и рамки. На большинстве объектов близко к Large за две "
       "трети времени, около 320 мс на кадр. Apache-2.0."),
    TR("Yalnızca tıklama ve kutu. Çoğu öznede Large'a yakın, sürenin üçte "
       "ikisinde, kare başına ~320 ms. Apache-2.0."));

SS_MSG(model_sam21_small_label,
    EN("SAM 2.1 Small"),  JA("SAM 2.1 Small"), ZH_HANS("SAM 2.1 Small"),
    ZH_HANT("SAM 2.1 Small"), KO("SAM 2.1 Small"), DE("SAM 2.1 Small"),
    FR("SAM 2.1 Small"),  ES("SAM 2.1 Small"), PT("SAM 2.1 Small"),
    IT("SAM 2.1 Small"),  NL("SAM 2.1 Small"), RU("SAM 2.1 Small"),
    TR("SAM 2.1 Small"));

SS_MSG(model_sam21_small_blurb,
    EN("Clicks and boxes only, ~255 ms per frame. The best speed-for-quality "
       "of the four: below Large the frame is mostly tracking, which does not "
       "care how big the backbone is. Apache-2.0."),
    JA("クリックと矩形のみ、1 フレームおよそ 255 ms。4つの中で速度対品質が"
       "いちばん良い選択です。Large 未満ではフレームの処理はほとんど追跡で、"
       "バックボーンの大きさはあまり効きません。Apache-2.0。"),
    ZH_HANS("只支持点击和拉框，每帧约 255 毫秒。四者中速度与质量的平衡最好："
            "在 Large 以下，每帧的工作主要是跟踪，而跟踪并不在乎主干有多大。"
            "Apache-2.0。"),
    ZH_HANT("只支援點擊和拉框，每格約 255 毫秒。四者中速度與品質的平衡最好："
            "在 Large 以下，每格的工作主要是追蹤，而追蹤並不在乎骨幹有多大。"
            "Apache-2.0。"),
    KO("클릭과 상자만 되며 프레임당 약 255 ms. 넷 중 속도 대비 품질이 가장 "
       "좋습니다. Large 아래에서는 프레임 처리의 대부분이 추적이고, 추적은 "
       "백본 크기에 크게 좌우되지 않습니다. Apache-2.0."),
    DE("Nur Klicks und Rahmen, etwa 255 ms je Bild. Das beste Verhältnis von "
       "Tempo zu Qualität der vier: unterhalb von Large ist ein Bild vor allem "
       "Nachverfolgung, und der ist die Größe des Rückgrats fast egal. "
       "Apache-2.0."),
    FR("Clics et cadres seulement, environ 255 ms par image. Le meilleur "
       "rapport vitesse/qualité des quatre : en dessous de Large, le travail "
       "par image est surtout du suivi, qui se moque de la taille du réseau. "
       "Apache-2.0."),
    ES("Solo clics y recuadros, unos 255 ms por fotograma. La mejor relación "
       "velocidad-calidad de los cuatro: por debajo de Large, el trabajo por "
       "fotograma es sobre todo seguimiento, al que le da igual el tamaño de "
       "la red. Apache-2.0."),
    PT("Só cliques e caixas, cerca de 255 ms por quadro. A melhor relação "
       "velocidade/qualidade dos quatro: abaixo do Large, o trabalho por "
       "quadro é sobretudo rastreamento, que não liga para o tamanho da rede. "
       "Apache-2.0."),
    IT("Solo clic e rettangoli, circa 255 ms per fotogramma. Il miglior "
       "rapporto velocità/qualità dei quattro: sotto Large il lavoro per "
       "fotogramma è soprattutto inseguimento, a cui la dimensione della rete "
       "importa poco. Apache-2.0."),
    NL("Alleen klikken en kaders, ongeveer 255 ms per beeld. De beste "
       "verhouding snelheid/kwaliteit van de vier: onder Large is het werk per "
       "beeld vooral volgen, en dat maalt niet om de grootte van het netwerk. "
       "Apache-2.0."),
    RU("Только щелчки и рамки, около 255 мс на кадр. Лучшее соотношение "
       "скорости и качества из четырёх: ниже Large работа над кадром — это в "
       "основном отслеживание, которому размер сети почти безразличен. "
       "Apache-2.0."),
    TR("Yalnızca tıklama ve kutu, kare başına ~255 ms. Dördü arasında hız/"
       "kalite dengesi en iyi olanı: Large'ın altında kare başına iş çoğunlukla "
       "izlemedir ve izleme omurganın büyüklüğüne pek aldırmaz. Apache-2.0."));

SS_MSG(model_sam21_tiny_label,
    EN("SAM 2.1 Tiny (fastest)"),
    JA("SAM 2.1 Tiny（最速）"), ZH_HANS("SAM 2.1 Tiny（最快）"),
    ZH_HANT("SAM 2.1 Tiny（最快）"), KO("SAM 2.1 Tiny(가장 빠름)"),
    DE("SAM 2.1 Tiny (am schnellsten)"), FR("SAM 2.1 Tiny (le plus rapide)"),
    ES("SAM 2.1 Tiny (el más rápido)"), PT("SAM 2.1 Tiny (o mais rápido)"),
    IT("SAM 2.1 Tiny (il più veloce)"), NL("SAM 2.1 Tiny (snelste)"),
    RU("SAM 2.1 Tiny (самая быстрая)"), TR("SAM 2.1 Tiny (en hızlı)"));

SS_MSG(model_sam21_tiny_blurb,
    EN("Clicks and boxes only, 76 MB. Only 4% quicker than Small at ~245 ms "
       "per frame, and it loses thin structure first -- take it for the "
       "download size, not the speed. Apache-2.0."),
    JA("クリックと矩形のみ、76 MB。1 フレームおよそ 245 ms で Small より 4% "
       "速いだけですし、細い構造から先に失われます。速度ではなく"
       "ダウンロードサイズのために選んでください。Apache-2.0。"),
    ZH_HANS("只支持点击和拉框，76 MB。每帧约 245 毫秒，只比 Small 快 4%，"
            "而且最先丢失细结构——选它是为了下载体积，不是速度。Apache-2.0。"),
    ZH_HANT("只支援點擊和拉框，76 MB。每格約 245 毫秒，只比 Small 快 4%，"
            "而且最先丟失細結構——選它是為了下載體積，不是速度。Apache-2.0。"),
    KO("클릭과 상자만 되며 76 MB. 프레임당 약 245 ms로 Small보다 4% 빠를 뿐이고 "
       "가느다란 구조를 가장 먼저 잃습니다. 속도가 아니라 내려받기 크기 때문에 "
       "고르세요. Apache-2.0."),
    DE("Nur Klicks und Rahmen, 76 MB. Mit etwa 245 ms je Bild nur 4 % "
       "schneller als Small, und feine Strukturen gehen zuerst verloren -- "
       "wegen der Downloadgröße nehmen, nicht wegen des Tempos. Apache-2.0."),
    FR("Clics et cadres seulement, 76 Mo. À environ 245 ms par image, à peine "
       "4 % plus rapide que Small, et c'est lui qui perd les structures fines "
       "en premier -- à prendre pour la taille du téléchargement, pas pour la "
       "vitesse. Apache-2.0."),
    ES("Solo clics y recuadros, 76 MB. Con unos 245 ms por fotograma, apenas "
       "un 4 % más rápido que Small, y es el primero en perder las "
       "estructuras finas: tómelo por el tamaño de descarga, no por la "
       "velocidad. Apache-2.0."),
    PT("Só cliques e caixas, 76 MB. Com cerca de 245 ms por quadro, apenas 4% "
       "mais rápido que o Small, e é o primeiro a perder estruturas finas -- "
       "escolha pelo tamanho do download, não pela velocidade. Apache-2.0."),
    IT("Solo clic e rettangoli, 76 MB. A circa 245 ms per fotogramma è appena "
       "il 4% più veloce di Small, ed è il primo a perdere le strutture "
       "sottili: lo prenda per la dimensione del download, non per la "
       "velocità. Apache-2.0."),
    NL("Alleen klikken en kaders, 76 MB. Met ongeveer 245 ms per beeld maar "
       "4% sneller dan Small, en het verliest fijne structuur het eerst -- "
       "neem het om de downloadgrootte, niet om de snelheid. Apache-2.0."),
    RU("Только щелчки и рамки, 76 МБ. При примерно 245 мс на кадр она быстрее "
       "Small всего на 4 % и первой теряет тонкие структуры — берите её ради "
       "размера загрузки, а не скорости. Apache-2.0."),
    TR("Yalnızca tıklama ve kutu, 76 MB. Kare başına ~245 ms ile Small'dan "
       "yalnızca %4 hızlı ve ince yapıyı ilk kaybeden o -- indirme boyutu "
       "için alın, hız için değil. Apache-2.0."));

// LEGAL -- human review in every language, see the block below.
SS_MSG(license_sam3_title,
    EN("SAM 3 License (Meta)"),
    JA("SAM 3 ライセンス（Meta）"),
    ZH_HANS("SAM 3 许可协议（Meta）"),
    ZH_HANT("SAM 3 授權條款（Meta）"),
    KO("SAM 3 라이선스(Meta)"),
    DE("SAM-3-Lizenz (Meta)"),
    FR("Licence SAM 3 (Meta)"),
    ES("Licencia de SAM 3 (Meta)"),
    PT("Licença do SAM 3 (Meta)"),
    IT("Licenza SAM 3 (Meta)"),
    NL("SAM 3-licentie (Meta)"),
    RU("Лицензия SAM 3 (Meta)"),
    TR("SAM 3 Lisansı (Meta)"));

SS_MSG(license_sam3_summary,
    EN("SAM 3 is Meta's model, not part of Spirula Studio, and it comes with "
       "its own licence -- which is not a standard one. It is free to use, "
       "including commercially, but only on Meta's terms, so we cannot ship "
       "it with the app or accept them for you.\n\n"
       "Please read it before continuing -- it is short, and it is the actual "
       "agreement, not this summary of it."),
    JA("SAM 3 は Meta のモデルで、Spirula Studio の一部ではなく、独自の"
       "ライセンスが付いています。それは標準的なライセンスではありません。"
       "商用を含めて無償で使えますが、あくまで Meta の条件のもとでです。"
       "そのため、当アプリに同梱することも、条件への同意を代行することも"
       "できません。\n\n"
       "続ける前にお読みください。短い文書ですし、実際の契約はこの要約では"
       "なくそちらです。"),
    ZH_HANS("SAM 3 是 Meta 的模型，不属于 Spirula Studio，并且带有它自己的"
            "许可协议——那不是一份标准协议。它可以免费使用，包括商业用途，"
            "但只在 Meta 的条件之下。因此我们既不能随应用一起分发它，"
            "也不能代你接受这些条件。\n\n"
            "请在继续之前阅读它——它很短，而且真正的协议是它，不是这段摘要。"),
    ZH_HANT("SAM 3 是 Meta 的模型，不屬於 Spirula Studio，並且帶有它自己的"
            "授權條款——那不是一份標準條款。它可以免費使用，包括商業用途，"
            "但只在 Meta 的條件之下。因此我們既不能隨應用一起散布它，"
            "也不能代你接受這些條件。\n\n"
            "請在繼續之前閱讀它——它很短，而且真正的協議是它，不是這段摘要。"),
    KO("SAM 3는 Meta의 모델로 Spirula Studio의 일부가 아니며, 자체 라이선스가 "
       "딸려 있습니다. 그것은 표준 라이선스가 아닙니다. 상업적 사용을 포함해 "
       "무료로 쓸 수 있지만 어디까지나 Meta의 조건 아래에서입니다. 그래서 저희는 "
       "이 모델을 앱과 함께 배포할 수도, 조건을 대신 수락할 수도 없습니다.\n\n"
       "계속하기 전에 읽어 주세요. 길지 않으며, 실제 계약은 이 요약이 아니라 "
       "그 문서입니다."),
    DE("SAM 3 ist Metas Modell, nicht Teil von Spirula Studio, und bringt "
       "eine eigene Lizenz mit -- keine übliche. Es ist kostenlos nutzbar, "
       "auch kommerziell, aber nur zu Metas Bedingungen; wir dürfen es daher "
       "weder mit der Anwendung ausliefern noch die Bedingungen für Sie "
       "annehmen.\n\n"
       "Bitte lesen Sie sie, bevor Sie fortfahren -- sie ist kurz, und sie "
       "ist die eigentliche Vereinbarung, nicht diese Zusammenfassung."),
    FR("SAM 3 est le modèle de Meta, il ne fait pas partie de Spirula Studio "
       "et il vient avec sa propre licence -- qui n'est pas une licence "
       "standard. Il est gratuit à utiliser, y compris commercialement, mais "
       "uniquement aux conditions de Meta ; nous ne pouvons donc ni le livrer "
       "avec l'application ni les accepter à votre place.\n\n"
       "Merci de la lire avant de continuer : elle est courte, et c'est elle "
       "l'accord véritable, pas ce résumé."),
    ES("SAM 3 es el modelo de Meta, no forma parte de Spirula Studio y viene "
       "con su propia licencia, que no es una licencia estándar. Su uso es "
       "gratuito, también comercial, pero solo en los términos de Meta, así "
       "que no podemos distribuirlo con la aplicación ni aceptarlos por "
       "usted.\n\n"
       "Léala antes de continuar: es breve, y es ella el acuerdo real, no "
       "este resumen."),
    PT("O SAM 3 é o modelo da Meta, não faz parte do Spirula Studio e vem com "
       "a própria licença -- que não é uma licença padrão. É gratuito, "
       "inclusive para uso comercial, mas só nos termos da Meta, então não "
       "podemos distribuí-lo com o aplicativo nem aceitá-los por você.\n\n"
       "Leia-a antes de continuar: é curta, e é ela o acordo de verdade, não "
       "este resumo."),
    IT("SAM 3 è il modello di Meta, non fa parte di Spirula Studio e ha una "
       "licenza propria, che non è una licenza standard. È gratuito, anche "
       "per uso commerciale, ma solo alle condizioni di Meta: non possiamo "
       "quindi distribuirlo con l'applicazione né accettarle al posto suo.\n\n"
       "La legga prima di proseguire: è breve, ed è lei l'accordo vero, non "
       "questo riassunto."),
    NL("SAM 3 is het model van Meta, hoort niet bij Spirula Studio en komt met "
       "een eigen licentie -- geen standaardlicentie. Het is gratis te "
       "gebruiken, ook commercieel, maar alleen op Meta's voorwaarden; we "
       "mogen het dus niet met de toepassing meeleveren en ze ook niet voor u "
       "aanvaarden.\n\n"
       "Lees ze voordat u doorgaat: ze zijn kort, en zij vormen de "
       "werkelijke overeenkomst, niet deze samenvatting."),
    RU("SAM 3 — модель Meta, она не входит в Spirula Studio и поставляется со "
       "своей лицензией, а она не стандартная. Пользоваться моделью можно "
       "бесплатно, в том числе коммерчески, но только на условиях Meta, "
       "поэтому мы не вправе ни поставлять её вместе с программой, ни "
       "принимать эти условия за вас.\n\n"
       "Прочитайте её, прежде чем продолжить: она короткая, и настоящее "
       "соглашение — это она, а не данная выжимка."),
    TR("SAM 3 Meta'nın modelidir, Spirula Studio'nun parçası değildir ve kendi "
       "lisansıyla gelir -- bu standart bir lisans değildir. Ticari kullanım "
       "dâhil ücretsizdir, ama yalnızca Meta'nın koşullarıyla; bu yüzden onu "
       "uygulamayla birlikte dağıtamayız ve koşulları sizin adınıza kabul "
       "edemeyiz.\n\n"
       "Devam etmeden önce lütfen okuyun -- kısadır ve gerçek sözleşme bu "
       "özet değil, o metindir."));

SS_MSG(license_sam2_title,
    EN("SAM 2.1 License (Meta, Apache-2.0)"),
    JA("SAM 2.1 ライセンス（Meta、Apache-2.0）"),
    ZH_HANS("SAM 2.1 许可协议（Meta，Apache-2.0）"),
    ZH_HANT("SAM 2.1 授權條款（Meta，Apache-2.0）"),
    KO("SAM 2.1 라이선스(Meta, Apache-2.0)"),
    DE("SAM-2.1-Lizenz (Meta, Apache-2.0)"),
    FR("Licence SAM 2.1 (Meta, Apache-2.0)"),
    ES("Licencia de SAM 2.1 (Meta, Apache-2.0)"),
    PT("Licença do SAM 2.1 (Meta, Apache-2.0)"),
    IT("Licenza SAM 2.1 (Meta, Apache-2.0)"),
    NL("SAM 2.1-licentie (Meta, Apache-2.0)"),
    RU("Лицензия SAM 2.1 (Meta, Apache-2.0)"),
    TR("SAM 2.1 Lisansı (Meta, Apache-2.0)"));

SS_MSG(license_sam2_summary,
    EN("SAM 2.1 is Meta's model, released under the Apache 2.0 licence. "
       "Nothing unusual to agree to; it is downloaded rather than bundled "
       "only to keep the app small."),
    JA("SAM 2.1 は Meta のモデルで、Apache 2.0 ライセンスで公開されています。"
       "特別に同意が必要なことはありません。同梱せずダウンロードにしているのは、"
       "アプリを小さく保つためだけです。"),
    ZH_HANS("SAM 2.1 是 Meta 的模型，以 Apache 2.0 许可协议发布。没有什么"
            "特别需要同意的；之所以下载而不是打包，只是为了让应用保持小巧。"),
    ZH_HANT("SAM 2.1 是 Meta 的模型，以 Apache 2.0 授權條款發布。沒有什麼"
            "特別需要同意的；之所以下載而不是打包，只是為了讓應用保持小巧。"),
    KO("SAM 2.1은 Meta의 모델이며 Apache 2.0 라이선스로 공개되어 있습니다. "
       "특별히 동의할 것은 없습니다. 함께 담지 않고 내려받게 한 것은 앱을 작게 "
       "유지하기 위해서일 뿐입니다."),
    DE("SAM 2.1 ist Metas Modell, veröffentlicht unter der Apache-2.0-Lizenz. "
       "Es ist nichts Ungewöhnliches zuzustimmen; heruntergeladen statt "
       "mitgeliefert wird es nur, damit die Anwendung klein bleibt."),
    FR("SAM 2.1 est le modèle de Meta, publié sous licence Apache 2.0. Rien "
       "d'inhabituel à accepter ; il est téléchargé plutôt que livré avec "
       "l'application uniquement pour la garder légère."),
    ES("SAM 2.1 es el modelo de Meta, publicado bajo la licencia Apache 2.0. "
       "No hay nada inusual que aceptar; se descarga en lugar de incluirse "
       "solo para que la aplicación siga siendo pequeña."),
    PT("O SAM 2.1 é o modelo da Meta, publicado sob a licença Apache 2.0. Não "
       "há nada de incomum a aceitar; ele é baixado em vez de embutido apenas "
       "para manter o aplicativo pequeno."),
    IT("SAM 2.1 è il modello di Meta, pubblicato con licenza Apache 2.0. Non "
       "c'è nulla di insolito da accettare; viene scaricato anziché incluso "
       "solo per tenere piccola l'applicazione."),
    NL("SAM 2.1 is het model van Meta, uitgebracht onder de Apache 2.0-"
       "licentie. Er valt niets ongebruikelijks te aanvaarden; het wordt "
       "gedownload in plaats van meegeleverd, alleen om de toepassing klein "
       "te houden."),
    RU("SAM 2.1 — модель Meta, выпущенная под лицензией Apache 2.0. Ничего "
       "необычного принимать не нужно; она загружается, а не поставляется в "
       "комплекте, лишь чтобы программа оставалась небольшой."),
    TR("SAM 2.1 Meta'nın modelidir ve Apache 2.0 lisansıyla yayımlanmıştır. "
       "Kabul edilecek olağandışı bir şey yok; uygulamanın küçük kalması için "
       "birlikte verilmek yerine indiriliyor."));

// ===========================================================================
// Model licence consent
//
// IRREVERSIBLE / LEGAL -- every message in this block gets human review in
// every language before shipping. A user is being asked to accept somebody
// else's terms; a translation that softens or overstates them is worse than
// no translation at all.
// ===========================================================================

SS_MSG(license_modal_title,
    EN("Model licence"), JA("モデルのライセンス"), ZH_HANS("模型许可协议"),
    ZH_HANT("模型授權條款"), KO("모델 라이선스"), DE("Modelllizenz"),
    FR("Licence du modèle"), ES("Licencia del modelo"),
    PT("Licença do modelo"), IT("Licenza del modello"), NL("Modellicentie"),
    RU("Лицензия модели"), TR("Model lisansı"));

SS_MSG(license_read,
    EN("Read the licence"),
    JA("ライセンスを読む"), ZH_HANS("阅读许可协议"), ZH_HANT("閱讀授權條款"),
    KO("라이선스 읽기"),   DE("Lizenz lesen"), FR("Lire la licence"),
    ES("Leer la licencia"), PT("Ler a licença"), IT("Leggi la licenza"),
    NL("Licentie lezen"),  RU("Прочитать лицензию"), TR("Lisansı oku"));

SS_MSG(license_copy_link,
    EN("Copy link"),     JA("リンクをコピー"), ZH_HANS("复制链接"),  ZH_HANT("複製連結"),
    KO("링크 복사"),      DE("Link kopieren"), FR("Copier le lien"),
    ES("Copiar el enlace"), PT("Copiar o link"), IT("Copia il link"),
    NL("Link kopiëren"), RU("Скопировать ссылку"), TR("Bağlantıyı kopyala"));

// {0} is a size like "707 MB".
SS_MSG(license_download_size,
    EN("Download: about {0}, kept for next time."),
    JA("ダウンロード: 約 {0}。次回以降は再利用されます。"),
    ZH_HANS("下载：约 {0}，之后会保留下来。"),
    ZH_HANT("下載：約 {0}，之後會保留下來。"),
    KO("내려받기: 약 {0}. 다음부터는 그대로 씁니다."),
    DE("Download: etwa {0}, bleibt für das nächste Mal erhalten."),
    FR("Téléchargement : environ {0}, conservé pour la prochaine fois."),
    ES("Descarga: unos {0}, se conserva para la próxima vez."),
    PT("Download: cerca de {0}, guardado para a próxima vez."),
    IT("Scaricamento: circa {0}, resta per la prossima volta."),
    NL("Download: ongeveer {0}, blijft bewaard voor de volgende keer."),
    RU("Загрузка: около {0}, сохраняется на будущее."),
    TR("İndirme: yaklaşık {0}, bir dahaki sefere saklanır."));

SS_MSG(license_accept_tick,
    EN("I have read and accept these terms"),
    JA("これらの条件を読み、同意します"),
    ZH_HANS("我已阅读并接受这些条款"),
    ZH_HANT("我已閱讀並接受這些條款"),
    KO("이 조건을 읽었고 이에 동의합니다"),
    DE("Ich habe diese Bedingungen gelesen und nehme sie an"),
    FR("J'ai lu et j'accepte ces conditions"),
    ES("He leído y acepto estos términos"),
    PT("Li e aceito estes termos"),
    IT("Ho letto e accetto queste condizioni"),
    NL("Ik heb deze voorwaarden gelezen en aanvaard ze"),
    RU("Я прочитал эти условия и принимаю их"),
    TR("Bu koşulları okudum ve kabul ediyorum"));

SS_MSG(license_download,
    EN("Download"),      JA("ダウンロード"),  ZH_HANS("下载"),     ZH_HANT("下載"),
    KO("내려받기"),       DE("Herunterladen"), FR("Télécharger"), ES("Descargar"),
    PT("Baixar"),        IT("Scarica"),      NL("Downloaden"),   RU("Загрузить"),
    TR("İndir"));

// {0} is the licence URL.
SS_MSG(license_no_browser,
    EN("Could not open a browser. The licence is at {0} (copied to the "
       "clipboard)."),
    JA("ブラウザを開けませんでした。ライセンスは {0} にあります"
       "（クリップボードにコピーしました）。"),
    ZH_HANS("无法打开浏览器。许可协议在 {0}（已复制到剪贴板）。"),
    ZH_HANT("無法開啟瀏覽器。授權條款在 {0}（已複製到剪貼簿）。"),
    KO("브라우저를 열지 못했습니다. 라이선스는 {0}에 있습니다(클립보드에 "
       "복사했습니다)."),
    DE("Es ließ sich kein Browser öffnen. Die Lizenz steht unter {0} (in die "
       "Zwischenablage kopiert)."),
    FR("Impossible d'ouvrir un navigateur. La licence est à l'adresse {0} "
       "(copiée dans le presse-papiers)."),
    ES("No se pudo abrir un navegador. La licencia está en {0} (copiada al "
       "portapapeles)."),
    PT("Não foi possível abrir um navegador. A licença está em {0} (copiada "
       "para a área de transferência)."),
    IT("Non è stato possibile aprire un browser. La licenza si trova in {0} "
       "(copiata negli appunti)."),
    NL("Er kon geen browser worden geopend. De licentie staat op {0} "
       "(gekopieerd naar het klembord)."),
    RU("Не удалось открыть браузер. Лицензия находится по адресу {0} "
       "(скопирован в буфер обмена)."),
    TR("Bir tarayıcı açılamadı. Lisans şu adreste: {0} (panoya kopyalandı)."));

// ===========================================================================
// Log lines this screen writes
// ===========================================================================

SS_MSG(log_masks_attached,
    EN("Using {0} as the masks for the images beside it."),
    JA("{0} を、その隣にある画像のマスクとして使います。"),
    ZH_HANS("把 {0} 用作它旁边那些图像的蒙版。"),
    ZH_HANT("把 {0} 用作它旁邊那些影像的遮罩。"),
    KO("{0}을(를) 그 옆 이미지들의 마스크로 사용합니다."),
    DE("{0} wird als Maskenordner für die daneben liegenden Bilder benutzt."),
    FR("{0} est utilisé comme masques pour les images qui se trouvent à côté."),
    ES("Se usa {0} como máscaras de las imágenes que están junto a ella."),
    PT("Usando {0} como as máscaras das imagens ao lado."),
    IT("Si usa {0} come maschere per le immagini che le stanno accanto."),
    NL("{0} wordt gebruikt als maskers voor de beelden ernaast."),
    RU("{0} используется как маски для соседних снимков."),
    TR("{0}, yanındaki görüntülerin maskeleri olarak kullanılıyor."));

SS_MSG(log_masks_orphaned,
    EN("Ignored {0}: that is a folder of masks, and the images they belong to "
       "were not picked. Add the images folder -- its masks are found on their "
       "own."),
    JA("{0} は無視しました。マスクのフォルダですが、対応する画像が選ばれて"
       "いません。画像のフォルダを追加してください。マスクは自動で見つかります。"),
    ZH_HANS("已忽略 {0}：那是一个蒙版文件夹，但它对应的图像没有被选中。"
            "请添加图像文件夹——它的蒙版会被自动找到。"),
    ZH_HANT("已忽略 {0}：那是一個遮罩資料夾，但它對應的影像沒有被選取。"
            "請新增影像資料夾——它的遮罩會被自動找到。"),
    KO("{0}은(는) 무시했습니다. 마스크 폴더인데 그것이 딸린 이미지가 선택되지 "
       "않았습니다. 이미지 폴더를 추가하세요. 마스크는 알아서 찾습니다."),
    DE("{0} wurde übergangen: das ist ein Maskenordner, und die zugehörigen "
       "Bilder wurden nicht gewählt. Fügen Sie den Bilderordner hinzu -- seine "
       "Masken werden von selbst gefunden."),
    FR("{0} a été ignoré : c'est un dossier de masques, et les images "
       "auxquelles ils appartiennent n'ont pas été choisies. Ajoutez le "
       "dossier d'images -- ses masques sont trouvés tout seuls."),
    ES("Se ignoró {0}: es una carpeta de máscaras y no se eligieron las "
       "imágenes a las que pertenecen. Añada la carpeta de imágenes: sus "
       "máscaras se encuentran solas."),
    PT("{0} foi ignorada: é uma pasta de máscaras, e as imagens a que "
       "pertencem não foram escolhidas. Adicione a pasta de imagens -- as "
       "máscaras dela são encontradas sozinhas."),
    IT("{0} è stata ignorata: è una cartella di maschere e le immagini a cui "
       "appartengono non sono state scelte. Aggiunga la cartella delle "
       "immagini: le sue maschere vengono trovate da sole."),
    NL("{0} is genegeerd: dat is een map met maskers, en de beelden waar ze "
       "bij horen zijn niet gekozen. Voeg de beeldmap toe -- de bijbehorende "
       "maskers worden vanzelf gevonden."),
    RU("{0} пропущена: это папка масок, а снимки, к которым они относятся, не "
       "выбраны. Добавьте папку со снимками — их маски находятся сами."),
    TR("{0} yok sayıldı: burası bir maske klasörü ve ait oldukları görüntüler "
       "seçilmedi. Görüntü klasörünü ekleyin -- maskeleri kendiliğinden "
       "bulunur."));

SS_MSG(log_clicks_dropped_input_gone,
    EN("Clicked objects dropped: {0}, the input they were drawn on, is no "
       "longer in the list."),
    JA("クリックした物体を破棄しました。それらを指定した入力 {0} が"
       "リストから外れたためです。"),
    ZH_HANS("已丢弃点选的物体：它们所依附的输入 {0} 已经不在列表里了。"),
    ZH_HANT("已丟棄點選的物體：它們所依附的輸入 {0} 已經不在清單裡了。"),
    KO("클릭한 물체를 버렸습니다. 그것들을 지정한 입력 {0}이(가) 목록에서 "
       "빠졌습니다."),
    DE("Angeklickte Objekte verworfen: {0}, die Eingabe, auf der sie "
       "eingezeichnet wurden, steht nicht mehr auf der Liste."),
    FR("Objets cliqués abandonnés : {0}, l'entrée sur laquelle ils avaient été "
       "tracés, ne figure plus dans la liste."),
    ES("Se descartaron los objetos marcados: {0}, la entrada sobre la que se "
       "marcaron, ya no está en la lista."),
    PT("Objetos clicados descartados: {0}, a entrada em que foram marcados, "
       "não está mais na lista."),
    IT("Oggetti cliccati scartati: {0}, l'ingresso su cui erano stati "
       "tracciati, non è più nell'elenco."),
    NL("Aangeklikte objecten vervallen: {0}, de invoer waarop ze waren gezet, "
       "staat niet meer in de lijst."),
    RU("Отмеченные объекты сброшены: вход {0}, на котором они были указаны, "
       "больше не в списке."),
    TR("Tıklanan nesneler bırakıldı: üzerlerinde işaretlendikleri {0} girdisi "
       "artık listede değil."));

SS_MSG(log_clicks_dropped_switched,
    EN("Clicked objects dropped: they were drawn on {0}, which is not the "
       "input being previewed now."),
    JA("クリックした物体を破棄しました。{0} で指定されたもので、いま"
       "プレビューしている入力とは異なるためです。"),
    ZH_HANS("已丢弃点选的物体：它们是在 {0} 上标注的，而现在预览的不是那个输入。"),
    ZH_HANT("已丟棄點選的物體：它們是在 {0} 上標註的，而現在預覽的不是那個輸入。"),
    KO("클릭한 물체를 버렸습니다. {0}에서 지정한 것인데, 지금 미리보기 중인 "
       "입력이 아닙니다."),
    DE("Angeklickte Objekte verworfen: sie wurden auf {0} eingezeichnet, und "
       "das ist nicht die jetzt betrachtete Eingabe."),
    FR("Objets cliqués abandonnés : ils avaient été tracés sur {0}, qui n'est "
       "pas l'entrée actuellement prévisualisée."),
    ES("Se descartaron los objetos marcados: se marcaron sobre {0}, que no es "
       "la entrada que se está previsualizando ahora."),
    PT("Objetos clicados descartados: foram marcados em {0}, que não é a "
       "entrada em prévia agora."),
    IT("Oggetti cliccati scartati: erano stati tracciati su {0}, che non è "
       "l'ingresso ora in anteprima."),
    NL("Aangeklikte objecten vervallen: ze waren gezet op {0}, en dat is niet "
       "de invoer die nu wordt bekeken."),
    RU("Отмеченные объекты сброшены: они указаны на входе {0}, а сейчас "
       "просматривается другой."),
    TR("Tıklanan nesneler bırakıldı: {0} üzerinde işaretlenmişlerdi, oysa şu "
       "an önizlenen girdi o değil."));

SS_MSG(log_drop_no_images,
    EN("Dropped folder contains no dataset or images: {0}"),
    JA("ドロップされたフォルダにデータセットも画像もありません: {0}"),
    ZH_HANS("拖入的文件夹里既没有数据集也没有图像：{0}"),
    ZH_HANT("拖入的資料夾裡既沒有資料集也沒有影像：{0}"),
    KO("끌어다 놓은 폴더에 데이터셋도 이미지도 없습니다: {0}"),
    DE("Der abgelegte Ordner enthält weder Datensatz noch Bilder: {0}"),
    FR("Le dossier déposé ne contient ni jeu de données ni images : {0}"),
    ES("La carpeta arrastrada no contiene ni conjunto de datos ni imágenes: {0}"),
    PT("A pasta arrastada não contém conjunto de dados nem imagens: {0}"),
    IT("La cartella trascinata non contiene né set di dati né immagini: {0}"),
    NL("De neergezette map bevat geen dataset of beelden: {0}"),
    RU("В перетащенной папке нет ни набора данных, ни изображений: {0}"),
    TR("Bırakılan klasörde ne veri kümesi ne de görüntü var: {0}"));

SS_MSG(log_drop_unsupported,
    EN("Unsupported dropped file: {0}"),
    JA("対応していないファイルがドロップされました: {0}"),
    ZH_HANS("拖入了不支持的文件：{0}"),
    ZH_HANT("拖入了不支援的檔案：{0}"),
    KO("지원하지 않는 파일을 끌어다 놓았습니다: {0}"),
    DE("Abgelegte Datei wird nicht unterstützt: {0}"),
    FR("Fichier déposé non pris en charge : {0}"),
    ES("Archivo arrastrado no compatible: {0}"),
    PT("Arquivo arrastado sem suporte: {0}"),
    IT("File trascinato non supportato: {0}"),
    NL("Niet-ondersteund bestand neergezet: {0}"),
    RU("Перетащенный файл не поддерживается: {0}"),
    TR("Desteklenmeyen dosya bırakıldı: {0}"));

SS_MSG(log_drop_while_training,
    EN("Dropped input ignored: stop training first"),
    JA("ドロップされた入力を無視しました。先に学習を停止してください"),
    ZH_HANS("已忽略拖入的输入：请先停止训练"),
    ZH_HANT("已忽略拖入的輸入：請先停止訓練"),
    KO("끌어다 놓은 입력을 무시했습니다. 먼저 학습을 멈추세요"),
    DE("Abgelegte Eingabe übergangen: erst das Training anhalten"),
    FR("Entrée déposée ignorée : arrêtez d'abord l'entraînement"),
    ES("Entrada arrastrada ignorada: detenga primero el entrenamiento"),
    PT("Entrada arrastada ignorada: pare o treinamento primeiro"),
    IT("Ingresso trascinato ignorato: fermi prima l'addestramento"),
    NL("Neergezette invoer genegeerd: stop eerst de training"),
    RU("Перетащенный вход пропущен: сначала остановите обучение"),
    TR("Bırakılan girdi yok sayıldı: önce eğitimi durdurun"));

SS_MSG(log_dataset_settings_changed,
    EN("Dataset settings changed; reloading dataset"),
    JA("データセットの設定が変わったため、読み込み直します"),
    ZH_HANS("数据集设置已更改，正在重新加载数据集"),
    ZH_HANT("資料集設定已變更，正在重新載入資料集"),
    KO("데이터셋 설정이 바뀌어 데이터셋을 다시 불러옵니다"),
    DE("Datensatzeinstellungen geändert; Datensatz wird neu geladen"),
    FR("Réglages du jeu de données modifiés ; rechargement en cours"),
    ES("Cambiaron los ajustes del conjunto de datos; recargándolo"),
    PT("As configurações do conjunto de dados mudaram; recarregando"),
    IT("Impostazioni del set di dati cambiate; ricaricamento in corso"),
    NL("Datasetinstellingen gewijzigd; dataset wordt opnieuw geladen"),
    RU("Параметры набора данных изменились; перезагрузка"),
    TR("Veri kümesi ayarları değişti; veri kümesi yeniden yükleniyor"));

}  // namespace dataset
}  // namespace msg
}  // namespace i18n
}  // namespace spirula

#include "i18n/EndCatalog.h"
