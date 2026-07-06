export camera_model=OPENCV  # or OPENCV_FISHEYE, THIN_PRISM_FISHEYE
export max_num_features=8192  # O(n^2) time complexity, less than 1024 can work
export image_path="images"  # or images_2, images_4
#export vocab_tree_path="../vocab_tree_flickr100K_words32K.bin"
export vocab_tree_path="../vocab_tree_faiss_flickr100K_words256K.bin"  # download from https://github.com/colmap/colmap/releases/download/3.11.1/vocab_tree_faiss_flickr100K_words256K.bin

# feature extraction; options to try:
# - replace `--ImageReader.single_camera_per_folder 1` with `--ImageReader.single_camera 0`, `--ImageReader.single_camera 1`, etc.
# - add `--ImageReader.mask_path ./masks`
# - specify `--ImageReader.camera_params "FX,FY,CX,CY,..."`, and optionally add `--Mapper.ba_refine_focal_length 0 --Mapper.ba_refine_extra_params 0` to mapper
# - `--SiftExtraction.estimate_affine_shape 1`, and use `--FeatureMatching.guided_matching 1` in feature matching
# - limit `--FeatureExtraction.max_image_size`
colmap feature_extractor --database_path database.db --image_path ./$image_path --ImageReader.single_camera_per_folder 1 --ImageReader.camera_model $camera_model --SiftExtraction.max_num_features $max_num_features
#colmap feature_extractor --database_path database.db --image_path ./$image_path --ImageReader.single_camera_per_folder 1 --ImageReader.camera_model $camera_model --FeatureExtraction.type ALIKED_N16ROT --AlikedExtraction.max_num_features $max_num_features --FeatureExtraction.max_image_size 750
#colmap feature_extractor --database_path database.db --image_path ./$image_path --ImageReader.single_camera_per_folder 1 --ImageReader.camera_model $camera_model --FeatureExtraction.type ALIKED_N32 --AlikedExtraction.max_num_features $max_num_features --FeatureExtraction.max_image_size 750 --AlikedExtraction.min_score 0.1

# pairwise feature matching
colmap vocab_tree_matcher --database_path database.db --VocabTreeMatching.vocab_tree_path $vocab_tree_path
#colmap sequential_matcher --database_path database.db --SequentialMatching.overlap 10 --SequentialMatching.quadratic_overlap 0
#colmap exhaustive_matcher --database_path database.db
#colmap exhaustive_matcher --database_path database.db --FeatureMatching.type ALIKED_BRUTEFORCE

mkdir -p sparse

# solve bundle adjustment
# optionally add `--Mapper.ba_global_backend CASPAR` for speed
colmap mapper --database_path database.db --image_path ./$image_path --output_path sparse --Mapper.ba_use_gpu 1
#colmap mapper --database_path database.db --image_path ./$image_path --output_path sparse --Mapper.ba_use_gpu 1 --Mapper.ba_refine_extra_params 0 --Mapper.init_min_tri_angle 4
#colmap global_mapper --database_path database.db --image_path ./$image_path --output_path sparse --GlobalMapper.ba_ceres_use_gpu 1 --GlobalMapper.ba_refine_extra_params 0

colmap bundle_adjuster --input_path sparse/0 --output_path sparse/0 --BundleAdjustment.refine_focal_length 1 --BundleAdjustment.refine_principal_point 1 --BundleAdjustment.refine_extra_params 1
