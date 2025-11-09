export camera_model=OPENCV  # or OPENCV_FISHEYE
export max_num_features=8192  # O(n^2) time complexity, less than 1024 can work
export image_path="images"  # or images_2, images_4
export vocab_tree_path="../vocab_tree_flickr100K_words32K.bin"

# feature extraction
# replace `--ImageReader.single_camera 1` with `--ImageReader.single_camera 0`, `--ImageReader.single_camera_per_folder 1`, etc.
colmap feature_extractor --database_path database.db --image_path ./$image_path --ImageReader.single_camera 1 --ImageReader.camera_model $camera_model --SiftExtraction.max_num_features $max_num_features

# pairwise feature matching
if true; then
    colmap vocab_tree_matcher --database_path database.db --VocabTreeMatching.vocab_tree_path $vocab_tree_path
elif false; then
    colmap sequential_matcher --database_path database.db --SequentialMatching.overlap 10 --SequentialMatching.quadratic_overlap 0
else
    colmap matches_importer --database_path database.db --match_list_path ./matches.txt --match_type pairs
fi

# solve bundle adjustment
mkdir sparse
if false; then
    glomap mapper --database_path database.db --output_path sparse
else
    colmap mapper --database_path database.db --image_path ./$image_path --output_path sparse
    colmap bundle_adjuster --input_path sparse/0 --output_path sparse/0 --BundleAdjustment.refine_principal_point 1
fi

# export poses and point cloud using nerfstudio
ns-process-data images --data ./$image_path --output-dir . --skip-image-processing --skip-colmap --colmap-model-path sparse/0/

if [[ "$image_path" != "images" ]]; then
cp transforms.json transforms_images.json
sed -E "s/\"images\//\"$image_path\//g" transforms_images.json > transforms.json
fi

# dense reconstruction
if false; then
mkdir dense
colmap image_undistorter --image_path $image_path --input_path ./sparse/0 --output_path dense --output_type COLMAP --max_image_size 512
if false; then
    colmap patch_match_stereo --workspace_path ./dense --workspace_format COLMAP --PatchMatchStereo.geom_consistency true
    colmap stereo_fusion --workspace_path ./dense --workspace_format COLMAP --input_type geometric --output_path ./dense/fused.ply
else
    colmap patch_match_stereo --workspace_path ./dense --workspace_format COLMAP --PatchMatchStereo.geom_consistency false
    colmap stereo_fusion --workspace_path ./dense --workspace_format COLMAP --input_type photometric --output_path ./dense/fused.ply
fi
#colmap poisson_mesher --input_path ./dense/fused.ply --output_path ./dense/meshed-poisson.ply --PoissonMeshing.depth 8
colmap delaunay_mesher --input_path ./dense --output_path ./dense/meshed-delaunay.ply --DelaunayMeshing.visibility_sigma 3 --DelaunayMeshing.quality_regularization 5 --DelaunayMeshing.max_side_length_factor inf
fi
