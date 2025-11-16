"""
Nerfstudio Template Config

Define your custom method here that registers with Nerfstudio CLI.
"""

from __future__ import annotations

from spirulae_splat.ns_datamanager import (
    SpirulaeDataManagerConfig, SpirulaeDataManager
)
from spirulae_splat.ns_model import SpirulaeModelConfig
from spirulae_splat.ns_pipeline import (
    SpirulaePipelineConfig, VanillaPipelineConfig
)
from spirulae_splat.ns_dataset import SpirulaeDataset
from spirulae_splat.ns_dataparser import Nerfstudio2DataParserConfig
from nerfstudio.configs.base_config import ViewerConfig
from nerfstudio.engine.optimizers import AdamOptimizerConfig, RAdamOptimizerConfig
from nerfstudio.engine.schedulers import (
    ExponentialDecaySchedulerConfig,
)
from nerfstudio.engine.trainer import TrainerConfig
from nerfstudio.plugins.types import MethodSpecification


_DEFAULT_DATAMANAGER_CONFIG = {
    '_target': SpirulaeDataManager,
    'dataparser': Nerfstudio2DataParserConfig(
        eval_mode="fraction", train_split_fraction=1.0,  # use all for training
        # eval_mode="interval", eval_interval=8,  # consistent with academic benchmarks
    ),
    'camera_res_scale_factor': 1.0,
    'eval_num_images_to_sample_from': -1,
    'eval_num_times_to_repeat_images': -1,
    'max_thread_workers': None
}

_DEFAULT_OPTIMIZERS = {
    "means": {
        "optimizer": AdamOptimizerConfig(lr=1.0e-4, eps=1e-15),
        "scheduler": ExponentialDecaySchedulerConfig(
            lr_final=1.0e-6, max_steps=30000,
        ),
    },
    "scales": {
        "optimizer": AdamOptimizerConfig(lr=0.005, eps=1e-15),
        "scheduler": None,
    },
    "quats": {
        "optimizer": AdamOptimizerConfig(lr=0.0005, eps=1e-15),
        "scheduler": None,
    },
    "features_dc": {
        "optimizer": AdamOptimizerConfig(lr=0.0025, eps=1e-15),
        "scheduler": None,
    },
    "features_sh": {
        "optimizer": AdamOptimizerConfig(lr=0.0025 / 20, eps=1e-15),
        "scheduler": None,
    },
    "features_ch": {
        "optimizer": AdamOptimizerConfig(lr=0.0025 / 5, eps=1e-15),
        "scheduler": None,
    },
    "opacities": {
        "optimizer": AdamOptimizerConfig(lr=0.05, eps=1e-15),
        "scheduler": None,
    },
    "background_color": {
        "optimizer": AdamOptimizerConfig(lr=0.0025, eps=1e-15),
        "scheduler": None
    },
    "background_sh": {
        "optimizer": AdamOptimizerConfig(lr=0.0025 / 5, eps=1e-15),
        "scheduler": None
    },
    "bilateral_grid": {
        "optimizer": AdamOptimizerConfig(lr=2e-3, eps=1e-15),
        "scheduler": ExponentialDecaySchedulerConfig(
            lr_final=1e-4, max_steps=30000, warmup_steps=1000, lr_pre_warmup=0
        ),
    },
    "bilateral_grid_geometry": {
        # "optimizer": AdamOptimizerConfig(lr=1e-2, eps=1e-15),
        "optimizer": AdamOptimizerConfig(lr=2e-3, eps=1e-15),
        "scheduler": ExponentialDecaySchedulerConfig(
            lr_final=1e-4, max_steps=30000, warmup_steps=1000, lr_pre_warmup=0
        ),
    },
    "camera_opt": {
        "optimizer": AdamOptimizerConfig(lr=1e-4, eps=1e-15),  # 1e-4
        "scheduler": ExponentialDecaySchedulerConfig(
            lr_final=5e-7, max_steps=30000, warmup_steps=1000, lr_pre_warmup=0
        ),
    },
}

_TRIANGLE_OPTIMIZERS = {**_DEFAULT_OPTIMIZERS}
_TRIANGLE_OPTIMIZERS["means"] = {
    "optimizer": AdamOptimizerConfig(lr=1.0e-4, eps=1e-15),
    "scheduler": ExponentialDecaySchedulerConfig(
        lr_final=1.0e-6, max_steps=30000,
    ),
}
# _TRIANGLE_OPTIMIZERS["scales"] = {
#     "optimizer": AdamOptimizerConfig(lr=0.005, eps=1e-15),
#     "scheduler": ExponentialDecaySchedulerConfig(
#         lr_final=0.0002, max_steps=30000,
#     ),
# }
# _TRIANGLE_OPTIMIZERS["quats"] = {
#     "optimizer": AdamOptimizerConfig(lr=0.0005, eps=1e-15),
#     "scheduler": ExponentialDecaySchedulerConfig(
#         lr_final=0.0001, max_steps=30000
#     ),
# }
_TRIANGLE_OPTIMIZERS["bilateral_grid"] = {
    "optimizer": AdamOptimizerConfig(lr=5e-4, eps=1e-15),
    "scheduler": ExponentialDecaySchedulerConfig(
        lr_final=1e-6, max_steps=30000, warmup_steps=1000, lr_pre_warmup=0
    ),
}


spirulae = MethodSpecification(
    config=TrainerConfig(
        method_name="spirulae",
        steps_per_eval_batch=0,
        steps_per_save=2000,
        max_num_iterations=30000,
        mixed_precision=False,
        pipeline=SpirulaePipelineConfig(
            datamanager=SpirulaeDataManagerConfig(
                **_DEFAULT_DATAMANAGER_CONFIG
            ),
            model=SpirulaeModelConfig(
                background_color="black",
                train_background_color=False,
            ),
            
        ),
        optimizers={
            **_DEFAULT_OPTIMIZERS
        },
        viewer=ViewerConfig(),
        vis="viewer",
    ),
    description="Spirulae 3DGS Default.",
)

spirulae_patched = MethodSpecification(
    config=TrainerConfig(
        method_name="spirulae-patched",
        steps_per_eval_batch=0,
        steps_per_save=2000,
        max_num_iterations=30000,
        mixed_precision=False,
        pipeline=SpirulaePipelineConfig(
            datamanager=SpirulaeDataManagerConfig(
                **_DEFAULT_DATAMANAGER_CONFIG,
                patch_batch_size=-1,
                patch_size=64,
                max_batch_per_epoch=512,
            ),
            model=SpirulaeModelConfig(
                packed=True,
                use_bvh=True,
                use_camera_optimizer=False,
                use_bilateral_grid=False,
                use_bilateral_grid_for_geometry=False,
                alpha_reg_weight=0.0,
                mcmc_max_screen_size=0.15,
                mcmc_max_screen_size_clip_hardness=1.5,
                # mcmc_max_world_size=1.0,
                background_color="black",
                train_background_color=False,
            ),
        ),
        optimizers={
            **_DEFAULT_OPTIMIZERS
        },
        viewer=ViewerConfig(),
        vis="viewer",
    ),
    description="Spirulae 3DGS with patch batching.",
)

spirulae_triangle = MethodSpecification(
    config=TrainerConfig(
        method_name="spirulae-triangle",
        steps_per_eval_batch=0,
        steps_per_save=2000,
        max_num_iterations=30000,
        mixed_precision=False,
        pipeline=SpirulaePipelineConfig(
            datamanager=SpirulaeDataManagerConfig(
                **_DEFAULT_DATAMANAGER_CONFIG
            ),
            model=SpirulaeModelConfig(
                primitive="opaque_triangle",
                kernel_radius=0.5,
                compute_depth_normal=True,
                sh_degree=0,
                bilagrid_shape=(16, 16, 8),
                stop_refine_at=30000,
                background_color="black",
                train_background_color=False,
                # alpha_reg_weight=0.0,
                mcmc_scale_reg=0.02,
                # erank_reg=1.0,
                # supersampling=2,
                mcmc_min_opacity=0.005,
                mcmc_noise_lr=5e5,  # or 0.0
                mcmc_max_screen_size=0.15,
                supervision_warmup=0,
                depth_supervision_weight=0.25,
                normal_supervision_weight=0.04,
                rgb_distortion_reg_weight=0.01,
                ssim_lambda=0.4,
            ),
            
        ),
        optimizers={
            **_TRIANGLE_OPTIMIZERS
        },
        viewer=ViewerConfig(),
        vis="viewer",
    ),
    description="Spirulae opaque triangle splatting.",
)

spirulae_triangle_patched = MethodSpecification(
    config=TrainerConfig(
        method_name="spirulae-triangle-patched",
        steps_per_eval_batch=0,
        steps_per_save=2000,
        max_num_iterations=30000,
        mixed_precision=False,
        pipeline=SpirulaePipelineConfig(
            datamanager=SpirulaeDataManagerConfig(
                **_DEFAULT_DATAMANAGER_CONFIG,
                patch_batch_size=-1,
                patch_size=64,
                max_batch_per_epoch=512,
            ),
            model=SpirulaeModelConfig(
                primitive="opaque_triangle",
                kernel_radius=0.5,
                compute_depth_normal=True,
                sh_degree=0,
                bilagrid_shape=(16, 16, 8),
                stop_refine_at=30000,
                background_color="black",
                train_background_color=False,
                # alpha_reg_weight=0.0,
                mcmc_scale_reg=0.02,
                # erank_reg=1.0,
                # supersampling=2,
                mcmc_min_opacity=0.005,
                mcmc_noise_lr=5e5,  # or 0.0
                supervision_warmup=0,
                depth_supervision_weight=0.25,
                normal_supervision_weight=0.04,
                rgb_distortion_reg_weight=0.01,
                ssim_lambda=0.4,
                reg_warmup_length=1000,

                packed=True,
                use_bvh=True,
                use_camera_optimizer=False,
                use_bilateral_grid=False,
                use_bilateral_grid_for_geometry=False,
                alpha_reg_weight=0.0,
                mcmc_max_screen_size=0.15,
                mcmc_max_screen_size_clip_hardness=1.5,
            ),
            
        ),
        optimizers={
            **_TRIANGLE_OPTIMIZERS
        },
        viewer=ViewerConfig(),
        vis="viewer",
    ),
    description="Spirulae opaque triangle splatting with patched batching.",
)
