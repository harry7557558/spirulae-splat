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
from spirulae_splat.modules.optimizer import FusedAdamOptimizerConfig, FusedNewtonOptimizerConfig

from spirulae_splat.ns_dataset import SpirulaeDataset
from spirulae_splat.ns_dataparser import Nerfstudio2DataParserConfig
from nerfstudio.configs.base_config import ViewerConfig
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
    "_dummy": {
        "optimizer": FusedAdamOptimizerConfig(lr=1.0, eps=0.0),
        "scheduler": None,
    },
    "means": {
        "optimizer": FusedAdamOptimizerConfig(lr=1.6e-4, eps=1e-15),
        "scheduler": ExponentialDecaySchedulerConfig(
            lr_final=1.6e-6, max_steps=30000,
        ),
    },
    "scales": {
        "optimizer": FusedAdamOptimizerConfig(lr=0.005, eps=1e-15),
        "scheduler": None,
    },
    "quats": {
        "optimizer": FusedAdamOptimizerConfig(lr=0.0005, eps=1e-15),
        "scheduler": None,
    },
    "features_dc": {
        "optimizer": FusedAdamOptimizerConfig(lr=0.0025, eps=1e-15),
        "scheduler": None,
    },
    "features_sh": {
        "optimizer": FusedAdamOptimizerConfig(
            lr=0.0025 / 20, eps=1e-15,
            tr=1.0e-6 / 20, tr_final=1.0e-8 / 20, tr_max_steps=30000
        ),
        "scheduler": None,
    },
    "features_ch": {
        "optimizer": FusedAdamOptimizerConfig(lr=0.0025 / 5, eps=1e-15),
        "scheduler": None,
    },
    "sv_sites": {
        "optimizer": FusedAdamOptimizerConfig(lr=0.01, eps=1e-15),
        "scheduler": None,
    },
    "sv_colors": {
        "optimizer": FusedAdamOptimizerConfig(lr=0.0005, eps=1e-15),
        "scheduler": None,
    },
    "opacities": {
        "optimizer": FusedAdamOptimizerConfig(lr=0.05, eps=1e-15),
        "scheduler": None,
    },
    "densities": {
        "optimizer": FusedAdamOptimizerConfig(lr=0.05, eps=1e-15),
        "scheduler": ExponentialDecaySchedulerConfig(
            lr_final=0.0005, max_steps=30000,
        ),
    },
    "background_color": {
        "optimizer": FusedAdamOptimizerConfig(lr=0.0025, eps=1e-15),
        "scheduler": None
    },
    "background_sh": {
        "optimizer": FusedAdamOptimizerConfig(lr=0.0025 / 5, eps=1e-15),
        "scheduler": None
    },
    "bilateral_grid": {
        "optimizer": FusedAdamOptimizerConfig(lr=2e-3, eps=1e-15),
        "scheduler": ExponentialDecaySchedulerConfig(
            lr_final=1e-4, max_steps=30000, warmup_steps=1000, lr_pre_warmup=0
        ),
    },
    "bilateral_grid_depth": {
        # "optimizer": FusedAdamOptimizerConfig(lr=1e-2, eps=1e-15),
        "optimizer": FusedAdamOptimizerConfig(lr=2e-3, eps=1e-15),
        "scheduler": ExponentialDecaySchedulerConfig(
            lr_final=1e-4, max_steps=30000, warmup_steps=2000, lr_pre_warmup=0, ramp="linear"
        ),
    },
    "bilateral_grid_normal": {
        "optimizer": FusedAdamOptimizerConfig(lr=5e-4, eps=1e-15),
        "scheduler": ExponentialDecaySchedulerConfig(
            lr_final=4e-5, max_steps=30000, warmup_steps=2000, lr_pre_warmup=0, ramp="linear"
        ),
    },
    "ppisp": {
        "optimizer": FusedAdamOptimizerConfig(lr=2e-3, eps=1e-15),
        "scheduler": ExponentialDecaySchedulerConfig(
            lr_final=2e-5, max_steps=30000, warmup_steps=500, lr_pre_warmup=2e-5
        ),
    },
    "camera_opt": {
        "optimizer": FusedAdamOptimizerConfig(lr=1e-4, eps=1e-15),  # 1e-4
        "scheduler": ExponentialDecaySchedulerConfig(
            lr_final=5e-7, max_steps=30000, warmup_steps=1000, lr_pre_warmup=0
        ),
    },
}

_DEFAULT_OPTIMIZERS_WITH_SCALE_SCHEDULER = {**_DEFAULT_OPTIMIZERS}
_DEFAULT_OPTIMIZERS_WITH_SCALE_SCHEDULER["scales"] = {
    # https://arxiv.org/abs/2603.08661
    "optimizer": FusedAdamOptimizerConfig(lr=0.02, eps=1e-15),
    # "optimizer": FusedAdamOptimizerConfig(lr=0.005, eps=1e-15),
    "scheduler": ExponentialDecaySchedulerConfig(
        lr_final=0.005, max_steps=10000, warmup_steps=1000, ramp="linear", lr_pre_warmup=0.005
    ),
}

_TRIANGLE_OPTIMIZERS = {**_DEFAULT_OPTIMIZERS_WITH_SCALE_SCHEDULER}
_TRIANGLE_OPTIMIZERS["means"] = {
    "optimizer": FusedAdamOptimizerConfig(lr=1.6e-4, eps=1e-15),
    "scheduler": ExponentialDecaySchedulerConfig(
        lr_final=1.6e-6, max_steps=30000,
    ),
}
# _TRIANGLE_OPTIMIZERS["scales"] = {
#     "optimizer": FusedAdamOptimizerConfig(lr=0.005, eps=1e-15),
#     "scheduler": ExponentialDecaySchedulerConfig(
#         lr_final=0.0002, max_steps=30000,
#     ),
# }
# _TRIANGLE_OPTIMIZERS["quats"] = {
#     "optimizer": FusedAdamOptimizerConfig(lr=0.0005, eps=1e-15),
#     "scheduler": ExponentialDecaySchedulerConfig(
#         lr_final=0.0001, max_steps=30000
#     ),
# }
_TRIANGLE_OPTIMIZERS["bilateral_grid"] = {
    "optimizer": FusedAdamOptimizerConfig(lr=5e-4, eps=1e-15),
    "scheduler": ExponentialDecaySchedulerConfig(
        lr_final=1e-6, max_steps=30000, warmup_steps=1000, lr_pre_warmup=0
    ),
}

_SECOND_ORDER_POSITION_OPTIMIZERS = {**_DEFAULT_OPTIMIZERS_WITH_SCALE_SCHEDULER}
_SECOND_ORDER_POSITION_OPTIMIZERS["means"] = {
    "optimizer": FusedNewtonOptimizerConfig(mode="mean", lr=1.0e-6, eps=1e-15),
    "scheduler": ExponentialDecaySchedulerConfig(
        lr_final=1.0e-8, max_steps=30000, #warmup_steps=1000, lr_pre_warmup=0
    ),
}

_SECOND_ORDER_OPTIMIZERS = {**_SECOND_ORDER_POSITION_OPTIMIZERS}
# _SECOND_ORDER_OPTIMIZERS["quats"] = {
#     "optimizer": FusedNewtonOptimizerConfig(mode="quat", lr=1.0e-6, eps=1e-15),
#     "scheduler": ExponentialDecaySchedulerConfig(
#         lr_final=1.0e-8, max_steps=30000, #warmup_steps=1000, lr_pre_warmup=0
#     ),
# }
_SECOND_ORDER_OPTIMIZERS["scales"] = {
    "optimizer": FusedNewtonOptimizerConfig(mode="scale", lr=1.0e-6, eps=1e-15),
    "scheduler": ExponentialDecaySchedulerConfig(
        lr_final=1.0e-8, max_steps=30000, #warmup_steps=1000, lr_pre_warmup=0
    ),
}
# TODO: investigate whether this messes up MCMC densification
# _SECOND_ORDER_OPTIMIZERS["opacities"] = {
#     "optimizer": FusedNewtonOptimizerConfig(mode="opacity", lr=1.0e-6, eps=1e-15),
#     "scheduler": ExponentialDecaySchedulerConfig(
#         lr_final=1.0e-8, max_steps=30000, #warmup_steps=3000, lr_pre_warmup=0
#     ),
# }
# _SECOND_ORDER_OPTIMIZERS["features_dc"] = {
#     "optimizer": FusedNewtonOptimizerConfig(mode="color", lr=1.0e-6, eps=1e-15),
#     "scheduler": ExponentialDecaySchedulerConfig(
#         lr_final=1.0e-8, max_steps=30000, #warmup_steps=1000, lr_pre_warmup=0
#     ),
# }


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
            ),
            
        ),
        optimizers={
            **_DEFAULT_OPTIMIZERS_WITH_SCALE_SCHEDULER
        },
        viewer=ViewerConfig(),
        vis="viewer",
    ),
    description="Spirulae 3DGS Default.",
)

spirulae_squared = MethodSpecification(
    config=TrainerConfig(
        method_name="spirulae^2",
        steps_per_eval_batch=0,
        steps_per_save=2000,
        max_num_iterations=30000,
        mixed_precision=False,
        pipeline=SpirulaePipelineConfig(
            datamanager=SpirulaeDataManagerConfig(
                **_DEFAULT_DATAMANAGER_CONFIG
            ),
            model=SpirulaeModelConfig(
                compute_hessian_diagonal="all",
                mcmc_noise_lr=5e5 * (1.6e-4 / 1.0e-6),
                # erank_reg=0.0,
                # erank_reg_s3=0.0,
            ),
        ),
        optimizers={
            **_SECOND_ORDER_OPTIMIZERS
        },
        viewer=ViewerConfig(),
        vis="viewer",
    ),
    description="Spirulae 3DGS Default.",
)

spirulae_squared_pos = MethodSpecification(
    config=TrainerConfig(
        method_name="spirulae^2-pos",
        steps_per_eval_batch=0,
        steps_per_save=2000,
        max_num_iterations=30000,
        mixed_precision=False,
        pipeline=SpirulaePipelineConfig(
            datamanager=SpirulaeDataManagerConfig(
                **_DEFAULT_DATAMANAGER_CONFIG
            ),
            model=SpirulaeModelConfig(
                compute_hessian_diagonal="position",
                mcmc_noise_lr=5e5 * (1.6e-4 / 1.0e-6),
            ),
        ),
        optimizers={
            **_SECOND_ORDER_POSITION_OPTIMIZERS
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
                max_batch_per_epoch=800,
            ),
            model=SpirulaeModelConfig(
                # packed=True,
                use_bvh=True,
                use_camera_optimizer=False,
                # use_bilateral_grid=False,
                # use_bilateral_grid_for_geometry=False,
                alpha_reg_weight=0.0,
                mcmc_max_screen_size=0.15,
                mcmc_max_screen_size_clip_hardness=1.5,
                # mcmc_max_world_size=1.0,
            ),
        ),
        optimizers={
            **_DEFAULT_OPTIMIZERS_WITH_SCALE_SCHEDULER
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
                **_DEFAULT_DATAMANAGER_CONFIG,
                compute_visibility_masks=True,
            ),
            model=SpirulaeModelConfig(
                primitive="opaque_triangle",
                kernel_radius=0.5,
                compute_depth_normal=True,
                sh_degree=0,
                bilagrid_shape=(16, 16, 8),
                stop_refine_at=30000,
                # alpha_reg_weight=0.0,
                mcmc_scale_reg=0.001,
                # erank_reg=1.0,
                # supersampling=2,
                mcmc_min_opacity=0.005,
                mcmc_noise_lr=5e5,  # or 0.0
                mcmc_max_screen_size=0.15,
                supervision_warmup=0,
                depth_supervision_weight=0.0,
                normal_supervision_weight=0.04,
                depth_reg_weight=0.01,
                normal_distortion_reg_weight=0.005,
                rgb_distortion_reg_weight=0.01,
                ssim_lambda=0.4,
                preallocate_splat_tensors=False,  # TODO
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
                compute_visibility_masks=True,
                patch_batch_size=-1,
                patch_size=64,
                max_batch_per_epoch=800,
            ),
            model=SpirulaeModelConfig(
                primitive="opaque_triangle",
                kernel_radius=0.5,
                compute_depth_normal=True,
                sh_degree=0,
                bilagrid_shape=(16, 16, 8),
                stop_refine_at=30000,
                # alpha_reg_weight=0.0,
                mcmc_scale_reg=0.001,
                # erank_reg=1.0,
                # supersampling=2,
                mcmc_min_opacity=0.005,
                mcmc_noise_lr=5e5,  # or 0.0
                supervision_warmup=0,
                depth_supervision_weight=0.0,
                normal_supervision_weight=0.04,
                depth_reg_weight=0.01,
                normal_distortion_reg_weight=0.005,
                rgb_distortion_reg_weight=0.01,
                ssim_lambda=0.4,
                reg_warmup_length=1000,
                preallocate_splat_tensors=False,  # TODO

                # packed=True,
                use_bvh=True,
                use_camera_optimizer=False,
                # use_bilateral_grid=False,
                # use_bilateral_grid_for_geometry=False,
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

spirulae_voxel = MethodSpecification(
    config=TrainerConfig(
        method_name="spirulae-voxel",
        steps_per_eval_batch=0,
        steps_per_save=2000,
        max_num_iterations=30000,
        mixed_precision=False,
        pipeline=SpirulaePipelineConfig(
            datamanager=SpirulaeDataManagerConfig(
                **_DEFAULT_DATAMANAGER_CONFIG,
            ),
            model=SpirulaeModelConfig(
                primitive="voxel",
                use_mcmc=False,
                sh_degree=0,
                background_color="black",
                train_background_color=False,
                # use_bilateral_grid=False,
                use_bilateral_grid_for_geometry=False,
                alpha_reg_weight=0.0,
                # alpha_loss_weight=0.0,
                # alpha_loss_weight_under=0.0,
                mcmc_opacity_reg=0.0001,
                mcmc_scale_reg=0.0,
                erank_reg=0.0,
                erank_reg_s3=0.0,
                depth_supervision_weight=0.0,
                supervision_warmup=0,
                depth_reg_weight=0.01,
                reg_warmup_length=100,
                distortion_reg_warmup=100,
                preallocate_splat_tensors=False,  # TODO
            ),
        ),
        optimizers={
            **_DEFAULT_OPTIMIZERS
        },
        viewer=ViewerConfig(),
        vis="viewer",
    ),
    description="Spirulae sparse voxel grid with patch batching.",
)



_MODEL_PRESET_CONFINED = dict(
    randomize_background=True,
    train_background_color=False,
)
_MODEL_PRESET_3DGS2TR_POS = dict(
    compute_hessian_diagonal="position",
    mcmc_noise_lr=5e5 * (1.6e-4 / 1.0e-6),
)
_MODEL_PRESET_3DGS2TR = dict(
    compute_hessian_diagonal="all",
    mcmc_noise_lr=5e5 * (1.6e-4 / 1.0e-6),
)
_MODEL_PRESET_RICH_TEXTURE = dict(
    mcmc_prob_grad_weight=1.0,
    mcmc_use_long_axis_split=True,
    mcmc_max_screen_size=0.10,
)
_MODEL_PRESET_NO_COLOR_SHIFT = dict(
    use_bilateral_grid=True,
    bilagrid_shape=(8, 8, 4),
    bilagrid_type="ppisp",
)


spirulae_preset_confined_low_texture = MethodSpecification(
    config=TrainerConfig(
        method_name="spirulae-preset-confined-low-texture",
        steps_per_eval_batch=0,
        steps_per_save=2000,
        max_num_iterations=30000,
        mixed_precision=False,
        pipeline=SpirulaePipelineConfig(
            datamanager=SpirulaeDataManagerConfig(
                **_DEFAULT_DATAMANAGER_CONFIG
            ),
            model=SpirulaeModelConfig(
                **_MODEL_PRESET_CONFINED,
                **_MODEL_PRESET_NO_COLOR_SHIFT,
            ),
        ),
        optimizers={
            **_DEFAULT_OPTIMIZERS_WITH_SCALE_SCHEDULER
        },
        viewer=ViewerConfig(),
        vis="viewer",
    ),
    description="Spirulae-Splat 3DGS preset for visually confined low-texture environments (e.g. indoor or ground-level city capture with largely textureless surfaces).",
)

spirulae_preset_confined = MethodSpecification(
    config=TrainerConfig(
        method_name="spirulae-preset-confined",
        steps_per_eval_batch=0,
        steps_per_save=2000,
        max_num_iterations=30000,
        mixed_precision=False,
        pipeline=SpirulaePipelineConfig(
            datamanager=SpirulaeDataManagerConfig(
                **_DEFAULT_DATAMANAGER_CONFIG
            ),
            model=SpirulaeModelConfig(
                **_MODEL_PRESET_3DGS2TR_POS,
                **_MODEL_PRESET_CONFINED,
                **_MODEL_PRESET_RICH_TEXTURE,
            ),
        ),
        optimizers={
            **_SECOND_ORDER_POSITION_OPTIMIZERS
        },
        viewer=ViewerConfig(),
        vis="viewer",
    ),
    description="Spirulae-Splat 3DGS preset for visually confined environments with moderate to rich texture (e.g. rich-texture indoor environments, ground-level city or foliage captures without much sky visible in input images).",
)

spirulae_squared_preset_confined = MethodSpecification(
    config=TrainerConfig(
        method_name="spirulae^2-preset-confined",
        steps_per_eval_batch=0,
        steps_per_save=2000,
        max_num_iterations=30000,
        mixed_precision=False,
        pipeline=SpirulaePipelineConfig(
            datamanager=SpirulaeDataManagerConfig(
                **_DEFAULT_DATAMANAGER_CONFIG
            ),
            model=SpirulaeModelConfig(
                **_MODEL_PRESET_3DGS2TR,
                **_MODEL_PRESET_CONFINED,
                **_MODEL_PRESET_RICH_TEXTURE,
            ),
        ),
        optimizers={
            **_SECOND_ORDER_OPTIMIZERS
        },
        viewer=ViewerConfig(),
        vis="viewer",
    ),
    description="Spirulae-Splat 3DGS preset for visually confined environments with moderate to rich texture (e.g. rich-texture indoor environments, ground-level city or foliage captures without much sky visible in input images).",
)

spirulae_preset_open_low_texture = MethodSpecification(
    config=TrainerConfig(
        method_name="spirulae-preset-open-low-texture",
        steps_per_eval_batch=0,
        steps_per_save=2000,
        max_num_iterations=30000,
        mixed_precision=False,
        pipeline=SpirulaePipelineConfig(
            datamanager=SpirulaeDataManagerConfig(
                **_DEFAULT_DATAMANAGER_CONFIG
            ),
            model=SpirulaeModelConfig(
                **_MODEL_PRESET_NO_COLOR_SHIFT,
                relative_scale=10.0,
            ),
        ),
        optimizers={
            **_DEFAULT_OPTIMIZERS_WITH_SCALE_SCHEDULER
        },
        viewer=ViewerConfig(),
        vis="viewer",
    ),
    description="Spirulae-Splat 3DGS preset for open low-texture environments (e.g. outdoor scenes containing large textureless concrete, water, sand, or snow surface).",
)

spirulae_preset_open = MethodSpecification(
    config=TrainerConfig(
        method_name="spirulae-preset-open",
        steps_per_eval_batch=0,
        steps_per_save=2000,
        max_num_iterations=30000,
        mixed_precision=False,
        pipeline=SpirulaePipelineConfig(
            datamanager=SpirulaeDataManagerConfig(
                **_DEFAULT_DATAMANAGER_CONFIG
            ),
            model=SpirulaeModelConfig(
                **_MODEL_PRESET_3DGS2TR_POS,
                **_MODEL_PRESET_RICH_TEXTURE,
            ),
        ),
        optimizers={
            **_SECOND_ORDER_POSITION_OPTIMIZERS
        },
        viewer=ViewerConfig(),
        vis="viewer",
    ),
    description="Spirulae-Splat 3DGS preset for open environments with moderate to rich texture (e.g. outdoor landscapes, drone capture).",
)

spirulae_squared_preset_open = MethodSpecification(
    config=TrainerConfig(
        method_name="spirulae^2-preset-open",
        steps_per_eval_batch=0,
        steps_per_save=2000,
        max_num_iterations=30000,
        mixed_precision=False,
        pipeline=SpirulaePipelineConfig(
            datamanager=SpirulaeDataManagerConfig(
                **_DEFAULT_DATAMANAGER_CONFIG
            ),
            model=SpirulaeModelConfig(
                **_MODEL_PRESET_3DGS2TR,
                **_MODEL_PRESET_RICH_TEXTURE,
            ),
        ),
        optimizers={
            **_SECOND_ORDER_OPTIMIZERS
        },
        viewer=ViewerConfig(),
        vis="viewer",
    ),
    description="Spirulae-Splat 3DGS preset for open environments with moderate to rich texture (e.g. outdoor landscapes, drone capture).",
)

spirulae_preset_centered_object = MethodSpecification(
    config=TrainerConfig(
        method_name="spirulae-preset-centered-object",
        steps_per_eval_batch=0,
        steps_per_save=2000,
        max_num_iterations=30000,
        mixed_precision=False,
        pipeline=SpirulaePipelineConfig(
            datamanager=SpirulaeDataManagerConfig(
                dataparser=Nerfstudio2DataParserConfig(
                    eval_mode="fraction",
                    train_split_fraction=1.0,
                    center_method="focus",
                ),
            ),
            model=SpirulaeModelConfig(
                **_MODEL_PRESET_RICH_TEXTURE,
                **_MODEL_PRESET_NO_COLOR_SHIFT,
                apply_loss_for_mask=True,
                mcmc_cap_max=100000,
            ),
        ),
        optimizers={
            # **_DEFAULT_OPTIMIZERS_WITH_SCALE_SCHEDULER
            **_SECOND_ORDER_POSITION_OPTIMIZERS
        },
        viewer=ViewerConfig(),
        vis="viewer",
    ),
    description="Spirulae-Splat 3DGS preset for small centered objects (e.g. small item on turntable, human avatar).",
)

spirulae_preset_academic_baseline = MethodSpecification(
    config=TrainerConfig(
        method_name="spirulae-preset-academic-baseline",
        steps_per_eval_batch=0,
        steps_per_save=2000,
        max_num_iterations=30000,
        mixed_precision=False,
        pipeline=SpirulaePipelineConfig(
            datamanager=SpirulaeDataManagerConfig(
                dataparser=Nerfstudio2DataParserConfig(
                    eval_mode="interval",
                    eval_interval=8,
                ),
                max_batch_per_epoch=9**9,
                compute_edge_maps=False,
            ),
            model=SpirulaeModelConfig(
                primitive="3dgs",
                background_color="black",
                train_background_color=False,
                compute_depth_normal=False,
                use_bilateral_grid=False,
                use_bilateral_grid_for_geometry=False,
                use_ppisp=False,
                mcmc_max_screen_size=float('inf'),
                suppress_initial_scales=False,
                depth_reg_weight=0.0,
                normal_reg_weight=0.0,
                alpha_reg_weight=0.0,
                alpha_loss_weight=0.0,
                alpha_loss_weight_under=0.0,
                erank_reg=0.0,
                erank_reg_s3=0.0,
                quat_norm_reg=0.0,
                randomize_background=0.0,
                normal_supervision_weight=0.0,
                mcmc_opacity_reg=0.01,
                mcmc_scale_reg=0.01,
                use_edge_aware_score=False,
            ),
        ),
        optimizers={
            **_DEFAULT_OPTIMIZERS
        },
        viewer=ViewerConfig(),
        vis="viewer",
    ),
    description="Spirulae-Splat 3DGS preset that replicates academic baseline (3DGS MCMC) as faithful as possible. Use this to reproduce metrics on academic benchmark datasets (e.g. Mip-NeRF 360).",
)
