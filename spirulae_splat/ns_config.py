"""
Nerfstudio Template Config

Define your custom method here that registers with Nerfstudio CLI.
"""

from __future__ import annotations

from spirulae_splat.ns_datamanager import (
    SpirulaeDataManagerConfig,
)
from spirulae_splat.ns_model import SpirulaeModelConfig
from spirulae_splat.ns_pipeline import (
    SpirulaePipelineConfig,
)
from nerfstudio.configs.base_config import ViewerConfig
from nerfstudio.data.dataparsers.nerfstudio_dataparser import NerfstudioDataParserConfig
from nerfstudio.engine.optimizers import AdamOptimizerConfig, RAdamOptimizerConfig
from nerfstudio.engine.schedulers import (
    ExponentialDecaySchedulerConfig,
)
from nerfstudio.engine.trainer import TrainerConfig
from nerfstudio.plugins.types import MethodSpecification


spirulae = MethodSpecification(
    config=TrainerConfig(
        method_name="spirulae",
        steps_per_eval_batch=0,
        steps_per_save=2500,
        max_num_iterations=30000,
        mixed_precision=False,
        pipeline=SpirulaePipelineConfig(
            datamanager=SpirulaeDataManagerConfig(
                dataparser=NerfstudioDataParserConfig(
                    train_split_fraction=1.0,
                    load_3D_points=True,
                ),
                camera_res_scale_factor=1.0,
                eval_num_images_to_sample_from=-1,
                eval_num_times_to_repeat_images=-1,
                eval_image_indices=(0,),
                cache_images='cpu',
                cache_images_type='uint8',
                max_thread_workers=None
            ),
            model=SpirulaeModelConfig(
                eval_num_rays_per_chunk=4096,
            ),
        ),
        optimizers={
            "means": {
                "optimizer": AdamOptimizerConfig(lr=1.6e-4, eps=1e-15),
                "scheduler": ExponentialDecaySchedulerConfig(
                    lr_final=1.6e-6,
                    max_steps=30000,
                ),
            },
            "scales": {
                "optimizer": AdamOptimizerConfig(lr=0.005, eps=1e-15),
                "scheduler": None,
            },
            "quats": {
                "optimizer": AdamOptimizerConfig(lr=0.001, eps=1e-15),
                "scheduler": None
            },
            "features_dc": {
                "optimizer": AdamOptimizerConfig(lr=0.0025, eps=1e-15),
                "scheduler": None,
            },
            "features_rest": {
                "optimizer": AdamOptimizerConfig(lr=0.0025 / 20, eps=1e-15),
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
            "camera_opt": {
                "optimizer": AdamOptimizerConfig(lr=1e-4, eps=1e-15),
                "scheduler": ExponentialDecaySchedulerConfig(
                    lr_final=5e-7, max_steps=30000, warmup_steps=1000, lr_pre_warmup=0
                ),
            },
        },
        viewer=ViewerConfig(num_rays_per_chunk=1 << 15),
        vis="viewer",
    ),
    description="Nerfstudio method template.",
)
