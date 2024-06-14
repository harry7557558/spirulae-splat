import torch
from torch import Tensor, nn

from nerfstudio.cameras.rays import RayBundle
from nerfstudio.field_components import MLP
from nerfstudio.field_components.encodings import SHEncoding
from nerfstudio.fields.base_field import Field, get_normalized_directions


class SineActivation(nn.Module):
    def __init__(self):
        super(SineActivation, self).__init__()
    
    def forward(self, x):
        return torch.sin(x)


class BG_Field(Field):
    def __init__(self, appearance_embedding_dim=0, implementation="torch"):
        super(BG_Field, self).__init__()

        self.direction_encoding = SHEncoding(
            levels=4,
            implementation=implementation,
        )

        self.mlp_background_color = MLP(
            in_dim=self.direction_encoding.get_out_dim()+appearance_embedding_dim,
            num_layers=4,
            layer_width=32,
            out_dim=3,
            activation=SineActivation(),
            out_activation=nn.Sigmoid(),
            implementation=implementation,
        )

    def get_background_rgb(self, ray_bundle: RayBundle, appearance_embedding=None) -> Tensor:
        """Predicts background colors at infinity."""
        directions = get_normalized_directions(ray_bundle.directions)
        directions_flat = self.direction_encoding(directions.view(-1, 3))
        if appearance_embedding is not None:
            x = torch.cat([directions_flat, appearance_embedding], dim=-1)
        else:
            x = directions_flat

        background_rgb = self.mlp_background_color(x).to(directions)

        return background_rgb
