"""U2Net with bilinear upsample replaced by pixel shuffle
Based on https://github.com/xuebinqin/U-2-Net/blob/master/model/u2net_refactor.py
License: Apache 2.0 (https://github.com/xuebinqin/U-2-Net?tab=Apache-2.0-1-ov-file)
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torchvision import transforms

import math

from typing import Optional


def _upsample_like(x, size):
    return nn.Upsample(size=size, mode='bilinear', align_corners=False)(x)

def _size_map(x, height):
    # {height: size} for Upsample
    size = list(x.shape[-2:])
    sizes = {}
    for h in range(1, height):
        sizes[h] = size
        size = [math.ceil(w / 2) for w in size]
    return sizes


class ConvBnRelu(nn.Module):
    def __init__(self, c_in, c_out, dilate=1):
        super(ConvBnRelu, self).__init__()

        self.main = nn.Sequential(
            nn.Conv2d(c_in, c_out, 3, padding=dilate, dilation=dilate),
            nn.BatchNorm2d(c_out),
            nn.ReLU(inplace=True),
        )

    def forward(self, x):
        return self.main(x)


class RSU(nn.Module):
    def __init__(self, name, height, in_ch, mid_ch, out_ch, dilated=False):
        super(RSU, self).__init__()
        self.name = name
        self.height = height
        self.dilated = dilated
        self._make_layers(height, in_ch, mid_ch, out_ch, dilated)

    def forward(self, x):
        x = self.rebnconvin(x)

        def unet(x, height=1):
            if height < self.height:
                x1 = getattr(self, f'rebnconv{height}')(x)
                if not self.dilated and height < self.height - 1:
                    x2 = unet(self.downsample(x1), height + 1)
                else:
                    x2 = unet(x1, height + 1)

                x = getattr(self, f'rebnconv{height}d')(torch.cat((x2, x1), 1))
                del x1
                del x2
                if self.dilated or height == 1:
                    return x
                return F.pixel_shuffle(x, 2)
            else:
                return getattr(self, f'rebnconv{height}')(x)

        return x + unet(x)

    def _make_layers(self, height, in_ch, mid_ch, out_ch, dilated=False):
        self.add_module('rebnconvin', ConvBnRelu(in_ch, out_ch))
        self.add_module('downsample', nn.MaxPool2d(2, stride=2, ceil_mode=True))

        self.add_module(f'rebnconv1', ConvBnRelu(out_ch, mid_ch))
        self.add_module(f'rebnconv1d', ConvBnRelu(mid_ch*2, out_ch))

        for i in range(2, height):
            dilate = 1 if not dilated else 2 ** (i - 1)
            self.add_module(f'rebnconv{i}', ConvBnRelu(mid_ch, mid_ch, dilate=dilate))
            self.add_module(f'rebnconv{i}d', ConvBnRelu(mid_ch*2, mid_ch if dilated else mid_ch*4, dilate=dilate))

        dilate = 2 if not dilated else 2 ** (height - 1)
        self.add_module(f'rebnconv{height}', ConvBnRelu(mid_ch, mid_ch, dilate=dilate))


class U2Net(nn.Module):
    def __init__(self, cfgs, out_ch):
        super(U2Net, self).__init__()
        self.out_ch = out_ch
        self._make_layers(cfgs)

    def forward(self, x):
        sizes = _size_map(x, self.height)
        maps = []

        def unet(x, height=1):
            if height < self.height:
                x1 = getattr(self, f'stage{height}')(x)
                x2 = unet(getattr(self, 'downsample')(x1), height + 1)
                x = getattr(self, f'stage{height}d')(torch.cat((x2, x1), 1))
            else:
                x = getattr(self, f'stage{height}')(x)
            x = F.pixel_shuffle(x, 2)
            side(x, height)
            return x

        def side(x, h):
            if True: x = F.avg_pool2d(x, 2)
            x = getattr(self, f'side{h}')(x)
            x = _upsample_like(x, sizes[1])
            maps.append(x)

        def fuse():
            maps.reverse()
            x = torch.cat(maps, 1)
            x = getattr(self, 'outconv')(x)
            return x

        unet(x)
        return fuse()

    def _make_layers(self, cfgs):
        self.height = int((len(cfgs) + 1) / 2)
        self.add_module('downsample', nn.MaxPool2d(2, stride=2, ceil_mode=True))
        for k, v in cfgs.items():
            self.add_module(k, RSU(v[0], *v[1]))
            if v[2] > 0:
                self.add_module(f'side{v[0][-1]}', nn.Conv2d(v[2], self.out_ch, 3, padding=1))
        self.add_module('outconv', nn.Conv2d(int(self.height * self.out_ch), self.out_ch, 1))


class Model(torch.nn.Module):
    def __init__(self, m1, m2):
        super(Model, self).__init__()

        self.mean = [0.485, 0.456, 0.406]
        self.std = [0.229, 0.224, 0.225]
        self.normalize = transforms.Normalize(mean=self.mean, std=self.std)

        self.u2net: Optional[U2Net] = U2Net({
            'stage1': ['En_1', (6, 3, m1, m2), -1],
            'stage2': ['En_2', (5, m2, m1, 2*m2), -1],
            'stage3': ['En_3', (4, 2*m2, 2*m1, 4*m2), -1],
            'stage4': ['En_4', (3, 4*m2, 4*m1, 4*m2, True), -1],
            'stage5': ['En_5', (3, 4*m2, 4*m1, 4*4*m2, True), 4*m2],
            'stage4d': ['De_4', (3, 8*m2, 4*m1, 4*4*m2, True), 4*m2],
            'stage3d': ['De_3', (4, 8*m2, 2*m1, 4*2*m2), 2*m2],
            'stage2d': ['De_2', (5, 4*m2, m1, 4*m2), m2],
            'stage1d': ['De_1', (6, 2*m2, m1, 4*m2), m2],
        }, out_ch=3)

    def forward(self, x0):
        x = x0
        x = self.normalize(x)

        r = self.u2net(x)

        y = torch.sigmoid(torch.logit(torch.clip(x0.float(),0.5/255,1-0.5/255))+torch.sinh(math.sqrt(0.5/math.pi)*r.float()))

        return y


torch.backends.cudnn.benchmark = False

model = None
dtype = torch.bfloat16

def infer(batch):
    batch = batch.contiguous(memory_format=torch.channels_last)

    batch_dtype = batch.dtype
    if batch_dtype == torch.uint8:
        batch = batch.float() / 255.0

    B, H, W, C = batch.shape
    s = 32
    batch = batch.permute(0, 3, 1, 2)[:, :, :s*(H//s), :s*(W//s)]

    global model
    if model is None:
        # TODO: load from URL when shipping this code out
        path = "/mnt/d/gs/ssplat-image-enhance/checkpoint.pth"
        model = Model(32, 64)
        model.load_state_dict(torch.load(path, weights_only=False, map_location="cpu"))
        model = model.to(memory_format=torch.channels_last)
        model = model.cuda().to(dtype)
        torch.cuda.empty_cache()

    with torch.no_grad():
        with torch.autocast(device_type="cuda", dtype=dtype):
            pred = model(batch.cuda().to(dtype)).to(device=batch.device, dtype=batch.dtype)

    if batch_dtype == torch.uint8:
        pred = (255.0 * pred).to(batch_dtype)
    return pred.permute(0, 2, 3, 1)
