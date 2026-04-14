import torch
import torch.nn.functional as F
import numpy as np


def render_equirectangular(base_poses, render_fn, W_out, H_out):
    """
    Args:
        base_poses: (B, 4, 4) World-to-Camera tensors.
        render_fn: function(poses, intrinsics, W, H) -> (B, H, W, 3)
        W_out, H_out: Output equirectangular dimensions.
    """
    device = base_poses.device
    B = base_poses.shape[0]
    
    # 1. Determine sufficient resolution for cube faces
    # To match the resolution at the equator: S = W_out / pi
    # To match the resolution at the poles: S = H_out
    # We'll use H_out as a high-quality heuristic.
    S = H_out 
    
    # 2. Define Intrinsics for 90-degree FOV
    # fx = S / (2 * tan(fov/2)) = S / (2 * tan(45deg)) = S / 2
    intrinsics = torch.tensor([[S/2, S/2, S/2, S/2]], device=device).repeat(B, 1)

    # 3. Define rotations for the 6 cube faces (Camera coordinates: +z forward, +x right, +y down)
    def rotation_matrix(axis, angle_deg):
        angle = np.radians(angle_deg)
        c, s = np.cos(angle), np.sin(angle)
        if axis == 'x':
            return torch.tensor([[1, 0, 0], [0, c, -s], [0, s, c]], dtype=torch.float32)
        if axis == 'y':
            return torch.tensor([[c, 0, s], [0, 1, 0], [-s, 0, c]], dtype=torch.float32)
        return torch.eye(3)

    # Rotations from base camera orientation to face orientation
    face_rots = {
        'front': torch.eye(3),
        'back':  rotation_matrix('y', 180),
        'left':  rotation_matrix('y', -90),
        'right': rotation_matrix('y', 90),
        'up':    rotation_matrix('x', -90),
        'down':  rotation_matrix('x', 90)
    }

    # 4. Render the 6 faces
    cube_faces = {}
    for name, R_f in face_rots.items():
        # Adjust World-to-Camera: P_new = R_f^-1 * P_base
        # Since R_f is a rotation matrix, R_f^-1 = R_f^T
        m = torch.eye(4, device=device)
        m[:3, :3] = R_f.T
        face_pose = m @ base_poses
        cube_faces[name] = render_fn(face_pose, intrinsics, S, S) # (B, S, S, 3)

    # 5. Create Equirectangular grid
    # Create normalized coordinates (phi: latitude, theta: longitude)
    theta = torch.linspace(-np.pi, np.pi, W_out, device=device)
    phi = torch.linspace(-np.pi/2, np.pi/2, H_out, device=device)
    phi_grid, theta_grid = torch.meshgrid(phi, theta, indexing='ij')

    # Convert spherical to Cartesian unit vectors (Camera Space)
    x = torch.cos(phi_grid) * torch.sin(theta_grid)
    y = torch.sin(phi_grid)  # In many CV conventions, y is down, so check your axis
    z = torch.cos(phi_grid) * torch.cos(theta_grid)
    
    # 6. Sample from cube faces using the unit vectors
    # We determine which face each vector hits based on the maximum absolute component
    vectors = torch.stack([x, y, z], dim=-1) # (H, W, 3)
    abs_vec = torch.abs(vectors)
    max_axis = torch.argmax(abs_vec, dim=-1)
    
    out_img = torch.zeros((B, H_out, W_out, 3), device=device)

    # Mapping logic for each face
    # We convert the 3D vector into 2D (u, v) coordinates for the hit face
    face_configs = [
        (2,  1, 'front', lambda v: ( v[...,0]/v[...,2],  v[...,1]/v[...,2])), # +z
        (2, -1, 'back',  lambda v: ( v[...,0]/v[...,2], -v[...,1]/v[...,2])), # -z
        (0,  1, 'right', lambda v: (-v[...,2]/v[...,0],  v[...,1]/v[...,0])), # +x
        (0, -1, 'left',  lambda v: (-v[...,2]/v[...,0], -v[...,1]/v[...,0])), # -x
        (1, -1, 'down',  lambda v: (-v[...,0]/v[...,1], -v[...,2]/v[...,1])), # +y
        (1,  1, 'up',    lambda v: ( v[...,0]/v[...,1], -v[...,2]/v[...,1])), # -y
    ]

    for axis_idx, sign, name, uv_fn in face_configs:
        mask = (max_axis == axis_idx) & (torch.sign(vectors[..., axis_idx]) == sign)
        if not mask.any(): continue
        
        # Get UVs for this face, range [-1, 1]
        u, v = uv_fn(vectors[mask])
        
        # grid_sample expects (N, H_out, W_out, 2) in range [-1, 1]
        # We process each image in the batch
        face_data = cube_faces[name].permute(0, 3, 1, 2) # (B, 3, S, S)
        grid = torch.stack([u, v], dim=-1).unsqueeze(0).repeat(B, 1, 1) # (B, num_pixels, 2)
        grid = grid.view(B, 1, -1, 2) 
        
        sampled = F.grid_sample(face_data, grid, align_corners=True) # (B, 3, 1, num_pixels)
        out_img[:, mask] = sampled.squeeze(2).permute(0, 2, 1)

    return out_img
