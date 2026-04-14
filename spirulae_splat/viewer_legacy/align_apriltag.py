import numpy as np
import cv2
import os
import json
import yaml

from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed


def load_frames(dataset_dir):
    transforms_path = os.path.join(dataset_dir, "transforms.json")
    if not os.path.exists(transforms_path):
        raise ValueError(f"transforms.json not found in {dataset_dir}")

    with open(transforms_path) as fp:
        transforms = json.load(fp)
    frames = transforms['frames']
    for frame in frames:
        frame['file_path'] = os.path.join(dataset_dir, frame['file_path'])
    if 'w' in transforms:
        for key in 'w h fl_x fl_y cx cy k1 k2 k3 k4 p1 p2 camera_model'.split():
            if key not in transforms:
                continue
            for frame in frames:
                frame[key] = transforms[key]
    return frames


def load_markers(dataset_dir):
    yaml_path = os.path.join(dataset_dir, "markers.yaml")
    if not os.path.exists(yaml_path):
        raise ValueError(f"markers.yaml not found in {dataset_dir}")
    with open(yaml_path, 'r') as fp:
        content = yaml.safe_load(fp)

    apriltags = {}
    for item in content.values():
        assert 'type' in item, "Missing attribute `type`"
        assert item['type'] == 'apriltag', "Only support apriltag marker at this time"
        assert 'family' in item and 'id' in item, "Missing apriltag attribute `family` and `id`"
        assert 'center' in item and isinstance(item['center'], list) and len(item['center']) == 3, 'Missing or incorrect attribute `center` (must be list of 3 numbers)'
        tag_id = (item['family'], item['id'])
        center = np.array(item['center'])
        assert tag_id not in apriltags, f"Apriltag (family, id) pair must be unique; You have more than one {tag_id}"
        apriltags[tag_id] = center

    assert len(apriltags) >= 3, "At least 3 apriltags needed"

    centers = np.stack([*apriltags.values()])
    U, S, Vt = np.linalg.svd(centers)
    S = sorted(S)
    assert S[1] > 0.02 * S[2], f"Apriltags must not be placed on the same line"

    return apriltags


def _detect_apriltags_one_image(image_path, families):

    # warning: apriltag package sometimes gives segfault
    try:
        import apriltag
    except ImportError:
        print("Failed to import apriltag. Please run `pip install apriltag` first.")
        exit(0)

    image = cv2.imread(image_path)
    if image is None:
        return []

    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    options = apriltag.DetectorOptions(families=families)
    detector = apriltag.Detector(options)
    results = detector.detect(gray)
    # results = [r for r in results if r.decision_margin > 50]
    return results


def detect_apriltags(frames, families):
    results = {}
    prev_num_results = len(results)
    while len(results) < len(frames):
        with ProcessPoolExecutor() as executor:
            futures = {executor.submit(_detect_apriltags_one_image, frame['file_path'], families): frame['file_path']
                       for frame in frames if frame['file_path'] not in results}
            
            for future in as_completed(futures):
                file_path = futures[future]
                try:
                    data = future.result()
                    results[file_path] = data
                except Exception as e:
                    # print(e)  # likely segfault in apriltag package
                    pass

        if len(results) == prev_num_results:
            break
        prev_num_results = len(results)

    results = {key: result for key, result in results.items() if len(result) > 0}
    return results


def closest_point_to_weighted_rays(ray_origin, ray_dir, weight):
    
    n = ray_origin.shape[0]
    I = np.eye(3)
    
    M = np.zeros((3, 3))
    b = np.zeros(3)
    
    for i in range(n):
        d = ray_dir[i]
        norm_d2 = np.dot(d, d)
        A = I - np.outer(d, d) / norm_d2
        wA = weight[i, 0] * A
        M += wA
        b += wA @ ray_origin[i]
    
    p = np.linalg.solve(M, b)
    
    t = np.sum((p - ray_origin) * ray_dir, axis=1) / np.sum(ray_dir**2, axis=1)
    
    return t, p.flatten()


def find_points_from_detections(frames, detections):
    frames = {frame['file_path']: frame for frame in frames}

    # exact detections as `id: [list of camera rays]`
    apriltags = {}
    for key, results in detections.items():
        frame = frames[key]
        c2w = np.array(frame['transform_matrix'])
        camera_model = frame['camera_model']
        # TODO: support other distortion free ones compatible with OPENCV
        assert camera_model in ["OPENCV", "OPENCV_FISHEYE"], "Only OPENCV and OPENCV_FISHEYE camera models are supported at this time"
        K = np.zeros((3, 3), dtype=np.float32)
        K[0,0], K[1,1], K[0,2], K[1,2] = [frame[key] for key in 'fl_x fl_y cx cy'.split()]
        dist_coeffs = np.array([
            frame[key] for key in
            ('k1 k2 p1 p2' if camera_model == "OPENCV" else 'k1 k2 k3 k4').split()
        ])
        cam_pos = c2w[:3, 3]
        for res in results:
            tag_id = (res.tag_family.decode('ascii'), res.tag_id)
            weight = res.decision_margin
            center = res.center
            try:
                undistort = cv2.undistortPoints if camera_model == "OPENCV" else cv2.fisheye.undistortPoints
                center = undistort(center.reshape((1, 1, 2)), K, dist_coeffs).squeeze()
                assert np.isfinite(center).all()
            except:
                continue
            ray = np.concatenate((center * [1, -1], [-1.0]))
            ray = c2w[:3, :3] @ ray
            ray /= np.linalg.norm(ray)
            if tag_id not in apriltags:
                apriltags[tag_id] = []
            apriltags[tag_id].append(np.concatenate((cam_pos, ray, [weight])))

    # find points
    results = {}
    for key, rays in apriltags.items():
        rays = np.stack(rays)
        ro, rd, w = rays[:, 0:3], rays[:, 3:6], rays[:, 6:7]
        t, pos = closest_point_to_weighted_rays(ro, rd, w)
        results[key] = pos
    return results


def align_points(points, markers):
    p1, p0 = [], []
    for key in points.keys():
        if key in markers:
            p1.append(points[key])
            p0.append(markers[key])
    p1, p0 = np.array(p1), np.array(p0)

    # Fit p0 ~ s * p1 @ R.T + t
    c1 = np.mean(p1, axis=0)
    c0 = np.mean(p0, axis=0)
    p1c, p0c = p1 - c1, p0 - c0
    H = np.dot(p1c.T, p0c)
    U, S, Vt = np.linalg.svd(H)
    V = Vt.T
    R = np.dot(V, U.T)
    if np.linalg.det(R) < 0:
        V[:, -1] *= -1
        R = np.dot(V, U.T)
    s = np.sqrt(np.sum(np.square(p0c)) / np.sum(np.square(p1c)))
    t = c0 - s * np.dot(R, c1)

    return s, R, t


def get_alignment(dataset_dir, verbose=False):
    frames = load_frames(dataset_dir)
    markers = load_markers(dataset_dir)

    if verbose:
        print("Detecting apriltags...")
    families = [*set(family for (family, id) in markers.keys())]
    detections = detect_apriltags(frames, families)
    if verbose:
        num_tags = sum([len(v) for v in detections.values()])
        unique_tags = len(set([(v.tag_family, v.tag_id) for v in sum([*detections.values()], [])]))
        print(f"{num_tags} apriltags ({unique_tags} unique) detected in {len(detections)} (out of {len(frames)}) images")
        print()

    points = find_points_from_detections(frames, detections)
    if verbose:
        print("Detected apriltag centers:")
        for (family, id), value in sorted(points.items()):
            print(f"family={family}, id={id}:", value)
        print()
        print("Given apriltag centers:")
        for (family, id), value in sorted(markers.items()):
            print(f"family={family}, id={id}:", value)
        print()

    s, R, t = align_points(points, markers)
    if verbose:
        print('Recovered rotation:', R, sep='\n', end='\n\n')
        print('Recovered translation:', t, sep='\n', end='\n\n')
        print('Recovered scale:', s, end='\n\n')
    return R, t, s


if __name__ == "__main__":
    dataset_dir = "/media/harry/d/gs/data/apriltag_ba_room"
    frames = load_frames(dataset_dir)
    markers = load_markers(dataset_dir)

    families = [*set(family for (family, id) in markers.keys())]
    detections = detect_apriltags(frames, families)
    points = find_points_from_detections(frames, detections)
    print(points)

    print(points)
    print(markers)
    s, R, t = align_points(points, markers)
    print(s)
    print(R)
    print(t)

