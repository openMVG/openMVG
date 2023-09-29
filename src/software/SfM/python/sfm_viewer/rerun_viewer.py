#!/usr/bin/env python3

#This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

# Copyright (c) 2023 Pierre Moulon

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""Example of using Rerun to log and visualize the output of openMVG's sparse reconstruction."""

# An openMVG sfm_data.bin reconstruction can be exported to sfm_data.json by calling:
#./openMVG_main_ConvertSfM_DataFormat -i 'sfm_data_in.bin' -o 'sfm_data.json'
#

import collections
import json
import numpy as np
import rerun as rr  # pip install rerun-sdk
from tqdm import tqdm

from argparse import ArgumentParser
from collections import defaultdict
from pathlib import Path
from scipy.spatial.transform import Rotation

Intrinsics = collections.namedtuple("Intrinsics", ["id", "model", "width", "height", "f", "cx", "cy"])

def read_and_log_sparse_reconstruction(args) -> None:
    print("Reading sparse OpenMVG reconstruction")
    print("Building visualization by logging to Rerun")

    print(f'Reading sfm_data JSON file: {args.sfm_data}')
    f = open(args.sfm_data)
    data = json.load(f)
    root_path = data["root_path"]

    views = data["views"]
    extrinsics = data['extrinsics']
    intrinsics = data['intrinsics']
    structure = data['structure']

    extrinsics_ids = [it['key'] for it in extrinsics]

    # Collect views
    print("Setup Views")
    for v in tqdm(views):
        value = v['value']['ptr_wrapper']['data']
        filename = value['filename']
        id_view = value['id_view']
        id_intrinsic = value['id_intrinsic']
        id_pose = value['id_pose']

        # Collect poses
        if id_pose not in extrinsics_ids:
            continue
        extrinsic = [it for it in extrinsics if it['key'] == id_pose]
        extrinsic = extrinsic[0]
        rot = extrinsic['value']['rotation']
        rot = Rotation.from_matrix(rot).as_matrix()
        center = np.array(extrinsic['value']['center'])
        t = - rot @ center
        rr.log_transform3d(f"world/camera/{filename}", rr.TranslationAndMat3(translation=t, matrix=rot), from_parent=True)
        rr.log_view_coordinates(f"world/camera/{filename}", xyz="RDF")  # X=Right, Y=Down, Z=Forward

        # Collect intrinsic data (show every image as if they were a pinhole camera)
        # We ignore camera distortion
        intrinsic = intrinsics[id_intrinsic]['value']['ptr_wrapper']['data']
        rr.log_pinhole(
            f"world/camera/{filename}",
                width=intrinsic['width'],
                height=intrinsic['height'],
                focal_length_px=intrinsic['focal_length'],
                principal_point_px=intrinsic['principal_point'],
        )
        if args.show_images and filename in args.show_images:
            rr.log_image_file(f"world/camera/{filename}", img_path=root_path + "/" + filename)

    # Collect 3D points
    print("Setup Structure")
    points = [obs['value']['X'] for obs in structure]
    tracks_ids = [obs['key'] for obs in structure]
    rr.log_points("world/points", points, keypoint_ids=tracks_ids)

    if args.show_keypoints or args.show_rays:
        print("Setup observations")

        # Collect 3d points observations
        points_per_image = defaultdict(list)
        XPoints_per_images = defaultdict(list)
        tracks_ids = defaultdict(list)
        for obs in tqdm(structure):
            id = obs['key']
            X = obs['value']['X']
            for x_it in obs['value']['observations']:
                x = x_it['value']['x']
                image_id = x_it['key']
                points_per_image[image_id].append(x)
                XPoints_per_images[image_id].append(X)

        for i in points_per_image.keys():
            value = views[i]['value']['ptr_wrapper']['data']
            filename = value['filename']
            id_pose = value['id_pose']
            id_view = value['id_view']

            if args.show_keypoints and filename in args.show_keypoints:
                # Print keypoint corresponding to the 3D points used in the SfM reconstruction
                rr.log_points(f"world/camera/{filename}", points_per_image[i])

            if args.show_rays and filename in args.show_rays:
                # Print camera rays
                if id_pose not in extrinsics_ids:
                    continue

                X = XPoints_per_images[id_view]
                extrinsic = [it for it in extrinsics if it['key'] == id_pose]
                extrinsic = extrinsic[0]
                center = np.array(extrinsic['value']['center'])
                centers = np.tile(center, (len(X), 1))
                interleaved_points = np.ravel(np.column_stack((centers, X)))
                interleaved_points = np.reshape(interleaved_points, [len(X)*2, 3])

                rr.log_line_segments(f"world/camera/rays/{filename}", interleaved_points, stroke_width=0.01)

    f.close()


def main() -> None:
    parser = ArgumentParser(description="Visualize the output of openMVG's sparse reconstruction.")
    parser.add_argument(
        "--sfm_data",
        type=Path,
        help="Which sfm_data.json scene to load",
        required=True
    )
    parser.add_argument(
        "--show_images",
        action="extend", nargs="+", type=str,
        help="Show the given image filename list (i.e --show-images 100_7100.JPG 100_7101.JPG"
    )
    parser.add_argument(
        "--show_rays",
        action="extend", nargs="+", type=str,
        help="Show the rays for the given image filename list (i.e --show-rays 100_7100.JPG 100_7101.JPG"
    )
    parser.add_argument(
        "--show_keypoints",
        action="extend", nargs="+", type=str,
        help="Show the keypoints used for the SfM scene for the given image filename list (i.e --show-keypoints 100_7100.JPG 100_7101.JPG"
    )
    rr.init("openMVG ReRun sfm_data viewer")
    rr.script_add_args(parser)
    args = parser.parse_args()

    rr.script_setup(args, "OpenMVG structure_from_motion ReRun Viewer")
    read_and_log_sparse_reconstruction(args)
    rr.script_teardown(args)


if __name__ == "__main__":
    main()
