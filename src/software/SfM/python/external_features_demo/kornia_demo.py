from argparse import ArgumentParser
from itertools import combinations
import json
import kornia as K
import numpy as np
import os
from pyvips import Image
import torch
import torchvision.transforms as transforms
import threading
from tqdm import tqdm

# // U T I L S ///////////////////////////////////////////////////////
def loadJSON():
  with open(args.input) as file:
    sfm_data = json.load(file)
  view_ids = {view['key']:view['value']['ptr_wrapper']['data']['filename'] for view in sfm_data['views']}
  image_paths = [os.path.join(sfm_data['root_path'], view['value']['ptr_wrapper']['data']['filename']) for view in sfm_data['views']]
  return view_ids, image_paths

def saveFeatures(keypoints, descriptors, scores, basename):
  torch.save(keypoints, os.path.join(args.matches, f'{basename}.keyp.pt'))
  torch.save(descriptors, os.path.join(args.matches, f'{basename}.desc.pt'))
  torch.save(scores, os.path.join(args.matches, f'{basename}.scor.pt'))

def loadFeatures(filename):
  keyp_filename = os.path.join(args.matches, f'{os.path.splitext(filename)[0]}.keyp.pt')
  desc_filename = os.path.join(args.matches, f'{os.path.splitext(filename)[0]}.desc.pt')
  scor_filename = os.path.join(args.matches, f'{os.path.splitext(filename)[0]}.scor.pt')
  keypoints = torch.load(os.path.join(args.matches, keyp_filename)).to(device)
  descriptors = torch.load(os.path.join(args.matches, desc_filename)).to(device)
  scores = torch.load(os.path.join(args.matches, scor_filename)).to(device)
  return keypoints, descriptors, scores

def saveFeaturesOpenMVG(basename, keypoints):
  with open(os.path.join(args.matches, f'{basename}.feat'), 'w') as feat:
    for x, y in keypoints.numpy():
      feat.write(f'{x} {y} 1.0 0.0\n')

def saveDescriptorsOpenMVG(basename, descriptors):
  with open(os.path.join(args.matches, f'{basename}.desc'), 'wb') as desc:
    desc.write(len(descriptors).to_bytes(8, byteorder='little'))
    desc.write(((descriptors.numpy() + 1) * 0.5 * 255).round(0).astype(np.ubyte).tobytes())

def saveMatchesOpenMVG(matches):
  with open(args.output, 'wb') as bin:
    bin.write((1).to_bytes(1, byteorder='little'))
    bin.write(len(matches).to_bytes(8, byteorder='little'))
    for index1, index2, idxs in matches:
      bin.write(index1.tobytes())
      bin.write(index2.tobytes())
      bin.write(len(idxs).to_bytes(8, byteorder='little'))
      bin.write(idxs.tobytes())

# // F E A T U R E ///////////////////////////////////////////////////
def featureExtraction():
  print('Extracting DISK features...')
  for image_path in tqdm(image_paths):
    img = Image.new_from_file(image_path, access='sequential')
    basename = os.path.splitext(os.path.basename(image_path))[0]

    if img.width % 2 != 0 or img.height % 2 != 0:
      img = img.crop(0, 0, img.width if img.width % 2 == 0 else img.width - 1, img.height if img.height % 2 == 0 else img.height - 1)

    max_res = args.max_resolution
    img_max, img_ratio = max(img.width, img.height), img.width / img.height

    ratio = 0
    while not ratio:
      scale = max_res / img_max
      scaled_width, scaled_height = round(img.width * scale), round(img.height * scale)
      if img_ratio == scaled_width / scaled_height:
        ratio = 1
      else:
        max_res -= 1

    img = transforms.ToTensor()(img.resize(scale, kernel='linear').numpy())[None, ...].to(device)

    features = disk(img, n=args.max_features, window_size=args.window_size, score_threshold=args.score_threshold, pad_if_not_divisible=True)[0].to('cpu')
    keypoints = torch.div(features.keypoints, scale)

    threading.Thread(target=lambda: saveFeatures(keypoints, features.descriptors, features.detection_scores, basename)).start()

    threading.Thread(target=lambda: saveFeaturesOpenMVG(basename, keypoints)).start()
    threading.Thread(target=lambda: saveDescriptorsOpenMVG(basename, features.descriptors)).start()

# // M A T C H I N G /////////////////////////////////////////////////
def featureMatching():
  print('Matching DISK features with LightGlue...')
  putative_matches = []
  for image1_index, image2_index in tqdm((np.loadtxt(args.pair_list, dtype=np.int32) if args.pair_list != None else np.asarray([*combinations(view_ids, 2)], dtype=np.int32))):
    keyp1, desc1, scor1 = loadFeatures(view_ids[image1_index])
    keyp2, desc2, scor2 = loadFeatures(view_ids[image2_index])

    lafs1 = K.feature.laf_from_center_scale_ori(keyp1[None], 96 * torch.ones(1, len(keyp1), 1, 1, device=device))
    lafs2 = K.feature.laf_from_center_scale_ori(keyp2[None], 96 * torch.ones(1, len(keyp2), 1, 1, device=device))

    dists, idxs = lightglue(desc1, desc2, lafs1, lafs2)

    putative_matches.append([image1_index, image2_index, idxs.cpu().numpy().astype(np.int32)])

  print('Saving putative matches...')
  saveMatchesOpenMVG(putative_matches)

if __name__ == '__main__':

  parser = ArgumentParser()
  parser.add_argument('-i', '--input', type=str, required=True, help='Path to the sfm_data file')
  parser.add_argument('--max_resolution', type=int, default=1024, help='Max image resolution')
  parser.add_argument('--max_features', type=int, default=4096, help='Max number of features to extract')
  parser.add_argument('-m', '--matches', type=str, required=True, help='Path to the matches directory')
  parser.add_argument('-p', '--pair_list', type=str, help='Path to the pair file')
  parser.add_argument('-o', '--output', type=str, help='Path to the output matches file')
  parser.add_argument('--preset', choices=['BOTH','EXTRACT','MATCH'], default='BOTH', help='Preset to run')
  parser.add_argument('--force_cpu', action='store_true', help='Force device to CPU')
  # DISK
  parser.add_argument('--window_size', type=int, default=5, help='DISK Non Maximum Suppression (NMS) radius (Must be odd)')
  parser.add_argument('--score_threshold', type=float, default=0.01, help='DISK keypoint detector confidence threshold')
  # LightGlue
  parser.add_argument('--n_layers', type=int, default=9, help='LightGlue layers')
  parser.add_argument('--depth_confidence', type=float, default=0.95, help='LightGlue early stopping (-1 - disable)')
  parser.add_argument('--width_confidence', type=float, default=0.99, help='LightGlue point pruning (-1 - disable)')
  parser.add_argument('--filter_threshold', type=float, default=0.99, help='LightGlue match threshold')
  args = parser.parse_args()

  view_ids, image_paths = loadJSON()
  if args.output == None:
    args.output = os.path.join(args.matches, 'matches.putative.bin')

  device = torch.device('cpu') if args.force_cpu else K.utils.get_cuda_device_if_available()

  config = {
    'lightglue': {
      'n_layers': args.n_layers,
      'depth_confidence': args.depth_confidence,
      'width_confidence': args.width_confidence,
      'filter_threshold': args.filter_threshold
    }
  }

  disk = K.feature.DISK().from_pretrained('depth').to(device)
  print('Loaded DISK model')
  lightglue = K.feature.LightGlueMatcher(params=config['lightglue']).to(device)

  with torch.inference_mode():
    if args.preset == 'EXTRACT' or args.preset == 'BOTH':
      featureExtraction()
    if args.preset == 'MATCH' or args.preset == 'BOTH':
      featureMatching()
