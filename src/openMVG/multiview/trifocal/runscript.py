import os
import subprocess
import sys
from itertools import combinations
import glob

OPENMVG_SFM_BIN = "/home/gabriel/openMVG-bin-official/Linux-x86_64-Release"
#path = "/home/gabriel/pavilion-multiview-3d-dataset/chair-sequence/images/midday/"
#path = "/home/gabriel/openMVG-bin-official/software/SfM/ImageDataset_SceauxCastle/images/"
#matches_dir = "/home/gabriel/openMVG-bin-official/software/SfM/ImageDataset_SceauxCastle/out/matches"
#reconstruction_dir = "/home/gabriel/openMVG-bin-official/software/SfM/ImageDataset_SceauxCastle/out/reconstruction_sequential"
#matches_dir = "pavillion/out/matches/"
#reconstruction_dir = "pavillion/out/reconstruction_sequential"
path = "/home/gabriel/flower-multiview-3d-dataset/Flower"
matches_dir = "/home/gabriel/flower-multiview-3d-dataset/Flower_out/matches"
reconstruction_dir = "/home/gabriel/flower-multiview-3d-dataset/Flower_out/reconstruction_sequential"
log = open("/home/gabriel/trifocalpy/flowerlogtmp","a")
#log = open("/home/gabriel/openMVG-bin-official/software/SfM/tee", "a")
file_list = glob.glob(os.path.join(os.getcwd(), path, "*.jpg"))
file_list = [i.strip(path)+"g" for i in file_list]
#file_list = glob.glob(os.path.join(os.getcwd(), path, "*.JPG"))
#file_list = [i.strip(path)+"G" for i in file_list]
corpus = combinations(file_list,2)
corpus2 = combinations(file_list,3)
listcomb = []
#print(file_list)
for j in corpus2:
    listcomb.append(j)
# filter repeateables
for i in listcomb:
    if i[0] == i[1] or i[1] == i[2] or i[0] == i[2]:
        print("equal")
        listcomb.remove(i)
print("triplets: ", len(listcomb))
print("Do Sequential/Incremental reconstruction 3-view")
log.write("triplet tests!\n")
final_list = [i for i in range(3,31)]
fails = 0
partials = 0
full = 0
for i in listcomb:
    # print(file_list.index(i[0]),file_list.index(i[1]),file_list.index(i[2]))
    precons = subprocess.Popen([os.path.join(OPENMVG_SFM_BIN,
        "openMVG_main_SfM"), "--sfm_engine", "INCREMENTAL", "--input_file",
        matches_dir+"/sfm_data.json", "--match_dir", matches_dir,
        "--output_dir", reconstruction_dir, "-x", i[0], "-y", i[1], "-z",
        i[2]], stderr=subprocess.STDOUT, stdout=subprocess.PIPE,shell=True)
    result = precons.stdout.read()
    result = result.decode()
    print("result:\n ",result)
    with open("/home/gabriel/openMVG-bin-official/software/SfM/tee","a") as out:
        if ("Robust estimation failed to compute calibrated trifocal tensor for this triplet" in result) or "failed" in result:
            out.write(f"{file_list.index(i[0])},{file_list.index(i[1])},{file_list.index(i[2])}"+
                    "-> failed: no coverage!\n")
            fails += 1
        else:
            if f"Camera calibrated: {len(file_list)} from {len(file_list)} input images" in result:
                out.write(f"{file_list.index(i[0])},{file_list.index(i[1])},{file_list.index(i[2])}-> "+
                        "Worked!\n")
                full += 1
            else:
                for j in final_list:
                    if f"Camera calibrated: {j} from {len(file_list)} input images" in result:
                        out.write(f"{file_list.index(i[0])},{file_list.index(i[1])},{file_list.index(i[2])}-> "
                                +f"Worked partially!#Camera calibrated: {j} from 31 input images\n")
                        partials += 1

print("Fail %: ", fails/len(listcomb) *100, "%")
print("Full %: ", full/len(listcomb) *100, "%")
print("Partial %: ", partials/len(listcomb) *100, "%")



