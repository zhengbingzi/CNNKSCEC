#合并Dnase+Hic的正负样本 训练验证-测试
#训练验证
import numpy as np
import os
import re
def extract_chromosome_number(filename):
    match = re.search(r'chr(\d+)', filename)
    if match:
        return int(match.group(1))
    return None

def process_npz_files(folder_path, label, chromosomes):
    data_list = []
    for filename in os.listdir(folder_path):
        if filename.endswith('.npz'):
            chr_num = extract_chromosome_number(filename)
            if chr_num in chromosomes:
                file_path = os.path.join(folder_path, filename)
                with np.load(file_path) as data:
                    dnase = data['dnase_data']  # Use actual key
                    hic = data['hic_data']  # Use actual key
                    label_col = data['labels']  # Use actual key

                    if len(label_col) == len(dnase):  # Ensuring consistency
                        for i in range(len(label_col)):
                            data_list.append((dnase[i], hic[i], label_col[i]))
    return data_list

def save_combined_data(positive_folder, negative_folder, output_path, chromosomes):
    pos_data = process_npz_files(positive_folder, 1, chromosomes)
    neg_data = process_npz_files(negative_folder, 0, chromosomes)

    all_data = pos_data + neg_data
    np.savez(output_path, dnase_data=np.array([x[0] for x in all_data]), 
                        hic_data=np.array([x[1] for x in all_data]),
                        labels=np.array([x[2] for x in all_data]))

# Folder paths
positive_folder = '/path/of/DNase-HiC/Positive/'
negative_folder = '/path/of/DNase-HiC/Negative/'
output_path = '/path/of/DNase-HiC/Train-Val-Test/chr1-19_trainval.npz'
chromosomes = list(range(1, 20))
save_combined_data(positive_folder, negative_folder, output_path, chromosomes)
print("trainval success!")
############################################################################################
#合并Dnase+Hic的正负样本 训练验证测试
#测试
positive_folder = '/path/of/DNase-HiC/Positive/'
negative_folder = '/path/of/DNase-HiC/Negative/'
output_path = '/path/of/DNase-HiC/Train-Val-Test/chr20-22_test.npz'
chromosomes = list(range(20, 23))
save_combined_data(positive_folder, negative_folder, output_path, chromosomes)
print("test success!")