import os
import glob
import re
from collections import defaultdict

# === 设置输入目录（请根据实际修改）===
input_dir = '/path/of/Negative_name_sort/'  # ← 修改为你的实际文件夹路径

os.chdir(input_dir)

# 匹配所有拆分文件
all_files = glob.glob("negative_chr*-5KB-sort_*.txt")

# 按染色体分类文件，例如 chr1、chr2、chrX 等
chr_files = defaultdict(list)
pattern = re.compile(r'negative_(chr[\w\d]+)-5KB-sort_(\d+)\.txt')

for f in all_files:
    match = pattern.search(f)
    if match:
        chr_name = match.group(1)
        chr_files[chr_name].append(f)

# 遍历每条染色体，按序号排序并合并文件
for chr_name, files in chr_files.items():
    # 排序按文件末尾数字顺序
    sorted_files = sorted(files, key=lambda x: int(pattern.search(x).group(2)))
    output_filename = f'negative_{chr_name}-5KB-sort.txt'

    with open(output_filename, 'w') as outfile:
        for f in sorted_files:
            with open(f, 'r') as infile:
                outfile.write(infile.read())

    print(f"合并完成: {output_filename}")
