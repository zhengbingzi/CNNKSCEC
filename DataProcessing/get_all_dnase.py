# generate_dnase_matrix_from_chr.py (optimized)
import numpy as np
import pyBigWig
from datetime import datetime


def getbigwig(bw, chrom, start, end, window_size=21, res=5000):
    """从 bigWig 中提取信号并 reshape"""
    chrom_length = bw.chroms(chrom)

    # 边界修正
    start = max(0, start)
    end = min(end, chrom_length)
    if start >= end:
        return np.zeros((window_size, res), dtype=np.float32)

    # 获取数据 (直接 numpy 格式, 避免多余拷贝)
    sample = bw.values(chrom, start, end, numpy=True)
    sample[np.isnan(sample)] = 0

    expected_size = window_size * res
    if sample.size < expected_size:
        sample = np.pad(sample, (0, expected_size - sample.size), 'constant')
    else:
        sample = sample[:expected_size]

    return sample.reshape(window_size, res)


def generate_dnase_matrix(anchor_file, bigwig_file, res=5000):
    """
    输入整条染色体的中心点文件，生成对应的 DNase 矩阵
    每一行 -> 一个 21x21 矩阵
    """
    all_matrix = []
    with pyBigWig.open(bigwig_file) as bw, open(anchor_file) as f:
        for n, line in enumerate(f, 1):
            chrom, x, y = line.split()[:3]
            x, y = int(x), int(y)

            if n % 1000 == 0:  # 调整打印频率
                print(f"[{datetime.now()}] {chrom} 已处理 {n} 个点", end="\r")

            # x 窗口
            win_x = getbigwig(bw, chrom, (x - 11) * res, (x + 10) * res, res=res)
            win_x = win_x.mean(axis=1, keepdims=True)  # (21,1)

            # y 窗口
            win_y = getbigwig(bw, chrom, (y - 11) * res, (y + 10) * res, res=res)
            win_y = win_y.mean(axis=1, keepdims=True)  # (21,1)

            # 外积 (直接广播，避免 np.dot 的函数调用开销)
            window_atac = win_x * win_y.T  # (21,21)
            all_matrix.append(window_atac)

    return np.array(all_matrix, dtype=np.float32)


if __name__ == '__main__':
    chrom = "20"  # 举例 chr20
    anchor_file = f"/path/of/chr_all_centerpoint_txt/chr{chrom}_matrixsize21_centerpoints.txt"
    bigwig_file = "/path/of/merged.bigWig"
    save_file = f"/path/of/DNase-ALL/dnase_chr{chrom}.npy"

    res = 5000
    dnase_matrix = generate_dnase_matrix(anchor_file, bigwig_file, res=res)
    np.save(save_file, dnase_matrix)

    print(f"\nchr{chrom} 已完成，矩阵 shape = {dnase_matrix.shape}, 保存到 {save_file}")

