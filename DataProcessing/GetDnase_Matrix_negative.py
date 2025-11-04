#生成Dnase矩阵 
#negative
from datetime import datetime
print(datetime.now())
#生成Dnase矩阵
import numpy as np
import pyBigWig
from tqdm import tqdm
import numpy as np
import datetime
import time  
def getbigwig(file, chrome, start, end):
    bw = pyBigWig.open(file)
    
    # 获取染色体的最大位置范围
    chrom_length = bw.chroms(chrome)
    
    # 确保 start 和 end 在有效范围内
    if start < 0:
        start = 0
    if end > chrom_length:
        end = chrom_length
    if start >= end:
        # return np.zeros((21, 10000))  # 返回全零的矩阵
        return np.zeros((21, 5000)) 
    # 获取实际范围内的数据
    sample = np.array(bw.values(chrome, start, end))
    
    # 将 NaN 值替换为 0
    sample[np.isnan(sample)] = 0
    
    # 期望的数组大小为 21 * 5000 = 105000
    # expected_size = 21 * 10000
    expected_size = 21 * 5000
    # 如果数据不够，使用 0 填充
    if sample.size < expected_size:
        sample = np.pad(sample, (0, expected_size - sample.size), 'constant', constant_values=0)
    
    # 如果数据超出，进行截取
    else:
        sample = sample[:expected_size]
    
    bw.close()
    
    return sample

def generateATAC(files_center_point, files_dnase, resou):
    """
    Generate training set
    :param coords: List of tuples containing coord bins
    :param width: Distance added to center. width=11 makes 21x21 windows
    :return: yields paired positive/negative samples for training
    """
    all_matrix = []
    with open(files_center_point) as f:
        n = 0
        for line in f:
            n += 1

            line = line.split()
            x, y = int(line[1]), int(line[2])
            chromname = line[0]
            print('当前染色体为', chromname, ' 进行了', n, "行", end='\r')
            ##########求x的atac
            window_x = getbigwig(files_dnase, chromname, (x - 10 - 1) * resou, (x + 10 + 1 - 1) * resou)
            window_x = window_x.reshape(21, 5000)

            # window_x = window_x.reshape(21, 10000)
            # window_x = window_x.reshape(21, 25000)
            window_x = [window_x.mean(axis=1)]  # 行
            #####求y的atac
            window_y = getbigwig(files_dnase, chromname, (y - 10 - 1) * resou, (y + 10 + 1 - 1) * resou)
            # window_y = window_y.reshape(21, 10000)
            window_y = window_y.reshape(21, 5000)
            window_y = [window_y.mean(axis=1)]
            window_atac = np.dot(np.transpose(window_x), window_y)
            all_matrix.append(window_atac)
    return all_matrix

if __name__ == '__main__':
    a = ['1','2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19','20', '21', '22']
    # a = ['20', '21', '22']
    for n in range(len(a)):                                               
        files_center_point = '/path/of/Negative_name_sort/negative_chr'+str(a[n])+'-5KB-sort.txt'

        files_dnase = "/path/of/merged.bigWig" #基于 DNase-seq 的 BigWig 格式文件
        resou = 5000
        all_matrix1 = generateATAC(files_center_point, files_dnase, resou)
        all_matrix_array = np.array(all_matrix1)
        np.save('基于 DNase-seq 的 BigWig 格式文件/DNase/dnase_chr' + str(a[n]) + '_0.npy', all_matrix_array)
        print('chr' + str(a[n]) + '已完成')