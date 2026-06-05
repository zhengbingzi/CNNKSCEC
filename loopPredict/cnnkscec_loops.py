#predicte
import warnings
# 忽略 FutureWarning
import warnings
import numpy as np
import torch
from torch.utils.data import DataLoader, TensorDataset
import numpy as np
import torch
from torch.utils.data import DataLoader, TensorDataset
from model import DeepCNN
from utils import KL
warnings.filterwarnings("ignore", category=FutureWarning)
# 加载数据
Data = np.load("chr20-22_predict.npz", allow_pickle=True)



# 提取数据
post = Data['position_data']  
seqs = Data['hic_data'].astype(np.float32)
atac_info = Data['dnase_data'].astype(np.float32)

# 对 seqs 和 atac_info 进行预处理
for i in range(len(atac_info)):
    atac_info[i] = np.log10(1 + atac_info[i]*10)
    atac_info[i] = atac_info[i]/np.max(atac_info[i]+1)

    seqs[i] = np.log10(1 + seqs[i] * 10)
    seqs[i] = seqs[i]/np.max(seqs[i]+1)

# 提取 infrancy（Hi-C 矩阵的中心点值）
infrancy = seqs[:, 10, 10]

# 生成随机索引
shuffle_indices = np.random.permutation(len(seqs))

# 按相同顺序打乱 post, seqs, atac_info 和 infrancy
post_shuffled = post[shuffle_indices]
seqs_shuffled = seqs[shuffle_indices]
atac_info_shuffled = atac_info[shuffle_indices]
infrancy_shuffled = infrancy[shuffle_indices]

# 堆叠 seqs 和 atac_info 作为输入
inputs_shuffled = np.stack([seqs_shuffled, atac_info_shuffled], axis=1)


inputs_tensor = torch.tensor(inputs_shuffled, dtype=torch.float32)

# 创建 dataset 和 dataloader
dataset = TensorDataset(inputs_tensor)
batch_size = 50  
test_loader = DataLoader(dataset, batch_size=batch_size, shuffle=False)
# 加载最佳模型权重
model = DeepCNN()
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model.to(device)
best_model_path =  "model_cnnklscc+ECA+CBAM+p0.2t4l0.003_student12.pth"#保存好的模型参数，同上
model.load_state_dict(torch.load(best_model_path))
model.eval()

# 预测函数
def predict(model, data_loader):
    model.eval()
    all_preds = []
    with torch.no_grad():  # 禁用梯度计算
        for inputs in data_loader:
            inputs = inputs[0].to(device)  # 仅使用输入张量
            outputs, _ = model(inputs)
            outputs = outputs.squeeze()
            #preds = (outputs > 0.5).float()  # 二分类预测，阈值0.5
            all_preds.append(outputs.cpu().numpy())  # 转为 numpy 数组
    all_preds = np.concatenate(all_preds)
    return all_preds
# 执行预测
predictions = predict(model, test_loader)

# 打乱后的 position, predictions 和 infrancy 组合
results = np.column_stack((post_shuffled, predictions, infrancy_shuffled))

# 提取并格式化为 BEDPE 格式
chrom1 = results[:, 0].astype(str)  
chrom2 = results[:, 0].astype(str)  
start1 = (results[:, 1].astype(int)-1) * 5000
end1 = results[:, 1].astype(int)* 5000  
start2 = (results[:, 2].astype(int)-1) * 5000  
end2 = results[:, 2].astype(int) * 5000 

# 合并 BEDPE 所需的6列，加上预测分数和 Hi-C 矩阵中心点值
bedpe_data = np.column_stack((chrom1, start1, end1, chrom2, start2, end2, predictions.astype(float), infrancy_shuffled.astype(float)))
print(bedpe_data.shape)
# 保存为 BEDPE 文件
output_bedpe_file = "predictions_chr20-22_cnnklsccecabm+p0.2t4l0.003_student12.bedpe"#保存路径
np.savetxt(output_bedpe_file, bedpe_data, delimiter='\t', fmt='%s')
print(f"Saved !")