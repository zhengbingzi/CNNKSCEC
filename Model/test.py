# test.py

import torch
import numpy as np
from torch.utils.data import DataLoader, TensorDataset
from sklearn.metrics import precision_score, recall_score, accuracy_score, f1_score, precision_recall_curve, auc
import warnings
from model import DeepCNN 

# 忽略警告信息
warnings.filterwarnings("ignore", category=FutureWarning)

# 设置测试数据路径和模型路径
test_data_path = "/mnt/sdf/zhengbingzi/Othertools/CNNKSCEC/5kb/newDNase-HiC1/Train-Val-Test/chr20-22_test.npz"
best_model_path = "/mnt/sdf/zhengbingzi/CNNKSCEC/Model/ModelWeight/model_cnnklscc+ECA+CBAM+p0.2t4l0.003_student12.pth"

# 加载测试数据
def load_test_data(test_data_path):
    Data = np.load(test_data_path, allow_pickle=True)
    seqs = []
    labels = []
    atac_info = []
    seqs.extend(Data['hic_data'].astype(np.float32))
    atac_info.extend(Data['dnase_data'].astype(np.float32))
    labels.extend(Data['labels'].astype(np.float32))
    seqs = np.array(seqs)
    labels = np.array(labels)
    atac_info = np.array(atac_info)

    # 预处理数据
    for i in range(len(atac_info)):
        atac_info[i] = np.log10(1 + atac_info[i]*10)
        atac_info[i] = atac_info[i] / np.max(atac_info[i] + 1)

        seqs[i] = np.log10(1 + seqs[i] * 10)
        seqs[i] = seqs[i] / np.max(seqs[i] + 1)

    # 拼接序列数据
    inputs = np.stack([seqs, atac_info], axis=1)  # shape: (340363, 2, 21, 21)

    # 转换为 PyTorch 张量
    inputs_tensor = torch.tensor(inputs, dtype=torch.float32)
    labels_tensor = torch.tensor(labels, dtype=torch.float32)

    # 创建数据集和数据加载器
    dataset = TensorDataset(inputs_tensor, labels_tensor)
    test_loader = DataLoader(dataset, batch_size=50, shuffle=False)
    
    return test_loader


# 加载模型
def load_model(model, best_model_path):
    model.load_state_dict(torch.load(best_model_path))
    model.eval()  # 设置模型为评估模式
    return model


# 预测函数
def predict(model, data_loader, device):
    model.eval()
    all_preds, all_labels = [], []
    with torch.no_grad():  # 禁用梯度计算
        for inputs, labels in data_loader:
            inputs = inputs.to(device)
            outputs, _ = model(inputs)  # Assuming model outputs a single value per input
            outputs = outputs.squeeze()
            preds = (outputs > 0.5).float()  # 二分类预测，阈值0.5
            all_preds.append(preds.cpu().numpy())  # 转为 numpy 数组
            all_labels.append(labels.cpu().numpy())  # 真实标签
    all_preds = np.concatenate(all_preds)
    all_labels = np.concatenate(all_labels)
    return all_preds, all_labels


# 计算并输出评估指标
def evaluate_model(true_labels, predictions):
    precision = precision_score(true_labels, predictions)
    recall = recall_score(true_labels, predictions)
    accuracy = accuracy_score(true_labels, predictions)
    f1 = f1_score(true_labels, predictions)

    precision_vals, recall_vals, _ = precision_recall_curve(true_labels, predictions)
    pr_auc = auc(recall_vals, precision_vals)

    print(f"Precision: {precision:.4f}")
    print(f"Recall: {recall:.4f}")
    print(f"Accuracy: {accuracy:.4f}")
    print(f"F1 Score: {f1:.4f}")
    print(f"PR AUC: {pr_auc:.4f}")


# 主函数
def main():
    # 定义设备
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    
    # 加载数据
    test_loader = load_test_data(test_data_path)
    
    # 加载最佳模型
    model = DeepCNN()  # Ensure your model is defined or imported from another file
    model = load_model(model, best_model_path)
    model.to(device)

    # 进行预测
    predictions, true_labels = predict(model, test_loader, device)

    # 评估模型
    evaluate_model(true_labels, predictions)


if __name__ == "__main__":
    main()
