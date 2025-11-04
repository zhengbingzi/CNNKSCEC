# train.py

import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import sys
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset, random_split
from sklearn.metrics import precision_score, recall_score, accuracy_score, f1_score, precision_recall_curve, auc
import warnings
import random
from tqdm import tqdm
from model import DeepCNN
from utils import KL,evaluate, save_best_model, set_random_seeds  # Import helper functions

#train-val
warnings.filterwarnings("ignore", category=FutureWarning)
Data = np.load("/mnt/sdf/zhengbingzi/Othertools/CNNKSCEC/5kb/newDNase-HiC1/Train-Val-Test/chr1-19_trainval.npz", allow_pickle=True)
seqs = []
labels = []
atac_info = []
seqs.extend(Data['hic_data'].astype(np.float32))
atac_info.extend(Data['dnase_data'].astype(np.float32))
labels.extend(Data['labels'].astype(np.float32))
seqs = np.array(seqs)
labels = np.array(labels)
atac_info = np.array(atac_info)
#labels = labels.reshape((labels.shape[0], 1))
print(seqs.shape)
print(atac_info.shape)
print(labels.shape)
# 生成随机索引
indices = np.random.permutation(len(labels))
# 根据随机索引同时打乱 hic_matrices, dnase_data 和 labels
seqs = seqs[indices]
atac_info = atac_info[indices]
labels = labels[indices]
for i in range(len(atac_info)):
    atac_info[i] = np.log10(1 + atac_info[i]*10)
    atac_info[i] = atac_info[i]/np.max(atac_info[i]+1)

    seqs[i] = np.log10(1 + seqs[i] * 10)
    seqs[i] = seqs[i]/np.max(seqs[i]+1)
inputs = np.stack([seqs, atac_info], axis=1)  

inputs_tensor = torch.tensor(inputs, dtype=torch.float32)
labels_tensor = torch.tensor(labels, dtype=torch.float32)


dataset = TensorDataset(inputs_tensor, labels_tensor)

# Train-validation split (4:1)
train_size = int(0.8 * len(dataset))
val_size = len(dataset) - train_size
train_dataset, val_dataset = random_split(dataset, [train_size, val_size])

# Set batch size to 50
batch_size = 50

train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
# Define model
model = DeepCNN()  
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
#device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
model.to(device)
criterion = nn.BCELoss() 
optimizer = optim.Adam(model.parameters(), lr=0.001)
# 设置随机种子
seed = 42
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)
torch.cuda.manual_seed_all(seed)  # 如果使用多 GPU
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
# Function to evaluate model
def evaluate(model, data_loader):
    model.eval()
    all_preds, all_labels = [], []
    with torch.no_grad():
        for inputs, labels in data_loader:
            inputs, labels = inputs.to(device), labels.to(device)
            outputs, _ = model(inputs)  # Assuming model outputs a single value per input
            outputs = outputs.squeeze()
            preds = (outputs > 0.5).float()
            all_preds.append(preds.cpu().numpy().reshape(-1))
            all_labels.append(labels.cpu().numpy().reshape(-1))
    
    all_preds = np.concatenate(all_preds)
    all_labels = np.concatenate(all_labels)
    return all_preds, all_labels

# Training loop
best_f1 = 0
num_epochs = 40
best_model_path = "/mnt/sdf/zhengbingzi/CNNKSCEC/Model/ModelWeight/model_cnnklscc+ECA+CBAM+p0.2t4l0.003_student12(1).pth"#保存参数的路径

for epoch in range(num_epochs):
    model.train()
    running_loss = 0.0
    print(f"Model is on: {next(model.parameters()).device}")
    # Add tqdm to the training loop
    for inputs, labels in tqdm(train_loader, desc=f"Epoch {epoch + 1}", leave=True):
        inputs, labels = inputs.to(device), labels.to(device)
        # Forward pass
        optimizer.zero_grad()
        outputs, fea = model(inputs)  # Assuming model outputs a single value per input
        outputs = outputs.squeeze()
        
        # Compute loss
        teacher = fea[0].detach()
        loss = criterion(outputs, labels) + 0.003*(KL(fea[1],teacher) + KL(fea[2],teacher))
        
        
        running_loss += loss.item()
        # Backward pass and optimization
        loss.backward()
        optimizer.step()
    
    epoch_loss = running_loss / len(train_loader)
    # Validation step
    val_preds, val_labels = evaluate(model, val_loader)
    
    # Calculate metrics
    val_f1 = f1_score(val_labels, val_preds)
    val_precision = precision_score(val_labels, val_preds)
    val_recall = recall_score(val_labels, val_preds)
    val_accuracy = accuracy_score(val_labels, val_preds)

    precision_vals, recall_vals, _ = precision_recall_curve(val_labels, val_preds)
    pr_auc = auc(recall_vals, precision_vals)
    
    print(f"Epoch {epoch + 1}: Loss={epoch_loss:.4f}, F1={val_f1:.4f}, Precision={val_precision:.4f}, Recall={val_recall:.4f}, Accuracy={val_accuracy:.4f}, PR AUC={pr_auc:.4f}")
    
    # Save the best model based on F1 score
    if val_f1 > best_f1:
        best_f1 = val_f1
        torch.save(model.state_dict(), best_model_path)
        print(f"Best model saved with F1={best_f1:.4f}")

# Load the best model
model.load_state_dict(torch.load(best_model_path,map_location=device))
model.to(device)
model.eval()