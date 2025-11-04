# utils.py

import torch
import numpy as np
import random
import datetime
from sklearn.metrics import precision_score, recall_score, accuracy_score, f1_score

# Evaluate the model on the validation set
def evaluate(model, data_loader, device):
    model.eval()
    all_preds, all_labels = [], []
    with torch.no_grad():
        for inputs, labels in data_loader:
            inputs, labels = inputs.to(device), labels.to(device)
            outputs, _ = model(inputs)  
            outputs = outputs.squeeze()
            preds = (outputs > 0.5).float()
            all_preds.append(preds.cpu().numpy().reshape(-1))
            all_labels.append(labels.cpu().numpy().reshape(-1))

    all_preds = np.concatenate(all_preds)
    all_labels = np.concatenate(all_labels)
    return all_preds, all_labels

# KL divergence function
def KL(outputs, targets):
    log_softmax_outputs = F.log_softmax(outputs / 4.0, dim=1)
    softmax_targets = F.softmax(targets / 4.0, dim=1)
    return -(log_softmax_outputs * softmax_targets).sum(dim=1).mean()

# Save the best model based on F1 score
def save_best_model(model, val_f1, best_f1, best_model_path):
    if val_f1 > best_f1:
        best_f1 = val_f1
        torch.save(model.state_dict(), best_model_path)
        print(f"Best model saved with F1={best_f1:.4f}")
    return best_f1

# Set random seeds for reproducibility
def set_random_seeds(seed=42):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)  # If using multiple GPUs
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
