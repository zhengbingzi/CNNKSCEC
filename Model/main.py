# main.py

from train import train
from utils import set_random_seeds
from torch.utils.data import DataLoader
from model import DeepCNN

def main():
    # Set random seed
    set_random_seeds(seed=42)

    # Initialize model, criterion, optimizer, etc.
    model = DeepCNN()
    criterion = nn.BCELoss()
    optimizer = optim.Adam(model.parameters(), lr=0.001)
    
    # Define DataLoader, train_loader, val_loader
    train_loader = DataLoader(...)  # Define your DataLoader here
    val_loader = DataLoader(...)    # Define your DataLoader here
    
    # Train model
    best_model_path = "best_model.pth"
    num_epochs = 40
    model = train(model, train_loader, val_loader, criterion, optimizer, num_epochs, device, best_model_path)

if __name__ == "__main__":
    main()
