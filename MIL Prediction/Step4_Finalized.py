# Import Package
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import h5py
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score, accuracy_score, roc_auc_score, classification_report
from sklearn.metrics import accuracy_score, roc_auc_score, classification_report
import warnings
warnings.filterwarnings("ignore")

# Base Setting
# Directory Path
H5AD_PATH   = r"C:\Users\p3ngu\OneDrive\桌面\Stats170AB\fulldata.h5ad"
OUTPUT_DIR  = "C:/Users/p3ngu/OneDrive/桌面/Stats170AB/new_Step4/step04_outputs/"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Column name
DONOR_COL    = "donor_id"
AGE_COL      = "age_or_mean_of_age_range"
CELLTYPE_COL = "ann_level_3"
EMBEDDING_KEY = "X_scANVI"

# Machine Learning Setting
MAX_CELLS_PER_DONOR = 300
MIN_CELLS_PER_DONOR = 30
LATENT_DIM = 30
MIN_AGE = 5
DISEASE_COL = "lung_condition"
HEALTHY_LABEL = "healthy"

ATTN_HIDDEN  = 64
MLP_HIDDEN = 32
BATCH_SIZE = 8
EPOCHS = 150
PATIENCE = 20
LR = 1e-3
TEST_SIZE = 0.2
SEED = 42

torch.manual_seed(SEED)
np.random.seed(SEED)
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Torch Set")


print("\n")
# 1. Read Data
print("1. Reading data")

def read_categorical(h5, field):
    data = h5[field]
    if isinstance(data, h5py.Dataset):
        return data[:]
    categories = data["categories"][:].astype(str)
    codes = data["codes"][:]
    result = np.where(codes >= 0, categories[codes], "NA")
    return result

with h5py.File(H5AD_PATH, "r") as f:
    obs = f["obs"]

    donor_ids = read_categorical(obs, DONOR_COL)
    age_raw = obs[AGE_COL][:]
    age = age_raw.astype(float)
    disease_labels = read_categorical(obs, DISEASE_COL)
    cell_types = read_categorical(obs,CELLTYPE_COL)

    available_obsm = list(f["obsm"].keys())
    preferred = ["X_scANVI", "X_scanvi_emb", "X_scvi", "X_pca"]
    EMBEDDING_KEY_ACTUAL = next(
        (k for k in preferred if k in available_obsm),
        next((k for k in available_obsm if "umap" not in k.lower()), available_obsm[0])
    )
    emb_raw = f["obsm"][EMBEDDING_KEY_ACTUAL][:]
    n_obs   = len(donor_ids)
    if emb_raw.shape[0] != n_obs and emb_raw.shape[1] == n_obs:
        emb_raw = emb_raw.T
    emb = emb_raw
    LATENT_DIM = emb.shape[1]
    print(f"  Embedding dim : {LATENT_DIM}")
    n_cells = len(donor_ids)


    # Checking
    print(f"  Total cells: {n_cells:,}")
    print(f"  Embedding shape: {emb.shape}")
    print(f"  Unique donors: {len(np.unique(donor_ids)):,}")



print("\n")
# 2. Build Donor Bags
print("2. Build Donor Bags")
rng = np.random.default_rng(SEED)
bags = {}

unique_donors = np.unique(donor_ids)
for donor in unique_donors:
    mask = donor_ids == donor
    idx = np.where(mask)[0]
    donor_age = age[idx[0]]
    if np.isnan(donor_age):
        continue
    if donor_age < MIN_AGE:
        continue

    # Build disease status
    donor_condition = disease_labels[idx[0]]
    disease_label   = 0.0 if str(donor_condition).lower() == HEALTHY_LABEL else 1.0
    if len(idx) < MIN_CELLS_PER_DONOR:
        continue
    if len(idx) > MAX_CELLS_PER_DONOR:
        idx = rng.choice(idx, size=MAX_CELLS_PER_DONOR, replace=False)

    bags[donor] = {
        "emb": emb[idx],
        "age": float(donor_age),
        "disease_label": disease_label,
        "condition": str(donor_condition),
        "cell_types": cell_types[idx].tolist(),
    }

# Donor Bags built:
donor_list = list(bags.keys())
ages_all = np.array([bags[d]["age"] for d in donor_list])
disease_all = np.array([bags[d]["disease_label"] for d in donor_list])

healthy_num  = int((disease_all == 0).sum())
disease_num = int((disease_all == 1).sum())
# Checking
print(f"  Valid donors: {len(donor_list)}")
print(f"  Median cells/donor: {np.median([len(bags[d]['emb']) for d in donor_list]):.0f}")
print(f"  Healthy numbers: {healthy_num}, Diseased numbers: {disease_num}")



print("\n")
# 3. Split Training Data
print("3. Split Training Data")

train_ids, test_ids = train_test_split(
    donor_list, test_size=TEST_SIZE, random_state=SEED
)
print(f"  Train donors : {len(train_ids)}")
print(f"  Test donors  : {len(test_ids)}")



print("\n")
# 4. Building Class of Bags
print("4. Building Class of Bags")

class DonorBagDataset(Dataset):
    def __init__(self, donor_ids, bags):
        self.donor_ids = donor_ids
        self.bags      = bags

    def __len__(self):
        return len(self.donor_ids)

    def __getitem__(self, idx):
        d   = self.donor_ids[idx]
        bag = self.bags[d]
        X   = torch.tensor(bag["emb"],      dtype=torch.float32)
        y   = torch.tensor(bag["disease_label"], dtype=torch.float32)
        return X, y, d


def collate_bags(batch):
    Xs, ys, ds = zip(*batch)
    max_n  = max(x.size(0) for x in Xs)
    n_feat = Xs[0].size(1)
    X_pad  = torch.zeros(len(Xs), max_n, n_feat)
    masks  = torch.zeros(len(Xs), max_n, dtype=torch.bool)
    for i, x in enumerate(Xs):
        n = x.size(0)
        X_pad[i, :n] = x
        masks[i, :n] = True
    return X_pad, masks, torch.stack(ys), list(ds)


class GatedAttentionMIL(nn.Module):
    def __init__(self, in_dim, attn_hidden, mlp_hidden):
        super().__init__()
        self.V = nn.Linear(in_dim, attn_hidden)
        self.U = nn.Linear(in_dim, attn_hidden)
        self.w = nn.Linear(attn_hidden, 1)
        self.mlp = nn.Sequential(
            nn.Linear(in_dim, mlp_hidden),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(mlp_hidden, 1),
        )

    def forward(self, X, mask):
        # (B, N, H)
        A = torch.tanh(self.V(X)) * torch.sigmoid(self.U(X))
        # (B, N)
        A = self.w(A).squeeze(-1)
        A = A.masked_fill(~mask, float("-inf"))
        # (B, N)
        alpha = torch.softmax(A, dim=1)
        # (B, in_dim)
        z_bag = (alpha.unsqueeze(-1) * X).sum(dim=1)
        # (B,) raw logit
        logit = self.mlp(z_bag).squeeze(-1)
        return logit, alpha

print("  Finish building class")



print("\n")
# 5. Training Model
print("5. Training Model")

train_ds = DonorBagDataset(train_ids, bags)
test_ds  = DonorBagDataset(test_ids,  bags)
train_dl = DataLoader(train_ds, batch_size=BATCH_SIZE, shuffle=True, collate_fn=collate_bags)
test_dl  = DataLoader(test_ds,  batch_size=BATCH_SIZE, shuffle=False, collate_fn=collate_bags)

model     = GatedAttentionMIL(LATENT_DIM, ATTN_HIDDEN, MLP_HIDDEN).to(DEVICE)
optimizer = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=1e-3)
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, patience=8, factor=0.5)
criterion = nn.BCEWithLogitsLoss()

best_val   = float("inf")
best_state = None
no_improve = 0
history    = {"train": [], "val": []}

for epoch in range(1, EPOCHS + 1):
    model.train()
    t_losses = []
    for X_b, mask_b, y_b, _ in train_dl:
        X_b, mask_b, y_b = X_b.to(DEVICE), mask_b.to(DEVICE), y_b.to(DEVICE)
        optimizer.zero_grad()
        y_hat, _ = model(X_b, mask_b)
        loss = criterion(y_hat, y_b)
        loss.backward()
        nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        t_losses.append(loss.item())

    model.eval()
    v_losses = []
    with torch.no_grad():
        for X_b, mask_b, y_b, _ in test_dl:
            X_b, mask_b, y_b = X_b.to(DEVICE), mask_b.to(DEVICE), y_b.to(DEVICE)
            y_hat, _ = model(X_b, mask_b)
            v_losses.append(criterion(y_hat, y_b).item())

    t_loss = np.mean(t_losses)
    v_loss = np.mean(v_losses)
    scheduler.step(v_loss)
    history["train"].append(t_loss)
    history["val"].append(v_loss)

    if epoch % 25 == 0 or epoch == 1:
        print(f"  Epoch {epoch:3d}  train={t_loss:.4f}  val={v_loss:.4f}")

    if v_loss < best_val - 1e-5:
        best_val   = v_loss
        best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
        no_improve = 0
    else:
        no_improve += 1
        if no_improve >= PATIENCE:
            print(f"  Early stopping at epoch {epoch}")
            break

model.load_state_dict(best_state)
print(f"  Best normalized MSE: {best_val:.3f}")



print("\n")
# 6. Evaluation
print("6. Evaluating")

model.eval()
pred_rows = []
attn_rows = []

with torch.no_grad():
    for X_b, mask_b, y_b, donors_b in test_dl:
        X_b, mask_b = X_b.to(DEVICE), mask_b.to(DEVICE)
        y_hat, alpha = model(X_b, mask_b)

        probs = torch.sigmoid(y_hat)
        for i, donor in enumerate(donors_b):
            n_real = mask_b[i].sum().item()
            a_vals = alpha[i, :n_real].cpu().numpy()
            cts = bags[donor]["cell_types"][:n_real]

            true_label = bags[donor]["disease_label"]
            pred_prob  = float(probs[i].cpu())
            pred_label = 1.0 if pred_prob >= 0.5 else 0.0
            condition  = bags[donor]["condition"]

            pred_rows.append({
                "donor_id": donor,
                "condition": condition,
                "true_label": int(true_label),
                "pred_prob": round(pred_prob, 4),
                "pred_label": int(pred_label),
                "correct": true_label == pred_label,
                "age": bags[donor]["age"],
            })

            for ct, av in zip(cts, a_vals):
                attn_rows.append({
                    "donor_id": donor,
                    "cell_type": ct,
                    "alpha": float(av),
                    "condition": condition,
                    "true_label": int(true_label),
                })

pred_df = pd.DataFrame(pred_rows)
attn_df = pd.DataFrame(attn_rows)

acc  = accuracy_score(pred_df["true_label"], pred_df["pred_label"])
try:
    auc = roc_auc_score(pred_df["true_label"], pred_df["pred_prob"])
except Exception:
    auc = float("nan")
print(f"  Test donors : {len(pred_df)}")
print(f"  Accuracy : {acc:.4f}")
print(f"  ROC-AUC : {auc:.4f}")
print(f"\n  Classification report:")
print(classification_report(pred_df["true_label"], pred_df["pred_label"],
      target_names=["Healthy", "Diseased"]))

# Mean attention of cells
celltype_attn = (
    attn_df.groupby("cell_type")["alpha"]
    .mean()
    .sort_values(ascending=False)
    .reset_index()
    .rename(columns={"alpha": "mean_alpha"})
)
print(f"\n  Top 10 highest-attention cell types (most predictive of age):")
print(celltype_attn.head(10).to_string(index=False))



print("\n")
print("Plotting")

pred_df.to_csv(os.path.join(OUTPUT_DIR, "04_donor_age_predictions.csv"),    index=False)
attn_df.to_csv(os.path.join(OUTPUT_DIR, "04_cell_attention_weights.csv"),   index=False)
celltype_attn.to_csv(os.path.join(OUTPUT_DIR, "04_celltype_attention.csv"), index=False)
torch.save(best_state, os.path.join(OUTPUT_DIR, "04_mil_best_model.pt"))

# ── Plot 1: Training curve ────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(history["train"], label="Train MSE", color="#2166ac")
ax.plot(history["val"],   label="Val MSE",   color="#d62728")
ax.set_xlabel("Epoch"); ax.set_ylabel("MSE Loss (normalised age)")
ax.set_title("MIL Training Curve")
ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "04_training_curve.png"), dpi=200)
plt.close()

# ── Plot 2: ROC curve ────────────────────────────────────────
from sklearn.metrics import roc_curve
if not np.isnan(auc):
    fpr, tpr, _ = roc_curve(pred_df["true_label"], pred_df["pred_prob"])
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(fpr, tpr, color="#d62728", lw=2, label=f"ROC (AUC={auc:.3f})")
    ax.plot([0,1],[0,1],"k--", lw=0.8, alpha=0.4)
    ax.set_xlabel("False Positive Rate"); ax.set_ylabel("True Positive Rate")
    ax.set_title("ROC Curve — Healthy vs Diseased")
    ax.legend(); ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "04_roc_curve.png"), dpi=200)
    plt.close()

# ── Plot 2b: Predicted probability by true label ─────────────
fig, ax = plt.subplots(figsize=(6, 4))
for label, name, color in [(0,"Healthy","#2166ac"),(1,"Diseased","#d62728")]:
    sub = pred_df[pred_df["true_label"]==label]["pred_prob"]
    ax.hist(sub, bins=15, alpha=0.6, label=name, color=color)
ax.axvline(0.5, color="black", lw=1, linestyle="--", label="Decision boundary")
ax.set_xlabel("Predicted disease probability")
ax.set_ylabel("Number of donors")
ax.set_title("Predicted probability distribution\nHealthy vs Diseased")
ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "04_pred_prob_distribution.png"), dpi=200)
plt.close()

# ── Plot 3: Cell type attention ranking ──────────────────────
top_ct = celltype_attn.head(20)
fig, ax = plt.subplots(figsize=(8, 6))
ax.barh(top_ct["cell_type"][::-1], top_ct["mean_alpha"][::-1],
        color="#d62728", edgecolor="white", height=0.7)
ax.set_xlabel("Mean attention weight (α)")
ax.set_title("Top 20 cell types by MIL attention weight\n"
             "(most predictive of healthy vs diseased)")
ax.grid(axis="x", alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "04_celltype_attention_ranking.png"), dpi=200)
plt.close()

# ── Plot 4: Attention by condition (healthy vs diseased) ─────
top6_cts = celltype_attn.head(6)["cell_type"].tolist()
fig, axes = plt.subplots(2, 3, figsize=(12, 7))
for ax, ct in zip(axes.flatten(), top6_cts):
    sub = attn_df[attn_df["cell_type"] == ct]
    for label, name, color in [(0,"Healthy","#2166ac"),(1,"Diseased","#d62728")]:
        grp = sub[sub["true_label"]==label]["alpha"]
        ax.hist(grp, bins=15, alpha=0.5, label=name, color=color, density=True)
    ax.set_title(ct.replace("_"," "), fontsize=9, fontweight="bold")
    ax.set_xlabel("α weight", fontsize=8)
    ax.set_ylabel("Density", fontsize=8)
    ax.legend(fontsize=7)
    ax.grid(alpha=0.2)
plt.suptitle("Attention weight distribution: healthy vs diseased — top 6 cell types",
             fontsize=10, fontweight="bold")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "04_attention_by_condition.png"), dpi=200)
plt.close()

print(f"\n✓ All outputs saved to: {OUTPUT_DIR}/")
print(f"  04_donor_age_predictions.csv")
print(f"  04_cell_attention_weights.csv")
print(f"  04_celltype_attention.csv")
print(f"  04_mil_best_model.pt")
print(f"  04_training_curve.png")
print(f"  04_true_vs_pred_age.png")
print(f"  04_celltype_attention_ranking.png")
print(f"  04_attention_vs_age.png")
print(f"\n  Accuracy={acc:.3f} | AUC={auc:.3f}")




