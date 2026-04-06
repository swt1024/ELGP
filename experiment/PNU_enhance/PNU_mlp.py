# nnPU_MLP.py
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import roc_auc_score, average_precision_score


class MLP(nn.Module):
    """
    Simple MLP.
    Output is logit (no sigmoid here).
    """
    def __init__(self, in_dim: int, hidden_dims=(128, 64), dropout: float = 0.0):
        super().__init__()
        layers = []
        prev = in_dim
        for h in hidden_dims:
            layers.append(nn.Linear(prev, h))
            layers.append(nn.ReLU())
            if dropout and dropout > 0:
                layers.append(nn.Dropout(dropout))
            prev = h
        layers.append(nn.Linear(prev, 1))
        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


def logistic_loss(logits: torch.Tensor, t: torch.Tensor) -> torch.Tensor:
    """
    Logistic loss: log(1 + exp(-t*logit)) = softplus(-t*logit)
    t in {+1, -1}
    """
    return F.softplus(-t * logits)


def _safe_mean(x: torch.Tensor, anchor: torch.Tensor) -> torch.Tensor:
    """
    If x is empty, return 0 * anchor.sum() to keep the graph (no backward crash).
    """
    if x.numel() == 0:
        return 0.0 * anchor.sum()
    return x.mean()

def pairwise_ranking_loss(
    pos_logits: torch.Tensor,
    unl_logits: torch.Tensor,
    num_u_per_p: int = 32,
) -> torch.Tensor:
    """
    Pairwise logistic ranking loss (RankNet/BPR style):
      E_{p~P, u~U} [ softplus(-(s_p - s_u)) ]  ==  -E log sigmoid(s_p - s_u)

    We *sample* unlabeled examples for each positive to avoid O(|P||U|).
    If either set is empty, returns 0 with a valid graph.
    """
    pos_logits = pos_logits.view(-1)
    unl_logits = unl_logits.view(-1)

    # Keep graph even if empty
    anchor = (pos_logits.sum() if pos_logits.numel() > 0 else 0.0) + (unl_logits.sum() if unl_logits.numel() > 0 else 0.0)

    if pos_logits.numel() == 0 or unl_logits.numel() == 0:
        return 0.0 * anchor

    # sample U indices for each P
    m = pos_logits.shape[0]
    n_u = unl_logits.shape[0]
    k = int(min(max(num_u_per_p, 1), n_u))

    # idx shape: (m, k)
    idx = torch.randint(low=0, high=n_u, size=(m, k), device=unl_logits.device)
    u_samp = unl_logits[idx]                 # (m, k)
    p_rep = pos_logits.unsqueeze(1).expand(-1, k)  # (m, k)

    # softplus(-(p - u)) = softplus(u - p)
    loss = F.softplus(u_samp - p_rep)
    return loss.mean()

# ----------------------------------------------------------------------
# PNU (Sakai et al., 2017) - NEW
# ----------------------------------------------------------------------
def pn_risk(logits: torch.Tensor, z: torch.Tensor) -> torch.Tensor:
    """
    Supervised PN risk on labeled samples only.
    z: +1 (P), -1 (N), 0 (U). Here we only use P and N.
    R_PN = E_P[loss(+1)] + E_N[loss(-1)]
    """
    logits = logits.view(-1)
    z = z.view(-1)

    p_mask = (z == 1)
    n_mask = (z == -1)

    lp = _safe_mean(logistic_loss(logits[p_mask], torch.ones_like(logits[p_mask])), logits)
    ln = _safe_mean(logistic_loss(logits[n_mask], -torch.ones_like(logits[n_mask])), logits)
    return lp + ln


def pu_risk_unbiased(logits: torch.Tensor, z: torch.Tensor, pi: float, num_u_per_p: int = 32) -> torch.Tensor:
    """
    Unbiased PU risk estimator (same structure as the unbiased part in nnPU):
    R_PU = pi * E_P[loss(+1)] + E_U[loss(-1)] - pi * E_P[loss(-1)]
    z: +1 (P), 0 (U), -1 (N) (N ignored here)
    """
    logits = logits.view(-1)
    z = z.view(-1)

    p_mask = (z == 1)
    u_mask = (z == 0)

    logits_p = logits[p_mask]
    logits_u = logits[u_mask]

    Rp_pos = pi * _safe_mean(logistic_loss(logits_p, torch.ones_like(logits_p)), logits)
    Lu_neg = _safe_mean(logistic_loss(logits_u, -torch.ones_like(logits_u)), logits)
    #Ru_rank = pairwise_ranking_loss(logits_p, logits_u, num_u_per_p=num_u_per_p)
    Lp_neg = _safe_mean(logistic_loss(logits_p, -torch.ones_like(logits_p)), logits)
    #return Rp_pos + (Lu_neg - pi * Lp_neg)
    return Rp_pos + (Lu_neg - float(pi) * Lp_neg)

def nnpu_risk(
    logits: torch.Tensor,
    z: torch.Tensor,
    pi: float,
    beta: float = 0.0,
    num_u_per_p: int = 32
) -> torch.Tensor:
    """
    nnPU risk (Kiryo et al., 2017).
    s: +1 for positive (P), 0 for unlabeled (U)
    logits: shape (B,) or (B,1)
    pi: class prior P(y=+1)
    beta/gamma: non-negative correction params

    IMPORTANT:
    - Never return a constant tensor that breaks backward.
    - Never use tensor directly in Python if-condition (use .item() only for branching).
    """
    logits = logits.view(-1)
    z = z.view(-1)

    p_mask = (z == 1)
    u_mask = (z == 0)

    # safe mean: if empty, return 0 * logits.sum() so grad graph is kept
    def safe_mean(x: torch.Tensor) -> torch.Tensor:
        if x.numel() == 0:
            return 0.0 * logits.sum()
        return x.mean()

    logits_p = logits[p_mask]
    logits_u = logits[u_mask]

    # positive risk: pi * E_p[l(f(x), +1)]
    Rp = pi * safe_mean(logistic_loss(logits_p, torch.ones_like(logits_p)))

    # negative risk estimator:
    # Rn_hat = E_u[l(f(x), -1)] - pi * E_p[l(f(x), -1)]
    Lu_neg = safe_mean(logistic_loss(logits_u, -torch.ones_like(logits_u)))
    #Ru_rank = pairwise_ranking_loss(logits_p, logits_u, num_u_per_p=num_u_per_p)
    Lp_neg = safe_mean(logistic_loss(logits_p, -torch.ones_like(logits_p)))
    Rn_hat = Lu_neg - pi * Lp_neg
    #Rn_hat = Lu_neg - pi * Ru_rank

    # non-negative correction (Kiryo):
    # Use .item() only for branching; the loss tensor still keeps grad.
    if Rn_hat.item() < -beta:
        risk = Rp
    else:
        risk = Rp + Rn_hat + beta

    return risk


def pnpu_risk(
    logits: torch.Tensor,
    z: torch.Tensor,
    pi: float,
    gamma: float,
    beta: float = 0.0,
) -> torch.Tensor:
    """
    PNU risk (Sakai et al., 2017):
    R^γ_PNPU = (1-γ) R_PN + γ R_PU

    Inputs:
      z: +1(P), -1(N), 0(U)
      pi: class prior P(y=+1) in the *population/unlabeled* mixture
      gamma
    """
    Rpn = pn_risk(logits, z)
    Rpu = nnpu_risk(logits, z, pi=pi, beta=beta)
    #Rpu = pu_risk_unbiased(logits, z, pi=pi)
    return (1.0 - gamma) * Rpn + gamma * Rpu


# ----------------------------------------------------------------------
# Prediction / proxy eval (unchanged)
# ----------------------------------------------------------------------
@torch.no_grad()
def predict_proba(model: nn.Module, X: np.ndarray, device: str, batch_size: int = 4096) -> np.ndarray:
    """
    Return probabilities sigmoid(logit).
    """
    model.eval()
    probs = []
    n = X.shape[0]
    for i in range(0, n, batch_size):
        xb = torch.tensor(X[i:i + batch_size], dtype=torch.float32, device=device)
        logits = model(xb).squeeze(-1)
        p = torch.sigmoid(logits).cpu().numpy()
        probs.append(p)
    return np.concatenate(probs, axis=0)


def eval_proxy_metrics(scores: np.ndarray, s_obs: np.ndarray, topk_list=(500, 1000, 2000, 7000)):
    """
    PU proxy eval:
    - s_obs=1: known positive
    - s_obs=0: unlabeled
    """
    y = s_obs.astype(int)
    if len(np.unique(y)) > 1:
        roc = roc_auc_score(y, scores)
        pr = average_precision_score(y, scores)
    else:
        roc, pr = float("nan"), float("nan")

    order = np.argsort(-scores)
    topk_hits = {}
    for K in topk_list:
        K2 = min(K, len(scores))
        idx = order[:K2]
        topk_hits[K] = float(np.sum(y[idx] == 1)) / float(K2)

    return roc, pr, topk_hits
