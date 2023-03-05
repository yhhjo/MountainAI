import torch
import torch.nn.functional as F
import torch.nn as nn

MASKS = {
    "none": None,
    "future": "future"
}

class DotProductAttention(nn.Module):
    """Scaled dot product attention from d2l.ai. Assumes query, ley, and value dimensions are all equal.
    Parameters
    ----------
    dropout:
        dropout rate, default = 0.1
    
    """
    def __init__(self, dropout=0):
        """Dropout rate"""
        super().__init__()
        self._dropout = nn.Dropout(dropout)


    def forward(self, q: torch.Tensor, k: torch.Tensor, v:torch.Tensor, mask=None) -> torch.Tensor:
        """Computes scaled dot product attention given queries, keys, and values.
        Parameters
        ----------
        q:
            query tensor, shape (batch_size, d_input, d_model)
        k:
            key tensor, shape (batch_size, d_input, d_model)
        v:
            value tensor, shape (batch_size, d_input, d_model)
        mask: 
            str, default: "None"
        Returns
        ----------
        attention scores: 
            shape (batch_size, d_input, d_model)
        """
        K = q.shape[1]  # Sequence length

        dk = q.shape[1]**.5 
        I = (q@k.transpose(1,2))

        # Syntax from (Cohen et. al 2022)
        if mask == MASKS.get("future"):
            future_mask = torch.triu(torch.ones((K, K)), diagonal=1).bool()
            I = I.masked_fill(future_mask, float('-inf'))

        self._attention_weights = F.softmax(I/dk, dim=-1)
        attention = self._dropout(self._attention_weights)@v
        return attention