import torch
import torch.nn.functional as F
import torch.nn as nn

class DotProductAttention(nn.Module):
    """Scaled dot product attention from d2l.ai. Assumes query, ley, and value dimensions are all equal.
    Parameters
    ----------
    dropout:
        dropout rate
    Returns
    ----------
    :
    """
    def __init__(self, dropout=0.1):
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
        Returns
        ----------
        attention scores: 
            shape (batch_size, d_input, d_model)
        """
        K = q.shape[1]

        rescale = q.shape[1]**.5 
        scores = (q@k.transpose(1,2)) / rescale

        # Cohen et. al (2022)
        ## TODO: Implement masking
        if mask == "next":
            future_mask = torch.triu(torch.ones((K, K)), diagonal=1).bool()
            self._scores = self._scores.masked_fill(future_mask, float('-inf'))

        self.attention_weights = F.softmax(scores, dim=-1)
        attention_scores = self._dropout(self.attention_weights)@v
        return attention_scores