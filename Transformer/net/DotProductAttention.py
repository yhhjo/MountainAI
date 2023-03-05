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
        I = (q@k.transpose(1,2))/dk

        # Syntax from (Cohen et. al 2022)
        """
        By setting K to a fixed value, we can use a 2D tensor to represent the attention mask, 
        where the element at position (i, j) represents the attention weight between the i-th and j-th tokens. 
        If we want to mask out the attention from future tokens for the i-th token, 
        we can simply set the values of the i-th row of the attention mask to -inf. 
        This will ensure that the attention weights from the i-th token to all future tokens 
        (i.e., those with an index greater than i) will be set to zero after applying the softmax function, 
        effectively masking out these tokens. 
        """
        if mask == MASKS.get("future"):
            future_mask = torch.triu(torch.ones((K, K)), diagonal=1).bool()
            I = I.masked_fill(future_mask, float('-inf'))

        self._attention_weights = F.softmax(I, dim=-1)
        attention = self._dropout(self._attention_weights)@v
        return attention